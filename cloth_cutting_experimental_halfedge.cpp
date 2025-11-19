#include <gl/freeglut.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <queue>
#include <set>
#include <map>

struct Vec2 { float x, y; Vec2() {} Vec2(float X, float Y) :x(X), y(Y) {} };
struct Vec3 { float x, y, z; Vec3() {} Vec3(float X, float Y, float Z) :x(X), y(Y), z(Z) {} };
struct Vertex { Vec3 p; Vec3 v; float mass; bool fixed; };
struct Triangle { int a, b, c; bool active; };
struct Edge { int v1, v2; bool active; };
struct Spring { int a, b; float rest; bool active; };

// Global state
int N = 20;
float spacing = 0.05f;
std::vector<Vertex> verts;
std::vector<Triangle> triangles;
std::vector<Spring> springs;
std::vector<Edge> edges;

// simulation params
Vec3 gravity = Vec3(0, -9.8f, 0);
float k_spring = 700.0f;
float damping = 0.2f;
float timeStep = 0.006f;
float separationOffset = 0.15f;

// mouse cut
bool mouseDown = false;
int mx0, my0, mx1, my1;
int winW = 900, winH = 700;

// operators
static inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline Vec3 operator*(const Vec3& a, float s) { return Vec3(a.x * s, a.y * s, a.z * s); }
static inline float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline float len(const Vec3& a) { return sqrtf(dot(a, a)); }
static inline Vec3 normalize(const Vec3& a) { float L = len(a); return L > 1e-6f ? a * (1.0f / L) : Vec3(0, 0, 0); }

Vec2 worldToScreen(const Vec3& p) {
    float scale = 800.0f;
    return Vec2(winW * 0.5f + p.x * scale, winH * 0.5f - p.y * scale);
}

Vec3 screenToWorld(int sx, int sy) {
    float scale = 800.0f;
    return Vec3((sx - winW * 0.5f) / scale, (winH - sy - winH * 0.5f) / scale, 0);
}

bool segSegIntersect2D(float x1, float y1, float x2, float y2,
    float x3, float y3, float x4, float y4, float* t = nullptr) {
    float denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (fabsf(denom) < 1e-6f) return false;

    float t1 = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
    float t2 = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / denom;

    if (t) *t = t1;
    return (t1 > 0.01f && t1 < 0.99f && t2 > 0.01f && t2 < 0.99f);
}

struct IntersectionPoint {
    int triIdx;
    int edgeInTri;
    int va, vb;
    float t;
    Vec3 pos;
};

Vec3 computeCutDirection(int sx0, int sy0, int sx1, int sy1) {
    Vec3 w0 = screenToWorld(sx0, sy0);
    Vec3 w1 = screenToWorld(sx1, sy1);
    Vec3 cutDir = w1 - w0;
    float L = len(cutDir);
    if (L > 1e-6f) cutDir = cutDir * (1.0f / L);
    Vec3 perpDir = Vec3(-cutDir.y, cutDir.x, 0);
    return normalize(perpDir);
}

std::vector<IntersectionPoint> detectIntersections(int sx0, int sy0, int sx1, int sy1) {
    std::vector<IntersectionPoint> intersections;

    for (int triIdx = 0; triIdx < (int)triangles.size(); triIdx++) {
        if (!triangles[triIdx].active) continue;

        Triangle& tri = triangles[triIdx];
        int verts_tri[3] = { tri.a, tri.b, tri.c };

        for (int e = 0; e < 3; e++) {
            int va = verts_tri[e];
            int vb = verts_tri[(e + 1) % 3];

            Vec2 sa = worldToScreen(verts[va].p);
            Vec2 sb = worldToScreen(verts[vb].p);

            float t;
            if (segSegIntersect2D(sa.x, sa.y, sb.x, sb.y,
                (float)sx0, (float)sy0, (float)sx1, (float)sy1, &t)) {
                IntersectionPoint ip;
                ip.triIdx = triIdx;
                ip.edgeInTri = e;
                ip.va = va;
                ip.vb = vb;
                ip.t = t;
                ip.pos = verts[va].p + (verts[vb].p - verts[va].p) * t;
                intersections.push_back(ip);
            }
        }
    }

    return intersections;
}

// FIX 1: Build adjacency EXCLUDING cut edges
void buildTriangleAdjacency(std::vector<std::vector<int>>& adjList,
                            const std::set<std::pair<int,int>>& cutEdges,
                            const std::set<int>& newVerts)
{
    int numTris = triangles.size();
    adjList.assign(numTris, std::vector<int>());

    std::map<std::pair<int,int>, std::vector<int>> edgeToTris;

    for (int i = 0; i < numTris; i++) {
        if (!triangles[i].active) continue;

        Triangle& tri = triangles[i];
        int v[3] = {tri.a, tri.b, tri.c};

        for (int e = 0; e < 3; e++) {
            int v1 = v[e];
            int v2 = v[(e+1)%3];

            int a = std::min(v1, v2);
            int b = std::max(v1, v2);
            std::pair<int,int> edge(a,b);

            // NEW FIX: If either vertex is newly created, skip this edge
            if (newVerts.count(v1) || newVerts.count(v2))
                continue;

            // OLD CHECK: Skip exact cut edges
            if (cutEdges.count(edge))
                continue;

            edgeToTris[edge].push_back(i);
        }
    }

    for (auto& entry : edgeToTris) {
        auto& tris = entry.second;
        for (int i = 0; i < (int)tris.size(); i++)
            for (int j = i+1; j < (int)tris.size(); j++) {
                adjList[tris[i]].push_back(tris[j]);
                adjList[tris[j]].push_back(tris[i]);
            }
    }
}

std::vector<int> findConnectedComponents(const std::vector<std::vector<int>>& adjList) {
    int numTris = triangles.size();
    std::vector<int> component(numTris, -1);
    int currentComponent = 0;

    for (int start = 0; start < numTris; start++) {
        if (!triangles[start].active || component[start] >= 0) continue;

        std::queue<int> q;
        q.push(start);
        component[start] = currentComponent;

        while (!q.empty()) {
            int curr = q.front();
            q.pop();

            for (int neighbor : adjList[curr]) {
                if (triangles[neighbor].active && component[neighbor] < 0) {
                    component[neighbor] = currentComponent;
                    q.push(neighbor);
                }
            }
        }

        currentComponent++;
    }

    return component;
}

void performCut(int sx0, int sy0, int sx1, int sy1) {
    printf("\n=== CUT OPERATION ===\n");

    auto intersections = detectIntersections(sx0, sy0, sx1, sy1);
    if (intersections.empty()) {
        printf("No intersections\n");
        return;
    }
    printf("Intersections: %d\n", (int)intersections.size());

    // Track cut edges and split vertices
    std::set<std::pair<int, int>> cutEdges;
    std::map<std::pair<int, int>, int> edgeToNewVert;
    std::set<int> newVerts;

    for (const auto& ip : intersections) {
        int va = std::min(ip.va, ip.vb);
        int vb = std::max(ip.va, ip.vb);
        std::pair<int, int> edge = std::make_pair(va, vb);

        if (edgeToNewVert.count(edge)) continue;

        // Create new vertex at intersection
        Vertex newV;
        newV.p = ip.pos;
        newV.v = verts[ip.va].v + (verts[ip.vb].v - verts[ip.va].v) * ip.t;
        newV.mass = 0.05f;
        newV.fixed = false;

        int newIdx = verts.size();
        verts.push_back(newV);
        edgeToNewVert[edge] = newIdx;
        cutEdges.insert(edge);
        newVerts.insert(newIdx);
    }

    printf("New vertices: %d\n", (int)edgeToNewVert.size());
    printf("Cut edges: %d\n", (int)cutEdges.size());

    // Re-triangulate affected triangles
    std::set<int> affectedTris;
    for (const auto& ip : intersections) {
        affectedTris.insert(ip.triIdx);
    }

    std::vector<Triangle> newTriangles;

    for (int triIdx : affectedTris) {
        Triangle& tri = triangles[triIdx];
        tri.active = false;

        int va = tri.a, vb = tri.b, vc = tri.c;

        auto getNewVert = [&](int v1, int v2) -> int {
            int a = std::min(v1, v2);
            int b = std::max(v1, v2);
            std::map<std::pair<int, int>, int>::iterator it = edgeToNewVert.find(std::make_pair(a, b));
            return (it != edgeToNewVert.end()) ? it->second : -1;
            };

        int vab = getNewVert(va, vb);
        int vbc = getNewVert(vb, vc);
        int vca = getNewVert(vc, va);

        int splitCount = (vab >= 0 ? 1 : 0) + (vbc >= 0 ? 1 : 0) + (vca >= 0 ? 1 : 0);

        if (splitCount == 0) {
            tri.active = true;
        }
        else if (splitCount == 1) {
            if (vab >= 0) {
                newTriangles.push_back({ va, vab, vc, true });
                newTriangles.push_back({ vab, vb, vc, true });
            }
            else if (vbc >= 0) {
                newTriangles.push_back({ va, vb, vbc, true });
                newTriangles.push_back({ va, vbc, vc, true });
            }
            else {
                newTriangles.push_back({ va, vb, vca, true });
                newTriangles.push_back({ vca, vb, vc, true });
            }
        }
        else if (splitCount == 2) {
            if (vab >= 0 && vbc >= 0) {
                newTriangles.push_back({ va, vab, vbc, true });
                newTriangles.push_back({ vab, vb, vbc, true });
                newTriangles.push_back({ va, vbc, vc, true });
            }
            else if (vbc >= 0 && vca >= 0) {
                newTriangles.push_back({ va, vb, vbc, true });
                newTriangles.push_back({ vbc, vc, vca, true });
                newTriangles.push_back({ va, vbc, vca, true });
            }
            else {
                newTriangles.push_back({ va, vab, vca, true });
                newTriangles.push_back({ vab, vb, vc, true });
                newTriangles.push_back({ vab, vc, vca, true });
            }
        }
        else if (splitCount == 3) {
            newTriangles.push_back({ va, vab, vca, true });
            newTriangles.push_back({ vb, vbc, vab, true });
            newTriangles.push_back({ vc, vca, vbc, true });
            newTriangles.push_back({ vab, vbc, vca, true });
        }
    }

    triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
    printf("New triangles: %d\n", (int)newTriangles.size());

    // Build adjacency EXCLUDING cut edges
    std::vector<std::vector<int>> adjList;
    buildTriangleAdjacency(adjList, cutEdges, newVerts);

    // Find connected components
    std::vector<int> component = findConnectedComponents(adjList);

    // Count actual components
    int maxComp = -1;
    for (size_t i = 0; i < component.size(); i++) {
        if (triangles[i].active && component[i] > maxComp) {
            maxComp = component[i];
        }
    }
    int numComponents = maxComp + 1;

    printf("Components found: %d\n", numComponents);

    // FIX 2: Duplicate vertices AFTER BFS, based on actual components
    if (numComponents > 1) {
        Vec3 cutDir = computeCutDirection(sx0, sy0, sx1, sy1);

        // Find which components each vertex belongs to
        std::map<int, std::set<int>> vertToComps;
        for (size_t i = 0; i < triangles.size(); i++) {
            if (!triangles[i].active || component[i] < 0) continue;

            vertToComps[triangles[i].a].insert(component[i]);
            vertToComps[triangles[i].b].insert(component[i]);
            vertToComps[triangles[i].c].insert(component[i]);
        }

        // Duplicate boundary vertices (vertices belonging to multiple components)
        std::map<std::pair<int, int>, int> vertCompToNew;

        for (std::map<int, std::set<int>>::iterator it = vertToComps.begin();
            it != vertToComps.end(); ++it) {
            int vertIdx = it->first;
            const std::set<int>& comps = it->second;

            if (comps.size() <= 1) continue;

            std::vector<int> compVec(comps.begin(), comps.end());

            // Original stays with first component
            vertCompToNew[std::make_pair(vertIdx, compVec[0])] = vertIdx;

            // Duplicate for other components with offset
            for (size_t c = 1; c < compVec.size(); c++) {
                Vertex newV = verts[vertIdx];
                float sign = (c % 2 == 0) ? 1.0f : -1.0f;
                newV.p = newV.p + cutDir * (separationOffset * sign);

                int newIdx = verts.size();
                verts.push_back(newV);
                vertCompToNew[std::make_pair(vertIdx, compVec[c])] = newIdx;
            }
        }

        printf("Duplicated %d boundary vertices\n", (int)(vertCompToNew.size() - vertToComps.size()));

        // Reassign triangle vertices based on their component
        for (size_t i = 0; i < triangles.size(); i++) {
            if (!triangles[i].active || component[i] < 0) continue;

            int comp = component[i];

            auto remap = [&](int v) -> int {
                std::map<std::pair<int, int>, int>::iterator it =
                    vertCompToNew.find(std::make_pair(v, comp));
                return (it != vertCompToNew.end()) ? it->second : v;
                };

            triangles[i].a = remap(triangles[i].a);
            triangles[i].b = remap(triangles[i].b);
            triangles[i].c = remap(triangles[i].c);
        }
    }

    // Rebuild edges
    edges.clear();
    std::set<std::pair<int, int>> edgeSet;

    for (const auto& tri : triangles) {
        if (!tri.active) continue;

        auto addEdge = [&](int v1, int v2) {
            int a = std::min(v1, v2);
            int b = std::max(v1, v2);
            edgeSet.insert(std::make_pair(a, b));
            };

        addEdge(tri.a, tri.b);
        addEdge(tri.b, tri.c);
        addEdge(tri.c, tri.a);
    }

    for (std::set<std::pair<int, int>>::iterator it = edgeSet.begin();
        it != edgeSet.end(); ++it) {
        edges.push_back({ it->first, it->second, true });
    }

    // Rebuild springs (springs are automatically excluded from cut edges 
    // because we rebuilt edges from triangles, which now use duplicated vertices)
    springs.clear();
    for (const auto& edge : edges) {
        if (!edge.active) continue;
        if (edge.v1 >= (int)verts.size() || edge.v2 >= (int)verts.size()) continue;

        Spring s;
        s.a = edge.v1;
        s.b = edge.v2;
        s.rest = len(verts[edge.v1].p - verts[edge.v2].p);
        s.active = true;
        springs.push_back(s);
    }

    printf("Final: Verts=%d, Edges=%d, Springs=%d\n",
        (int)verts.size(), (int)edges.size(), (int)springs.size());
    printf("=== CUT COMPLETE ===\n\n");
}

void initCloth() {
    verts.clear();
    triangles.clear();
    springs.clear();
    edges.clear();

    float half = (N - 1) * spacing * 0.5f;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            Vertex V;
            V.p = Vec3((i * spacing) - half, ((N - 1 - j) * spacing) - 0.6f, 0.0f);
            V.v = Vec3(0, 0, 0);
            V.mass = 0.05f;
            V.fixed = (j == 0);
            verts.push_back(V);
        }
    }

    auto idx = [](int i, int j) { return j * N + i; };

    for (int j = 0; j < N - 1; j++) {
        for (int i = 0; i < N - 1; i++) {
            int a = idx(i, j), b = idx(i + 1, j), c = idx(i, j + 1), d = idx(i + 1, j + 1);
            triangles.push_back({ a, c, b, true });
            triangles.push_back({ b, c, d, true });
        }
    }

    std::set<std::pair<int, int>> edgeSet;
    for (const auto& tri : triangles) {
        auto addEdge = [&](int v1, int v2) {
            int a = std::min(v1, v2), b = std::max(v1, v2);
            edgeSet.insert(std::make_pair(a, b));
            };
        addEdge(tri.a, tri.b);
        addEdge(tri.b, tri.c);
        addEdge(tri.c, tri.a);
    }

    for (std::set<std::pair<int, int>>::iterator it = edgeSet.begin();
        it != edgeSet.end(); ++it) {
        edges.push_back({ it->first, it->second, true });
    }

    for (const auto& edge : edges) {
        Spring s;
        s.a = edge.v1;
        s.b = edge.v2;
        s.rest = len(verts[edge.v1].p - verts[edge.v2].p);
        s.active = true;
        springs.push_back(s);
    }

    printf("\n=== CLOTH INITIALIZED ===\n");
    printf("Vertices: %d | Triangles: %d | Springs: %d\n\n",
        (int)verts.size(), (int)triangles.size(), (int)springs.size());
}

void stepPhysics(float dt) {
    int V = verts.size();
    std::vector<Vec3> forces(V, Vec3(0, 0, 0));

    for (int i = 0; i < V; i++) {
        if (!verts[i].fixed) {
            forces[i] = forces[i] + gravity * verts[i].mass;
        }
    }

    for (const auto& s : springs) {
        if (!s.active) continue;
        int a = s.a, b = s.b;
        if (a >= V || b >= V) continue;

        Vec3 diff = verts[a].p - verts[b].p;
        float L = len(diff);
        if (L > 1e-6f) {
            Vec3 dir = diff * (1.0f / L);
            float fs = -k_spring * (L - s.rest);
            Vec3 relVel = verts[a].v - verts[b].v;
            float fdamp = -damping * dot(relVel, dir);
            Vec3 F = dir * (fs + fdamp);
            if (!verts[a].fixed) forces[a] = forces[a] + F;
            if (!verts[b].fixed) forces[b] = forces[b] - F;
        }
    }

    for (int i = 0; i < V; i++) {
        if (verts[i].fixed) continue;
        Vec3 accel = forces[i] * (1.0f / verts[i].mass);
        verts[i].v = verts[i].v + accel * dt;
        verts[i].v = verts[i].v * 0.998f;
        verts[i].p = verts[i].p + verts[i].v * dt;
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0f, -0.05f, -1.6f);
    glRotatef(12.0f, 1, 0, 0);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glColor3f(0.3f, 0.45f, 0.65f);
    glBegin(GL_TRIANGLES);
    for (const auto& tri : triangles) {
        if (!tri.active) continue;
        if (tri.a >= (int)verts.size() || tri.b >= (int)verts.size() || tri.c >= (int)verts.size()) continue;
        Vec3 a = verts[tri.a].p;
        Vec3 b = verts[tri.b].p;
        Vec3 c = verts[tri.c].p;
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
        glVertex3f(c.x, c.y, c.z);
    }
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);

    glLineWidth(1.0f);
    glColor3f(0.6f, 0.7f, 0.8f);
    glBegin(GL_LINES);
    for (const auto& e : edges) {
        if (!e.active) continue;
        if (e.v1 >= (int)verts.size() || e.v2 >= (int)verts.size()) continue;
        Vec3 a = verts[e.v1].p;
        Vec3 b = verts[e.v2].p;
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
    }
    glEnd();

    glPointSize(3.0f);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < verts.size(); i++) {
        if (verts[i].fixed) glColor3f(1, 0.4f, 0.2f);
        else glColor3f(0.9f, 0.9f, 0.9f);
        glVertex3f(verts[i].p.x, verts[i].p.y, verts[i].p.z);
    }
    glEnd();

    if (mouseDown) {
        glLineWidth(3.0f);
        glColor3f(1, 0.2f, 0.2f);
        glBegin(GL_LINES);
        Vec3 w0 = screenToWorld(mx0, my0);
        Vec3 w1 = screenToWorld(mx1, my1);
        glVertex3f(w0.x, w0.y, w0.z);
        glVertex3f(w1.x, w1.y, w1.z);
        glEnd();
    }

    glutSwapBuffers();
}

void idle() {
    const int steps = 3;
    for (int i = 0; i < steps; i++) {
        stepPhysics(timeStep / steps);
    }
    glutPostRedisplay();
}

void mouseFunc(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            mouseDown = true;
            mx0 = x; my0 = y;
            mx1 = x; my1 = y;
        }
        else {
            mouseDown = false;
            mx1 = x; my1 = y;
            performCut(mx0, my0, mx1, my1);
        }
    }
}

void motionFunc(int x, int y) {
    if (mouseDown) {
        mx1 = x; my1 = y;
    }
}

void keyboard(unsigned char key, int x, int y) {
    if (key == 'r' || key == 'R') {
        initCloth();
    }
}

void setupGL() {
    glClearColor(0.08f, 0.08f, 0.08f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glViewport(0, 0, winW, winH);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = (float)winW / (float)winH;
    glOrtho(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect, 0.1, 10.0);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(winW, winH);
    glutCreateWindow("CORRECTLY FIXED Cloth Cutting");

    initCloth();
    setupGL();

    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);
    glutKeyboardFunc(keyboard);

    printf("Controls:\n");
    printf("  Left-drag = Cut cloth\n");
    printf("  R = Reset\n\n");
    glutMainLoop();
    return 0;
}p