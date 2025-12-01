#include <gl/freeglut.h>
#include <vector>
#include <cmath>
#include <cstdio>

struct Vec2 { float x, y; Vec2() {} Vec2(float X, float Y) :x(X), y(Y) {} };
struct Vec3 { float x, y, z; Vec3() {} Vec3(float X, float Y, float Z) :x(X), y(Y), z(Z) {} };
struct Vertex { Vec3 p; Vec3 v; float mass; bool fixed; };
struct Spring { int a, b; float rest; bool active; };

int N = 28;                 // grid resolution (NxN)
float spacing = 0.04f;      // spacing between particles
std::vector<Vertex> verts;
std::vector<Spring> springs;
std::vector<int> triIdx;    // not used for physics; optional for display as triangles

// simulation params
Vec3 gravity = Vec3(0, -9.8f, 0);
float k_spring = 350.0f;
float damping = 0.02f;
float timeStep = 0.008f;

// mouse cut
bool mouseDown = false;
int mx0, my0, mx1, my1;
int winW = 900, winH = 700;

// helpers
static inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline Vec3 operator*(const Vec3& a, float s) { return Vec3(a.x * s, a.y * s, a.z * s); }
static inline float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline float len(const Vec3& a) { return sqrtf(dot(a, a)); }

// project world coords to screen (orthographic) - simple mapping for interactions
Vec2 worldToScreen(const Vec3& p) {
    // our cloth lies in X-Y plane roughly centered near origin; map to window
    float scale = 800.0f;
    float cx = winW * 0.5f + p.x * scale;
    float cy = winH * 0.5f - p.y * scale;
    return Vec2(cx, cy);
}
Vec3 screenToWorld(int sx, int sy) {
    float scale = 800.0f;
    float x = (sx - winW * 0.5f) / scale;
    float y = (winH - sy - winH * 0.5f) / scale;
    return Vec3(x, y, 0);
}

// 2D segment-segment intersection (strict), returns true if intersects
bool segSegIntersect2D(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4) {
    auto orient = [](float ax, float ay, float bx, float by, float cx, float cy) {
        return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
        };
    float o1 = orient(x1, y1, x2, y2, x3, y3);
    float o2 = orient(x1, y1, x2, y2, x4, y4);
    float o3 = orient(x3, y3, x4, y4, x1, y1);
    float o4 = orient(x3, y3, x4, y4, x2, y2);
    return (o1 * o2 < 0.0f) && (o3 * o4 < 0.0f);
}

// build grid & springs (structural + shear + optional bend)
void initCloth() {
    verts.clear(); springs.clear(); triIdx.clear();
    // create NxN grid centered around origin
    float half = (N - 1) * spacing * 0.5f;
    for (int j = 0;j < N;j++) {
        for (int i = 0;i < N;i++) {
            Vertex V;
            V.p = Vec3((i * spacing) - half, ((N - 1 - j) * spacing) - 0.6f, 0.0f); // hang below origin
            V.v = Vec3(0, 0, 0);
            V.mass = 0.05f;
            V.fixed = (j == 0); // fix top corners
            verts.push_back(V);
        }
    }
    auto idx = [&](int i, int j) { return j * N + i; };
    // triangles for optional draw
    for (int j = 0;j < N - 1;j++) {
        for (int i = 0;i < N - 1;i++) {
            int a = idx(i, j), b = idx(i + 1, j), c = idx(i, j + 1), d = idx(i + 1, j + 1);
            triIdx.push_back(a); triIdx.push_back(c); triIdx.push_back(b);
            triIdx.push_back(b); triIdx.push_back(c); triIdx.push_back(d);
        }
    }
    // structural springs (grid edges)
    for (int j = 0;j < N;j++) {
        for (int i = 0;i < N;i++) {
            if (i < N - 1) {
                int a = idx(i, j), b = idx(i + 1, j);
                Spring s; s.a = a; s.b = b; s.rest = len(verts[a].p - verts[b].p); s.active = true; springs.push_back(s);
            }
            if (j < N - 1) {
                int a = idx(i, j), b = idx(i, j + 1);
                Spring s; s.a = a; s.b = b; s.rest = len(verts[a].p - verts[b].p); s.active = true; springs.push_back(s);
            }
            // shear springs
            if (i < N - 1 && j < N - 1) {
                { int a = idx(i, j), b = idx(i + 1, j + 1); Spring s = { a,b,len(verts[a].p - verts[b].p),true }; springs.push_back(s); }
                { int a = idx(i + 1, j), b = idx(i, j + 1); Spring s = { a,b,len(verts[a].p - verts[b].p),true }; springs.push_back(s); }
            }
            // optional bend springs (two-away) - helps stiffness; skip for speed
        }
    }
}

// physics: semi-implicit Euler
void stepPhysics(float dt) {
    int V = verts.size();
    std::vector<Vec3> forces(V, Vec3(0, 0, 0));
    // gravity
    for (int i = 0;i < V;i++) if (!verts[i].fixed) {
        forces[i] = forces[i] + gravity * verts[i].mass;
    }
    // springs
    for (auto& s : springs) if (s.active) {
        int a = s.a, b = s.b;
        Vec3 pa = verts[a].p, pb = verts[b].p;
        Vec3 diff = pa - pb;
        float L = len(diff);
        if (L > 1e-6f) {
            Vec3 dir = diff * (1.0f / L);
            float fs = -k_spring * (L - s.rest);
            // damping along spring
            Vec3 relVel = verts[a].v - verts[b].v;
            float fdamp = -damping * dot(relVel, dir);
            Vec3 F = dir * (fs + fdamp);
            if (!verts[a].fixed) forces[a] = forces[a] + F;
            if (!verts[b].fixed) forces[b] = forces[b] - F;
        }
    }
    // integrate velocities & positions
    for (int i = 0;i < V;i++) {
        if (verts[i].fixed) continue;
        Vec3 accel = forces[i] * (1.0f / verts[i].mass);
        // semi-implicit Euler
        verts[i].v = verts[i].v + accel * dt;
        // global damping
        verts[i].v = verts[i].v * 0.999f;
        verts[i].p = verts[i].p + verts[i].v * dt;
    }
}

// drawing helpers
void drawClothLines() {
    glLineWidth(1.2f);
    glBegin(GL_LINES);
    for (auto& s : springs) {
        if (!s.active) continue;
        Vec3 a = verts[s.a].p, b = verts[s.b].p;
        glColor3f(0.2f, 0.5f, 1.0f);
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
    }
    glEnd();
}
void drawVerts() {
    glPointSize(4.0f);
    glBegin(GL_POINTS);
    for (auto& v : verts) {
        if (v.fixed) glColor3f(1, 0.4f, 0.2f); else glColor3f(0.9f, 0.9f, 0.9f);
        glVertex3f(v.p.x, v.p.y, v.p.z);
    }
    glEnd();
}
void drawTrianglesWire() {
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0;i < triIdx.size();i += 3) {
        Vec3 a = verts[triIdx[i]].p;
        Vec3 b = verts[triIdx[i + 1]].p;
        Vec3 c = verts[triIdx[i + 2]].p;
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
        glVertex3f(c.x, c.y, c.z);
    }
    glEnd();
}

// process a cut segment from screen points (mx0,my0) -> (mx1,my1)
// removes any spring whose edge intersects the cut in screen-space projection
void processCutLine(int sx0, int sy0, int sx1, int sy1) {
    // map spring endpoints to screen coords
    for (auto& s : springs) {
        if (!s.active) continue;
        Vec2 A = worldToScreen(verts[s.a].p);
        Vec2 B = worldToScreen(verts[s.b].p);
        if (segSegIntersect2D(A.x, A.y, B.x, B.y, (float)sx0, (float)sy0, (float)sx1, (float)sy1)) {
            s.active = false; // cut the spring
        }
    }
}

// GLUT callbacks
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW); glLoadIdentity();
    // simple view transform
    glTranslatef(0.0f, -0.05f, -1.6f);
    glRotatef(12.0f, 1, 0, 0);

    // draw
    drawClothLines();
    drawVerts();
    // optionally triangle fill (visual artifact without remesh)
    // drawTrianglesWire();

    // draw cut line if dragging
    if (mouseDown) {
        glLineWidth(2.0f);
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
    // fixed-step substepping for stability
    const int steps = 2;
    for (int i = 0;i < steps;i++) stepPhysics(timeStep / steps);
    glutPostRedisplay();
}

// mouse handling: left click & drag draws cut line; on release process it
void mouseFunc(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            mouseDown = true; mx0 = x; my0 = y; mx1 = x; my1 = y;
        }
        else {
            mouseDown = false; mx1 = x; my1 = y;
            processCutLine(mx0, my0, mx1, my1);
        }
    }
}
void motionFunc(int x, int y) {
    if (mouseDown) { mx1 = x; my1 = y; }
}

// keyboard: r to reset
void keyboard(unsigned char key, int x, int y) {
    if (key == 'r' || key == 'R') { initCloth(); }
}

// setup GL
void setupGL() {
    glClearColor(0.08f, 0.08f, 0.08f, 1.0f);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glViewport(0, 0, winW, winH);
    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    float aspect = (float)winW / (float)winH;
    // orthographic-ish projection for easier mapping
    glOrtho(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect, 0.1, 10.0);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(winW, winH);
    glutCreateWindow("Cloth Cutting - MassSpring (left-drag to cut, R reset)");

    initCloth();
    setupGL();

    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);
    glutKeyboardFunc(keyboard);

    printf("Left-drag to cut springs. Press R to reset.\n");
    glutMainLoop();
    return 0;

}
