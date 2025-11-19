# Cloth_Cutting_Mesh_Manipulation

Real-time physics-based cloth cutting system with dynamic triangulation, topology updates, and mesh separation.

This project implements multiple **progressive versions** of cloth simulation and cutting.
Each version represents a step toward a fully functional, topologically correct, real-time cloth cutting engine.

---

# **Repository Structure**

```
/cloth-cutting-project
│
├── cloth_sim_basic_spring_cut.cpp
├── cloth_cutting_full_mesh_split.cpp
├── cloth_cutting_experimental_halfedge.cpp
└── cloth_cutting_debug_version_of_full_mesh_and_halfedge.cpp
└── README.md
```

---

# **Project Overview**

Modern cloth cutting requires more than deleting springs —
it must **dynamically update mesh topology** and **physically separate pieces**.

This project implements:

* Mass-spring cloth physics
* Screen-space cutting
* Edge intersection detection
* New vertex insertion
* Full re-triangulation
* BFS connected component detection
* Vertex duplication for separation
* Dynamic spring & edge rebuild

Each version improves or tests a part of this pipeline.

---

# **Version Overview**

---

## **1. cloth_sim_basic_spring_cut.cpp**

**Earliest version** — simple springs-only cutting.

### Features:

* Mass-spring cloth grid
* Basic physics simulation
* User cutting removes springs only
* NO mesh splitting
* NO retriangulation
* NO separation

### Use Case:

Baseline code for testing physics and rendering.

---

## **2. cloth_cutting_full_mesh_split.cpp**

**Main full-featured cutting engine.**

### Implements:

* Segment-edge intersection detection
* Insertion of new intersection vertices
* Triangle splitting for 0–3 edge cuts
* Adjacency graph build
* BFS connected-component detection
* Boundary vertex duplication
* Updating edges and springs
* Actual physical separation

This is the **correctly designed, production-ready** version.

---

## **3. cloth_cutting_experimental_halfedge.cpp**

Prototype of a **half-edge based** approach.

### Contains:

* Early half-edge structure
* Experiments with edge/vertex splitting
* Debug printing for topology tests
* Incomplete splitting logic

### Purpose:

Research version, not final.

---

## **4. cloth_cutting_debug_version_of_full_mesh_and_halfedge.cpp**

A **diagnostic build** used during development.

### Features:

* Verbose logging
* Checks for adjacency issues
* Prints components, vertices, triangles
* Tests cut-edge logic
* Helps identify separation bugs

### Purpose:

Debugging tool to validate correctness.
