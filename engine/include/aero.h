#ifndef AERO_H
#define AERO_H

#define MAX_PANELS  256
#define MAX_VERTS   64

typedef struct { double x, y, z; } Vec3;

/* ── Standard parametric wing ── */
typedef struct {
    double span, root_chord, tip_chord;
    double sweep_deg, dihedral_deg;
    int    nx, ny;
} WingGeometry;

/* ── Arbitrary polygon wing (free-draw mode) ── */
typedef struct {
    double verts_x[MAX_VERTS];  /* planform polygon x coords (chord dir) */
    double verts_y[MAX_VERTS];  /* planform polygon y coords (span dir)  */
    int    n_verts;             /* number of vertices                    */
    double dihedral_deg;        /* dihedral applied uniformly            */
    int    nx;                  /* chordwise panel divisions             */
    int    ny;                  /* spanwise panel divisions              */
} PolyWing;

/* ── Single VLM panel ── */
typedef struct {
    Vec3   corners[4];
    Vec3   control;
    Vec3   normal;
    double area;
} Panel;

/* ── Solver results ── */
typedef struct {
    double CL, CD, CM;
    double LD_ratio;
    double stall_angle;
    int    converged;
} AeroResult;

/* ── Function declarations ── */
void       generate_mesh        (const WingGeometry *wing,
                                  Panel *panels, int total_panels);
void       generate_poly_mesh   (const PolyWing *pw,
                                  Panel *panels, int *out_n);
AeroResult run_vlm_solver       (const WingGeometry *wing,
                                  double alpha_deg, double velocity);
AeroResult run_poly_vlm_solver  (const PolyWing *pw,
                                  double alpha_deg, double velocity);
double     compute_lift         (double CL, double rho, double V, double S);
double     compute_drag         (double CD, double rho, double V, double S);

#endif