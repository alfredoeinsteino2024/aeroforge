#ifndef AERO_H
#define AERO_H

/* ─── Basic 3D vector ─────────────────────────────── */
typedef struct {
    double x, y, z;
} Vec3;

/* ─── A single panel on the wing surface ─────────────
   VLM divides the wing into a grid of panels.
   Each panel has 4 corner points and a control point
   where we enforce the no-penetration condition.      */
typedef struct {
    Vec3 corners[4];   /* panel vertices (A, B, C, D)  */
    Vec3 control;      /* 3/4 chord control point       */
    Vec3 normal;       /* unit normal vector            */
    double area;       /* panel area                    */
} Panel;

/* ─── Wing geometry input ─────────────────────────── */
typedef struct {
    double span;        /* total wingspan (m)           */
    double root_chord;  /* chord at root (m)            */
    double tip_chord;   /* chord at tip (m)             */
    double sweep_deg;   /* leading edge sweep (degrees) */
    double dihedral_deg;/* dihedral angle (degrees)     */
    int    nx;          /* panels along chord           */
    int    ny;          /* panels along span            */
} WingGeometry;

/* ─── Aerodynamic results ─────────────────────────── */
typedef struct {
    double CL;          /* lift coefficient             */
    double CD;          /* induced drag coefficient     */
    double CM;          /* pitching moment coefficient  */
    double LD_ratio;    /* lift-to-drag ratio           */
    double stall_angle; /* estimated stall AoA (deg)    */
    int    converged;   /* 1 = solver converged, 0 = no */
} AeroResult;

/* ─── Function declarations ──────────────────────── */
void       generate_mesh   (const WingGeometry *wing, Panel *panels, int total_panels);
AeroResult run_vlm_solver  (const WingGeometry *wing, double alpha_deg, double velocity);
double     compute_lift    (double CL, double rho, double V, double S);
double     compute_drag    (double CD, double rho, double V, double S);

#endif /* AERO_H */