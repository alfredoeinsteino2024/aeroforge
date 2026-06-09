#include <math.h>
#include <string.h>
#include <stdio.h>
#include "aero.h"

#define PI         3.14159265358979323846
#define DEG2RAD    (PI / 180.0)
#define MAX_PANELS 256

static Panel  s_panels[MAX_PANELS];
static double s_AIC[MAX_PANELS][MAX_PANELS];
static double s_RHS[MAX_PANELS];
static double s_gamma[MAX_PANELS];

/* ── Mesh generator ─────────────────────────────────────────── */
void generate_mesh(const WingGeometry *wing, Panel *panels, int total_panels) {
    int    nx        = wing->nx;
    int    ny        = wing->ny;
    double hs        = wing->span / 2.0;
    double sweep_rad = wing->sweep_deg    * DEG2RAD;
    double dih_rad   = wing->dihedral_deg * DEG2RAD;

    for (int j = 0; j < ny; j++) {
        double eta0 = (double) j      / ny;
        double eta1 = (double)(j + 1) / ny;
        double c0   = wing->root_chord + (wing->tip_chord - wing->root_chord) * eta0;
        double c1   = wing->root_chord + (wing->tip_chord - wing->root_chord) * eta1;
        double xle0 = eta0 * hs * tan(sweep_rad);
        double xle1 = eta1 * hs * tan(sweep_rad);
        double y0   = eta0 * hs;
        double y1   = eta1 * hs;
        double z0   = eta0 * hs * tan(dih_rad);
        double z1   = eta1 * hs * tan(dih_rad);

        for (int i = 0; i < nx; i++) {
            double xi0 = (double) i      / nx;
            double xi1 = (double)(i + 1) / nx;
            int    idx = j * nx + i;
            Panel *p   = &panels[idx];

            /* corners: 0=root-LE, 1=root-TE, 2=tip-TE, 3=tip-LE */
            p->corners[0] = (Vec3){ xle0 + xi0*c0,  y0, z0 };
            p->corners[1] = (Vec3){ xle0 + xi1*c0,  y0, z0 };
            p->corners[2] = (Vec3){ xle1 + xi1*c1,  y1, z1 };
            p->corners[3] = (Vec3){ xle1 + xi0*c1,  y1, z1 };

            /* control point at 3/4 chord midspan */
            p->control = (Vec3){
                0.5*(xle0+xle1) + 0.75*0.5*(c0+c1),
                0.5*(y0+y1),
                0.5*(z0+z1)
            };
            p->normal = (Vec3){ 0.0, 0.0, 1.0 };
            p->area   = 0.5*(c0+c1) * (y1-y0) / nx;
        }
    }
}

/* ── Biot-Savart for ONE finite vortex segment A→B ─────────────
   Returns induced velocity at P with unit circulation.
   Formula: V = (1/4π) * (r1×r2)/|r1×r2|² * (r0·r1/|r1| - r0·r2/|r2|)
   where r0=B-A, r1=P-A, r2=P-B                                  */
static Vec3 biot_savart(Vec3 P, Vec3 A, Vec3 B) {
    Vec3 r1 = { P.x-A.x, P.y-A.y, P.z-A.z };
    Vec3 r2 = { P.x-B.x, P.y-B.y, P.z-B.z };
    Vec3 r0 = { B.x-A.x, B.y-A.y, B.z-A.z };

    /* cross product r1 × r2 */
    Vec3 cr = {
        r1.y*r2.z - r1.z*r2.y,
        r1.z*r2.x - r1.x*r2.z,
        r1.x*r2.y - r1.y*r2.x
    };
    double cr2 = cr.x*cr.x + cr.y*cr.y + cr.z*cr.z;

    double r1m = sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z);
    double r2m = sqrt(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);

    if (cr2 < 1e-12 || r1m < 1e-12 || r2m < 1e-12)
        return (Vec3){0,0,0};

    double d1 = (r0.x*r1.x + r0.y*r1.y + r0.z*r1.z) / r1m;
    double d2 = (r0.x*r2.x + r0.y*r2.y + r0.z*r2.z) / r2m;
    double K  = (d1 - d2) / (4.0 * PI * cr2);

    return (Vec3){ K*cr.x, K*cr.y, K*cr.z };
}

/* ── Semi-infinite vortex from A going in direction +x ─────────
   Trailing vortex goes from A downstream to +infinity.
   Induced velocity at P:
   V = (1/4π) * (e×r1) / (|r1| * (|r1| - e·r1))
   where e = unit vector in free-stream direction (+x here)       */
static Vec3 semi_infinite(Vec3 P, Vec3 A, double sign) {
    Vec3 r1  = { P.x-A.x, P.y-A.y, P.z-A.z };
    double r1m = sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z);

    /* e = (1,0,0) so e×r1 = (0*r1z - 0*r1y, 0*r1x - 1*r1z, 1*r1y - 0*r1x)
                            = (0, -r1z, r1y)                                 */
    double edotr = r1.x;   /* e·r1 = r1x since e=(1,0,0) */
    double denom = r1m * (r1m - edotr);

    if (fabs(denom) < 1e-12 || r1m < 1e-12)
        return (Vec3){0,0,0};

    double K = sign / (4.0 * PI * denom);
    return (Vec3){ 0.0, K*(-r1.z), K*(r1.y) };
}

/* ── Horseshoe vortex influence coefficient ─────────────────────
   Horseshoe = bound leg (A→B) + right trailing (B→+inf)
                                + left trailing  (+inf→A)
   Bound leg sits at quarter-chord of source panel.
   A = root end, B = tip end (increasing y).                      */
static double horseshoe_w(Vec3 cp, Panel *src) {
    /* quarter-chord line of source panel */
    Vec3 A = {
        src->corners[0].x + 0.25*(src->corners[1].x - src->corners[0].x),
        src->corners[0].y,
        src->corners[0].z
    };
    Vec3 B = {
        src->corners[3].x + 0.25*(src->corners[2].x - src->corners[3].x),
        src->corners[3].y,
        src->corners[3].z
    };

    /* bound vortex A→B */
    Vec3 vb = biot_savart(cp, A, B);

    /* right trailing: B → +x infinity  (sign = +1) */
    Vec3 vr = semi_infinite(cp, B, +1.0);

    /* left trailing: A → +x infinity, but this leg has opposite
       circulation so it runs from +inf→A, equivalent to -(A→+inf) */
    Vec3 vl = semi_infinite(cp, A, -1.0);

    /* normal component (normal = +z for flat wing) */
    double wz = vb.z + vr.z + vl.z;
    return wz;
}

/* ── Main VLM solver ────────────────────────────────────────── */
AeroResult run_vlm_solver(const WingGeometry *wing,
                           double alpha_deg, double velocity) {
    AeroResult result;
    memset(&result, 0, sizeof(result));

    int N = wing->nx * wing->ny;
    if (N > MAX_PANELS) {
        printf("[VLM] %d panels exceeds MAX %d\n", N, MAX_PANELS);
        return result;
    }

    generate_mesh(wing, s_panels, N);

    double alpha = alpha_deg * DEG2RAD;

    /* RHS = downwash needed to cancel free-stream normal component
       For flat wing (normal = +z): w_needed = -V*sin(alpha) ≈ -V*alpha */
    for (int i = 0; i < N; i++)
        s_RHS[i] = -velocity * sin(alpha);

    /* Build AIC matrix */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            s_AIC[i][j] = horseshoe_w(s_panels[i].control, &s_panels[j]);

    /* Gaussian elimination with partial pivoting */
    for (int col = 0; col < N; col++) {
        /* find pivot */
        int pivot = col;
        for (int row = col+1; row < N; row++)
            if (fabs(s_AIC[row][col]) > fabs(s_AIC[pivot][col]))
                pivot = row;

        /* swap rows */
        if (pivot != col) {
            for (int k = 0; k < N; k++) {
                double t = s_AIC[col][k];
                s_AIC[col][k]   = s_AIC[pivot][k];
                s_AIC[pivot][k] = t;
            }
            double t    = s_RHS[col];
            s_RHS[col]  = s_RHS[pivot];
            s_RHS[pivot]= t;
        }

        if (fabs(s_AIC[col][col]) < 1e-14) continue;

        for (int row = col+1; row < N; row++) {
            double f = s_AIC[row][col] / s_AIC[col][col];
            for (int k = col; k < N; k++)
                s_AIC[row][k] -= f * s_AIC[col][k];
            s_RHS[row] -= f * s_RHS[col];
        }
    }

    /* back substitution */
    for (int i = N-1; i >= 0; i--) {
        s_gamma[i] = s_RHS[i];
        for (int j = i+1; j < N; j++)
            s_gamma[i] -= s_AIC[i][j] * s_gamma[j];
        if (fabs(s_AIC[i][i]) > 1e-14)
            s_gamma[i] /= s_AIC[i][i];
    }

    /* Kutta-Joukowski: L = rho * V * sum(gamma * dy) */
    double rho = 1.225;
    double L   = 0.0;
    double S   = 0.0;

    for (int i = 0; i < N; i++) {
        double dy = fabs(s_panels[i].corners[3].y - s_panels[i].corners[0].y);
        L += rho * velocity * s_gamma[i] * dy;
        S += s_panels[i].area;
    }

    /* multiply by 2 for full wing (we solved half span) */
    L *= 2.0;
    S *= 2.0;

    double q  = 0.5 * rho * velocity * velocity;
    double AR = (wing->span * wing->span) / S;

    result.CL        = L / (q * S);
    result.CD        = (result.CL * result.CL) / (PI * AR);
    result.CM        = -0.25 * result.CL;
    result.LD_ratio  = (fabs(result.CD) > 1e-9)
                       ? result.CL / result.CD : 0.0;
    result.stall_angle = 15.0;
    result.converged = 1;

    return result;
}

double compute_lift(double CL, double rho, double V, double S) {
    return 0.5 * rho * V * V * S * CL;
}

double compute_drag(double CD, double rho, double V, double S) {
    return 0.5 * rho * V * V * S * CD;
}