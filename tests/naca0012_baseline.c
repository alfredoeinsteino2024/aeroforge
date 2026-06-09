#include <stdio.h>
#include <math.h>
#include "../engine/include/aero.h"

int main(void) {

    printf("Starting validation...\n");
    fflush(stdout);

    WingGeometry wing;
    wing.span         = 10.0;
    wing.root_chord   = 1.5;
    wing.tip_chord    = 1.5;
    wing.sweep_deg    = 0.0;
    wing.dihedral_deg = 0.0;
    wing.nx           = 1;   /* simplest possible: 1 chordwise panel */
    wing.ny           = 4;   /* 4 spanwise panels                    */

    int N = wing.nx * wing.ny;
    printf("Panel count: %d\n", N);
    fflush(stdout);

    /* ── single AoA debug run at 5 degrees ── */
    double alpha_deg = 5.0;
    double velocity  = 50.0;
    double alpha_rad = alpha_deg * 3.14159265 / 180.0;

    printf("\nDebug at AoA=5 deg, V=50 m/s\n");
    printf("Expected CL ~ %.4f  (2*pi*alpha for thin airfoil)\n",
           2.0 * 3.14159265 * alpha_rad);
    fflush(stdout);

    AeroResult r = run_vlm_solver(&wing, alpha_deg, velocity);

    printf("Got CL      = %.6f\n", r.CL);
    printf("Got CD      = %.6f\n", r.CD);
    printf("Got L/D     = %.6f\n", r.LD_ratio);
    printf("Converged   = %d\n",   r.converged);
    fflush(stdout);

    /* ── sweep table ── */
    printf("\nAoA(deg) |   CL    |   CD    |  L/D\n");
    printf("---------|---------|---------|------\n");

    for (int a = 0; a <= 12; a += 2) {
        AeroResult res = run_vlm_solver(&wing, (double)a, velocity);
        printf("  %3d deg  | %7.4f | %7.4f | %6.2f\n",
               a, res.CL, res.CD, res.LD_ratio);
        fflush(stdout);
    }

    printf("\nDone.\n");
    return 0;
}