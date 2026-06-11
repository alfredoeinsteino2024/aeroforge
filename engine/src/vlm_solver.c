#include <math.h>
#include <string.h>
#include <stdio.h>
#include "aero.h"

#define PI      3.14159265358979323846
#define D2R     (PI/180.0)

/* ── Static memory (never on the stack) ── */
static Panel  s_panels[MAX_PANELS];
static double s_AIC[MAX_PANELS][MAX_PANELS];
static double s_RHS[MAX_PANELS];
static double s_gamma[MAX_PANELS];

/* ════════════════════════════════════════
   MESH — parametric trapezoidal wing
   ════════════════════════════════════════ */
void generate_mesh(const WingGeometry *w, Panel *panels, int total) {
    int nx=w->nx, ny=w->ny;
    double hs=w->span/2.0;
    double sw=tan(w->sweep_deg*D2R), dh=tan(w->dihedral_deg*D2R);

    for (int j=0;j<ny;j++) {
        double e0=(double)j/ny, e1=(double)(j+1)/ny;
        double c0=w->root_chord+(w->tip_chord-w->root_chord)*e0;
        double c1=w->root_chord+(w->tip_chord-w->root_chord)*e1;
        double xle0=e0*hs*sw, xle1=e1*hs*sw;
        double y0=e0*hs, y1=e1*hs;
        double z0=e0*hs*dh, z1=e1*hs*dh;

        for (int i=0;i<nx;i++) {
            double xi0=(double)i/nx, xi1=(double)(i+1)/nx;
            Panel *p=&panels[j*nx+i];
            p->corners[0]=(Vec3){xle0+xi0*c0, z0, y0};
            p->corners[1]=(Vec3){xle0+xi1*c0, z0, y0};
            p->corners[2]=(Vec3){xle1+xi1*c1, z1, y1};
            p->corners[3]=(Vec3){xle1+xi0*c1, z1, y1};
            p->control   =(Vec3){0.5*(xle0+xle1)+0.75*0.5*(c0+c1),
                                  0.5*(z0+z1), 0.5*(y0+y1)};
            p->normal    =(Vec3){0,1,0};
            p->area      =0.5*(c0+c1)*(y1-y0)/nx;
        }
    }
}

/* ════════════════════════════════════════
   MESH — arbitrary polygon wing
   Takes a 2D polygon (x=chord, y=span)
   and meshes it into nx×ny panels by
   sampling spanwise slices through the
   polygon and finding chord intersections.
   ════════════════════════════════════════ */
void generate_poly_mesh(const PolyWing *pw, Panel *panels, int *out_n) {
    int nx=pw->nx, ny=pw->ny;
    double dh=tan(pw->dihedral_deg*D2R);
    int count=0;

    /* find spanwise extents of polygon */
    double ymin=pw->verts_y[0], ymax=pw->verts_y[0];
    for (int i=1;i<pw->n_verts;i++) {
        if (pw->verts_y[i]<ymin) ymin=pw->verts_y[i];
        if (pw->verts_y[i]>ymax) ymax=pw->verts_y[i];
    }
    double hs=(ymax-ymin)/2.0;

    for (int j=0;j<ny&&count+nx<=MAX_PANELS;j++) {
        double eta0=(double)j/ny, eta1=(double)(j+1)/ny;
        double y0=ymin+(ymax-ymin)*eta0;
        double y1=ymin+(ymax-ymin)*eta1;
        double ymid=0.5*(y0+y1);

        /* find chord extent at ymid by ray-casting through polygon */
        double xmin_s=1e9, xmax_s=-1e9;
        int nv=pw->n_verts;
        for (int i=0;i<nv;i++) {
            int i2=(i+1)%nv;
            double ya=pw->verts_y[i], yb=pw->verts_y[i2];
            double xa=pw->verts_x[i], xb=pw->verts_x[i2];
            if ((ya<=ymid&&yb>ymid)||(yb<=ymid&&ya>ymid)) {
                double t=(ymid-ya)/(yb-ya);
                double xi=xa+t*(xb-xa);
                if (xi<xmin_s) xmin_s=xi;
                if (xi>xmax_s) xmax_s=xi;
            }
        }
        if (xmin_s>=xmax_s) continue; /* degenerate slice */

        /* similarly for y0 and y1 */
        double xmin0=1e9,xmax0=-1e9,xmin1=1e9,xmax1=-1e9;
        for (int i=0;i<nv;i++) {
            int i2=(i+1)%nv;
            double ya=pw->verts_y[i],yb=pw->verts_y[i2];
            double xa=pw->verts_x[i],xb=pw->verts_x[i2];
            /* y0 */
            if ((ya<=y0&&yb>y0)||(yb<=y0&&ya>y0)) {
                double t=(y0-ya)/(yb-ya);
                double xi=xa+t*(xb-xa);
                if (xi<xmin0) xmin0=xi;
                if (xi>xmax0) xmax0=xi;
            }
            /* y1 */
            if ((ya<=y1&&yb>y1)||(yb<=y1&&ya>y1)) {
                double t=(y1-ya)/(yb-ya);
                double xi=xa+t*(xb-xa);
                if (xi<xmin1) xmin1=xi;
                if (xi>xmax1) xmax1=xi;
            }
        }
        if (xmin0>=xmax0||xmin1>=xmax1) continue;

        double c0=xmax0-xmin0, c1=xmax1-xmin1;
        double z0=(y0-ymin-hs)*dh, z1=(y1-ymin-hs)*dh;

        for (int i=0;i<nx&&count<MAX_PANELS;i++) {
            double xi0=(double)i/nx, xi1=(double)(i+1)/nx;
            Panel *p=&panels[count++];
            p->corners[0]=(Vec3){xmin0+xi0*c0, z0, y0};
            p->corners[1]=(Vec3){xmin0+xi1*c0, z0, y0};
            p->corners[2]=(Vec3){xmin1+xi1*c1, z1, y1};
            p->corners[3]=(Vec3){xmin1+xi0*c1, z1, y1};
            p->control   =(Vec3){0.5*(xmin0+xmin1)+0.75*0.5*(c0+c1),
                                  0.5*(z0+z1), 0.5*(y0+y1)};
            p->normal    =(Vec3){0,1,0};
            p->area      =0.5*(c0+c1)*(y1-y0)/nx;
        }
    }
    *out_n=count;
}

/* ════════════════════════════════════════
   BIOT-SAVART + HORSESHOE VORTEX
   ════════════════════════════════════════ */
static Vec3 biot_savart(Vec3 P, Vec3 A, Vec3 B) {
    Vec3 r1={P.x-A.x,P.y-A.y,P.z-A.z};
    Vec3 r2={P.x-B.x,P.y-B.y,P.z-B.z};
    Vec3 r0={B.x-A.x,B.y-A.y,B.z-A.z};
    Vec3 cr={r1.y*r2.z-r1.z*r2.y,
             r1.z*r2.x-r1.x*r2.z,
             r1.x*r2.y-r1.y*r2.x};
    double cr2=cr.x*cr.x+cr.y*cr.y+cr.z*cr.z;
    double r1m=sqrt(r1.x*r1.x+r1.y*r1.y+r1.z*r1.z);
    double r2m=sqrt(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z);
    if (cr2<1e-12||r1m<1e-12||r2m<1e-12) return (Vec3){0,0,0};
    double d1=(r0.x*r1.x+r0.y*r1.y+r0.z*r1.z)/r1m;
    double d2=(r0.x*r2.x+r0.y*r2.y+r0.z*r2.z)/r2m;
    double K=(d1-d2)/(4.0*PI*cr2);
    return (Vec3){K*cr.x,K*cr.y,K*cr.z};
}

static Vec3 semi_inf(Vec3 P, Vec3 A, double sign) {
    Vec3 r1={P.x-A.x,P.y-A.y,P.z-A.z};
    double r1m=sqrt(r1.x*r1.x+r1.y*r1.y+r1.z*r1.z);
    double edotr=r1.x;
    double denom=r1m*(r1m-edotr);
    if (fabs(denom)<1e-12||r1m<1e-12) return (Vec3){0,0,0};
    double K=sign/(4.0*PI*denom);
    return (Vec3){0,K*(-r1.z),K*r1.y};
}

static double horseshoe_w(Vec3 cp, Panel *src) {
    Vec3 A={src->corners[0].x+0.25*(src->corners[1].x-src->corners[0].x),
            src->corners[0].y, src->corners[0].z};
    Vec3 B={src->corners[3].x+0.25*(src->corners[2].x-src->corners[3].x),
            src->corners[3].y, src->corners[3].z};
    Vec3 vb=biot_savart(cp,A,B);
    Vec3 vr=semi_inf(cp,B,+1.0);
    Vec3 vl=semi_inf(cp,A,-1.0);
    return vb.y+vr.y+vl.y;  /* normal component (y=up for flat wing) */
}

/* ════════════════════════════════════════
   CORE VLM SOLVE  (shared by both modes)
   ════════════════════════════════════════ */
static AeroResult vlm_solve(Panel *panels, int N,
                             double alpha_deg, double velocity) {
    AeroResult res; memset(&res,0,sizeof(res));
    if (N<1||N>MAX_PANELS) return res;

    double alpha=alpha_deg*D2R;

    for (int i=0;i<N;i++) s_RHS[i]=-velocity*sin(alpha);

    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            s_AIC[i][j]=horseshoe_w(panels[i].control,&panels[j]);

    /* Gaussian elimination */
    for (int col=0;col<N;col++) {
        int piv=col;
        for (int r=col+1;r<N;r++)
            if (fabs(s_AIC[r][col])>fabs(s_AIC[piv][col])) piv=r;
        if (piv!=col) {
            double t;
            for (int k=0;k<N;k++){t=s_AIC[col][k];s_AIC[col][k]=s_AIC[piv][k];s_AIC[piv][k]=t;}
            t=s_RHS[col];s_RHS[col]=s_RHS[piv];s_RHS[piv]=t;
        }
        if (fabs(s_AIC[col][col])<1e-14) continue;
        for (int r=col+1;r<N;r++) {
            double f=s_AIC[r][col]/s_AIC[col][col];
            for (int k=col;k<N;k++) s_AIC[r][k]-=f*s_AIC[col][k];
            s_RHS[r]-=f*s_RHS[col];
        }
    }
    for (int i=N-1;i>=0;i--) {
        s_gamma[i]=s_RHS[i];
        for (int j=i+1;j<N;j++) s_gamma[i]-=s_AIC[i][j]*s_gamma[j];
        if (fabs(s_AIC[i][i])>1e-14) s_gamma[i]/=s_AIC[i][i];
    }

    double rho=1.225, L=0, S=0;
    for (int i=0;i<N;i++) {
        double dz=fabs(panels[i].corners[3].z-panels[i].corners[0].z);
        L+=rho*velocity*s_gamma[i]*dz;
        S+=panels[i].area;
    }
    L*=2.0; S*=2.0;

    double q=0.5*rho*velocity*velocity;
    double AR=(S>1e-9)?(4.0*(panels[N-1].corners[3].z*panels[N-1].corners[3].z)/S):6.0;

    res.CL       =(q*S>1e-9)?L/(q*S):0;
    res.CD       =(res.CL*res.CL)/(PI*AR);
    res.CM       =-0.25*res.CL;
    res.LD_ratio =(res.CD>1e-9)?res.CL/res.CD:0;
    res.stall_angle=15.0;
    res.converged=1;
    return res;
}

/* ════════════════════════════════════════
   PUBLIC API
   ════════════════════════════════════════ */
AeroResult run_vlm_solver(const WingGeometry *w,
                           double alpha_deg, double vel) {
    generate_mesh(w, s_panels, w->nx*w->ny);
    return vlm_solve(s_panels, w->nx*w->ny, alpha_deg, vel);
}

AeroResult run_poly_vlm_solver(const PolyWing *pw,
                                double alpha_deg, double vel) {
    int N=0;
    generate_poly_mesh(pw, s_panels, &N);
    return vlm_solve(s_panels, N, alpha_deg, vel);
}

double compute_lift(double CL,double rho,double V,double S){return 0.5*rho*V*V*S*CL;}
double compute_drag(double CD,double rho,double V,double S){return 0.5*rho*V*V*S*CD;}