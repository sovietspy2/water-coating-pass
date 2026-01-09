#include "pass.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef PASS_DEBUG
#define PASS_DEBUG 0
#endif

#if PASS_DEBUG
#define DBG_PRINTF(...) printf(__VA_ARGS__)
#else
#define DBG_PRINTF(...) ((void)0)
#endif

/* Hard cap a szomszédok számára: megfogja a k^2 robbanást sűrű régiókban */
#ifndef NEIGH_CAP
#define NEIGH_CAP 96
#endif

static vec3 v3(double x, double y, double z) { vec3 r = {x,y,z}; return r; }

static vec3 v3_add(vec3 a, vec3 b) { return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static vec3 v3_sub(vec3 a, vec3 b) { return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static vec3 v3_scale(vec3 a, double s) { return v3(a.x*s, a.y*s, a.z*s); }

static double v3_dot(vec3 a, vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

static vec3 v3_cross(vec3 a, vec3 b) {
    return v3(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

static double v3_norm2(vec3 a) { return v3_dot(a,a); }
static double v3_norm(vec3 a) { return sqrt(v3_norm2(a)); }

static double v3_dist(vec3 a, vec3 b) { return v3_norm(v3_sub(a,b)); }
static double v3_dist2(vec3 a, vec3 b) { return v3_norm2(v3_sub(a,b)); }

static int v3_normalize(vec3 in, vec3 *out) {
    double n = v3_norm(in);
    if (n < EPS) return 0;
    *out = v3_scale(in, 1.0/n);
    return 1;
}

/* --- Atom sugarak (σ_i) becslése ---
 * Megjegyzés: itt minimális elemtérképet használunk (H, C, N, O, S).
 * Ha ismeretlen, konzervatív alapértéket adunk vissza.
 */
static double sigma_from_atom(const ap *a) {
    /* A PDB atomnév (pdb_type) tipikusan pl. " C  ", "CA", " O  " stb.
       Itt egyszerűen megkeressük az első betűt, ami A-Z vagy a-z. */
    char e = '\0';
    for (size_t t = 0; t < 4; t++) {
        char c = a->pdb_type[t];
        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) { e = c; break; }
    }
    if (e >= 'a' && e <= 'z') e = (char)(e - 32);

    /* PASS táblázat szerinti gyakori sugarak (Å), egyszerűsítve. */
    switch (e) {
        case 'H': return 1.20;
        case 'O': return 1.52;
        case 'N': return 1.55;
        case 'C': return 1.70;
        case 'S': return 1.80;
        default:  return 1.70; /* fallback */
    }
}

/* Kinyerjük egy atom 3D koordinátáját vec3-ként. */
static vec3 pos_from_atom(const ap *a) {
    return v3(a->x_coord, a->y_coord, a->z_coord);
}

/* --- Uniform grid (cell list) --- */
typedef struct {
    vec3 min, max;
    double cell;
    int nx, ny, nz;
    int *head; /* grid cell -> first index */
    int *next; /* index -> next index in cell */

    const vec3 *pos;
    const double *sigma;
    size_t n;

    double max_sigma;
} grid3d;

static int clampi(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static int grid_index(const grid3d *g, int ix, int iy, int iz) {
    return (iz * g->ny + iy) * g->nx + ix;
}

static void grid_free(grid3d *g) {
    if (!g) return;
    free(g->head);
    free(g->next);
    g->head = NULL;
    g->next = NULL;
}

static int grid_build(grid3d *g, const vec3 *pos, const double *sigma, size_t n, double sigma_p) {
    if (!g || !pos || !sigma || n == 0) return 0;
    memset(g, 0, sizeof(*g));

    g->pos = pos;
    g->sigma = sigma;
    g->n = n;

    vec3 mn = pos[0];
    vec3 mx = pos[0];
    double max_sigma = sigma[0];

    for (size_t i = 0; i < n; i++) {
        vec3 p = pos[i];
        if (p.x < mn.x) mn.x = p.x; if (p.x > mx.x) mx.x = p.x;
        if (p.y < mn.y) mn.y = p.y; if (p.y > mx.y) mx.y = p.y;
        if (p.z < mn.z) mn.z = p.z; if (p.z > mx.z) mx.z = p.z;
        if (sigma[i] > max_sigma) max_sigma = sigma[i];
    }

    g->min = mn;
    g->max = mx;
    g->max_sigma = max_sigma;

    /* cell size: a tipikus “interakciós sugár” nagyságrendje */
    g->cell = (max_sigma + sigma_p);
    if (g->cell < 1e-6) g->cell = 1.0;

    g->nx = (int)floor((g->max.x - g->min.x) / g->cell) + 1;
    g->ny = (int)floor((g->max.y - g->min.y) / g->cell) + 1;
    g->nz = (int)floor((g->max.z - g->min.z) / g->cell) + 1;
    if (g->nx < 1) g->nx = 1;
    if (g->ny < 1) g->ny = 1;
    if (g->nz < 1) g->nz = 1;

    size_t grid_sz = (size_t)g->nx * (size_t)g->ny * (size_t)g->nz;

    g->head = (int*)malloc(grid_sz * sizeof(int));
    g->next = (int*)malloc(n * sizeof(int));
    if (!g->head || !g->next) { grid_free(g); return 0; }

    for (size_t i = 0; i < grid_sz; i++) g->head[i] = -1;

    for (size_t i = 0; i < n; i++) {
        vec3 p = pos[i];
        int ix = (int)floor((p.x - g->min.x) / g->cell);
        int iy = (int)floor((p.y - g->min.y) / g->cell);
        int iz = (int)floor((p.z - g->min.z) / g->cell);
        ix = clampi(ix, 0, g->nx - 1);
        iy = clampi(iy, 0, g->ny - 1);
        iz = clampi(iz, 0, g->nz - 1);

        int gi = grid_index(g, ix, iy, iz);
        g->next[i] = g->head[gi];
        g->head[gi] = (int)i;
    }

    return 1;
}

/* Bridgeable dist²-vel (sqrt nélkül) */
static int pair_bridgeable(vec3 A, double sigma_A, vec3 B, double sigma_B, double sigma_p) {
    double RA = sigma_A + sigma_p;
    double RB = sigma_B + sigma_p;

    double d2 = v3_dist2(A, B);
    double sum = RA + RB;
    double diff = fabs(RA - RB);

    if (d2 > (sum * sum)) return 0;
    if (d2 < (diff * diff)) return 0;
    return 1;
}

/* --- K-nearest jellegű szomszéd limitálás (egyszerű, O(n*cap)) --- */
static void neigh_try_add(
    int *neigh,
    double *neigh_d2,
    int *nnb_io,
    int cap,
    int cand,
    double cand_d2
) {
    int nnb = *nnb_io;
    if (nnb < cap) {
        neigh[nnb] = cand;
        neigh_d2[nnb] = cand_d2;
        *nnb_io = nnb + 1;
        return;
    }

    /* keressük a legrosszabbat (legtávolabbit) a tartott listában */
    int worst_i = 0;
    double worst_d2 = neigh_d2[0];
    for (int i = 1; i < cap; i++) {
        if (neigh_d2[i] > worst_d2) {
            worst_d2 = neigh_d2[i];
            worst_i = i;
        }
    }

    if (cand_d2 < worst_d2) {
        neigh[worst_i] = cand;
        neigh_d2[worst_i] = cand_d2;
    }
}

/* --- egyszerű dinamikus vec3 lista --- */
typedef struct {
    vec3 *data;
    size_t n;
    size_t cap;
} vec3_list;

static void v3list_free(vec3_list *l) {
    if (!l) return;
    free(l->data);
    l->data = NULL;
    l->n = l->cap = 0;
}

static int v3list_push(vec3_list *l, vec3 p) {
    if (l->n == l->cap) {
        size_t new_cap = (l->cap == 0) ? 1024 : l->cap * 2;
        vec3 *tmp = (vec3*)realloc(l->data, new_cap * sizeof(vec3));
        if (!tmp) return 0;
        l->data = tmp;
        l->cap = new_cap;
    }
    l->data[l->n++] = p;
    return 1;
}

/* Ütközés grid-del: csak közeli cellák */
static int clashes_with_any_grid(vec3 p, const grid3d *g, double sigma_p) {
    if (!g || !g->head || !g->next) return 0;

    int ix = (int)floor((p.x - g->min.x) / g->cell);
    int iy = (int)floor((p.y - g->min.y) / g->cell);
    int iz = (int)floor((p.z - g->min.z) / g->cell);
    ix = clampi(ix, 0, g->nx - 1);
    iy = clampi(iy, 0, g->ny - 1);
    iz = clampi(iz, 0, g->nz - 1);

    for (int dz = -1; dz <= 1; dz++) {
        int zz = iz + dz; if (zz < 0 || zz >= g->nz) continue;
        for (int dy = -1; dy <= 1; dy++) {
            int yy = iy + dy; if (yy < 0 || yy >= g->ny) continue;
            for (int dx = -1; dx <= 1; dx++) {
                int xx = ix + dx; if (xx < 0 || xx >= g->nx) continue;

                int gi = grid_index(g, xx, yy, zz);
                for (int ai = g->head[gi]; ai != -1; ai = g->next[ai]) {
                    double min_d = g->sigma[ai] + sigma_p;
                    double lim = (min_d - CLASH_TOL);
                    if (v3_dist2(p, g->pos[ai]) < (lim * lim)) return 1;
                }
            }
        }
    }
    return 0;
}

/* New-point weeding grid-del */
typedef struct {
    vec3 min;
    double cell;
    int nx, ny, nz;
    int *head;
    int *next;
    vec3 *pts;
    size_t n;
    size_t cap;
} grid_pts;

static void grid_pts_free(grid_pts *gp) {
    if (!gp) return;
    free(gp->head);
    free(gp->next);
    free(gp->pts);
    gp->head = NULL;
    gp->next = NULL;
    gp->pts  = NULL;
    gp->n = gp->cap = 0;
}

static int grid_pts_init(grid_pts *gp, vec3 mn, vec3 mx, double cell, size_t cap_pts) {
    if (!gp) return 0;
    memset(gp, 0, sizeof(*gp));
    gp->min = mn;
    gp->cell = (cell < 1e-6) ? 1.0 : cell;

    gp->nx = (int)floor((mx.x - mn.x) / gp->cell) + 1;
    gp->ny = (int)floor((mx.y - mn.y) / gp->cell) + 1;
    gp->nz = (int)floor((mx.z - mn.z) / gp->cell) + 1;
    if (gp->nx < 1) gp->nx = 1;
    if (gp->ny < 1) gp->ny = 1;
    if (gp->nz < 1) gp->nz = 1;

    size_t grid_sz = (size_t)gp->nx * (size_t)gp->ny * (size_t)gp->nz;
    gp->head = (int*)malloc(grid_sz * sizeof(int));
    gp->next = (int*)malloc(cap_pts * sizeof(int));
    gp->pts  = (vec3*)malloc(cap_pts * sizeof(vec3));
    if (!gp->head || !gp->next || !gp->pts) { grid_pts_free(gp); return 0; }

    for (size_t i = 0; i < grid_sz; i++) gp->head[i] = -1;
    gp->n = 0;
    gp->cap = cap_pts;
    return 1;
}

static int grid_pts_index(const grid_pts *gp, int ix, int iy, int iz) {
    return (iz * gp->ny + iy) * gp->nx + ix;
}

static int grid_pts_ensure_capacity(grid_pts *gp) {
    if (gp->n < gp->cap) return 1;

    size_t new_cap = (gp->cap == 0) ? 1024 : gp->cap * 2;
    int *new_next = (int*)realloc(gp->next, new_cap * sizeof(int));
    vec3 *new_pts = (vec3*)realloc(gp->pts,  new_cap * sizeof(vec3));
    if (!new_next || !new_pts) {
        /* ha egyik realloc fail, a másik lehet már sikerült; kezeljük óvatosan */
        if (new_next) gp->next = new_next;
        if (new_pts)  gp->pts  = new_pts;
        return 0;
    }
    gp->next = new_next;
    gp->pts  = new_pts;
    gp->cap  = new_cap;
    return 1;
}

static int too_close_to_new_grid(vec3 p, const grid_pts *gp, double weed_dist) {
    double wd2 = weed_dist * weed_dist;

    int ix = (int)floor((p.x - gp->min.x) / gp->cell);
    int iy = (int)floor((p.y - gp->min.y) / gp->cell);
    int iz = (int)floor((p.z - gp->min.z) / gp->cell);
    ix = clampi(ix, 0, gp->nx - 1);
    iy = clampi(iy, 0, gp->ny - 1);
    iz = clampi(iz, 0, gp->nz - 1);

    for (int dz = -1; dz <= 1; dz++) {
        int zz = iz + dz; if (zz < 0 || zz >= gp->nz) continue;
        for (int dy = -1; dy <= 1; dy++) {
            int yy = iy + dy; if (yy < 0 || yy >= gp->ny) continue;
            for (int dx = -1; dx <= 1; dx++) {
                int xx = ix + dx; if (xx < 0 || xx >= gp->nx) continue;

                int gi = grid_pts_index(gp, xx, yy, zz);
                for (int pi = gp->head[gi]; pi != -1; pi = gp->next[pi]) {
                    if (v3_dist2(p, gp->pts[pi]) < wd2) return 1;
                }
            }
        }
    }
    return 0;
}

static int grid_pts_add(grid_pts *gp, vec3 p) {
    if (!grid_pts_ensure_capacity(gp)) return 0;

    int ix = (int)floor((p.x - gp->min.x) / gp->cell);
    int iy = (int)floor((p.y - gp->min.y) / gp->cell);
    int iz = (int)floor((p.z - gp->min.z) / gp->cell);
    ix = clampi(ix, 0, gp->nx - 1);
    iy = clampi(iy, 0, gp->ny - 1);
    iz = clampi(iz, 0, gp->nz - 1);

    int gi = grid_pts_index(gp, ix, iy, iz);
    int idx = (int)gp->n;

    gp->pts[idx] = p;
    gp->next[idx] = gp->head[gi];
    gp->head[gi] = idx;
    gp->n++;
    return 1;
}

/* Új HOH oxigén atom kitöltése. */
static void make_hoh_oxygen(
    ap *out,
    int model_ser,
    int atom_ser,
    int res_ser,
    char chain_id,
    vec3 p
) {
    memset(out, 0, sizeof(*out));

    out->model_ser = model_ser;
    strncpy(out->pdb_token, "HETATM", sizeof(out->pdb_token) - 1);
    out->atom_ser = atom_ser;
    strncpy(out->pdb_type, "O", sizeof(out->pdb_type) - 1);
    strncpy(out->pdb_alt_loc, "", sizeof(out->pdb_alt_loc) - 1);
    strncpy(out->res_type, "H2O", sizeof(out->res_type) - 1);

    out->chain[0] = chain_id;
    out->chain[1] = '\0';

    out->res_ser = res_ser;

    strncpy(out->pdb_achar, "", sizeof(out->pdb_achar) - 1);

    out->x_coord = p.x;
    out->y_coord = p.y;
    out->z_coord = p.z;

    out->occ = 1.0;
    out->b_factor = 0.0;

    strncpy(out->mol2_type, "", sizeof(out->mol2_type) - 1);
    out->mol2_charge = 0.0;
    strncpy(out->pdbqt_type, "", sizeof(out->pdbqt_type) - 1);
}

/* -------- Appendix A szerinti számítás -------- */
int three_point_sphere_geometry(
    vec3 i, double sigma_i,
    vec3 j, double sigma_j,
    vec3 k, double sigma_k,
    double sigma_p,
    vec3 *p_plus,
    vec3 *p_minus,
    double *U_out,
    double *V_out,
    double *h_out
) {
    /* (1) Kiterjesztett sugarak: R_i = σ_i + σ_p stb. */
    double R_i = sigma_i + sigma_p;
    double R_j = sigma_j + sigma_p;
    double R_k = sigma_k + sigma_p;

    /* (2) Lokális bázis felépítése: x' i->j irány, y' i-j-k síkban, z' merőleges. */
    vec3 ij = v3_sub(j, i);
    double d_ij = v3_norm(ij);
    if (d_ij < EPS) return 0;

    vec3 ex;
    if (!v3_normalize(ij, &ex)) return 0; /* x' egységvektor */

    vec3 ik = v3_sub(k, i);

    /* k koordinátája az x' tengelyen: x_k = (ik · ex) */
    double x_k = v3_dot(ik, ex);

    /* y' irány: ik-ből kivesszük az x' komponensét */
    vec3 ik_perp = v3_sub(ik, v3_scale(ex, x_k));
    double y_k = v3_norm(ik_perp);
    if (y_k < EPS) return 0; /* i, j, k közel kollineáris -> instabil */

    vec3 ey;
    if (!v3_normalize(ik_perp, &ey)) return 0; /* y' egységvektor */

    vec3 ez = v3_cross(ex, ey); /* z' (jobbkéz-szabály) */

    /* (3) U és V számítása (trilateration az i-j-k síkon)
       Az Appendix A ábrájának megfelelően p = (U, V, h) a lokális bázisban. */

    /* U az i->j tengely menti hely: */
    double U = (d_ij*d_ij + R_i*R_i - R_j*R_j) / (2.0*d_ij);

    /* k lokális koordinátái: (x_k, y_k, 0)
       V a síkon belüli második koordináta: */
    double d_ik2 = v3_norm2(ik); /* x_k^2 + y_k^2 ugyanaz, csak stabilan */
    double V = (d_ik2 + R_i*R_i - R_k*R_k - 2.0*x_k*U) / (2.0*y_k);

    /* (4) h számítása Pithagorasz-szerűen az i középpontú gömbből:
       U^2 + V^2 + h^2 = R_i^2  ->  h^2 = R_i^2 - U^2 - V^2 */
    double h2 = R_i*R_i - U*U - V*V;
    if (h2 < -1e-8) return 0; /* nincs valós megoldás */
    if (h2 < 0.0) h2 = 0.0;   /* kis negatív numerikus hiba levágása */
    double h = sqrt(h2);

    /* (5) p visszaalakítása globális koordinátába:
       p = i + U*ex + V*ey ± h*ez */
    vec3 base = v3_add(i, v3_add(v3_scale(ex, U), v3_scale(ey, V)));

    if (p_plus)  *p_plus  = v3_add(base, v3_scale(ez,  h));
    if (p_minus) *p_minus = v3_add(base, v3_scale(ez, -h));

    if (U_out) *U_out = U;
    if (V_out) *V_out = V;
    if (h_out) *h_out = h;

    /* Ha h = 0, a két megoldás egybeesik (érintő eset). */
    return (h < EPS) ? 1 : 2;
}

int pass_like_coating(
    ap **atoms_io,
    int *n_atoms_io,
    int *cap_io,
    int model_ser,
    double sigma_p,
    double weed_dist,
    int n_layers,
    int *next_atom_ser,
    int *next_res_ser,
    char chain_id
) {
    if (!atoms_io || !n_atoms_io || !cap_io || !next_atom_ser || !next_res_ser) return 0;
    if (n_layers <= 0) return 0;

    ap *atoms = *atoms_io;
    int n_atoms = *n_atoms_io;
    int cap = *cap_io;

    int total_added = 0;

    printf("pass\n");

    for (int layer = 0; layer < n_layers; layer++) {
        size_t n_substrate = (size_t)n_atoms;
        printf("Current layer: %d \n", layer + 1);

        vec3 *pos = (vec3*)malloc(n_substrate * sizeof(vec3));
        double *sig = (double*)malloc(n_substrate * sizeof(double));
        if (!pos || !sig) { free(pos); free(sig); break; }

        for (size_t i = 0; i < n_substrate; i++) {
            pos[i] = pos_from_atom(&atoms[i]);
            sig[i] = sigma_from_atom(&atoms[i]);
        }

        grid3d g;
        if (!grid_build(&g, pos, sig, n_substrate, sigma_p)) {
            fprintf(stderr, "ERROR: grid_build failed\n");
            free(pos); free(sig);
            break;
        }

        /* Globális weeding grid (végső összefésüléshez) */
        grid_pts ng;
        if (!grid_pts_init(&ng, g.min, g.max, weed_dist, 1024)) {
            grid_free(&g);
            free(pos); free(sig);
            break;
        }

        /* A szálak által gyűjtött jelöltek egyesítve ide jönnek */
        vec3_list all_candidates = {0};

        /* ---------- OpenMP: párhuzamos a ciklus ---------- */
        #pragma omp parallel
        {
            /* Szál-lokális jelölt lista + szál-lokális weeding (csak saját jelölteken) */
            vec3_list local = {0};

            grid_pts local_ng;
            int local_ng_ok = grid_pts_init(&local_ng, g.min, g.max, weed_dist, 1024);

            int *neigh = (int*)malloc((size_t)NEIGH_CAP * sizeof(int));
            double *neigh_d2 = (double*)malloc((size_t)NEIGH_CAP * sizeof(double));

            #pragma omp for schedule(dynamic, 16)
            for (size_t a = 0; a + 2 < n_substrate; a++) {
                vec3 i = pos[a];
                double sigma_i = sig[a];

                /* cell-radius dinamikusan, de csak ahhoz, hogy a grid cellákat bejárjuk */
                double cutoff_max = (sigma_i + sigma_p) + (g.max_sigma + sigma_p);
                int r = (int)ceil(cutoff_max / g.cell);
                if (r < 1) r = 1;

                int ix = (int)floor((i.x - g.min.x) / g.cell);
                int iy = (int)floor((i.y - g.min.y) / g.cell);
                int iz = (int)floor((i.z - g.min.z) / g.cell);
                ix = clampi(ix, 0, g.nx - 1);
                iy = clampi(iy, 0, g.ny - 1);
                iz = clampi(iz, 0, g.nz - 1);

                int nnb = 0;

                /* Szomszédok gyűjtése: már itt pontos bridgeable szűrés + K-cap */
                for (int dz = -r; dz <= r; dz++) {
                    int zz = iz + dz; if (zz < 0 || zz >= g.nz) continue;
                    for (int dy = -r; dy <= r; dy++) {
                        int yy = iy + dy; if (yy < 0 || yy >= g.ny) continue;
                        for (int dx = -r; dx <= r; dx++) {
                            int xx = ix + dx; if (xx < 0 || xx >= g.nx) continue;

                            int gi = grid_index(&g, xx, yy, zz);
                            for (int bi = g.head[gi]; bi != -1; bi = g.next[bi]) {
                                if ((size_t)bi <= a) continue;

                                vec3 jpos = pos[(size_t)bi];
                                double sigma_j = sig[(size_t)bi];

                                /* pontos páros szűrés (ez sok szemetet kidob) */
                                double d2 = v3_dist2(i, jpos);
                                double RA = sigma_i + sigma_p;
                                double RB = sigma_j + sigma_p;
                                double sum = RA + RB;
                                double diff = fabs(RA - RB);
                                if (d2 > (sum * sum)) continue;
                                if (d2 < (diff * diff)) continue;

                                if (neigh && neigh_d2) {
                                    neigh_try_add(neigh, neigh_d2, &nnb, NEIGH_CAP, bi, d2);
                                }
                            }
                        }
                    }
                }

                /* Tripletezés a capped szomszédokból */
                for (int nb1 = 0; nb1 + 1 < nnb; nb1++) {
                    size_t b = (size_t)neigh[nb1];
                    vec3 j = pos[b];
                    double sigma_j = sig[b];

                    if (!pair_bridgeable(i, sigma_i, j, sigma_j, sigma_p)) continue;

                    for (int nb2 = nb1 + 1; nb2 < nnb; nb2++) {
                        size_t c = (size_t)neigh[nb2];
                        vec3 k = pos[c];
                        double sigma_k = sig[c];

                        if (!pair_bridgeable(i, sigma_i, k, sigma_k, sigma_p)) continue;
                        if (!pair_bridgeable(j, sigma_j, k, sigma_k, sigma_p)) continue;

                        vec3 p_plus, p_minus;
                        double U, V, h;
                        int nsol = three_point_sphere_geometry(
                            i, sigma_i,
                            j, sigma_j,
                            k, sigma_k,
                            sigma_p,
                            &p_plus, &p_minus,
                            &U, &V, &h
                        );
                        if (nsol == 0) continue;

                        for (int s = 0; s < nsol; s++) {
                            vec3 p = (s == 0) ? p_plus : p_minus;

                            /* Ütközés: globális szubsztráthoz (read-only, OK párhuzamosan) */
                            if (clashes_with_any_grid(p, &g, sigma_p)) continue;

                            /* Szál-lokális weed (csak a saját jelöltek között) */
                            if (local_ng_ok) {
                                if (too_close_to_new_grid(p, &local_ng, weed_dist)) continue;
                                (void)grid_pts_add(&local_ng, p);
                                (void)v3list_push(&local, p);
                            } else {
                                (void)v3list_push(&local, p);
                            }
                        }
                    }
                }
            }

            free(neigh);
            free(neigh_d2);
            grid_pts_free(&local_ng);

            /* Szál-lokális listák összefésülése a globáliszba */
            #pragma omp critical
            {
                for (size_t i = 0; i < local.n; i++) {
                    (void)v3list_push(&all_candidates, local.data[i]);
                }
            }

            v3list_free(&local);
        } /* omp parallel */

        /* ---------- Globális weeding + append (determinista, szálbiztos) ---------- */
        for (size_t t = 0; t < all_candidates.n; t++) {
            vec3 p = all_candidates.data[t];

            if (too_close_to_new_grid(p, &ng, weed_dist)) continue;
            if (!grid_pts_add(&ng, p)) break;
        }

        for (size_t t = 0; t < ng.n; t++) {
            if (n_atoms == cap) {
                size_t new_cap_atoms = (cap == 0) ? 1024 : (size_t)cap * 2;
                ap *tmp = (ap*)realloc(atoms, new_cap_atoms * sizeof(ap));
                if (!tmp) break;
                atoms = tmp;
                cap = (int)new_cap_atoms;
            }

            make_hoh_oxygen(
                &atoms[n_atoms],
                model_ser,
                (*next_atom_ser)++,
                (*next_res_ser)++,
                chain_id,
                ng.pts[t]
            );

            n_atoms++;
            total_added++;
        }

        v3list_free(&all_candidates);
        grid_pts_free(&ng);
        grid_free(&g);
        free(pos);
        free(sig);
    }

    *atoms_io = atoms;
    *n_atoms_io = n_atoms;
    *cap_io = cap;

    return total_added;
}
