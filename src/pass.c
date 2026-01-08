#include "pass.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

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

/* Két gömb "bridge-elhetősége" adott σ_p mellett (gyors szűrés tripletekhez).
   Ha két kiterjesztett sugarú gömb (R_a = σ_a+σ_p, R_b = σ_b+σ_p) nem metsz,
   akkor a háromgömb-metszés sem fog megoldást adni. */
static int pair_bridgeable(vec3 A, double sigma_A, vec3 B, double sigma_B, double sigma_p) {
    double d = v3_dist(A,B);
    double RA = sigma_A + sigma_p;
    double RB = sigma_B + sigma_p;
    if (d > (RA + RB)) return 0;                 /* túl messze */
    if (d < fabs(RA - RB)) return 0;             /* egyik "benne lenne" a másikban metszés nélkül */
    return 1;
}

/* Ütközésvizsgálat: a jelölt p pont ne lógjon bele egyetlen meglévő atomba sem.
   Feltétel: dist(p, atom) >= σ_atom + σ_p (kis toleranciával). */
static int clashes_with_any(
    vec3 p,
    const ap *atoms,
    size_t n_atoms,
    double sigma_p
) {
    for (size_t idx = 0; idx < n_atoms; idx++) {
        vec3 c = pos_from_atom(&atoms[idx]);
        double sigma_c = sigma_from_atom(&atoms[idx]);
        double min_d = sigma_c + sigma_p;
        if (v3_dist(p, c) < (min_d - CLASH_TOL)) {
            return 1; /* clash */
        }
    }
    return 0; /* ok */
}

/* Weeding: ne vegyünk fel túl közeli új pontot ugyanabban a layerben. */
static int too_close_to_new(
    vec3 p,
    const vec3 *new_pts,
    size_t n_new,
    double weed_dist
) {
    for (size_t i = 0; i < n_new; i++) {
        if (v3_dist(p, new_pts[i]) < weed_dist) return 1;
    }
    return 0;
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

    /* PDB rekord típusa: a víz tipikusan HETATM. */
    strncpy(out->pdb_token, "HETATM", sizeof(out->pdb_token) - 1);

    out->atom_ser = atom_ser;

    /* Atomnév/típus: oxigén */
    strncpy(out->pdb_type, "O", sizeof(out->pdb_type) - 1);

    /* Alternatív hely jelölés: üres */
    strncpy(out->pdb_alt_loc, "", sizeof(out->pdb_alt_loc) - 1);

    /* Maradék: HOH */
    strncpy(out->res_type, "HOH", sizeof(out->res_type) - 1);

    /* Lánc */
    out->chain[0] = chain_id;
    out->chain[1] = '\0';

    out->res_ser = res_ser;

    /* Insertion code */
    strncpy(out->pdb_achar, "", sizeof(out->pdb_achar) - 1);

    /* Koordináták */
    out->x_coord = p.x;
    out->y_coord = p.y;
    out->z_coord = p.z;

    /* Tipikus PDB mezők */
    out->occ = 1.0;
    out->b_factor = 0.0;

    /* Opcionális típusok üresen maradnak */
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

/* -------- "Bevonás" iteratív HOH hozzáadással -------- */
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
        /* A layer elején lefagyasztjuk a "szubsztrát" atomokat:
           csak ezekből képzünk tripleteket. Az előző layerek HOH-jai már benne vannak,
           de az ebben a layerben frissen hozzáadottak még nem. */
        size_t n_substrate = n_atoms;
        printf("Current layer: %d \n", layer+1);

        /* Dinamikus lista az új pontok (p) tárolására ebben a layerben */
        size_t new_cap = 1024;
        size_t n_new = 0;
        vec3 *new_pts = (vec3*)malloc(new_cap * sizeof(vec3));
        if (!new_pts) break;

        /* Triplet bejárás: O(n^3) - nagy rendszernél lassú lehet, de oktatási célra OK. */
        for (size_t a = 0; a + 2 < n_substrate; a++) {
            vec3 i = pos_from_atom(&atoms[a]);
            double sigma_i = sigma_from_atom(&atoms[a]);

            ap* atom = &atoms[a];

            printf("Atom for %d \n", atom->atom_ser);

            for (size_t b = a + 1; b + 1 < n_substrate; b++) {
                vec3 j = pos_from_atom(&atoms[b]);
                double sigma_j = sigma_from_atom(&atoms[b]);

                /* Gyors páros "bridge" szűrés i-j-re */
                if (!pair_bridgeable(i, sigma_i, j, sigma_j, sigma_p)) continue;

                for (size_t c = b + 1; c < n_substrate; c++) {

                    //printf("third for \n");
                    vec3 k = pos_from_atom(&atoms[c]);
                    double sigma_k = sigma_from_atom(&atoms[c]);

                    /* Gyors páros szűrések */
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

                    /* Mindkét megoldást (ha van) megpróbáljuk felvenni. */
                    for (int s = 0; s < nsol; s++) {
                        vec3 p = (s == 0) ? p_plus : p_minus;

                        if (clashes_with_new_points(p, new_pts, n_new, SIGMA_HOH_O)) continue;

                        /* Ütközés: p ne lógjon bele meglévő atomokba */
                        if (clashes_with_any(p, atoms, n_atoms, sigma_p)) continue;

                        /* Weeding: ugyanabban a layerben ne legyen túl közel más új ponthoz */
                        if (too_close_to_new(p, new_pts, n_new, weed_dist)) continue;

                        /* Hozzáadjuk az új ponthoz */
                        if (n_new == new_cap) {
                            new_cap *= 2;
                            vec3 *tmp = (vec3*)realloc(new_pts, new_cap * sizeof(vec3));
                            if (!tmp) { n_new = 0; break; }
                            new_pts = tmp;
                        }
                        new_pts[n_new++] = p;
                    }
                }
            }
        }

        /* Az összes új pontból HOH oxigéneket gyártunk és appendeljük az atomlistához. */
        for (size_t t = 0; t < n_new; t++) {
            if (n_atoms == cap) {
                /* Kapacitás növelés */
                size_t new_cap_atoms = (cap == 0) ? 1024 : cap * 2;
                ap *tmp = (ap*)realloc(atoms, new_cap_atoms * sizeof(ap));
                if (!tmp) break;
                atoms = tmp;
                cap = new_cap_atoms;
            }

            make_hoh_oxygen(
                &atoms[n_atoms],
                model_ser,
                (*next_atom_ser)++,
                (*next_res_ser)++,
                chain_id,
                new_pts[t]
            );

            n_atoms++;
            total_added++;
        }

        free(new_pts);
    }

    /* visszaírjuk a kimeneteket */
    *atoms_io = atoms;
    *n_atoms_io = n_atoms;
    *cap_io = cap;

    return total_added;
}

/* Segédfüggvény: megkeresi a legnagyobb atom_ser és res_ser értéket,
   hogy az új HOH atomok egyedi sorszámokat kapjanak. */
 void find_next_serials(const ap *atoms, int n_atoms, int *next_atom_ser, int *next_res_ser)
{
    int max_atom_ser = 0;
    int max_res_ser  = 0;

    for (int i = 0; i < n_atoms; i++) {
        if (atoms[i].atom_ser > max_atom_ser) max_atom_ser = atoms[i].atom_ser;
        if (atoms[i].res_ser  > max_res_ser)  max_res_ser  = atoms[i].res_ser;
    }

    *next_atom_ser = max_atom_ser + 1;
    *next_res_ser  = max_res_ser  + 1;
}

int clashes_with_new_points(vec3 p, const vec3 *new_pts, int n_new, double sigma_new)
{
    /* két oxigén ne fedje egymást: dist >= 2*sigma_new */
    double min_d = 2.0 * sigma_new;
    for (int i = 0; i < n_new; i++) {
        if (v3_dist(p, new_pts[i]) < (min_d - CLASH_TOL)) return 1;
    }
    return 0;
}

