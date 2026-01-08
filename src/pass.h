

#include "global.h"
#include "macros.h"

#ifndef PASS_COAT_H
#define PASS_COAT_H

#include <stddef.h>  // size_t

/* HOH oxigén sugár PyMOL vdW szerint kb. 1.52 Å (PASS táblázatból) */
static const double SIGMA_HOH_O = 1.52;
static const double EPS = 1e-9;
static const double CLASH_TOL = 1e-4; /* numerikus tolerancia ütközés ellenőrzéshez */

/* Egyszerű 3D vektor típus a számolásokhoz. */
typedef struct {
    double x, y, z;
} vec3;

/*
 * Appendix A: Three-Point Sphere Geometry (jelölések az ábrához igazítva)
 *
 * Bemenet:
 *   i, j, k  - a három atom középpontja (pontok)
 *   sigma_i, sigma_j, sigma_k - a három atom sugara (σ_i, σ_j, σ_k)
 *   sigma_p  - a próbagömb sugara (σ_p)
 *
 * Kimenet:
 *   p_plus, p_minus - a két lehetséges megoldás p (a háromszög síkja fölött / alatt)
 *   U_out, V_out, h_out - az ábrán szereplő U, V és h (opcionális; lehet NULL)
 *
 * Visszatérés:
 *   0: nincs valós megoldás
 *   1: egy megoldás (érintő eset, h = 0)
 *   2: két megoldás (h > 0)
 */
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
);

/*
 * PASS-szerű "bevonás" (egyszerűsített): új HOH oxigéneket adunk hozzá a listához.
 *
 * atoms_io / n_atoms_io / cap_io:
 *   dinamikus tömbkezelés: ha nincs elég kapacitás, realloc történik.
 *
 * model_ser:
 *   az új atomok model/frame sorszáma.
 *
 * sigma_p:
 *   próbagömb sugara (σ_p).
 *
 * weed_dist:
 *   minimális távolság két újonnan felvett pont között ugyanazon iterációban.
 *
 * n_layers:
 *   hány bevonási iteráció fusson le.
 *
 * next_atom_ser / next_res_ser:
 *   külső számlálók, hogy folyamatosan egyedi sorszámot adjunk.
 *
 * chain_id:
 *   pl. 'W' vagy 'A' – ide kerül a HOH lánc azonosítója.
 *
 * Visszatérés: összesen hány HOH oxigént adott hozzá.
 */
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
);


void find_next_serials(const ap *atoms, int n_atoms, int *next_atom_ser, int *next_res_ser);

int clashes_with_new_points(vec3 p, const vec3 *new_pts, int n_new, double sigma_new);
#endif /* PASS_COAT_H */
