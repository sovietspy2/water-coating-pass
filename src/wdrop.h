

#include "global.h"
#include "macros.h"

#ifndef PASS_COAT_H
#define PASS_COAT_H

#include <stddef.h>  // size_t

/* HOH oxygen radius per PyMOL vdW is about 1.52 Å (from the PASS table) */
static const double SIGMA_HOH_O = 1.52;
static const double EPS = 1e-9;
static const double CLASH_TOL = 1e-4; /* numerical tolerance for clash checking */

/* Simple 3D vector type for the calculations. */
typedef struct {
    double x, y, z;
} vec3;

/*
 * Appendix A: Three-Point Sphere Geometry (notation aligned with the figure)
 *
 * Input:
 *   i, j, k  - the centers of the three atoms (points)
 *   sigma_i, sigma_j, sigma_k - the radii of the three atoms (σ_i, σ_j, σ_k)
 *   sigma_p  - the probe sphere radius (σ_p)
 *
 * Output:
 *   p_plus, p_minus - the two possible solutions p (above / below the triangle plane)
 *   U_out, V_out, h_out - the U, V and h shown in the figure (optional; may be NULL)
 *
 * Returns:
 *   0: no real solution
 *   1: one solution (tangent case, h = 0)
 *   2: two solutions (h > 0)
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
 * PASS-like "coating" (simplified): add new HOH oxygens to the list.
 *
 * atoms_io / n_atoms_io / cap_io:
 *   dynamic array handling: if capacity is insufficient, a realloc happens.
 *
 * model_ser:
 *   the model/frame serial number of the new atoms.
 *
 * sigma_p:
 *   probe sphere radius (σ_p).
 *
 * weed_dist:
 *   minimum distance between two newly added points within the same iteration.
 *
 * n_layers:
 *   how many coating iterations to run.
 *
 * next_atom_ser / next_res_ser:
 *   external counters so that we keep assigning unique serial numbers.
 *
 * chain_id:
 *   e.g. 'W' or 'A' – the HOH chain identifier goes here.
 *
 * Returns: the total number of HOH oxygens added.
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

int clashes_with_new_points(vec3 p, const vec3 *new_pts, int n_new, double sigma_new);

char *pdb_to_edited(const char *in, int layer);

void find_next_serials(const ap *atoms, int n_atoms, int *next_atom_ser, int *next_res_ser);

#endif /* PASS_COAT_H */
