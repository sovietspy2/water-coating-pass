#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "input.h"
#include "../src/print.h"
#include "macros.h"
#include "../src/pass.h"

char g_v;
char g_pli_name [MAX_FILENAME_LENGTH];
double g_p_top;
char g_tpy;
char g_c_scheme;
char g_ref_name [MAX_FILENAME_LENGTH];
char g_r;
char g_top_file_name [MAX_FILENAME_LENGTH];
char g_pli;

int g_pdb_ref_no = 0;   // atom counter
ap *g_pdb_ref = NULL;   // atoms array

void find_next_serials(const ap *atoms, int n_atoms, int *next_atom_ser, int *next_res_ser)
{
    int max_atom_ser = 0;
    int max_res_ser  = 0;

    if (!next_atom_ser || !next_res_ser) return;

    /* Üres lista esetén kezdjük 1-ről */
    if (!atoms || n_atoms <= 0) {
        *next_atom_ser = 1;
        *next_res_ser  = 1;
        return;
    }

    for (int i = 0; i < n_atoms; i++) {
        if (atoms[i].atom_ser > max_atom_ser) max_atom_ser = atoms[i].atom_ser;
        if (atoms[i].res_ser  > max_res_ser)  max_res_ser  = atoms[i].res_ser;
    }

    *next_atom_ser = max_atom_ser + 1;
    *next_res_ser  = max_res_ser  + 1;
}

int main(int argc, char *argv[])
{


    printf("Starting app");

    strncpy(g_ref_name, "1PSV.pdb", MAX_FILENAME_LENGTH - 1);
    g_ref_name[MAX_FILENAME_LENGTH - 1] = '\0';

    g_pdb_ref = read_in_pdb(g_ref_name, &g_pdb_ref_no, -2);
    if (!g_pdb_ref || g_pdb_ref_no <= 0) {
        fprintf(stderr, "Hiba: PDB beolvasas sikertelen vagy ures.\n");
        return 1;
    }

    /* --- 2) PASS-szerű paraméterek --- */
    const double sigma_p   = 1.8;  /* próbagömb sugara (σ_p) */
    const double weed_dist = 3.5;  /* minimális távolság két új pont között */
    const int n_layers     = 1;    /* hány "bevonási" iterációt futtassunk */

    const int n_runs = 1;

    /* --- 3) Dinamikus tömbkapacitás kezelése ---
       read_in_pdb() után nem tudjuk biztosan a kapacitást; induljunk abból,
       hogy pont akkora, mint a beolvasott atomszám, és a pass_like_coating()
       majd realloc-ol, ha kell. */
    int g_pdb_ref_cap = g_pdb_ref_no;

    /* --- 4) Következő sorszámok az új HOH atomokhoz --- */
    int next_atom_ser = 1;
    int next_res_ser  = 1;
    find_next_serials(g_pdb_ref, g_pdb_ref_no, &next_atom_ser, &next_res_ser);

    /* --- 5) Bevonás futtatása --- */
    for (int run = 0; run < n_runs; run++) {
        int added = pass_like_coating(
            &g_pdb_ref,
            &g_pdb_ref_no,
            &g_pdb_ref_cap,
            /* model_ser */ 1,
            /* sigma_p   */ sigma_p,
            /* weed_dist */ weed_dist,
            /* n_layers  */ n_layers,
            /* serialok  */ &next_atom_ser, &next_res_ser,
            /* chain_id  */ 'W'    /* HOH lánc azonosító */
        );

        printf("RUN %d: Added HOH oxygens: %d\n", run + 1, added);
        if (added <= 0) {
            printf("Nincs uj pont.\n");
            break;
        }
    }

    /* --- 6) Kiírás --- */
    print_pdb_file(g_pdb_ref, g_pdb_ref_no, "1PSV_edited.pdb");

    free(g_pdb_ref);

    return 0;
}
