#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "input.h"
#include "print.h"
#include "pass.h"

char g_v;
char g_pli_name[MAX_FILENAME_LENGTH];
double g_p_top;
char g_tpy;
char g_c_scheme;
char g_ref_name[MAX_FILENAME_LENGTH];
char g_r;
char g_top_file_name[MAX_FILENAME_LENGTH];
char g_pli;

int g_pdb_ref_no = 0; // atom counter
ap *g_pdb_ref = NULL; // atoms array


int main(int argc, char *argv[]) {

    char file[200];
    // pass params
    double sigma_p = 1.8; /* próbagömb sugara */
    double weed_dist = 3.5; /* minimális távolság két új pont között */
    int n_layers = 1;

    if (argc >= 2) strcpy(file, argv[1]);
    if (argc >= 3) sigma_p = atof(argv[2]);
    if (argc >= 4) weed_dist = atof(argv[3]);
    if (argc >= 5) n_layers = atoi(argv[4]);

    printf("Parameters:\n");
    printf("  file      : %s\n", file);
    printf("  sigma_p   : %.3f\n", sigma_p);
    printf("  weed_dist : %.3f\n", weed_dist);
    printf("  n_layers  : %d\n", n_layers);

    strncpy(g_ref_name, file, MAX_FILENAME_LENGTH - 1);
    g_ref_name[MAX_FILENAME_LENGTH - 1] = '\0';

    g_pdb_ref = read_in_pdb(g_ref_name, &g_pdb_ref_no, 1);
    if (!g_pdb_ref || g_pdb_ref_no <= 0) {
        fprintf(stderr, "Hiba: PDB beolvasas sikertelen vagy ures.\n");
        return 1;
    }

    int g_pdb_ref_cap = g_pdb_ref_no;

    int next_atom_ser = 1;
    int next_res_ser = 1;
    find_next_serials(g_pdb_ref, g_pdb_ref_no, &next_atom_ser, &next_res_ser);

    int added = pass_like_coating(
        &g_pdb_ref,
        &g_pdb_ref_no,
        &g_pdb_ref_cap,
        1,
        sigma_p,
        weed_dist,
        n_layers,
        &next_atom_ser, &next_res_ser,
        'W'
    );

    printf("Added HOH oxygens: %d\n", added);

    char *output_file = pdb_to_edited(file);
    printf("%s created \n", output_file);
    print_pdb_file(g_pdb_ref, g_pdb_ref_no, output_file);

    free(g_pdb_ref);

    return 0;
}
