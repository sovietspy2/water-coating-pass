#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "input.h"
#include "print.h"
#include "wdrop.h"

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
    struct timespec start_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    WdropProgramOptions opt;
    if (!parse_program_args(argc, argv, &opt)) {
        print_program_usage(argv[0]);
        return EXIT_FAILURE;
    }

    printf("Parameters:\n");
    printf("  file      : %s\n", opt.file);
    printf("  sigma_p   : %.3f\n", opt.sigma_p);
    printf("  weed_dist : %.3f\n", opt.weed_dist);
    printf("  n_layers  : %d\n", opt.n_layers);

    snprintf(g_ref_name, sizeof(g_ref_name), "%s", opt.file);

    g_pdb_ref = read_in_pdb(g_ref_name, &g_pdb_ref_no, -1);
    if (!g_pdb_ref || g_pdb_ref_no <= 0) {
        fprintf(stderr, "Error: PDB reading was unsuccessful, or PDB file is empty.\n");
        return EXIT_FAILURE;
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
        opt.sigma_p,
        opt.weed_dist,
        opt.n_layers,
        &next_atom_ser, &next_res_ser,
        'W'
    );

    printf("Added HOH oxygens: %d\n", added);

    char *output_file = pdb_to_edited(opt.file, opt.n_layers);
    printf("%s created\n", output_file);

    print_pdb_file(g_pdb_ref, g_pdb_ref_no, output_file);
    print_logs(output_file, added, start_time);

    free(g_pdb_ref);
    return EXIT_SUCCESS;
}