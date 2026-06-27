#include "../src/print.h"

char *create_log_file_name(const char *input, const char *sep)
{
    if (!input) input = "";
    if (!sep)   sep   = " ";

    time_t now = time(NULL);
    struct tm tm_local;
    if (localtime_r(&now, &tm_local) == NULL) return NULL;

    char dt[32];
    if (strftime(dt, sizeof dt, "wdrop-%Y-%m-%d %H:%M:%S", &tm_local) == 0) return NULL;

    int needed = snprintf(NULL, 0, "%s%s%s", input, sep, dt);
    if (needed < 0) return NULL;

    char *out = malloc((size_t)needed + 1);
    if (!out) return NULL;

    int written = snprintf(out, (size_t)needed + 1, "%s%s%s", input, sep, dt);
    if (written < 0 || written > needed) { free(out); return NULL; }

    return out;
}

char *print_pdb_file (ap *pdb, int atom_num, char file_out [MAX_FILENAME_LENGTH]) {

    FILE *outfile;
    int j;
    char *line = NULL;
    char ter [5] = "TER\n";
    char end [5] = "END\n";
    char header [MAX_LINE_LENGTH];

    char *print_pdb_line (ap *pdb, int index, char pdbqt_ind, int pdbqt_rank);

    outfile = fopen(file_out, "w");
    if (!outfile) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", file_out);
        return NULL;
    }

    sprintf(header,"%s","REMARK Input coordinates.\nMODEL        1\n");
    fputs(header,outfile);

    for (j = 0; j < atom_num; j++) {
        line = print_pdb_line(pdb, j, ' ',j);
        fputs(line,outfile);
        free(line);
        if (j == atom_num-1) fputs(end,outfile);
        else if (strcmp ((pdb+j)->chain,(pdb+j+1)->chain) != 0) fputs(ter,outfile);
    }
    fclose(outfile);
    return(file_out);
}

int print_logs(const char *filepath, int added, struct timespec start_time)
{
    char *file_name = create_log_file_name(filepath, "-");
    if (!file_name) return -1;

    FILE *f = fopen(file_name, "a");
    if (!f) { free(file_name); return -1; }

    // adding water molecule count
    fprintf(f, "added water O atoms: %d\n", added);

    // adding execution time
    struct timespec end_time;
    clock_gettime(CLOCK_MONOTONIC, &end_time);

    double elapsed_ms =
        (double)(end_time.tv_sec - start_time.tv_sec) * 1000.0 +
        (double)(end_time.tv_nsec - start_time.tv_nsec) / 1000000.0;

    fprintf(f, "runtime: %.3f ms\n", elapsed_ms);

    fclose(f);
    free(file_name);
    return 0;
}

char *print_pdb_line (ap *pdb, int index, char pdbqt_ind, int pdbqt_rank) {

    char *line = NULL;
    line = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));

    if ( pdbqt_ind == 't' ) {
        sprintf(line,"%-6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%10.3f%s%-4s%s",

        "ATOM",
        (pdb+index)->atom_ser,
        " ", /* just space */
         (pdb+index)->pdb_type,
        (pdb+index)->pdb_alt_loc,
        (pdb+index)->res_type,
        " ", /* just space */
        (pdb+index)->chain,
        (pdb+index)->res_ser,		// pdbqt_rank was here
        (pdb+index)->pdb_achar,
        "   ", /* just space */
        (pdb+index)->x_coord,
        (pdb+index)->y_coord,
        (pdb+index)->z_coord,
        (pdb+index)->occ,
        (pdb+index)->b_factor,
        (pdb+index)->mol2_charge,
        " ", /* just space */
        (pdb+index)->pdbqt_type,
        "\n");
    }

    else {
        sprintf(line,"%-6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%s",

        "ATOM",
        (pdb+index)->atom_ser,
        " ", /* just space */
         (pdb+index)->pdb_type,
        (pdb+index)->pdb_alt_loc,
        (pdb+index)->res_type,
        " ", /* just space */
        (pdb+index)->chain,
        (pdb+index)->res_ser,
        (pdb+index)->pdb_achar,
        "   ", /* just space */
        (pdb+index)->x_coord,
        (pdb+index)->y_coord,
        (pdb+index)->z_coord,
        (pdb+index)->occ,
        (pdb+index)->b_factor,

        "\n");
    }
    return(line);
}