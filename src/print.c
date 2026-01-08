#include "../src/print.h"

char *print_pdb_file (ap *pdb, int atom_num, char file_out [MAX_FILENAME_LENGTH]) {

    FILE *outfile;
    int j;
    char *line = NULL;
    char ter [5] = "TER\n";
    char end [5] = "END\n";
    char header [MAX_LINE_LENGTH];

    char *print_pdb_line (ap *pdb, int index, char pdbqt_ind, int pdbqt_rank);

    outfile=fopen(file_out,"w");

    sprintf(header,"%s","REMARK Input coordinates.\n");
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