#ifndef INPUT_H
#define INPUT_H

#include "global.h"
#include <stdbool.h>

extern char g_v;

void str0 (char string [], int length);

ap *read_in_pdb (char *file_name, int *atom_num, int mdl_serno);

typedef struct {
    char file[MAX_FILENAME_LENGTH];
    double sigma_p;
    double weed_dist;
    int n_layers;
} WdropProgramOptions;

void init_program_options(WdropProgramOptions *opt);
bool parse_program_args(int argc, char *argv[], WdropProgramOptions *opt);
void print_program_usage(const char *progname);

#endif /* INPUT_H */