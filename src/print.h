#include "global.h"
#include <time.h>

char *print_pdb_file (ap *pdb, int atom_num, char file_out [MAX_FILENAME_LENGTH]);

int print_logs(const char *filepath, int added, struct timespec start_time);