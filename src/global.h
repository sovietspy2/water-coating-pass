#ifndef GLOBAL_H
#define GLOBAL_H

#include "../src/macros.h"

struct pdb_atoms {
	int model_ser;
	char pdb_token [PDB_TOKEN_LENGTH+1];
	int atom_ser;
	char pdb_type [PDB_TYPE_LENGTH+1];
	char pdb_alt_loc [PDB_ALT_LOC_LENGTH+1];
	char res_type [RES_TYPE_LENGTH+1];
	char chain [CHAIN_LENGTH+1];
	int res_ser;
	char pdb_achar [PDB_ACHAR_LENGTH+1];
	double x_coord;
	double y_coord;
	double z_coord;
	double occ;
	double b_factor;
	char mol2_type [MOL2_TYPE_LENGTH+1];
	double mol2_charge;
	char pdbqt_type [PDBQT_TYPE_LENGTH+1];
};

typedef struct pdb_atoms ap;

#endif