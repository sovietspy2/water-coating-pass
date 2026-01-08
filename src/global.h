#ifndef GLOBAL_H
#define GLOBAL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <memory.h>
#include "../src/macros.h"

struct pdb_atoms {
	int model_ser;                         // NMR vagy trajektória fájlok esetén a modell/frame sorszáma
	char pdb_token [PDB_TOKEN_LENGTH+1];   // A PDB sor eleji azonosítója (pl. "ATOM  " vagy "HETATM")
	int atom_ser;                          // Az atom egyedi sorszáma a fájlon belül
	char pdb_type [PDB_TYPE_LENGTH+1];     // Az atom neve/típusa (pl. "CA", "N", "O")
	char pdb_alt_loc [PDB_ALT_LOC_LENGTH+1]; // Alternatív konformáció jelölése (ha egy atom több helyen is lehet)
	char res_type [RES_TYPE_LENGTH+1];     // Az aminosav vagy maradék hárombetűs kódja (pl. "THR", "CYS")
	char chain [CHAIN_LENGTH+1];           // A fehérjelánc azonosító betűje (pl. "A", "B")
	int res_ser;                           // Az aminosav sorszáma a láncon belül
	char pdb_achar [PDB_ACHAR_LENGTH+1];   // "Insertion code" - extra karakter az aminosav sorszámához
	double x_coord;                        // Az atom X koordinátája (Angströmben)
	double y_coord;                        // Az atom Y koordinátája (Angströmben)
	double z_coord;                        // Az atom Z koordinátája (Angströmben)
	double occ;                            // Occupancy (kitöltöttség) - az atom előfordulási valószínűsége
	double b_factor;                       // B-faktor - az atom hőmozgását/bizonytalanságát jelző érték
	char mol2_type [MOL2_TYPE_LENGTH+1];   // Mol2 formátum szerinti atomtípus (erőtér számításokhoz)
	double mol2_charge;                    // Az atom parciális töltése (elektrosztatikai számításokhoz)
	char pdbqt_type [PDBQT_TYPE_LENGTH+1]; // PDBQT formátum szerinti atomtípus (pl. AutoDock-hoz)
};

typedef struct pdb_atoms ap;

struct pdb_topology {
	int atom_ser;
	char pdb_type [PDB_TYPE_LENGTH+1];
	char pdb_alt_loc [PDB_ALT_LOC_LENGTH+1];
	char res_type [RES_TYPE_LENGTH+1];
	char chain [CHAIN_LENGTH+1];
	int res_ser;
	double occ; 
	double b_factor;
	char pdb_achar [PDB_ACHAR_LENGTH+1];	
	};

typedef struct pdb_topology tp;

struct pdb_distance {			// aimed for ligand-target interactions to select the closest target residues; ligand = waters sometimes
	int ligand_ser_no;
	int ligand_atom_ser_no;
	char ligand_pdb_type [PDB_TYPE_LENGTH+1];	
	char ligand_res_type [RES_TYPE_LENGTH+1];
	char ligand_chain [CHAIN_LENGTH+1];			 
	int ligand_res_ser;
	double ligand_x_coord;
	double ligand_y_coord;
	double ligand_z_coord;
	double ligand_b_factor;
	int target_atom_ser_no;
	char target_pdb_type [PDB_TYPE_LENGTH+1];	
	char target_res_type [RES_TYPE_LENGTH+1];
	char target_chain [CHAIN_LENGTH+1];			 
	int target_res_ser;
	double target_b_factor;
	double target_dist; 
	};

typedef struct pdb_distance pd;

struct hnet_edges {
	int wat_atm_ser;
	int wat_res_ser;
	int ptr_atm_ser;
	int ptr_res_ser;
	char wat_res_typ [RES_TYPE_LENGTH+1];
	char ptr_atm_typ [PDB_TYPE_LENGTH+1];
	char ptr_res_typ [RES_TYPE_LENGTH+1];
	char wat_chn [CHAIN_LENGTH+1];
	char ptr_chn [CHAIN_LENGTH+1];
	char ptr_type;		// l,t,b
	double wat_b_factor;
	double ptr_b_factor;
	double dist;
	char used;		// y=used
	char wat_dyn;
	char ptr_dyn;
	char edge_dyn;		// for drawing static/dynamic subnnets
	};

typedef struct hnet_edges hs;

struct hnet_node {
	int atm_ser;
	int res_ser;
	char res_typ [RES_TYPE_LENGTH+1];
	char chain [CHAIN_LENGTH+1];
	char type;		// w,l,t,b
	double b_factor;	// mobility for w otherwise not used
	int degree;
	int sol_degree;
	char node_dyn;		// for drawing static/dynamic subnnets
	};

typedef struct hnet_node hd;


struct res_list {
	int ligand_ser_no;		// also frame_no
	int target_atom_ser_no;		// word "target" is irrelevant, just left in the code from the beginning
	char target_res_type [RES_TYPE_LENGTH+1];
	char target_chain [CHAIN_LENGTH+1];			 
	int target_res_ser;
	double x_coord;
	double y_coord;
	double z_coord;	
	double b_factor;
	int label_mind_dist;		// =1 if the residue (water) is in close contact (or if clustered in cluster_by_olist); or if used up in plist
	int water_no;			// water count in plist
	double avg_water_no;		// for merged and cons plists
	int cluster_no;			// # rows no in olist
	int cluster2_no;		// # clusters in olist
	char plist_type;
	int plist_ser_no;
	} ;

typedef struct res_list re;


struct pool_list {
	float x_coord;
	float y_coord;
	float z_coord;	
	int aser;		// atom serial no
	int rser;		// res serial no
	} ;

typedef struct pool_list pl;

struct pool_info {
	double dmax;
	int nmin;
	int nmax;
	char t;
	char l;
	char w;
	} ;

typedef struct pool_info pf;

struct water {
	int atom_ser;
	int res_ser;
	int filt_atom_ser;
	int filt_res_ser;
	double dist;
	int counter;
	};

typedef struct water wr;

struct matching_ref_water {
	int ref_water_ser_no;
	int ref_water_atom_ser_no;
	int ref_water_res_ser_no;
	double distance;
	double b_factor;
	int indicator;				// 0 if no match found for the prediction
	};

typedef struct matching_ref_water mr;

struct thresholds {
	double min_ctol;
	double min_ptol;
	double min_stol;
	int num_ctol;
	int num_ptol;
	int num_stol;
	double step_ctol;
	double step_ptol;
	double step_stol;
	double lowb_final_ref_wats;		// to be removed
	};

typedef struct thresholds ts;

struct input_report {
	int target_traj_atom_no;
	int ligand_traj_atom_no;
	int waters_traj_atom_no;
	int target_traj_ca_no;
	};

typedef struct input_report ir;

#endif