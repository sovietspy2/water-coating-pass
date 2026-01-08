// Constants
#define MAX_FILENAME_LENGTH 100		
#define MAX_LINE_LENGTH 200
#define MIN_LINE_LENGTH 3
#define MIN_LINE_LENGTH2 30
#define MIN_LINE_LENGTH3 10
#define VERY_BIG_DISTANCE 1000000
#define MAX_OUT_TEXT_LENGTH 5000
#define DEF_MODEL_NO 1
#define PDB_TOKEN_LENGTH 6     		// standard pdb definitions 
#define PDB_ATOM_NUM_LENGTH 5 
#define PDB_TYPE_LENGTH 4
#define PDB_ALT_LOC_LENGTH 1
#define RES_TYPE_LENGTH 3
#define CHAIN_LENGTH 1
#define PDB_RES_NUM_LENGTH 4
#define PDB_ACHAR_LENGTH 1
#define PDB_COORD_LENGTH 8
#define PDB_OCC_B_LENGTH 6
#define MOL2_TYPE_LENGTH 5
#define PDBQT_CHARGE_LENGTH 10
#define PDBQT_TYPE_LENGTH 4
#define MDL_NUM_LENGTH 4
#define NUM_LENGTH 20
#define REMARK_TOKEN_LENGTH 6
#define DEF_MAX_B_FACTOR 30.000		// defaults
#define DEF_CLUST_TOLERANCE 1.000			// -c
#define DEF_DIST_TOLERANCE 3.500			// -d
#define DEF_SEP_TOLERANCE 1.500				// match tolerance -s
#define DEF_PRED_LIST_LENGTH 50.000
#define DEF_IDENTICAL_PREDICTION_TOLERANCE 2.500	// -p
#define CONS_SCALE 0.5
#define	DEF_MIN_CTOL   0.250				// MIN + NUM*STEP = MAX
#define	DEF_MIN_PTOL   1.000
#define	DEF_MIN_STOL   2.000
#define	DEF_NUM_CTOL   8
#define	DEF_NUM_PTOL   4
#define	DEF_NUM_STOL   0
#define	DEF_STEP_CTOL  0.250
#define	DEF_STEP_PTOL  0.500
#define	DEF_STEP_STOL  0
#define EVAL_STEPC 0.5
#define EVAL_STEPI 0.5
#define EVAL_STEPS 0.5
#define DEF_P_TOP 50.000
#define DEF_MIN_DIST_TOLERANCE 1.750
#define MIN_B_FACTOR 0.00
#define MAX_B_FACTOR 100.00
#define DEF_PRINT_RANK_NUM 100
#define MAX_OCC 100.00
#define DEF_MOB_TOLERANCE 50.000

// Default char
#define DEF_MOBYREF_FILENAME "system_ref.pdb"
#define DEF_PLI_FILENAME "system.pli"
#define DEF_TOPOL_FILENAME "system_tpy.pdb"
#define DEF_TARGET_FILENAME "target.pdb"
#define DEF_REF_FILENAME "reference.pdb"
#define DEF_WATERS_FILENAME "waters.pdb"
#define DEF_E_WISH_NAME "energy_types.prm"
#define DEF_V_MODE 's'
#define DEF_R_MODE ' '
#define AUTO_KEYW "Auto"
#define WAT_KEYW "WAT"
#define DEF_E 'i'					// other mode s
#define DEF_FRAME_FILENAME "system.xtc"
#define DEF_PROG_MODE "Prediction"	// other mode(s): Analysis & more will come
#define ALT1_PROG_MODE "Analysis"
#define ALT2_PROG_MODE "Editing"
#define ALT3_PROG_MODE "NetDraw"
#define DEF_M 'p'
#define DEF_FRAME_RANGE "0-10"
#define DEF_PRINT_RANKS "30"
#define DEF_F_O_FILENAME "output.pdb"
#define DEF_DLG_FILENAME "ligand_target.dlg"
#define DEF_I_TYPE 'n'
#define MODEL_TAG "mdl"
#define ND_STYLE "filled"
#define WAT_ND_COL "#486FB4"
#define BLK_ND_COL "#FF3300"
#define LIG_ND_COL "#AFF1F5"
#define TAR_ND_COL "#C2C4E8"
#define WAT_ND_COL_RGB "\'72,111,180\'"
#define BLK_ND_COL_RGB "\'255,51,0\'"
#define LIG_ND_COL_RGB "\'175,241,245\'"
#define TAR_ND_COL_RGB "\'194,196,232\'"
#define ED_STYLE "solid"
#define DYN_COL "#000000"
#define STA_COL "#ff0000"
#define DYN_COL_RGB "\'0,0,0\'"
#define STA_COL_RGB "\'255,0,0\'"

// Header & error char
#define PROGRAM "MobyWat"
#define MOBYWAT_NAME "mobywat"
#define PANTERG_NAME "pantherg"
#define PANTOOL_NAME "pantools"
#define NAME_LENGTH 7
#define INTRO "\
__________________________________________________________________________________\n\
\n\
MobyWat\n\
Calculation of hydration structures of molecular surfaces and interfaces\n\
__________________________________________________________________________________\n\
                                                          ===Ver=1.1=10=04=2016===\n"
#define EPILOGUE2 "\
__________________________________________________________________________________\n"
#define EPILOGUE "\
__________________________________________________________________________________\n\
                                                                ===by=C=Hetenyi===\n"
#define WARN_BOXATN printf("%s","WARNING Number of atoms in box is not equal to number of atoms assigned to ligand+target+waters.\n");
#define MEMERR { printf("ERROR: No Memory Available.\n"); exit(0); }
#define PROGTER { printf("%s%s",PROGRAM," terminates.\n"); exit(0); }
#define TRAJ_NOIFACE_FILE_HEADER "Number of selected interface waters in each frame\n"
#define TRAJ_NOSFACE_FILE_HEADER "Number of selected surface waters in each frame\n"
#define OLIST_IF_FILE_HEADER "Occupancy list of representative interface waters\n"
#define OLIST_SF_FILE_HEADER "Occupancy list of representative surface waters\n"
#define ALL_CLUST_WAT_FILE_I_HEADER "REMARK PDB list of all clustered interface water OW atoms\n"
#define ALL_CLUST_WAT_FILE_S_HEADER "REMARK PDB list of all clustered surface water OW atoms\n"
#define ALL_CLUST_WAT_FILE_HEADER "REMARK Residue type column holds predicted Water ser. numbers\n\
REMARK Occupancy column holds original Atom ser. numbers\n\
REMARK B-factor column holds Frame ser. numbers\n"
#define CLIST_IF_FILE_HEADER "REMARK PDB list of average representatives of clustered interface water OW atoms\n"
#define CLIST_SF_FILE_HEADER "REMARK PDB list of average representatives of clustered surface water OW atoms\n"
#define CLIST_WAT_FILE_HEADER "REMARK Residue type column holds predicted Water ser. numbers\n\
REMARK Occupancy column holds Cluster ser. numbers\n\
REMARK B-factor column holds Occupancy counts of clusters\n"
#define TRAJ_ALLIFACE_FILE_HEADER "REMARK PDB List of interface water OW atoms ordered by frames\n\
REMARK Occupancy column holds original Atom ser. numbers\n\
REMARK B-factor column holds original Residue ser. numbers\n"
#define TRAJ_ALLSFACE_FILE_HEADER "REMARK PDB List of surface water OW atoms ordered by frames\n\
REMARK Occupancy column holds original Atom ser. numbers\n\
REMARK B-factor column holds original Residue ser. numbers\n"
#define REF_IF_WAT_FILE_HEADER "REMARK PDB List of reference interface water OW atoms \n"
#define REF_SF_WAT_FILE_HEADER "REMARK PDB List of reference surface water OW atoms \n"

// Functions
#define INP_M(X) printf("%s%s%s","   ",X," Mode.\n");   
#define INP_TY(X) printf("%s%s%s","       >>> ",X," Analysis.\n"); 
#define INP_P(X,Y,Z) printf("%s%c%s%-20.3f%s%-s%s","   -",X,"  ",Y,"  ",Z,"\n");
#define INP_R(X,Y,Z) printf("%s%c%s%-20d%s%-s%s","   -",X,"  ",Y,"  ",Z,"\n");
#define INP_S(X,Y) printf("%s%c%s%-40.40s%s","   -",X,"  ",Y,"\n");
#define NL printf("%s","\n");
#define HELP printf("%s","Type \"mobywat\" for quick help!\n");
#define INPUT(X,Y,Z) printf("%s%s%s%-20.20s%s%-s%s","   -",X,"  ",Y,"  ",Z," file)\n");
#define OUT_T(X) printf("%s%s%s","\n   Output ",X,"\n");
#define PDB_CANNOT_OPEN(X) printf("%s%s%s","ERROR! The ",X," file cannot be opened or missing model(s).\n");
#define OUT(X) printf("%s%-40.40s%s","   ",X,"\n");
#define WARN_ATN(X) printf("%s%s%s","WARNING The are more atoms in ",X," section than in the entire box.\n");
#define PRINT_LINE(X,Y,Z,I) { for (I = 0; I < Y; I++) fprintf(X,"%s",Z); fprintf(X,"%s","\n"); }
#define PRINT_PDB_FILES(A,B,C) { str0 (g_temp, MAX_FILENAME_LENGTH); \
			sprintf(g_temp,"%s%s%s%s","O_",g_root_name,A,"\0"); \
			print_pdb_file (B, C, g_temp); \
			print_log_file ( g_log_file_name, g_temp, 'a' ); \
			print_log_file ( g_log_file_name, "\n", 'a' ); \
			str0 (g_temp, MAX_FILENAME_LENGTH); }
#define PRINT_PLIST(A,B,C,D) { str0 (out_file_name, MAX_FILENAME_LENGTH); \
		sprintf(out_file_name,"%s%s%s%s","O_",root_name,A,"\0"); \
		print_plist_waters(B,C, out_file_name, no_frames, D); \
		print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL }
#define PRINT_MLIST_PLIST(A,B,C,D,E) { str0 (out_file_name, MAX_FILENAME_LENGTH); \
		sprintf(out_file_name,"%s%s%s%s","O_",root_name,A,"\0"); \
		print_mlist_plist(B, no_cl_members, C, D, s_tol, b_factor, out_file_name, E); \
		print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL }
#define PRINT_SUCC_MATR(A,B) { str0 (out_file_name, MAX_FILENAME_LENGTH); \
		sprintf(out_file_name,"%s%s%s%d%s%s","O_",root_name,"_",i+1,A,"\0"); \
		print_succ_matr(B, i, no_lowb_ref_wats, out_file_name); \
		print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL }
#define NLL { print_log_file ( g_log_file_name, "\n", 'a' ); }
#define PAR2LOGFILE_DC(A,B) { str0 (g_temp, MAX_FILENAME_LENGTH); \
			sprintf(g_temp,"%-37s%3s%42.3f%s",A," : ",B,"\0"); \
			print_log_file ( g_log_file_name, g_temp, 'a' ); NLL }
#define PAR2LOGFILE_CC(A,B) { str0 (g_temp, MAX_FILENAME_LENGTH); \
			sprintf(g_temp,"%-37s%3s%42s%s",A," : ",B,"\0"); \
			print_log_file ( g_log_file_name, g_temp, 'a' ); NLL }
#define PAR2LOGFILE_IC(A,B) { str0 (g_temp, MAX_FILENAME_LENGTH); \
			sprintf(g_temp,"%-37s%3s%42d%s",A," : ",B,"\0"); \
			print_log_file ( g_log_file_name, g_temp, 'a' ); NLL }
#define PAR2LOGFILE_LC(A,B) { str0 (g_temp, MAX_FILENAME_LENGTH); \
			sprintf(g_temp,"%-37s%3s%42ld%s",A," : ",B,"\0"); \
			print_log_file ( g_log_file_name, g_temp, 'a' ); NLL }
#define PAR2LOGFILE_EC(A,B) { str0 (g_temp, MAX_FILENAME_LENGTH); \
			sprintf(g_temp,"%-37s%3s%42.3e%s",A," : ",B,"\0"); \
			print_log_file ( g_log_file_name, g_temp, 'a' ); NLL }
#define FILENAME2GTEMP(A) { str0 (g_temp, MAX_FILENAME_LENGTH); \
			strcpy(g_temp, A); }
#define READ_SWITCH_TEXT(A,B,C,D,E,F) { if ( i+1 <= argc )  A = B; else A = ' ';  \
			if ( A == B && i+1 < argc ) { if ( argv[i+1][0] != '-' ) strcpy(C, argv[i+1]); else if (A != ' ') A = 'e'; } else if (A != ' ') A = 'e';  \
			if ( A == B && i+1 == argc ) A = 'e'; \
			if ( A == 'e' ) { printf("%s%s", D, E); F } } 
#define READ_SWITCH_NUM(A,B,C,D,E,F,G) { if ( i+1 < argc && strncmp(*(argv+i)+2,A,B) == 0 && argv[i+1][0] != '-' && C(argv[i+1]) > 0 ) D = C(argv[i+1]); \
			if ( strncmp(*(argv+i)+2,A,B) == 0 && (i+1 == argc || argv[i+1][0] == '-') ) \
			{ D = E; printf(F,G,D,"\n"); } }
#define READ_SWITCH_NUM2(A,B,C,D,E,F,G) { if ( i+1 < argc && strncmp(*(argv+i)+2,A,B) == 0 && argv[i+1][0] != '-' && C(argv[i+1]) >= 0 ) D = C(argv[i+1]); \
			if ( strncmp(*(argv+i)+2,A,B) == 0 && (i+1 == argc || argv[i+1][0] == '-') ) \
			{ D = E; printf(F,G,D,"\n"); } }
#define READ_SWITCH_TEXT_2(A,B,C,D,E,F,G,H) { if (i+1 <= argc && strncmp(*(argv+i)+2,A,B) == 0 ) C = D; else C = ' '; \
			if ( C == D && i+1 < argc ) { if ( argv[i+1][0] != '-' ) strcpy(E, argv[i+1]); else if (C != ' ') C = 'e'; } else if (C != ' ') C = 'e'; \
			if ( C == D && i+1 == argc ) C = 'e'; \
			if ( C == 'e' ) { printf("%s%s", F, G); H } }

#define TERMINATE { HELP exit(1); }

