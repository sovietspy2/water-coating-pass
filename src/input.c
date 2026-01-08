#include "../src/input.h"

ap *read_in_pdb (char *file_name, int *atom_num, int mdl_serno) {
	
	FILE *f_pdb;
	char line [MAX_LINE_LENGTH]="";
	char atom_token []="ATOM  ";
	char hetatm_token []="HETATM";
	char model_token []="MODEL";
	char endmdl_token []="ENDMDL";
	char temp [NUM_LENGTH]="";
	char temp2 [NUM_LENGTH]="";	
	char *delims = " \0";
	int counter = 0;
	int line_indicator = 0;
	int temp_mdl_no = -10000;
	int i;
	ap *molecule = NULL;
	ap *x = NULL;

	str0 (temp, NUM_LENGTH);
	str0 (temp2, NUM_LENGTH);
	str0 (line, MAX_LINE_LENGTH);

	f_pdb = fopen(file_name, "r");
	if (f_pdb != NULL) { 
		while (fgets(line, sizeof(line), f_pdb)) {
			if ( mdl_serno > -1 ) {					// mdl_serno > -1 if mdl file to be read
				if ( strncmp (model_token, line, 5) == 0 ) {
					strncpy (temp, line+10, MDL_NUM_LENGTH);
					temp [MDL_NUM_LENGTH] = '\0';
 					temp_mdl_no = atoi(temp);
					str0 (temp, NUM_LENGTH);
					if (temp_mdl_no == mdl_serno) line_indicator = 1; else line_indicator = 0;
					}
				if ( line_indicator == 1 && strncmp (endmdl_token, line, 6) == 0 ) break;
				}
			if ( mdl_serno <= -1 ) line_indicator = 1;		// mdl_serno <= -1 if separate pdb file to be read
			if ( line_indicator == 1 && ( strncmp (atom_token, line, PDB_TOKEN_LENGTH) == 0 || strncmp (hetatm_token, line, PDB_TOKEN_LENGTH) == 0 ) )  {
				if ( counter == 0) molecule = (ap *) malloc(sizeof(ap));
					else {
						x = realloc( molecule, (counter+1) * sizeof(ap) );
						if (x) molecule = x;
						x = NULL;
						}
				(molecule+counter)->model_ser = DEF_MODEL_NO;
				strncpy ((molecule+counter)->pdb_token, line, PDB_TOKEN_LENGTH);	
				(molecule+counter)->pdb_token[PDB_TOKEN_LENGTH] = '\0';				
				strncpy (temp, line+6, PDB_ATOM_NUM_LENGTH);
				temp [PDB_ATOM_NUM_LENGTH] = '\0';
 				(molecule+counter)->atom_ser = atoi(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+12, PDB_TYPE_LENGTH);
				strncpy (temp2, strtok(temp, delims), PDB_TYPE_LENGTH) ;			// strtok involved to delete empty first columns and spaces
				if ( temp2[0] == '0' || temp2[0] == '1' || temp2[0] == '2' || temp2[0] == '3' || temp2[0] == '4' || temp2[0] == '5' || temp2[0] == '6' || temp2[0] == '7' || temp2[0] == '8' || temp2[0] == '9' )
					strncpy ((molecule+counter)->pdb_type, temp2+1, PDB_TYPE_LENGTH); // delete numerical first columns
					else strncpy ((molecule+counter)->pdb_type, temp2, PDB_TYPE_LENGTH) ;	
				str0 (temp, NUM_LENGTH);
				str0 (temp2, NUM_LENGTH);
				(molecule+counter)->pdb_type[PDB_TYPE_LENGTH] = '\0';
				(molecule+counter)->pdb_alt_loc [0] = line[16];
				(molecule+counter)->pdb_alt_loc [PDB_ALT_LOC_LENGTH] = '\0';				
				strncpy ((molecule+counter)->res_type, line+17, RES_TYPE_LENGTH);	
				(molecule+counter)->res_type[RES_TYPE_LENGTH] = '\0';
				(molecule+counter)->chain [0] = line[21];
				(molecule+counter)->chain [CHAIN_LENGTH] = '\0';
				strncpy (temp, line+22, PDB_RES_NUM_LENGTH);			
				temp [PDB_RES_NUM_LENGTH] = '\0';
				(molecule+counter)->res_ser = atoi(temp);
				str0 (temp, NUM_LENGTH);
				(molecule+counter)->pdb_achar [0] = line[26];
				(molecule+counter)->pdb_achar [PDB_ACHAR_LENGTH] = '\0';				
				strncpy (temp, line+30, PDB_COORD_LENGTH);
				temp [PDB_COORD_LENGTH] = '\0';					
 				(molecule+counter)->x_coord = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+38, PDB_COORD_LENGTH);
				temp [PDB_COORD_LENGTH] = '\0';					
 				(molecule+counter)->y_coord = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+46, PDB_COORD_LENGTH);
				temp [PDB_COORD_LENGTH] = '\0';					
 				(molecule+counter)->z_coord = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+54, PDB_OCC_B_LENGTH);
				temp [PDB_OCC_B_LENGTH] = '\0';					
 				(molecule+counter)->occ = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+60, PDB_OCC_B_LENGTH);
				temp [PDB_OCC_B_LENGTH] = '\0';					
 				(molecule+counter)->b_factor = atof(temp);
				str0 (temp, NUM_LENGTH);
				counter++;
				}
			str0 (line, MAX_LINE_LENGTH); 
			}	
		*(int *)atom_num = counter;


		if ( (g_v == 'v' || g_v == 'd') && (molecule+counter-1)->atom_ser != *(int *)atom_num ) {
			printf("%s%s\n","\nWARNING! Number of ATOM/HETATM entries does not match with maximum atom serial number in file ",file_name); 
			printf("%s%d%s","WARNING! Number of ATOM/HETATM entries:   ",*(int *)atom_num,"\n");	
			printf("%s%d%s","WARNING! Maximum atom serial number:      ",(molecule+counter-1)->atom_ser,"\n");			
			printf("WARNING! This may happen e.g. if a unique atom serial number is assigned to TER tokens and/or >= 100000 atoms are listed.\n");
			}	
				
		// Check for alternate locations

		counter = 0;
		for (i = 0; i < *(int *)atom_num; i++)
			if ( (molecule+i)->pdb_alt_loc [0] != ' ' ) {
				if ( counter == 0 ) printf("%s%s%s", "ERROR! Problems occured during processing file ",file_name," !\n"); 
				printf("%s%d%s%d%s", "ERROR! Alternate location # ",counter+1," was detected at atom # ",(molecule+i)->atom_ser," !\n");
				counter++;
				}
		if (counter > 0) {
			printf("%s","Select only one location per atom and re-run the program!\n"); exit (1);
			}				
		fclose(f_pdb);
		return (molecule); 

	} else return (NULL); 
	}

char find_token_in_name (char name[MAX_FILENAME_LENGTH], char token [MAX_FILENAME_LENGTH], char *delims) {

	char *temp;
	int token_length = strlen(token);
	char token_indicator = 'n';	// not found

	temp = strtok(name, delims);

	while ( temp != NULL ) {
		if ( strncmp(temp, token, token_length) == 0 ) token_indicator = 'y'; 	// found
		temp = strtok(NULL, delims);
		}

	return token_indicator;
	}

void str0 (char string [], int length) {
	int i;
	for (i = 0; i < length; i++)
		string[i] = '\0';
	} 
