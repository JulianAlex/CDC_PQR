#include "cdcProt_head.h"
 
void initResLeng(int *res_leng, double *res_q){ 

  res_leng[0]  = 10;  res_q[0]  =  0;  // ALA 0  
  res_leng[1]  = 24;  res_q[1]  = +1;  // ARG 1  
  res_leng[2]  = 14;  res_q[2]  =  0;  // ASN 2  
  res_leng[3]  = 12;  res_q[3]  = -1;  // ASP 3   
  res_leng[4]  = 11;  res_q[4]  =  0;  // CYS 4  
  res_leng[5]  = 17;  res_q[5]  =  0;  // GLN 5  
  res_leng[6]  = 15;  res_q[6]  = -1;  // GLU 6  
  res_leng[7]  =  7;  res_q[7]  =  0;  // GLY 7  
  res_leng[8]  = 17;  res_q[8]  =  0;  // HSD 8  
  res_leng[9]  = 18;  res_q[9]  = +1;  // HSP 9  
  res_leng[10] = 19;  res_q[10] =  0;  // ILE 10 
  res_leng[11] = 19;  res_q[11] =  0;  // LEU 11 
  res_leng[12] = 22;  res_q[12] = +1;  // LYS 12 
  res_leng[13] = 17;  res_q[13] =  0;  // MET 13 
  res_leng[14] = 20;  res_q[14] =  0;  // PHE 14 
  res_leng[15] = 14;  res_q[15] =  0;  // PRO 15 
  res_leng[16] = 11;  res_q[16] =  0;  // SER 16 
  res_leng[17] = 14;  res_q[17] =  0;  // THR 17 
  res_leng[18] = 24;  res_q[18] =  0;  // TRP 18 
  res_leng[19] = 21;  res_q[19] =  0;  // TYR 19 
  res_leng[20] = 16;  res_q[20] =  0;  // VAL 20 

  res_leng[21] = 98;  res_q[21] = -0;  // LUT 21 
  res_leng[22] = 100; res_q[22] = -0;  // NEX 22 
  res_leng[23] = 102; res_q[23] = -1;  // PGD 23 
  res_leng[24] = 100; res_q[24] =  0;  // XAT 24 

  res_leng[25] = 30;  res_q[25] = +1;  // ARG_2, 25 
  res_leng[26] = 17;  res_q[26] =  0;  // SER_2, 26
 
}

//========================================================================================

void assignProtResID(int n_atom, int n_res, int *atmsPerRes, int *res_indx, int *resid, 
		     int *resid_at, char resi_type[n_atom][5], int *res_leng){

  int i = 0;
  int count = 0;    

  while(i < n_atom)
  {
    if( strcmp( resi_type[i], "ALA") == 0) 
      resid[count] = 0;
    else if( strcmp( resi_type[i], "ARG") == 0)
      resid[count] = 1;
    else if( strcmp( resi_type[i], "ASN") == 0)  
      resid[count] = 2;
    else if( strcmp( resi_type[i], "ASP") == 0) 
      resid[count] = 3;
    else if( strcmp( resi_type[i], "CYS") == 0)
      resid[count] = 4;
    else if( strcmp( resi_type[i], "GLN") == 0)
      resid[count] = 5;
    else if( strcmp( resi_type[i], "GLU") == 0)
      resid[count] = 6;
    else if( strcmp( resi_type[i], "GLY") == 0)
      resid[count] = 7;
    else if( strcmp( resi_type[i], "HSD") == 0)
      resid[count] = 8;
    else if( strcmp( resi_type[i], "HSP") == 0)
      resid[count] = 9;
    else if( strcmp( resi_type[i], "ILE") == 0)
      resid[count] = 10;
    else if( strcmp( resi_type[i], "LEU") == 0)
      resid[count] = 11;
    else if( strcmp( resi_type[i], "LYS") == 0)
      resid[count] = 12;
    else if( strcmp( resi_type[i], "MET") == 0)
      resid[count] = 13;
    else if( strcmp( resi_type[i], "PHE") == 0)
      resid[count] = 14;
    else if( strcmp( resi_type[i], "PRO") == 0)
      resid[count] = 15;
    else if( strcmp( resi_type[i], "SER") == 0)
      resid[count] = 16;
    else if( strcmp( resi_type[i], "THR") == 0)
      resid[count] = 17;
    else if( strcmp( resi_type[i], "TRP") == 0)
      resid[count] = 18;
    else if( strcmp( resi_type[i], "TYR") == 0)
      resid[count] = 19;
    else if( strcmp( resi_type[i], "VAL") == 0)
      resid[count] = 20;
    else if( strcmp( resi_type[i], "LUT") == 0)
      resid[count] = 21;
    else if( strcmp( resi_type[i], "NEX") == 0)
      resid[count] = 22;
    else if( strcmp( resi_type[i], "PGD") == 0)  
      resid[count] = 23;                         
    else if( strcmp( resi_type[i], "XAT") == 0)  
      resid[count] = 24;             
    else{ 
      printf("Unknown Residue  %s ! Exit.\n", resi_type[i]);
      exit(EXIT_FAILURE);
    }        
 
    //    printf("%8d %8d %8s %8d %8d %8d\n", 
    //	   i, count, resi_type[i], resid[count], atmsPerRes[count], res_leng[resid[count]]);
    
    i += atmsPerRes[count];   
    count++;  
   
  } 


  // check length of residues. 
  
  count = 0; i=0;
 
  while(i < n_atom)
  {
    // printf("%8d %8d %8s %8d %8d %8d\n", 
    //	   i, count, resi_type[i], resid[count], atmsPerRes[count], res_leng[resid[count]]);

    if(atmsPerRes[count] != res_leng[resid[count]]){
    {
      printf("New Res-ID for Residue  %6d  %4s\n", i, resi_type[i] );
      if ( strcmp( resi_type[i], "ARG") == 0) 
      {
	resid[count] = 25;
      }
      else if ( strcmp( resi_type[i], "SER") == 0) 
      {
      	resid[count] = 26;
      }
      else 
      {
	printf("\nWarning: A Residue has wrong length!!  %6d  %4s\n\n", i, resi_type[i]); 
      }
    }

    }
    i += atmsPerRes[count];   
    count++;  
  }
  
}

//========================================================================================================

void initPigmentLeng(int *res_leng){ 

  res_leng[0]  = 137;  // CLA 0  
  res_leng[1]  = 125;  // CLB 1  
  res_leng[2]  = 125;  // CLC 2  
  res_leng[3]  =  92;  // CLD 3   
  res_leng[4]  = 110;  // CLE 4  
  res_leng[5]  =  78;  // CLF 5  
  res_leng[6]  = 122;  // CLG 6  
  res_leng[7]  = 116;  // CLH 7  
  res_leng[8]  =  86;  // CLI 8  
  res_leng[9]  =  89;  // CLJ 9  
  res_leng[10] =  77;  // CHA 10 
  res_leng[11] = 127;  // CHB 11 
  res_leng[12] = 124;  // CHC 12 
  res_leng[13] = 121;  // CHD 13 
  res_leng[14] = 130;  // CHE 14 
  res_leng[15] =  97;  // CHF 15 
  res_leng[16] =  73;  // CHG 16 
  res_leng[17] = 106;  // CHH 17 
  res_leng[18] =  72;  // CHI 18 
  res_leng[19] = 112;  // CHJ 19 
  res_leng[20] = 136;  // CHL 20 
}

//========================================================================================

void assignPigmResID(int n_atom, int n_res, int *atmsPerRes, int *res_indx, int *resid, 
		     int *resid_at, char resi_type[n_atom][5], int *res_leng){

  int i=0, count=0;    

  while(i < n_atom)
  {
    if(      strcmp( resi_type[i], "CLA") == 0) 
      resid[count] = 0;
    else if( strcmp( resi_type[i], "CLB") == 0)
      resid[count] = 1;
    else if( strcmp( resi_type[i], "CLC") == 0)  
      resid[count] = 2;
    else if( strcmp( resi_type[i], "CLD") == 0) 
      resid[count] = 3;
    else if( strcmp( resi_type[i], "CLE") == 0)
      resid[count] = 4;
    else if( strcmp( resi_type[i], "CLF") == 0)
      resid[count] = 5;
    else if( strcmp( resi_type[i], "CLG") == 0)
      resid[count] = 6;
    else if( strcmp( resi_type[i], "CLH") == 0)
      resid[count] = 7;
    else if( strcmp( resi_type[i], "CLI") == 0)
      resid[count] = 8;
    else if( strcmp( resi_type[i], "CLJ") == 0)
      resid[count] = 9;
    else if( strcmp( resi_type[i], "CHA") == 0)
      resid[count] = 10;
    else if( strcmp( resi_type[i], "CHB") == 0)
      resid[count] = 11;
    else if( strcmp( resi_type[i], "CHC") == 0)
      resid[count] = 12;
    else if( strcmp( resi_type[i], "CHD") == 0)
      resid[count] = 13;
    else if( strcmp( resi_type[i], "CHE") == 0)
      resid[count] = 14;
    else if( strcmp( resi_type[i], "CHF") == 0)
      resid[count] = 15;
    else if( strcmp( resi_type[i], "CHG") == 0)
      resid[count] = 16;
    else if( strcmp( resi_type[i], "CHH") == 0)
      resid[count] = 17;
    else if( strcmp( resi_type[i], "CHI") == 0)
      resid[count] = 18;
    else if( strcmp( resi_type[i], "CHJ") == 0)
      resid[count] = 19;
    else if( strcmp( resi_type[i], "CHL") == 0)
      resid[count] = 20;
    else{ 
      printf("Unknown Residue  %s ! Exit.\n", resi_type[i]);
      exit(EXIT_FAILURE);
    }        
 
    //printf("%8d %8d %8s %8d %8d %8d\n", 
    //	   i, count, resi_type[i], resid[count], atmsPerRes[count], res_leng[resid[count]]);

    i += atmsPerRes[count];   
    count++;     
  } 
  
  // check length of residues: 
  
  count = 0; i=0;
 
  while(i < n_atom)
  {
    // printf("%8d %8d %8s %8d %8d %8d\n", 
    //	   i, count, resi_type[i], resid[count], atmsPerRes[count], res_leng[resid[count]]);

    if(atmsPerRes[count] != res_leng[resid[count]])
    {
      printf("\nWarning: A Residue has wrong length!!  %6d  %4s\n\n", i, resi_type[i]); 
    }
    
    i += atmsPerRes[count];   
    count++;  
  }

}

//================================================================================================

void assignPigmCharges(int n_char, int n_atom, 
		       double *q_gr_in, double *q_ex_in, double *q_gr_out, double *q_ex_out, 
		       char pigm_type_char[n_char][5], char atom_type_char[n_char][5], 
		       char atom_type[n_atom][5], char pigm_type[n_atom][5]){
  int i=0, j=0;

  for(i=0; i < n_atom; i++) 
  {
    for(j=0; j < n_char; j++)
    {
      if( (strcmp( pigm_type[i], pigm_type_char[j]) == 0) && (strcmp( atom_type[i], atom_type_char[j]) == 0) )
      {
	q_gr_out[i] = q_gr_in[j];
	q_ex_out[i] = q_ex_in[j];
      }	
    } 
  }   
}    

//=========================================================================================================

// Check if pigments got correct charges 

void outputPigmKoorChar(char *out_dat, int nn, double **xyz, double *q_gr, double *q_ex, 
		       char atom_type[nn][5], char pigm_type[nn][5])
{
  FILE *fp;
  int i;

  fp = fopen(out_dat, "w");
  if (fp == NULL) {
    (void)printf("Can not open %s\n", out_dat);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<nn; i++){
    fprintf(fp, "%5s %5s   %8.3lf %8.3lf %8.3lf   %8.3lf  %8.3lf\n", 
	   atom_type[i], pigm_type[i], xyz[i][0], xyz[i][1], xyz[i][2], q_gr[i], q_ex[i]); 
  }
  fclose(fp);

}

//================================================================================================
