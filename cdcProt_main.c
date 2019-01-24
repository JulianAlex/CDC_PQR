/****************************************************************************************************
   
  CDC-Program for pqr-files, Julian Adolphs, 2018

  Compile with:  "make" 

  Execute with:  "./cdcProt"  or  "./cdcProt > output.txt"

  --------------------------------------------------------------------------------------------------

  Program consists of: 

  Makefile 
  cdcProt_main.c, cdcProt_input.c, cdcProt_attach.c, cdcProt_calculate.c
  cdcProt_head.h, cdcProt_functions.h 


  Input-Files: 

  1) pqr-file, for instance "cp29.pqr", contains protein, carotinoids, lipids, etc and pigments. 
  2) "pigment_charges.dat", contains all pigment charges 

  Output-Files:

  3) "pigment_koord.dat" contains all pigment coordinats,   
  4) "protein_koord.dat" contains all protein coordinats. 

  5) "pigm_koor_char.dat" control-output. 
      check, if pigments got assigned charges correctly 

  -------------------------------------------------------------------------------------------------

  The program computes site energy shifts of pigment transition energies   
  caused by the protein environment (including carotenoids, lipids, other pigments). 

  The program uses the CDC-method (charge density coupling), describend in:
 
  Calculation of pigment transition energies in the FMO protein, Photosynth. Res. (2008) 95:197-209
  J. Adolphs, F. Müh, M.E. Madjet, T. Renger. File: Adol2008.pdf 

****************************************************************************************************/

 
#include "cdcProt_head.h"   
      
int main(int argc, char *argv[ ]){       
     
  int k=0, nrow=0, n_lines=0, n_ges_pigm_atms=0, n_pigm=0, n_ges_prot_atms=0;  
  int n_res=0, n_char=0;

  char *dat_1 = "cp29.pqr";             // pqr-Input-Files
  char *dat_2 = "pigment_charges.dat"; 
      
  char *dat_3 = "pigment_koord.dat";    // output
  char *dat_4 = "protein_koord.dat";  

  char *dat_5 = "pigm_koor_char.dat";   // control output
  
  n_lines = eingabeCountInputLines(dat_1);  

  eingabePqrFile(dat_1, dat_3, dat_4, n_lines);     // input pqr-File
  
  n_char = eingabeCountInputLines(dat_2);           // number of all pigment charges
  n_ges_pigm_atms = eingabeCountInputLines(dat_3);  // number of all pigment atoms
  n_ges_prot_atms = eingabeCountInputLines(dat_4);  // number of protein atoms    
 
  if( n_ges_pigm_atms + n_ges_prot_atms != n_lines)
  {
    printf("Warning: Atoms from input-file lost!! %d + %d != %d\n\n", 
	   n_ges_pigm_atms, n_ges_prot_atms, n_lines); 
    printf("Check in files: %s\n %s\n %s\n\n", dat_1, dat_3, dat_4); 
  }
  

  // allocate arrays (with initialisation):

  double **koor_pigm = calloc(n_ges_pigm_atms, sizeof(double*));
  nrow =  n_ges_pigm_atms;
  while(nrow--)
    koor_pigm[nrow] = calloc(3, sizeof(double));

  double **koor_prot = calloc( n_ges_prot_atms, sizeof(double*));
  nrow =  n_ges_prot_atms;
  while(nrow--)
    koor_prot[nrow] = calloc(3, sizeof(double));

  double *prot_char_ges = calloc( n_ges_prot_atms, sizeof(double));

  double *q_gr_in = calloc( n_char, sizeof(double) );  // from pigment_charges.dat
  double *q_ex_in = calloc( n_char, sizeof(double) );  

  double *q_gr = calloc( n_ges_pigm_atms, sizeof(double)); // assign pigments from pqr-file 
  double *q_ex = calloc( n_ges_pigm_atms, sizeof(double));

  double *resi_q = calloc( N_PROT_RESI, sizeof(double));
  int *resi_leng = calloc( N_PROT_RESI, sizeof(int));

  int *resi_indx = calloc( n_ges_prot_atms, sizeof(int) ); // Res-indx atom-wise
  int *resid_at  = calloc( n_ges_prot_atms, sizeof(int) ); // Res-ID atom-wise

  int *pigm_indx  = calloc( n_ges_pigm_atms, sizeof(int) ); 
  int *pigm_id_at = calloc( n_ges_pigm_atms, sizeof(int) ); // Pigm-ID atom-wise

  int *pigm_leng = calloc( N_PIGM_CHAR, sizeof(int)); 

  char atom_type_prot[n_ges_prot_atms][5];           
  char  res_type_prot[n_ges_prot_atms][5];           

  char   pigm_type_at[n_ges_pigm_atms][5]; // pigment type atom-wise coordinate-file
  char atom_type_pigm[n_ges_pigm_atms][5]; // all pigments from pqr    

  char pigm_type_char[n_char][5];         // pigment type charge-file
  char atom_type_char[n_char][5];         // pigment type charge-file

  for(k=0; k<n_ges_prot_atms; k++)
  {
    strcpy( atom_type_prot[k], "Err");   // init. string-arrays with "Err" 
    strcpy( res_type_prot[k], "Err");    // for finding errors in output easily 
  }                                      
  for(k=0; k<n_ges_pigm_atms; k++) 
  {
    strcpy( pigm_type_at[k], "Err");
    strcpy( atom_type_pigm[k], "Err");
  }
  for(k=0; k<n_char; k++) 
  {
    strcpy( pigm_type_char[k], "Err");
    strcpy( pigm_type_char[k], "Err");
  }


  //--------------------------------- Pigmente: --------------------------------------------------------------------------------

  initPigmentLeng(pigm_leng);  // init. der pigment lengths

  n_pigm = eingabePigmKoord(dat_3, n_ges_pigm_atms, koor_pigm, pigm_indx, atom_type_pigm, pigm_type_at);

  int *atmsPerPigm = calloc( n_pigm, sizeof(int) ); // pigment-lengths  
  int *pigm_id = calloc( n_pigm, sizeof(int) );     // pigment-ID, pigment-wise  
  char pigm_type_pw[n_pigm][5];                     // pigment-type, pigment-wise
 
  for(k=0; k<n_pigm; k++)
  {
    strcpy(pigm_type_pw[k], "Err");   
  }

  countAtmsPerResidue(n_ges_pigm_atms, n_pigm, atmsPerPigm, koor_pigm, pigm_indx, atom_type_pigm, pigm_type_at, pigm_type_pw); 

  eingabePigmentCharges(dat_2, n_char, q_gr_in, q_ex_in, atom_type_char, pigm_type_char); 

  assignPigmResID(n_ges_pigm_atms, n_pigm, atmsPerPigm, pigm_indx, pigm_id, pigm_id_at, pigm_type_at, pigm_leng);

  assignPigmCharges(n_char, n_ges_pigm_atms, q_gr_in, q_ex_in, q_gr, q_ex, 
  		    pigm_type_char, atom_type_char, atom_type_pigm, pigm_type_at);
  
  outputPigmKoorChar(dat_5, n_ges_pigm_atms, koor_pigm, q_gr, q_ex, atom_type_pigm, pigm_type_at);


  // -------------------------------  Protein: ------------------------------------------------------------------------------

  initResLeng(resi_leng, resi_q);  // init. amino acid lengths

  n_res = eingabeProtKoord(dat_4, n_ges_prot_atms, koor_prot, prot_char_ges, resi_indx, atom_type_prot, res_type_prot);
  
  int *atmsPerRes = calloc( n_res, sizeof(int) );   // Residue lengths
  int *resid = calloc( n_res, sizeof(int) );        // residue-id, residue-wise  
  char res_type[n_res][5];                    
 
  for(k=0; k<n_res; k++)
  {
    strcpy(res_type[k], "Err");   
  }

  countAtmsPerResidue(n_ges_prot_atms, n_res, atmsPerRes, koor_prot, resi_indx, atom_type_prot, res_type_prot, res_type); 

  assignProtResID(n_ges_prot_atms, n_res, atmsPerRes, resi_indx, resid, resid_at, res_type_prot, resi_leng);

  calcChargeSum(n_ges_prot_atms, n_res, prot_char_ges, resid, resi_leng, resi_q, res_type);

  
  // ------------- Interaction: -------------------------------------------------------------------------------- 

  calcEnergShift(n_ges_prot_atms, n_res, n_ges_pigm_atms, n_pigm, atmsPerPigm, atmsPerRes, koor_prot, koor_pigm, 
    		 prot_char_ges, q_gr, q_ex, resid_at, resid, resi_leng, res_type_prot, atom_type_prot);


  //============================================================================================================

  printf("\n");
  printf("Number of AminoAcid Atoms, Number of Residues: %6d  %5d\n", 
   n_ges_prot_atms, n_res); 
  printf("Number of Pigment Atoms, Number of Pigments:   %6d  %5d\n", 
	 n_ges_pigm_atms, n_pigm); 
  printf("\n");

  
  //============================================================================================================


  // --- free allocated memory:
  
  free(prot_char_ges); prot_char_ges = NULL;
  free(q_gr_in); q_gr_in = NULL;
  free(q_ex_in); q_ex_in = NULL;
  free(q_gr); q_gr = NULL;
  free(q_ex); q_ex = NULL;
  free(resi_q); resi_q = NULL;
  free(resi_leng); resi_leng = NULL;
  free(resi_indx); resi_indx = NULL;
  free(resid_at); resid_at = NULL;
  free(pigm_indx); pigm_indx = NULL;
  free(pigm_id_at); pigm_id_at = NULL;
  free(pigm_leng); pigm_leng = NULL;
  free(atmsPerPigm); atmsPerPigm = NULL;
  free(pigm_id); pigm_id = NULL;
  free(atmsPerRes); atmsPerRes = NULL;
  free(resid); resid = NULL;

  nrow =  n_ges_pigm_atms;
  while(nrow--)
    free(koor_pigm[nrow]); 
  free(koor_pigm); koor_pigm = NULL;

  nrow =  n_ges_prot_atms;
  while(nrow--)
    free(koor_prot[nrow]);
  free(koor_prot); koor_prot = NULL;

 
  return (0);
  
}


