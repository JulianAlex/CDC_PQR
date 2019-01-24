#include "cdcProt_head.h"


int eingabeCountInputLines(char *in_dat)
{

  FILE *in_file;

  int ch, n=0;
 
  in_file = fopen(in_dat, "r");
  if (in_file == NULL) 
  {
    (void)printf("Can not open %s\n", in_dat);
    exit(EXIT_FAILURE);
  }

  while (1) {            
    ch = fgetc(in_file);
    if (ch == '\n')
      ++n;
    if (ch == EOF)
      break;
  }
  
  fclose(in_file);

  printf("Input-File '%s' has %6d Lines\n\n", in_dat, n); 

  return n;
}


// ===========================================================================================================


void eingabePqrFile(char *in_dat, char *out_dat_1, char *out_dat_2, int n_lines)
{

  FILE *fp, *fp1, *fp2;   // file_pointer 

  double f1;
  int i, atom_num; 
  int res_indx[n_lines];
  double xyz[n_lines][3], q[n_lines]; 
  char str1[5], str2[5], chain[2];
  char atom_type[n_lines][5];           
  char res_type[n_lines][5];           


  for(i=0; i<n_lines; i++)
    res_indx[i] = 0;

  fp = fopen(in_dat, "r");
  if (fp == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<n_lines; i++)
  {
    fscanf(fp, "%s %d %s %s %s %d   %lf %lf %lf %lf %lf %s\n", 
	   str1, &atom_num, atom_type[i], res_type[i], chain, &res_indx[i], 
	   &xyz[i][0], &xyz[i][1], &xyz[i][2], &q[i], &f1, str2); 
  }
  fclose(fp);


  fp1 = fopen(out_dat_1, "w");
  if (fp1 == NULL) {
    (void)printf("Can not open %s\n", out_dat_1);
    exit(EXIT_FAILURE);
  }
  fp2 = fopen(out_dat_2, "w");
  if (fp2 == NULL) {
    (void)printf("Can not open %s\n", out_dat_2);
    exit(EXIT_FAILURE);
  }

  // Pigments and protein, lipids, carotinoids discriminate and output

 for(i=0; i<n_lines; i++)
  {
    if((strcmp(res_type[i],"CLA")==0)||(strcmp(res_type[i],"CLB")==0)||(strcmp(res_type[i],"CLC")==0)||
       (strcmp(res_type[i],"CLD")==0)||(strcmp(res_type[i],"CLE")==0)||(strcmp(res_type[i],"CLF")==0)||
       (strcmp(res_type[i],"CLG")==0)||(strcmp(res_type[i],"CLH")==0)||(strcmp(res_type[i],"CLI")==0)||
       (strcmp(res_type[i],"CLJ")==0)||(strcmp(res_type[i],"CHA")==0)||(strcmp(res_type[i],"CHB")==0)||
       (strcmp(res_type[i],"CHC")==0)||(strcmp(res_type[i],"CHD")==0)||(strcmp(res_type[i],"CHE")==0)||
       (strcmp(res_type[i],"CHF")==0)||(strcmp(res_type[i],"CHG")==0)||(strcmp(res_type[i],"CHH")==0)||
       (strcmp(res_type[i],"CHI")==0)||(strcmp(res_type[i],"CHJ")==0)||(strcmp(res_type[i],"CHL")==0) )
    {
      fprintf(fp1, "%5s %5s %5d   %8.3lf %8.3lf %8.3lf\n", 
	      atom_type[i], res_type[i], res_indx[i], xyz[i][0], xyz[i][1], xyz[i][2]); 
    }
    else 
      fprintf(fp2, "%5s %5s %5d   %8.3lf %8.3lf %8.3lf  %8.3lf\n", 
	      atom_type[i], res_type[i], res_indx[i], xyz[i][0], xyz[i][1], xyz[i][2], q[i]); 
  }

  fclose(fp1); 
  fclose(fp2); 


}


//=========================================================================================================


int eingabeProtKoord(char *in_dat, int n_atom, double **xyz, double *q, int *res_indx, 
		     char atom_type[n_atom][5], char amac_type[n_atom][5])
{
  FILE *fp_in;

  int i, k, nres, new_indx[n_atom];

  for(i=0; i<n_atom; i++)
    new_indx[i] = 0;


  fp_in = fopen(in_dat, "r");
  if (fp_in == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<n_atom; i++){
    fscanf(fp_in, "%s %s %d   %lf %lf %lf %lf \n", 
	   atom_type[i], amac_type[i], &res_indx[i], &xyz[i][0], &xyz[i][1], &xyz[i][2], &q[i]); 
  }
  fclose(fp_in);



  // Residue numbers need to be corrected
  // (in output file proteins and pigments are mixed, interrupting the running numbers) 

  k = 1;
  while (k < n_atom)
  {
    if (res_indx[k] > res_indx[k-1]) 
    {
      for(i=k; i<n_atom; i++)
	new_indx[i] = new_indx[k-1] + 1;
    }
    k++;
  }  
 

  //  for(i=0; i<n_atom; i++)
  //    printf("%5d %5d %5d\n", i, res_indx[i], new_indx[i]); 
 

  for(k=0; k<n_atom; k++)
  {
    res_indx[k] = new_indx[k];
  }  
  

  if(res_indx[0] == 1)
  {
    nres = res_indx[n_atom-1];

    for(i=0; i<n_atom; i++)
      res_indx[i] = res_indx[i] - 1;
  }
  else if( res_indx[0] == 0 ) 
  {
    nres = res_indx[n_atom-1] + 1;
  }
  else
  {
    printf("Something wrong with numbering in file %s", in_dat);
    exit(EXIT_FAILURE);
  }


  fp_in = fopen(in_dat, "w");
  if (fp_in == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<n_atom; i++){
    fprintf(fp_in, "%5s %5s %5d   %8.3lf %8.3lf %8.3lf  %8.3lf\n", 
	    atom_type[i], amac_type[i], res_indx[i], xyz[i][0], xyz[i][1], xyz[i][2], q[i]); 
  }
  fclose(fp_in);

  printf("Number of Protein-Residues: %6d\n\n", nres); 

  return(nres);

}

//=========================================================================================================


int eingabePigmKoord(char *in_dat, int n_atom, double **xyz, int *res_indx, 
		     char atom_type[n_atom][5], char amac_type[n_atom][5])
{
  FILE *fp_in;

  int i, k, nres, new_indx[n_atom];

  for(i=0; i<n_atom; i++)
    new_indx[i] = 0;

  fp_in = fopen(in_dat, "r");
  if (fp_in == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<n_atom; i++){
    fscanf(fp_in, "%s %s %d   %lf %lf %lf \n", 
	   atom_type[i], amac_type[i], &res_indx[i], &xyz[i][0], &xyz[i][1], &xyz[i][2]); 
  }
  fclose(fp_in);


  // Residue-Nummers need correction
  // (in output file proteins and pigments are mixed, interrupting the running numbers) 

  k = 1;
  while (k < n_atom)
  {
    if (res_indx[k] > res_indx[k-1]) 
    {
      for(i=k; i<n_atom; i++)
	new_indx[i] = new_indx[k-1] + 1;
    }
    k++;
  }  
 

  //  for(i=0; i<n_atom; i++)
  //    printf("%5d %5d %5d\n", i, res_indx[i], new_indx[i]); 
 

  for(k=0; k<n_atom; k++)
  {
    res_indx[k] = new_indx[k];
  }  
  

  if(res_indx[0] == 1)
  {
    nres = res_indx[n_atom-1];

    for(i=0; i<n_atom; i++)
      res_indx[i] = res_indx[i] - 1;
  }
  else if( res_indx[0] == 0 ) 
  {
    nres = res_indx[n_atom-1] + 1;
  }
  else
  {
    printf("Something wrong with numbering in file %s", in_dat);
    exit(EXIT_FAILURE);
  }


  fp_in = fopen(in_dat, "w");
  if (fp_in == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<n_atom; i++){
    fprintf(fp_in, "%5s %5s %5d   %8.3lf %8.3lf %8.3lf \n", 
	    atom_type[i], amac_type[i], res_indx[i], xyz[i][0], xyz[i][1], xyz[i][2]); 
  }
  fclose(fp_in);

  printf("Number of Pigments: %6d\n\n", nres); 

  return(nres);

}


// =====================================================================================================


void countAtmsPerResidue(int n_atom, int n_res, int *atmsPerRes, double **xyz, int *res_indx, char atom_type[n_atom][5], 
			  char amac_type[n_atom][5], char resi_type_resi[n_res][5])
{

  double res_idx;
  int k, cnt_atm=0, cnt_res=0;

  res_idx = res_indx[0];

  strcpy(resi_type_resi[0], amac_type[0]);  //printf("%s %s\n", resi_type_resi[0], amac_type[0]); 

  for(k=0; k<n_atom; k++)
  { 
    if(res_idx == res_indx[k]) 
      cnt_atm++;
    else
    {
      res_idx = res_indx[k];
      atmsPerRes[cnt_res] = cnt_atm;
      cnt_atm = 1;  
      cnt_res++; 
      strcpy(resi_type_resi[cnt_res], amac_type[k]);
    }
  }
  atmsPerRes[cnt_res] = cnt_atm;  // otherwise last residue gets atomnumber 0


  //  for(k=0; k<n_res; k++)
  //    printf("%3d %5d   %3s\n", k, atmsPerRes[k], resi_type_resi[k] );  printf("\n");

}

// ============================================================================================

void eingabePigmentCharges(char *in_dat, int n_char, double *q_gr, double *q_ex, 
			   char atom_type[n_char][5], char pigm_type[n_char][5])
{ 

  FILE *fp_in;

  double q_temp_diff[n_char];
  int i, num1, num2, temp_id[N_PIGM_CHAR];
    
  for(i=0; i < N_PIGM_CHAR; i++)
    temp_id[i] = 0;

  fp_in = fopen(in_dat, "r");
  if (fp_in == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(8);
  }

  for(i=0; i < n_char; i++)
  {
    fscanf(fp_in, "%d %s %s %d  %lf %lf %lf\n", 
	   &num1, atom_type[i], pigm_type[i], &num2, &q_gr[i], &q_ex[i], &q_temp_diff[i]); 
  }
  fclose(fp_in);


}


