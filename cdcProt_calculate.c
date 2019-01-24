#include "cdcProt_head.h"
 

void calcChargeSum(int n_atom, int n_res, double *q_prot, int *resid, 
		   int *res_leng, double *res_q, char resi_type_resi[n_res][5]){

  int i = 0; 
  int j = 0; 
  int sum_over = 0;

  double q_sum[n_res];
  for(i=0; i<n_res; i++)
    q_sum[i] = 0.0;
  

  int count = 0;
  for(i = 0; i < n_res; i++){       // 1050
    sum_over = res_leng[resid[i]];  // <=24
    for(j = 0; j < sum_over; j++){
      q_sum[i] += q_prot[count]; 
      count++;
    }
    //    printf("%8d %8d %8d %10.3lf  %10.3lf\n", i, j, resid[i], res_q[resid[i]], q_sum[i]);
  }
   
  
  printf("\n");
  for(i=0; i<n_res; i++){
    if( ABS(q_sum[i] - res_q[resid[i]]) > 0.0001 )
          printf("error in residue:  %s, %3d, type = %2d  charge_sum = %6.3lf, %6.3lf\n", 
		 resi_type_resi[i], i+1, resid[i], q_sum[i], res_q[resid[i]]);
  }
  printf("\n");
  
}

//==============================================================================================


void calcEnergShift(int n_atom_prot, int n_res, int n_atom_pigm, int n_pig, int *atmsPerPigm,  
		    int *atmsPerRes, double **koor_prot, double **koor_pigm, 
		    double *q_prot, double *q_gr, double *q_exc, int *typ_at,
		    int *resid, int *res_leng, char amac_type[n_atom_prot][5], 
		    char atom_type[n_atom_prot][5]){

  // Pigments = PartialCharges  &  AminoAcids = PartialCharges 

  int i=0, j=0, k=0, l=0, m=0, n=0, nrow=0, pigm_leng=0, count=0, count_m=0, count_j=0;

  double rrij=0.0, rij=0.0, qi=0.0, qj=0.0, sum=0.0; 

  double energy[n_pig], ener_amac[n_pig], ener_pigm[n_pig]; 
  double ener_pigpig[n_pig][n_pig];

  double **sum_array = calloc(n_pig, sizeof(double*));
  double **sum_array_2 = calloc(n_pig, sizeof(double*));
  nrow = n_pig;
  while(nrow--){
    sum_array[nrow]   = calloc(n_atom_prot, sizeof(double));
    sum_array_2[nrow] = calloc(n_res, sizeof(double));
  }

  for(m=0; m < n_pig; m++){
    energy[m]     = 0.0; 
    ener_amac[m] = 0.0; 
    ener_pigm[m]  = 0.0;
    for(n=0; n < n_pig; n++)
      ener_pigpig[m][n] = 0.0;
  }    


  // ---- interaction of pigment (Delta q) with charges of proteins: ---------------------------------------------

  count = 0; 
  
  for(m=0; m < n_pig; m++){  // loop over all pigments 

    pigm_leng = atmsPerPigm[m];

    for(i=0; i < pigm_leng; i++){  // loop over all pigments-partial-charges of pigment m
 
      qi = (q_exc[count+i] - q_gr[count+i])*EL_CHAR;  // delt_mu = mu_exc - mu_gr   
   
      for(j=0; j < n_atom_prot; j++){  // loop over all protein-partial-charges
      
	rrij = 0.0;
	for(l=0; l<3; l++){
	  //rij   = koor_pigm[i + m*pigm_leng][l] - koor_prot[j][l];
          rij   = koor_pigm[i + count][l] - koor_prot[j][l];
	  rrij += SQ(rij); 
	}
	rrij = sqrt(rrij)*ANG;      // rij in Angstrom
      
	qj = q_prot[j]*EL_CHAR;         // qi, qj in Coulomb

	sum = 1/(4*PI*EPS_0)*qi*qj/rrij; // in Joule

	ener_amac[m] += EVWZ/EL_CHAR*sum;      // in 1/cm  ( Joule => eV => 1/cm )	

	//if( (m == 0))// && (i==0) )
	//  printf("%3d %3d %5d  %10.3lf %10.3lf %10.3lf  %5s   %6s  %10.3lf    %10.3lf\n", 
	// m, i, j, qi/EL_CHAR, koor_pigm[i + m*pigm_leng][0], koor_prot[j][0], amac_type[j], atom_type[j], q_prot[j], ener_amac[m]);

	sum_array[m][j] += EVWZ/EL_CHAR*sum; 

      } //end j                

    } //end i

    count += pigm_leng;

    //printf("%3d  %3d  %10.3lf\n", m, i, ener_amac[m]); 

  } //end m     
  
  printf("\n");




  //---- interaction of pigments m with groundstatecharges of all other pigments: ------------------------

  count_m=0;    //    pigm_leng = atmsPerPigm[m];
  count_j=0;

  for(m=0; m < n_pig; m++)  // loop over all pigments 
  {
    count_j=0;
    for(j=0; j < n_pig; j++)    // loop over all pigments (without m)
    { 
      if( j != m )
      {  
	for(i=0; i < atmsPerPigm[m]; i++)  // loop over all pigment-partial-charges of pigment m 
	{
	  for(k=0; k < atmsPerPigm[j]; k++)  // loop over all pigment-partial-charges of pigment j
	  {
	    rrij = 0.0;
	    for(l=0; l<3; l++)
	    {
	      rij   = koor_pigm[i+count_m][l] - koor_pigm[k+count_j][l];      
	      rrij += SQ(rij); 
	    } 
	    rrij = sqrt( rrij )*ANG;      // rij in meter
	  	
	    qi = (q_exc[i+count_m] - q_gr[i+count_m])*EL_CHAR;  // Pigm m,  delt_mu = mu_exc - mu_gr 
	    qj = q_gr[k+count_j]*EL_CHAR;                       // Pigm j,  qi, qj in Coulomb  
	    
	    sum = 1/(4*PI*EPS_0)*qi*qj/rrij; // in Joule

	    ener_pigm[m] += EVWZ/EL_CHAR*sum;      // in 1/cm  ( Joule => eV => 1/cm )

	    ener_pigpig[m][j] += EVWZ/EL_CHAR*sum;  

	    //printf("%d  %d %d %d  %9.3lf  %9.3lf  %9.3lf  %9.3lf  %9.3lf  %9.3lf\n", 
	    //	   m, j, i, k, koor_pigm[k+count_j][0], rrij/ANG, qi/EL_CHAR, qj/EL_CHAR, EVWZ/EL_CHAR*sum, ener_pigm[m]);

	  } // end_for k       
	} // end_for i
      }  // end_if
      count_j += atmsPerPigm[j];
    }  // end_for j
    count_m += atmsPerPigm[m]; 
  } // end_for m

  printf("\n");

  //-------------------------------------------------------------------------------------------------------------------

  printf("Large interactions (>= %2.0lf) between Pigments n and m  [cm-1]\n\n", LARGE_COUPLING);
  for(m=0; m < n_pig; m++)
    for(n=0; n < m; n++)
      if( ABS(ener_pigpig[m][n]) >= LARGE_COUPLING)
	printf("%3d  %3d  %10.1lf\n", n+1, m+1, ener_pigpig[m][n]); 
  printf("\n"); 

  
  //--- interaction of single amino acids with pigments -------------------------------------------------


  /*
  for(m=0; m < n_pig; m++)
    for(j=0; j < n_atom_prot; j++)
      printf("%5d %5d %8.3lf\n", m, j, sum_array[m][j]);
  printf("\n"); 
  */
  
  for(m=0; m < n_pig; m++)
  {
    count=0;
    for(j=0; j < n_res; j++)
    {
      for(k=0; k < atmsPerRes[j]; k++)
      {
	sum_array_2[m][j] += sum_array[m][count+k]; 
	//	printf("%3d  %3d  %3d  %3d  %10.3lf\n", m, j, k, count, sum_array[m][count+k]);
      }
      count += atmsPerRes[j];
    }
  }
   
  printf("\n\n");
  printf(" ---  Contribution of single amino acids [cm-1] > MIN_CONTRIB = %3.0lf  --------------\n", MIN_CONTRIB );  
  printf("\n");

  for(m=0; m < n_pig; m++)
    for(j=0; j < n_res; j++)
      if(ABS(sum_array_2[m][j]) > MIN_CONTRIB)
	printf("%3d  %3d  %10.1lf\n", m+1, j+1, sum_array_2[m][j]); 
  

  //-------------------------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------


  for(i=0; i < n_pig; i++)
    energy[i] = ener_amac[i] + ener_pigm[i];

  printf("\n\n");
  printf(" --- pigments-partial-charges & protein-partial-charges [cm-1] --------------\n");  
  printf("\n");
  printf("        energy_prot  enery_pigm   energ(prot+pigm)  [cm-1]\n"); printf("\n");
  for(i=0; i < n_pig; i++)
    printf("%3d  %10.1lf  %10.1lf  %10.1lf\n", 
	   i+1, ener_amac[i], ener_pigm[i], energy[i]);
  printf("\n\n");  



  //-------------------------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------
   

  nrow = n_pig;
  while(nrow--)    
  {    
    free(sum_array[nrow]);
    free(sum_array_2[nrow]);
  }
  free(sum_array); sum_array = NULL;
  free(sum_array_2); sum_array_2 = NULL;



} // END CALCULATE


