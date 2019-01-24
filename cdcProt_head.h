#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "cdcProt_functions.h"   // function declarations

#define ANG 1E-10
#define ANG3 1E-30
#define MILLI 0.001
#define PI 3.141592654
#define EPS_0 8.854187817E-12
#define EL_CHAR 1.60217733E-19

#define SQ(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define ZP(x) pow(10,x)
#define ABS(x) sqrt( ((x) * (x)) )
#define EVWZ 8065.54       // (eV 2 1/cm)

#define EPS 1.


#define N_PIGM_CHAR 21   // Number of Pigment Types in Charge-File
#define SIZE_PIGM  137   // Max Number of Atoms per Pigment

#define N_PROT_RESI     27  //25 // Number of Protein Residue Types 
#define SIZE_PROT_RESI 102  // Max Number of Atoms per Protein Residue 

#define MIN_CONTRIB 100.
#define LARGE_COUPLING 10.
