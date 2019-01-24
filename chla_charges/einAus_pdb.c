#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/****************************************************************

Programm zum formatierten ausgeben

Aufruf: "einAus in.dat"

*****************************************************************/

#define ECHARGE 1.60217733*pow(10,-19)

main(int argc, char *argv[ ]){       

  FILE *in_file;
  char *cmd, *in_dat_1;
  int i, ch, n1=0;

  if(argc<2)
  {
    printf( "Missing argument\n" );
    return 1;
  }

  /*
  if((cmd = malloc(strlen(argv[1])+1)) == NULL)
  {
    printf( "Out of memory\n" );
    return 1;
  }
  */

  in_dat_1 = argv[1];


  // === Zählen der Zeilen ===========================

  in_file = fopen(in_dat_1, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_1);
    exit(8);
  }
  while (1) {            //zaehle zeilenzahl = atomzahl
    ch = fgetc(in_file);
    if (ch == '\n')
      ++n1;
    if (ch == EOF)
      break;
  } 
  fclose(in_file);


  double x[n1], y[n1], dydx[n1];

  char atom_type[n1][5];
  char amac_type[n1][4];
  char atom[n1][5];
  char str[2]="N", str2[3];
  int atom_num[n1], res_indx[n1], l1[n1], l2[n1];
  double xyz[n1][3], q[n1], r[n1];


  // === Einlesen des Files =========================


  in_file = fopen(in_dat_1, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_1);
    exit(8);
  }

  for(i=0; i<n1; i++){
    fscanf(in_file, "%d %s  %lf %lf %lf\n", 
	   &atom_num[i], atom_type[i], &xyz[i][0], &xyz[i][1], &xyz[i][2]); 
  }

  fclose(in_file);  


 // === Formatierte Ausgabe =========================

  for(i=0; i<n1; i++)
    printf("%5d %5s %5s %5d   %8.3lf %8.3lf %8.3lf \n", 
  	   i+1081, atom_type[i], "CLI", 10, xyz[i][0], xyz[i][1], xyz[i][2]); 


}
