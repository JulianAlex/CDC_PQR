# CDC_PQR
charge density coupling, here applied to CP29 (part of photosystem PSII, cp29.pqr)

 --------------------------------------------------------------------------------------------------
 
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
  


  ==================================================================================================



  CDC-Programm für pqr-files, Julian Adolphs, 2018


  Kompilieren mit:  "make" 

  Ausführen mit:  "./cdcProt"  bzw.  "./cdcProt > ausgabe.out"

  --------------------------------------------------------------------------------------------------

  Programm besteht aus: 

  Makefile 
  cdcProt_main.c, cdcProt_input.c, cdcProt_attach.c, cdcProt_calculate.c
  cdcProt_head.h, cdcProt_functions.h 


  Eingabe-Files: 

  1) pqr-File, der das Protein, Carotinoide, Lipide, etc UND die Pigmente enthaelt. (dat_1, s.u.)
  2) ein File, der alle Pigment-Ladungen enthaelt (dat_2, s.u.)


  Ausgabe-Files:

  3) "pigment_koord.dat" enthaelt alle Pigment-Koordinaten. Wird wieder eingelesen. 
  4) "protein_koord.dat" enthaelt alle Protein-Koordinaten. Wird wieder eingelesen.

  5) "pigm_koor_char.dat" dient nur als Kontroll-Ausgabe. Hier kann man checken, ob den Pigmenten 
      die Ladungen korrekt zugeordnet werden. 
