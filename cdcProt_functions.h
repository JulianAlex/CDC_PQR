int eingabeCountInputLines(char *in_dat);

void eingabePqrFile(char *in_dat, char *out_dat_1, char *out_dat_2, int n_lines);

int eingabeProtKoord(char *in_dat, int n_atom, double **xyz, double *q, int *res_indx,char atom_type[n_atom][5], char amac_type[n_atom][5]);

int eingabePigmKoord(char *in_dat, int n_atom, double **xyz, int *res_indx, char atom_type[n_atom][5], char amac_type[n_atom][5]);

void countAtmsPerResidue(int n_atom, int n_res, int *atmsPerRes, double **xyz, int *res_indx, char atom_type[n_atom][5],char amac_type[n_atom][5], char resi_type_resi[n_res][5]);

void eingabePigmentCharges(char *in_dat, int n_char, double *q_gr, double *q_ex,char atom_type[n_char][5], char pigm_type[n_char][5]);

void initResLeng(int *res_leng, double *res_q);

void assignProtResID(int n_atom, int n_res, int *atmsPerRes, int *res_indx, int *resid,int *resid_at, char resi_type[n_atom][5], int *res_leng);

void initPigmentLeng(int *res_leng); 

void assignPigmResID(int n_atom, int n_res, int *atmsPerRes, int *res_indx, int *resid,int *resid_at, char resi_type[n_atom][5], int *res_leng);

void assignPigmCharges(int n_char, int n_atom,double *q_gr_in, double *q_ex_in, double *q_gr_out, double *q_ex_out,char pigm_type_char[n_char][5], char atom_type_char[n_char][5],char atom_type[n_atom][5], char pigm_type[n_atom][5]);

void outputPigmKoorChar(char *out_dat, int nn, double **xyz, double *q_gr, double *q_ex,char atom_type[nn][5], char pigm_type[nn][5]);

void assignProtCharges(int n_prot_char, int n_res, int n_atom_amac, int *resid, int *res_leng,double **q_amac, double *q_prot, char atom_type_amac[n_prot_char][5],char atom_type[n_atom_amac][5], char resi_type[n_atom_amac][5]);

void calcChargeSum(int n_atom, int n_res, double *q_prot, int *resid,int *res_leng, double *res_q, char resi_type_resi[n_res][5]);

void calcEnergShift(int n_atom_prot, int n_res, int n_atom_pigm, int n_pig, int *atmsPerPigm,int *atmsPerRes, double **koor_prot, double **koor_pigm,double *q_prot, double *q_gr, double *q_exc, int *typ_at,int *resid, int *res_leng, char amac_type[n_atom_prot][5],char atom_type[n_atom_prot][5]);
