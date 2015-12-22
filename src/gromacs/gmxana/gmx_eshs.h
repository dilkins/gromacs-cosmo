typedef struct {
  float beta;
  float num;
} t_Corr;

typedef struct {
  char  MapFName[256];
  char  WaterName[256];
  char  CationName[256];
  char  AnionName[256];
  char  Polarization[256];
  char  Method[256];
  char  GPU[256];
  char  Fade[256];
  float cutoff; // cutoff for electic field calculation
  int   nbinq;  // number of bins in q
} t_Input;

typedef struct {
  float beta_gas[27];
  float gamma_gas[81];
  float D[27][212];
} t_Map;

typedef struct {
  int   atom[4];
  float q[4];
  float dip_x[3];   // only for H
  float dip_z[3];
  float qdp_xx[3];
  float qdp_yy[3];
  float qdp_zz[3];
  float qdp_xz[3];  // only for H
  float beta_mol[3][3][3];
  float beta_lab[3][3][3];
  float mu_ind;
} t_Water;

typedef struct {
  int    atom[1];
  float  q[1];
  float  beta[3];
} t_Ion;


void readInput(char const inputFName[256],t_Input *Input);
void readMap(t_Input *Input,t_Map *Map);
int keyWordS(char const *keyWord,char *Buffer,char  *ivalue,size_t LabelLength);
int keyWordF(char const *keyWord,char *Buffer,float *ivalue,size_t LabelLength);
int keyWordI(char const *keyWord,char *Buffer,int   *ivalue,size_t LabelLength);

int check_ion(t_topology *top,char name[256]);
int check_water(t_topology *top);
void identifyIon(t_topology *top,t_Ion *Ion,char name[256]);
void identifyTIP4P(t_topology *top,t_Water *Water);
void identifyAMOEBA(t_topology *top,t_Water *Water);

void calc_electronic_correlation(t_pbc *pbc,rvec *x,t_Water *Water,int nwaters,float R2,float inv_binwidth,int nbins,t_Corr *Correlation);

void calc_eshs(const char *fnTRX,t_topology *top,output_env_t oenv,char const input[256]);
void calc_water_beta(t_pbc *pbc,t_topology *top,t_Input *Input,t_Map *Map,rvec *x,t_Water *Water,t_Ion *Cation,t_Ion *Anion,int nwaters,int ncations,int nanions,float cutoff2);

void calc_ind_dipole_water(t_Input *Input,t_Water *Water);
void calc_beta_map(t_Map *Map,float E[9],float beta_liquid[27]);

void calc_efield_amoeba(t_pbc *pbc,t_topology *top,rvec *x,int i,rvec x_mol_i,rvec y_mol_i,rvec z_mol_i,t_Water *Water,t_Ion *Cation,t_Ion *Anion,int nwaters,int ncations,int nanions,float cutoff2,float E[27]);

void calc_efield_tip4p(t_pbc *pbc,t_topology *top,rvec *x,int i,rvec x_mol,rvec y_mol,rvec z_mol,t_Water *Water,t_Ion *Cation,t_Ion *Anion,int nwaters,int ncations,int nanions,float cutoff2,float E[27]);
