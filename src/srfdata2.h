struct srfdata {
   int nx, ny;
   int nt;
   float dt;
   char *fname;
   int itype, wstep, ntiskp;
   float **V;
};

struct srfdata srfdata_init(int nx, int ny, int nt, float dt, int itype, int wstep, 
        int ntiskp, char *fname);

void srfdata_loadts(struct srfdata *srf, int t);

void srfdata_getts(struct srfdata *srf, float *D, int t);

float *srfdata_gettrace(struct srfdata *srf, int k, int l);

void srfdata_getrow(struct srfdata *srf, int l, float **D);
void srfdata_getrow_wstep(struct srfdata *srf, int l, float **D);
void srfdata_getrows_wstep(struct srfdata *srf, int l, int nrows, float ***D);

void srfdata_getchunk(struct srfdata *srf, int l1, int l2, float ***D);
