typedef struct {
   int nbgx, nbgy, nbgz, nedx, nedy, nedz, nskpx, nskpy, nskpz, ntiskp; 
   float dt, tmax;
   char sxrgo[100], syrgo[100], szrgo[100];
   int itype;
   int readstep, writestep;
} FD3D_param;

void read_settings(FD3D_param *P, char*in3dfile);
