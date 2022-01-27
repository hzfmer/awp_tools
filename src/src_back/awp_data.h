void read_awp_timeseries(char *fname, int nx, int ny, int nz, int xi, int yi, int zi, 
    int nt, int nsites, float *buf);

void read_awp_timeseries_multi(char *fbase, int nx, int ny, int nz, int xi, int yi, int zi, 
    int nt, int wstep, int nskip, int nsites, float *buf);
