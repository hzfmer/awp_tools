#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
'''Add heterogeneities to mesh in the top block
Parameters
    std:  std of Von Karman
    nx, ny, nz:  dimensions of the output file
    mz:          number of layers of the ssh file
    fname_fr:    fractal ssh file
    fname_hom:   original "homogeneous" model
    fname_het:   output heterogeneous model
    nvar:        number of variables in the mesh (3/5/8)
    vsmin:       given from the homo mesh
    vpmax:       given by the curvilinear CFL
    vsmax:       vpmax / sqrt(2)
'''

import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# Modify these parameters when needed
std = 0.05
nx, ny, nz, mz = 3456, 3456, 108, 108
nvar = 3
vsmin = 100
fname_fr = "ssh_8m_h005l100.out"
fname_hom = "../la_habra_small_8m_gpu_dm_abc50/mesh_orig_p0.bin"
fname_het = "../la_habra_small_8m_gpu_dm_abc50/mesh_het_p0.bin"
sqrt2 = np.sqrt(2)

vpmax = np.genfromtxt('vpmax_depth.txt')
vsmax = vpmax / sqrt2
length_fr = nx * ny
length_hom = nx * ny * nvar 
length_het = nx * ny * nvar 

if rank >= nz:
    offset_fr = 4 * (nz - 1) * length_fr
else:
    offset_fr = 4 * rank * length_fr
offset_hom = 4 * rank * length_hom
offset_het = 4 * rank * length_het

for i in range(rank, nz, size):
    buf_fr = np.fromfile(fname_fr, dtype='float32', count=length_fr, offset=offset_fr)
    buf_hom = np.fromfile(fname_hom, dtype='float32', count=length_hom, offset=offset_hom).reshape(-1, nvar)
    buf_het = np.copy(buf_hom)
    if i < mz:
        vp = buf_hom[:, 0]
        vs = buf_hom[:, 1]
        rho = buf_hom[:, 2]
        vp = np.divide(1, vp, out=np.zeros_like(vp), where=vs>0)
        vs = np.divide(1, vs, out=np.zeros_like(vs), where=vs>0)
        vp = vp * (1 - std * buf_fr)
        vp[vp > vpmax[i]] = vpmax[i]
        vs = vs * (1 - std * buf_fr)
        vs[vs > vsmax[i]] = vsmax[i]
        vp[vp < vsmin * sqrt2] = vsmin * sqrt2
        buf_het[:, 0] = np.divide(1, vp, out=np.zeros_like(vp), where=vp>0)
        buf_het[:, 1] = np.divide(1, vs, out=np.zeros_like(vs), where=vs>0)
        buf_het[:, 2] = rho * (1 + std * buf_fr)


amode = MPI_MODE_WRONLY|MPI_MODE_CREATE
fh = MPI.File.Open(comm, fname_het, amode)
fh.Write_at_all(offset_het, buf_het)
