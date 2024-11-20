# ======================================================================
# © 2024. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
# is operated by Triad National Security, LLC for the U.S.
# Department of Energy/National Nuclear Security Administration. All
# rights in the program are reserved by Triad National Security, LLC,
# and the U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others
# acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
# ======================================================================

#--------------------------------------------------------------------------------------------------#
# Regrid (interpolate) Rosseland and Planck opacity data to a log-log rho-T grid.
#--------------------------------------------------------------------------------------------------#
import numpy as np
from scipy.interpolate import interp2d
import argparse

#-- parse command line arguments
parser = argparse.ArgumentParser(description="Re-grid opacity data to be log-log in density-temperature.")
parser.add_argument("fname", type=str, help="Opacity file to re-grid.")
parser.add_argument("--is_loglog", action="store_true", help="Avoid regridding if already log-log.")
parser.add_argument("--fname_new", type=str, default="kappa_new.txt", help="File name for new file created.")
parser.add_argument("--plot", action="store_true", help="If not log-log, plot colormap of original and regridded data.")
args = parser.parse_args()

NRho = -1
NT = -1
with open(args.fname, "r") as ff:
    svec = ff.readline().split(" ")
    NRho = int(svec[2])
    NT = int(svec[4])

print("NRho = ", NRho)
print("NT = ", NT)

opac = np.loadtxt(args.fname, skiprows=1)

print("opac.shape = ", opac.shape)
assert np.size(opac, 0) == NT * NRho, "np.size(opac, 0) != NT * NRho"

#-- density, temperature grids
Rho = np.unique(opac[:,0])
T = np.unique(opac[:,1])

assert np.size(Rho) == NRho, "np.size(Rho) != Rho"
assert np.size(T) == NT, "np.size(T) != NT"

r = np.logspace(np.log10(Rho[0]), np.log10(Rho[NRho - 1]), NRho)
tt = np.logspace(np.log10(T[0]), np.log10(T[NT - 1]), NT)

if (args.is_loglog):
    assert np.max(abs(r - Rho) / Rho) < 1e-6, "np.max(abs(r - Rho) / Rho) >= 1e-6"
    assert np.max(abs(tt - T) / T) < 1e-6, "np.max(abs(tt - T) / T) >= 1e-6"
else:
    print()
    print("Interpolating data to log-log rho-T grid...")
    #-- reshape to rho-T grid
    ross_old_opac = np.reshape(opac[:,2], (NRho, NT))
    plnk_old_opac = np.reshape(opac[:,3], (NRho, NT))
    #-- linearly interpolate in log-log rho-T
    ross_f2d = interp2d(np.log10(T), np.log10(Rho), ross_old_opac, kind='linear')
    plnk_f2d = interp2d(np.log10(T), np.log10(Rho), plnk_old_opac, kind='linear')
    ross_new_opac = ross_f2d(np.log10(tt), np.log10(r))
    plnk_new_opac = plnk_f2d(np.log10(tt), np.log10(r))
    #-- plot 4-panel of data
    if (args.plot):
        import matplotlib.pyplot as plt
        XOLD, YOLD = np.meshgrid(np.log10(T), np.log10(Rho))
        XNEW, YNEW = np.meshgrid(np.log10(tt), np.log10(r))
        fig, axes = plt.subplots(2, 2)
        im1 = axes[0, 0].pcolormesh(XOLD, YOLD, ross_old_opac)
        axes[0, 0].set_title('Old Rosseland')
        fig.colorbar(im1, ax=axes[0, 0])
        im2 = axes[0, 1].pcolormesh(XOLD, YOLD, plnk_old_opac)
        axes[0, 1].set_title('Old Planck')
        fig.colorbar(im2, ax=axes[0, 1])
        im3 = axes[1, 0].pcolormesh(XNEW, YNEW, ross_new_opac)
        axes[1, 0].set_title('New Rosseland')
        fig.colorbar(im3, ax=axes[1, 0])
        im4 = axes[1, 1].pcolormesh(XNEW, YNEW, plnk_new_opac)
        axes[1, 1].set_title('New Planck')
        fig.colorbar(im4, ax=axes[1, 1])
        plt.tight_layout()
        plt.show()
    #-- reset opacity array
    for i in range(NRho):
        for j in range(NT):
            k = j + NT * i
            opac[k,0] = r[i]
            opac[k,1] = tt[j]
            opac[k,2] = ross_new_opac[i,j]
            opac[k,3] = plnk_new_opac[i,j]

#-- save with new header containing min/max rho, T
hdr = "NRho "+str(NRho)+" NT "+str(NT) + "\n" \
    "rho_min {:.8e}".format(r[0])+" rho_max {:.8e}".format(r[NRho-1]) + "\n" \
    "T_min {:.8e}".format(tt[0]) + " T_max {:.8e}".format(tt[NT-1])
np.savetxt(args.fname_new, opac[:,2:], delimiter=" ", header=hdr, comments=" ", fmt="%.8e")
