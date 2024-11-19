#--------------------------------------------------------------------------------------------------#
# This script defaults to using the temperature and density bounds for opacity data found
# in Zhaohuan Zhu et al (2021), "Global 3D Radiation Hydrodynamic Simulations of Proto-Jupiter’s
# Convective EnvelopeConvective Envelope". Otherwise, it merely generates a power law for Rosseland
# and Planck-averaged opacities, and produces a table in an ASCII format equivalent to that of Zhu.
#--------------------------------------------------------------------------------------------------#
import numpy as np
import argparse

#-- parse command line arguments
parser = argparse.ArgumentParser(description="Create power-law grey opacity ASCII table (for testing).")
parser.add_argument("fname", type=str, help="base opacity file name to create (.txt is appended).")
parser.add_argument("--Nrho", type=int, default=4, help="Number of density points (assumes log space)")
parser.add_argument("--rho", type=float, nargs=2, default=[1e-14, 0.7943282347241912],
                    help="min, max densities [g/cc]")
parser.add_argument("--NT", type=int, default=4, help="Number of temperature points (assumes log space)")
parser.add_argument("--T", type=float, nargs=2, default=[1.0, 7943282.347242886],
                    help="min, max temperatures [K]")
parser.add_argument("--rho_ref", type=float, default=1e-6, help="Reference density [g/cc]")
parser.add_argument("--rho_exp", type=float, default=0.0, help="Exponent of density / ref. density")
parser.add_argument("--T_ref", type=float, default=11604, help="Reference density [K]")
parser.add_argument("--T_exp", type=float, default=0.0, help="Exponent of temperature / ref. temperature")
parser.add_argument("--Pkap_coef", type=float, default=0.01, help="Planck opacity coefficient [cm^2/g]")
parser.add_argument("--Rkap_coef", type=float, default=0.01, help="Rosseland opacity coefficient [cm^2/g]")
args = parser.parse_args()

#-- implement simple power law
def power_law_kappa(r, t):
    '''power_law_kappa (cgs): r=density (1st arg), t=temperature (2nd arg)'''
    kap_coef = [args.Rkap_coef, args.Pkap_coef]
    kap_pwrs = (r / args.rho_ref)**(args.rho_exp) * (t / args.T_ref)**(args.T_exp)
    return [kapc * kap_pwrs for kapc in kap_coef]

#-- set rho and T arrays
rho = np.logspace(np.log10(args.rho[0]), np.log10(args.rho[1]), args.Nrho)
T = np.logspace(np.log10(args.T[0]), np.log10(args.T[1]), args.NT)

#-- initialize the opacity array - columns: density, temperature, Rosseland, Planck
Nrow = args.Nrho * args.NT
opac = np.zeros((Nrow, 4))

#-- populate opacity array
for i in range(args.Nrho):
    for j in range(args.NT):
        k = j + args.NT * i
        opac[k,0] = rho[i]
        opac[k,1] = T[j]
        opac[k,2:] = power_law_kappa(rho[i], T[j])

#-- save opacity array
np.savetxt(args.fname+'.txt', opac, delimiter=" ",
           header="rho "+str(args.Nrho)+" T "+str(args.NT), comments=" ", fmt="%.8e")

#--------------------------------------------------------------------------------------------------#
# End of create_mean_opac_ascii.py
#--------------------------------------------------------------------------------------------------#
