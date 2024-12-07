import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root_scalar
from scipy.integrate import simpson

from constants import *
from functions import *

def Delta(T,V,N,m,s0):
    """
    Inputs:
    T = Kelvin
    V = nm^3
    N = nmols
    m = amu
    s0 = ps
    """
    delta = (2*s0)/(9*N) * (pi*kB*T/m)**(1/2) * (N/V)**(1/3) * (6/pi)**(2/3)
    # Delta should be dimensionless, but it's calculated as (ps * sqrt(eV / amu) / nm)
    # define the conversion factor:
    conv2 = ps * np.sqrt(eV / amu) / nm
    return delta*conv2

def equation_of_fluidicity(f, delta):
    if f < 0 or delta < 0:
        return 69420    # negative values of f or delta are not possible, return huge number to reset the root finger
    return 2*delta**(-9/2)*f**(15/2) - 6*delta**(-3)*f**(5) - delta**(-3/2)*f**(7/2) + 6*delta**(-3/2)*f**(5/2) + 2*f - 2

def calculate_fluidicity(delta):
    return root_scalar(lambda f: equation_of_fluidicity(f, delta), x0 = 0.01, x1 = 0.2, method='secant').root




    
    ##################################################################
    """
    FLUIDICITY CALCULATION AND DENSITY OF STATES DECOMPOSITION
    """
    ##################################################################

def TwoPhaseDecompose_Translational(nu, DOS_tr, T, V, N, m):
    s0_tr = DOS_tr[0]
    Delta_tr = Delta(T,V,N,m,s0_tr)
    f_tr = calculate_fluidicity(Delta_tr)

    DOS_tr_g = s0_tr/(1 + ((pi*s0_tr*nu)/(6*N*f_tr))**2)
    DOS_tr_s = DOS_tr - DOS_tr_g

    plt.plot(nu/c, c*DOS_tr_g,  label = 'tr_g')
    plt.plot(nu/c, c*DOS_tr_s,  label = 'tr_s')
    plt.plot(nu/c, c*DOS_tr,    label = 'tr')

    plt.gca().set(
        xlim = [0,1500],
        xlabel = r'Wavenumber $(cm^{-1})$',
        ylabel = r'$S(\nu)$ (cm)',
    )
    plt.legend()
    plt.savefig('DoS_tr.png', dpi = 200, bbox_inches = 'tight')
    plt.clf()

    return Delta_tr, f_tr, s0_tr, DOS_tr_g, DOS_tr_s


def TwoPhaseDecompose_Rotational(nu, DOS_rot, T, V, N, m):
    s0_rot = DOS_rot[0]
    Delta_rot = Delta(T,V,N,m,s0_rot)
    f_rot = calculate_fluidicity(Delta_rot)

    DOS_rot_g = s0_rot/(1 + ((pi*s0_rot*nu)/(6*N*f_rot))**2)
    DOS_rot_s = DOS_rot - DOS_rot_g

    # Check condition for solid DoS being negative because S0_rot[nu=0] is discontinuous due to bugs
    if np.any(DOS_rot_s < 0) == True:
        Log("""
        ! Rotational DOS (solid) appears to contain negative values. This might (but not necessarily) be a bug where S0_rot has an unusually large zero-frequency discontinuity.
        ! Interpolating s_0 to approximate correct value. This does not affect entropy calculations if there is no discontinuity.
        """, console=False)
        # Linear interpolate the correct value
        s0_rot_new = DOS_rot[1] + (DOS_rot[1]-DOS_rot[2])/(nu[1] - nu[2])*(nu[1]-nu[0])
        DOS_rot[0] = s0_rot_new

        Log("New value of s0_rot = {:.5f}".format(s0_rot), console=False)
    
        Delta_rot = Delta(T,V,N,m,s0_rot)
        f_rot = calculate_fluidicity(Delta_rot)
        DOS_rot_g = s0_rot/(1 + ((pi*s0_rot*nu)/(6*N*f_rot))**2)
        DOS_rot_s = DOS_rot - DOS_rot_g


    plt.plot(nu/c, c*DOS_rot_g,  label = 'rot_g')
    plt.plot(nu/c, c*DOS_rot_s,  label = 'rot_s')
    plt.plot(nu/c, c*DOS_rot,    label = 'rot')

    plt.gca().set(
        xlim = [0,1500],
        xlabel = r'Wavenumber $(cm^{-1})$',
        ylabel = r'$S(\nu)$ (cm)',
    )
    plt.legend()
    plt.savefig('DoS_rot.png', dpi = 200, bbox_inches = 'tight')
    plt.clf()

    return Delta_rot, f_rot, s0_rot, DOS_rot_g, DOS_rot_s


def TwoPhaseDecompose_Vibrational(nu, DOS_vib, T, V, N, m):
    # Fluidity 0 in case of vibrations, entire spectrum is solid
    s0_vib = DOS_vib[0]
    Delta_vib = 0
    f_vib = 0

    DOS_vib_g = np.zeros(nu.shape)
    DOS_vib_s = DOS_vib - DOS_vib_g

    plt.plot(nu/c, c*DOS_vib_g,  label = 'vib_g')
    plt.plot(nu/c, c*DOS_vib_s,  label = 'vib_s')
    plt.plot(nu/c, c*DOS_vib,    label = 'vib')

    plt.gca().set(
        xlim = [0,4000],
        xlabel = r'Wavenumber $(cm^{-1})$',
        ylabel = r'$S(\nu)$ (cm)',
    )
    plt.legend()
    plt.savefig('DoS_vib.png', dpi = 200, bbox_inches = 'tight')
    plt.clf()
    return Delta_vib, f_vib, s0_vib, DOS_vib_g, DOS_vib_s


    ##################################################################
    """
    ENTROPY CALCULATION
    """
    ##################################################################


def CalculateEntropy_Translational(nu, DOS_tr_s, DOS_tr_g, f_tr, Delta_tr, m, N, V, T):
    nu = nu.copy()
    kT = kB*T

    # Hard sphere packing fraction; Carnahan-Sterling hard sphere equation of state
    y = (f_tr**(5/2)) / (Delta_tr**(3/2))
    z = (1 + y + y**2 - y**3) / (1 - y)**3
    if z <= 0:
        Log("?! Hard sphere compressibility negative or zero! This should not happen. y, z = {}, {}".format(y,z))

    # Dimensionless hard sphere entropy
    # The quantity inside the SHS logarithm will be calculated in (amu / eV)^3/2 * nm^3 / ps^3 but it should be dimensionless.
    # define this conversion factor:
    conv3 = (amu / eV)**(3/2) * nm**3 / ps**3
    S_HS = 5/2 + np.log(conv3 * (2*pi*m*kB*T/h**2)**(3/2) * V/N * z/f_tr ) + (3*y*y - 4*y)/(1-y)**2

    # Weighing function
    nu = nu[np.nonzero(nu)]
    bhn = h/kT * nu
    W_gas = 1/3 * S_HS
    W_solid = bhn/(np.exp(bhn)-1) - np.log(1-np.exp(-bhn))
    DOS_tr_g = DOS_tr_g[np.nonzero(nu)]
    DOS_tr_s = DOS_tr_s[np.nonzero(nu)]

    # DEBUGGING
    # Sometimes S_translational turns out to be negative. Why is this happening? Which one of these things is negative: S_HS, W_solid
    Log("S_HS = {}".format(S_HS), console=False)
    plt.plot(nu, W_solid)
    plt.savefig('W_solid.png', dpi = 200, bbox_inches = 'tight')
    plt.clf()

    S_tr = simpson(DOS_tr_g*W_gas, x=nu) + simpson(DOS_tr_s*W_solid, x=nu)
    S_tr *= kB*eVtoJ
    return S_tr


def CalculateEntropy_Rotational(nu, DOS_rot_s, DOS_rot_g, I_1, I_2, I_3, Sigma, T):
    nu = nu.copy()
    kT = kB*T

    # Rotational entropy
    conv4 = (eV * ps**2) / (amu * nm**2)
    Theta1 = h**2/(8*pi**2*kB*I_1) * conv4
    Theta2 = h**2/(8*pi**2*kB*I_2) * conv4
    Theta3 = h**2/(8*pi**2*kB*I_3) * conv4
    S_R =  np.log( pi**(1/2)*e**(3/2)/Sigma * (T**3/(Theta1*Theta2*Theta3))**(1/2))

    # Weighing function
    nu = nu[np.nonzero(nu)]
    bhn = h/kT * nu
    W_gas = 1/3 * S_R
    W_solid = bhn/(np.exp(bhn)-1) - np.log(1-np.exp(-bhn))
    DOS_rot_g = DOS_rot_g[np.nonzero(nu)]
    DOS_rot_s = DOS_rot_s[np.nonzero(nu)]

    S_rot = simpson(DOS_rot_g*W_gas, x=nu) + simpson(DOS_rot_s*W_solid, x=nu)
    S_rot *= kB*eVtoJ
    return S_rot


def CalculateEntropy_Vibrational(nu, DOS_vib, T):
    nu = nu.copy()
    kT = kB*T

    # Weighing function
    nu = nu[np.nonzero(nu)]
    bhn = h/kT * nu
    W_solid = bhn/(np.exp(bhn)-1) - np.log(1-np.exp(-bhn))
    DOS_vib = DOS_vib[np.nonzero(nu)]

    S_vib = simpson(DOS_vib*W_solid, x=nu)
    S_vib *= kB*eVtoJ
    return S_vib