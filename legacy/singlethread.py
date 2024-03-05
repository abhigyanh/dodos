import numpy as np
import matplotlib.pyplot as plt
import time
import os
from scipy.integrate import simpson
from scipy.optimize import fsolve
from concurrent.futures import ProcessPoolExecutor
from functools import partial

from functions import *

pi = np.pi
c = 0.0299792   #cm/ps

nthreads = 1 #os.cpu_count()//2


# USER DEFINED INPUTS
velocities_fname = "../veloc.xvg"
positions_fname = "../pos.xvg"
temp = 300                      # Kelvin
kT = (8.314e-3)*(temp)          # amu nm^2 ps^-2

# mass array, usage: createMassArray(natoms, atomsPerMolecule, mass1, mass2, ...)
atomsPerMolecule = 3
masslist = [15.999, 1.001, 1.001]






if __name__ == '__main__':
    # clear up the log file
    cleanLogFile()

    # IMPORT VELOCITY FILES AND PRINT OUT SOME BASIC INFORMATION
    veloc = loadxvg(velocities_fname)
    pos = loadxvg(positions_fname)

    dt = veloc[1,0]
    nframes = veloc.shape[0]
    simtime = nframes*dt
    natoms = (veloc.shape[1]-1)//3
    nmols = natoms//atomsPerMolecule

    Log('nframes = {}, simtime = {} ps, dumpfreq = {} ps'.format(nframes, simtime, dt))
    Log('Atoms = {}, Molecules = {}, {} atoms per molecule'.format(natoms, nmols, atomsPerMolecule))
    Log("Velocity array shape = {}".format(veloc.shape))
    Log("Velocity array dtype = {}".format(veloc.dtype))
    Log("Coordinate array shape = {}".format(pos.shape))
    Log("Coordinate array dtype = {}".format(pos.dtype))
    Log("Temperature set to {} K".format(temp))

    # SETUP ARRAYS
    m = createMassArray(natoms, atomsPerMolecule, masslist)
    mtot = sum(masslist)
    
    x = pos[:,1::3]
    y = pos[:,2::3]
    z = pos[:,3::3]

    vx = veloc[:,1::3]
    vy = veloc[:,2::3]
    vz = veloc[:,3::3]

    ##################################################################
    """
    V E L O C I T Y       D E C O M P O S I T I O N
    """
    ##################################################################
    Log("- Beginning velocity decomposition")

    x_rel = np.empty((nframes,natoms))
    y_rel = np.empty((nframes,natoms))
    z_rel = np.empty((nframes,natoms))

    vx_rel = np.empty((nframes,natoms))
    vy_rel = np.empty((nframes,natoms))
    vz_rel = np.empty((nframes,natoms))

    vx_tr = np.empty((nframes,natoms))
    vy_tr = np.empty((nframes,natoms))
    vz_tr = np.empty((nframes,natoms))

    vx_rot = np.empty((nframes,natoms))
    vy_rot = np.empty((nframes,natoms))
    vz_rot = np.empty((nframes,natoms))

    vx_vib = np.empty((nframes,natoms))
    vy_vib = np.empty((nframes,natoms))
    vz_vib = np.empty((nframes,natoms))

    VelocDecompTime = 0

    # def VelocityDecomposition(molID):
    for molID in range(nmols):
        t1 = time.time()
        # Fetch atom IDs for specified molID
        atomIDs = []
        for i in range(atomsPerMolecule):
            atomIDs.append(atomsPerMolecule*molID + i)
        Log("molecule ID = {}, atomIDs = {}".format(molID, atomIDs), console=False)

        for frame in range(nframes):
            vx_com = 0
            vy_com = 0
            vz_com = 0

            x_com = 0
            y_com = 0
            z_com = 0

            Lx = 0
            Ly = 0
            Lz = 0

            wx = 0
            wy = 0
            wz = 0

            # inertia tensor: https://farside.ph.utexas.edu/teaching/336k/Newton/node64.html
            Ixx = 0
            Iyy = 0
            Izz = 0
            Ixy = 0
            Iyz = 0
            Izx = 0 
            
            # Get center of mass variables
            for atomID in atomIDs:
                vx_com += m[atomID]/mtot * vx[frame,atomID]
                vy_com += m[atomID]/mtot * vy[frame,atomID]
                vz_com += m[atomID]/mtot * vz[frame,atomID]

                x_com += m[atomID]/mtot * x[frame,atomID]
                y_com += m[atomID]/mtot * y[frame,atomID]
                z_com += m[atomID]/mtot * z[frame,atomID]


            # Set relative positions and velocities
            for atomID in atomIDs:
                # Relative position w.r.t. COM
                x_rel[frame,atomID] = x[frame,atomID] - x_com
                y_rel[frame,atomID] = y[frame,atomID] - y_com
                z_rel[frame,atomID] = z[frame,atomID] - z_com
                # Relative velocity w.r.t COM
                vx_rel[frame,atomID] = vx[frame,atomID] - vx_com
                vy_rel[frame,atomID] = vy[frame,atomID] - vy_com
                vz_rel[frame,atomID] = vz[frame,atomID] - vz_com
                # Set translational velocity as COM velocity
                vx_tr[frame,atomID] = vx_com
                vy_tr[frame,atomID] = vy_com
                vz_tr[frame,atomID] = vz_com

            # CALCULATE L AND I
            for atomID in atomIDs:
                Lx += m[atomID]*(y_rel[frame,atomID]*vz_rel[frame,atomID] - z_rel[frame,atomID]*vy_rel[frame,atomID])
                Ly += m[atomID]*(z_rel[frame,atomID]*vx_rel[frame,atomID] - x_rel[frame,atomID]*vz_rel[frame,atomID])
                Lz += m[atomID]*(x_rel[frame,atomID]*vy_rel[frame,atomID] - y_rel[frame,atomID]*vx_rel[frame,atomID])

                Ixx += m[atomID]*(y_rel[frame,atomID]**2 + z_rel[frame,atomID]**2)
                Iyy += m[atomID]*(x_rel[frame,atomID]**2 + z_rel[frame,atomID]**2)
                Izz += m[atomID]*(x_rel[frame,atomID]**2 + y_rel[frame,atomID]**2)

                Ixy -= m[atomID]*(x_rel[frame,atomID]*y_rel[frame,atomID])
                Iyz -= m[atomID]*(y_rel[frame,atomID]*z_rel[frame,atomID])
                Izx -= m[atomID]*(z_rel[frame,atomID]*x_rel[frame,atomID])

            # NOW CALCULATE w (omega) BY SOLVING Iw = L
            I = np.array([[Ixx,Ixy,Izx], 
                          [Ixy,Iyy,Iyz], 
                          [Izx,Iyz,Izz]])
            L = np.array([Lx, Ly, Lz])
            eigenvalues, eigenvectors = np.linalg.eig(I)
            I_1 = eigenvalues[0]
            I_2 = eigenvalues[1]
            I_3 = eigenvalues[2]
            p_1 = eigenvectors[0]
            p_2 = eigenvectors[1]
            p_3 = eigenvectors[2]

            omega_1 = (p_1/I_1)*(L[0]*p_1[0] + L[1]*p_1[1] + L[2]*p_1[2])
            omega_2 = (p_2/I_2)*(L[0]*p_2[0] + L[1]*p_2[1] + L[2]*p_2[2])
            omega_3 = (p_3/I_3)*(L[0]*p_3[0] + L[1]*p_3[1] + L[2]*p_3[2])
            

            omega = omega_1 + omega_2 + omega_3
            wx = omega[0]
            wy = omega[1]
            wz = omega[2]

            # HENCE CALCULATE THE ROTATIONAL AND VIBRATIONAL
            for atomID in atomIDs:
                vx_rot[frame,atomID] = wy*z_rel[frame,atomID] - wz*y_rel[frame,atomID]
                vy_rot[frame,atomID] = wz*x_rel[frame,atomID] - wx*z_rel[frame,atomID]
                vz_rot[frame,atomID] = wx*y_rel[frame,atomID] - wy*x_rel[frame,atomID]

                vx_vib[frame,atomID] = vx[frame,atomID] - vx_com - vx_rot[frame,atomID]
                vy_vib[frame,atomID] = vy[frame,atomID] - vy_com - vy_rot[frame,atomID]
                vz_vib[frame,atomID] = vz[frame,atomID] - vz_com - vz_rot[frame,atomID]

                # if np.sqrt(vx_vib[frame,atomID]**2 + vy_vib[frame,atomID]**2 + vz_vib[frame,atomID]**2) > 0.1:
                    # Log('Vibrational components of atomic velocity quite large, is this correct?', console=False)

                # Log(str([frame, atomID, vx_tr[frame,atomID], vy_rot[frame,atomID], vz_vib[frame,atomID]]), console=False)

        Log("- Velocity decomposition done for molecule ID = {}, atomIDs = {}".format(molID, atomIDs), console=False)
        t2 = time.time()
        if molID % 20 == 0:
            print('Approx. time left = {:.2f} s'.format((t2-t1)*(nmols-molID)/nthreads))
        VelocDecompTime += t2-t1
        # return t2-t1
    '''
    # Do velocity decomposition by parallellizing over molecules
    executor = ProcessPoolExecutor(max_workers=nthreads)
    for t in executor.map(VelocityDecomposition, range(nmols)):
        VelocDecompTime += t/nthreads
    '''
    Log("- Finished doing velocity decomposition (took {:.1f} s, approx. {:.1f} ms per molecule)".format(VelocDecompTime, 1000*VelocDecompTime/nmols))

    # Save the decomposed velocities if you want lol
    # np.savetxt('decomposed-velocities.dat', np.column_stack((np.arange(nframes)*dt,vx_tr, vy_tr, vz_tr, vx_rot, vy_rot, vz_rot, vx_vib, vy_vib, vz_vib)))
    ##################################################################
    """
    E N D    O F     V E L O C I T Y       D E C O M P O S I T I O N
    """
    ##################################################################
    """
    F F T   Z O N E
    """
    ##################################################################
    Log("Beginning DoS calculation(s) for translational, rotational and vibrational components...")

    nu = np.fft.rfftfreq(n = nframes, d = dt)
    dos_tr = np.zeros(nframes//2+1)
    dos_rot = np.zeros(nframes//2+1)
    dos_vib = np.zeros(nframes//2+1)
    dos_tot = np.zeros(nframes//2+1)
    fft_normalization = dt/nframes

    avgwalltime = 0

    for atomID in range(natoms):
        t1 = time.time()

        sx_tr  = powerspec(vx_tr[:,atomID])*fft_normalization
        sx_rot = powerspec(vx_rot[:,atomID])*fft_normalization
        sx_vib = powerspec(vx_vib[:,atomID])*fft_normalization
        sx_tot = powerspec(vx[:,atomID])*fft_normalization

        sy_tr  = powerspec(vy_tr[:,atomID])*fft_normalization
        sy_rot = powerspec(vy_rot[:,atomID])*fft_normalization
        sy_vib = powerspec(vy_vib[:,atomID])*fft_normalization
        sy_tot = powerspec(vy[:,atomID])*fft_normalization

        sz_tr  = powerspec(vz_tr[:,atomID])*fft_normalization
        sz_rot = powerspec(vz_rot[:,atomID])*fft_normalization
        sz_vib = powerspec(vz_vib[:,atomID])*fft_normalization
        sz_tot = powerspec(vz[:,atomID])*fft_normalization

        dos_tr  += m[atomID]*(sx_tr + sy_tr + sz_tr)
        dos_rot += m[atomID]*(sx_rot + sy_rot + sz_rot)
        dos_vib += m[atomID]*(sx_vib + sy_vib + sz_vib)
        dos_tot += m[atomID]*(sx_tot + sy_tot + sz_tot)
        t2 = time.time()
        avgwalltime += 1/natoms * (t2-t1)
    
    # Multiply with dimensional constant
    dos_tr *= 2/kT
    dos_rot *= 2/kT
    dos_vib *= 2/kT
    dos_tot *= 2/kT


    # Integrate all DoS to get number of states
    integral_dos_tr = simpson(dos_tr, nu)
    integral_dos_rot = simpson(dos_rot, nu)
    integral_dos_vib = simpson(dos_vib, nu)
    integral_dos_tot = simpson(dos_tot, nu)


    Log("FFT time per atom: {:.1f} ms".format(avgwalltime*1000))
    Log("FFT time (total): {:.1f} s".format(avgwalltime*natoms))
    Log("Integral of DoS_tr = {}".format(integral_dos_tr))
    Log("Integral of DoS_rot = {}".format(integral_dos_rot))
    Log("Integral of DoS_vib = {}".format(integral_dos_vib))
    Log("Integral of DoS_total = {}".format(integral_dos_tot))

    plt.plot(nu/c, c*dos_tr,  label = 'tr')
    plt.plot(nu/c, c*dos_rot, label = 'rot')
    plt.plot(nu/c, c*dos_vib, label = 'vib')

    plt.gca().set(
        xlim = [0,4000],
        title = 'Density of States',
        xlabel = r'Wavenumber $(cm^{-1})$',
        ylabel = r'$S(\nu)$ (cm)',
    )
    plt.legend()
    plt.savefig('DoS_cm.png', dpi = 200, bbox_inches = 'tight')
    plt.clf()

    plt.plot(nu/c, c*dos_tot, label = 'Total DoS')

    plt.gca().set(
        xlim = [0,4000],
        title = 'Density of States',
        xlabel = r'Wavenumber $(cm^{-1})$',
        ylabel = r'$S(\nu)$ (cm)',
    )
    plt.legend()
    plt.savefig('DoS_total.png', dpi = 200, bbox_inches = 'tight')
    plt.clf()

