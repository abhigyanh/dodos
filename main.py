import numpy as np
from scipy.integrate import simpson
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial

from constants import *
from functions import *
from entropy import *




if __name__ == '__main__':
    positions_fname, velocities_fname, box_fname, temp, masslist, nthreads = userInputHandler(sys.argv[1:])
    # clear up the log file
    cleanLogFile()

    # IMPORT VELOCITY FILES AND PRINT OUT SOME BASIC INFORMATION
    veloc = loadxvg(velocities_fname)
    pos = loadxvg(positions_fname)
    box = loadxvg(box_fname)
    m = masslist
    mtot = sum(masslist)
    atomsPerMolecule = len(m)
    kT = (kB)*(temp)                    # eV

    dt = veloc[1,0]
    nframes = veloc.shape[0]
    simtime = nframes*dt
    natoms = (veloc.shape[1]-1)//3
    nmols = natoms//atomsPerMolecule

    boxX = np.mean(box[:,1])    
    boxY = np.mean(box[:,2])
    boxZ = np.mean(box[:,3])
    Volume = boxX * boxY * boxZ         # nm^3

    """Special code block to recalculate the volume if I'm reducing the trajectory to liquid-only"""
    SUBTRACT_ICE_VOL = True
    if SUBTRACT_ICE_VOL == True:
        # Subtract the volume of the ice block. It has same x and y dimensions, but z is custom defined by me ofc.
        Log("## Correcting volume by subtracting estimated volume of ice block! ##")
        Volume -= boxX*boxY*(2.8)


    sysDensity = nmols*mtot/Volume * amu/(nm**3)


    Log("Velocity array shape = {}".format(veloc.shape), console=False)
    Log("Velocity array dtype = {}".format(veloc.dtype), console=False)
    Log("Coordinate array shape = {}".format(pos.shape), console=False)
    Log("Coordinate array dtype = {}".format(pos.dtype), console=False)
    Log("Box array shape = {}".format(box.shape), console=False)
    Log("Box array dtype = {}".format(box.dtype), console=False)

    Log("""
    # BEGINNING DODOS #
    Trajectory and system information:
    - Number of frames = {}          Simulation time = {} ps          Dump frequency = {} ps
    - Number of atoms = {}           Number of molecules = {}         Atoms per molecule = {}
    - Atomic masses = {}             
    - Box volume = {:06.3f} nm^3         System density = {:.1f} kg/m^3
    - Temperature set to {} K        
    - Number of CPU threads used = {}
    """.format(
        nframes, simtime, dt, 
        natoms, nmols, atomsPerMolecule, 
        masslist, 
        Volume, sysDensity,
        temp, 
        nthreads)
    )

    X = pos[:,1::3]
    Y = pos[:,2::3]
    Z = pos[:,3::3]

    Vx = veloc[:,1::3]
    Vy = veloc[:,2::3]
    Vz = veloc[:,3::3]

    ##################################################################
    """
    V E L O C I T Y       D E C O M P O S I T I O N
    """
    ##################################################################
    Log("> Beginning velocity decomposition and calculation of velocity spectra")

    VelocDecompTime = 0

    Inertia1 = 0
    Inertia2 = 0
    Inertia3 = 0

    def VelocityDecomposition(molID):
        t1 = time.time()
        # Fetch atom IDs for specified molID
        atomIDs = []
        fullatomIDs = []
        for i in range(atomsPerMolecule):
            atomIDs.append(i)
            fullatomIDs.append(atomsPerMolecule*molID + i)
        Log("- molecule ID = {}, fullatomIDs = {}".format(molID, fullatomIDs), console=False)

        """
            Find and splice local position and velocity arrays from the full data
        """
        x = np.empty((nframes,atomsPerMolecule))
        y = np.empty((nframes,atomsPerMolecule))
        z = np.empty((nframes,atomsPerMolecule))

        vx = np.empty((nframes,atomsPerMolecule))
        vy = np.empty((nframes,atomsPerMolecule))
        vz = np.empty((nframes,atomsPerMolecule))

        for i,fatomID in enumerate(fullatomIDs):
            x[:,i] = X[:,fatomID].copy()
            vx[:,i] = Vx[:,fatomID].copy()
            y[:,i] = Y[:,fatomID].copy()
            vy[:,i] = Vy[:,fatomID].copy()
            z[:,i] = Z[:,fatomID].copy()
            vz[:,i] = Vz[:,fatomID].copy()

        """
            Define decomposed velocity arrays inside this function and later print them if needed
        """
        x_rel = np.empty((nframes,atomsPerMolecule))
        y_rel = np.empty((nframes,atomsPerMolecule))
        z_rel = np.empty((nframes,atomsPerMolecule))

        vx_rel = np.empty((nframes,atomsPerMolecule))
        vy_rel = np.empty((nframes,atomsPerMolecule))
        vz_rel = np.empty((nframes,atomsPerMolecule))

        vx_tr = np.empty((nframes,atomsPerMolecule))
        vy_tr = np.empty((nframes,atomsPerMolecule))
        vz_tr = np.empty((nframes,atomsPerMolecule))

        vx_rot = np.empty((nframes,atomsPerMolecule))
        vy_rot = np.empty((nframes,atomsPerMolecule))
        vz_rot = np.empty((nframes,atomsPerMolecule))

        vx_vib = np.empty((nframes,atomsPerMolecule))
        vy_vib = np.empty((nframes,atomsPerMolecule))
        vz_vib = np.empty((nframes,atomsPerMolecule))

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

            inertias = [I_1, I_2, I_3]
            inertias.sort()

            # HENCE CALCULATE THE ROTATIONAL AND VIBRATIONAL
            for atomID in atomIDs:
                vx_rot[frame,atomID] = wy*z_rel[frame,atomID] - wz*y_rel[frame,atomID]
                vy_rot[frame,atomID] = wz*x_rel[frame,atomID] - wx*z_rel[frame,atomID]
                vz_rot[frame,atomID] = wx*y_rel[frame,atomID] - wy*x_rel[frame,atomID]

                vx_vib[frame,atomID] = vx[frame,atomID] - vx_com - vx_rot[frame,atomID]
                vy_vib[frame,atomID] = vy[frame,atomID] - vy_com - vy_rot[frame,atomID]
                vz_vib[frame,atomID] = vz[frame,atomID] - vz_com - vz_rot[frame,atomID]



        Log("> Velocity decomposition done for molecule ID = {}".format(molID), console=False)

        # Now get the velocity spectrum and add it up 
        fft_normalization = dt/nframes
        dos_tr = np.zeros(nframes//2+1)
        dos_rot = np.zeros(nframes//2+1)
        dos_vib = np.zeros(nframes//2+1)
        dos_tot = np.zeros(nframes//2+1)
        
        for atomID in atomIDs:
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
    
        # Multiply with dimensional constant
        conv1 = nm**2 * amu / eV / ps**2
        dos_tr *= 2/kT * conv1
        dos_rot *= 2/kT * conv1
        dos_vib *= 2/kT * conv1
        dos_tot *= 2/kT * conv1

        t2 = time.time()
        if molID % (nmols//20) == 0:
            print('Est. time left = {:.0f} s'.format((t2-t1)*(nmols-molID)/nthreads))
        return [dos_tr, dos_rot, dos_vib, dos_tot, inertias[0], inertias[1], inertias[2]]

    """
        Do velocity decomposition by parallellizing over molecules
    """

    DOS_tr = np.zeros(nframes//2+1)
    DOS_rot = np.zeros(nframes//2+1)
    DOS_vib = np.zeros(nframes//2+1)
    DOS_tot = np.zeros(nframes//2+1)

    VelocDecompTime = time.time()
    executor = ProcessPoolExecutor(max_workers=nthreads)
    for result in executor.map(VelocityDecomposition, range(nmols)):
        DOS_tr += result[0]
        DOS_rot += result[1]
        DOS_vib += result[2]
        DOS_tot += result[3]
        Inertia1 += result[4]/nmols
        Inertia2 += result[5]/nmols
        Inertia3 += result[6]/nmols

    VelocDecompTime = -(VelocDecompTime - time.time())
    Log("> Finished velocity decomposition + spectrum (took {:.1f} s | {:.1f} ms per atom)".format(VelocDecompTime, 1000*VelocDecompTime/natoms))
    Log("- Principal axes of inertia = {}, {}, {} amu nm^2".format(Inertia1, Inertia2, Inertia3), console=False)

    ##################################################################
    """
    E N D    O F     V E L O C I T Y       D E C O M P O S I T I O N
    """
    ##################################################################



    Log("> Finishing up DoS calculation(s) for translational, rotational and vibrational components...")

    # Integrate all DoS to get number of states
    nu = np.fft.rfftfreq(n = nframes, d = dt)
    integral_dos_tr = simpson(DOS_tr, nu)
    integral_dos_rot = simpson(DOS_rot, nu)
    integral_dos_vib = simpson(DOS_vib, nu)
    integral_dos_tot = simpson(DOS_tot, nu)


    Log("""
    # DENSITY OF STATES STATISTICS #
                        Translational       Rotational          Vibrational         Total
    Integral of DoS     {:.1f}              {:.1f}              {:.1f}              {:.1f}
    """.format(integral_dos_tr,integral_dos_rot,integral_dos_vib,integral_dos_tot)
    )

    # Plot the DoS plots and export the DoS to file
    plotDOS(nu, DOS_tr, DOS_rot, DOS_vib, DOS_tot)
    np.savetxt(
        fname = 'DoS.txt', 
        X = np.column_stack(
            (nu, DOS_tr, DOS_rot, DOS_vib, DOS_tot)
            ),
        header = 'Frequency (ps-1)     DoS_tr (ps)     DoS_rot (ps)     DoS_vib (ps)     DoS_total (ps)',
    )


    ##################################################################
    """
    2 - P H A S E   T H E R M O D Y N A M I C S   &   E N T R O P Y
    """
    ##################################################################
    """
    1. Decompose S_tr, S_rot, S_vib into S_solid and S_gas
    2. Use the respective S_solid and S_gas to calculate Tr, Rot and Vib entropy using appropriate weighing functions
    """
    # 2-P-T fluidicity determination
    Delta_tr, f_tr, s0_tr, DOS_tr_g, DOS_tr_s = TwoPhaseDecompose_Translational(nu, DOS_tr, temp, Volume, nmols, sum(masslist))
    Delta_rot, f_rot, s0_rot, DOS_rot_g, DOS_rot_s = TwoPhaseDecompose_Rotational(nu, DOS_rot, temp, Volume, nmols, sum(masslist))
    Delta_vib, f_vib, s0_vib, DOS_vib_g, DOS_vib_s = TwoPhaseDecompose_Vibrational(nu, DOS_vib, temp, Volume, nmols, sum(masslist))


    # Supply solid and gaseous DOS separately, calculate entropy using weighted integrals
    S_tr = CalculateEntropy_Translational(nu, DOS_tr_s, DOS_tr_g, f_tr, Delta_tr, sum(masslist), nmols, Volume, temp)
    S_rot = CalculateEntropy_Rotational(nu, DOS_rot_s, DOS_rot_g, Inertia1, Inertia2, Inertia3, 2, temp)
    S_vib = CalculateEntropy_Vibrational(nu, DOS_vib, temp)

    Log("""
    # TWO PHASE THERMODYNAMICS RESULTS #
                        Translational       Rotational          Vibrational
    S_0 (ps)            {:08.2f}            {:08.2f}            {:08.2f}
    Delta               {:1.5f}             {:1.5f}             {:1.5f}
    Fluidicity          {:1.5f}             {:1.5f}             {:1.5f}
    Entropy (J/molK)    {:08.2f}            {:08.2f}            {:08.2f}
    """.format(s0_tr, s0_rot, s0_vib,
            Delta_tr, Delta_rot, Delta_vib,
            f_tr, f_rot, f_vib,
            S_tr/nmols, S_rot/nmols, S_vib/nmols
    ))

