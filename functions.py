import time
import os
import sys, getopt

import numpy as np
import matplotlib.pyplot as plt

from constants import *

logfile = "log.txt"

def userInputHandler(argv):
    helptext = 'Usage: dodos -x <positions.xvg> -v <velocities.xvg> -b <box.xvg> -T <temp> -m <masses.txt> -t <number of CPU threads>'

    temp = 300
    xfile = 'pos.xvg'
    vfile = 'veloc.xvg'
    bfile = 'box.xvg'
    massfile = 'masses.txt'
    nt = 1

    try:
        opts, args = getopt.getopt(argv,"hx:v:b:m:T:t:",["position=", "velocity=", "box=", "masses=", "temperature=", "threads="])
    except getopt.GetoptError:
        print(helptext)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print (helptext)
            sys.exit()
        elif opt in ("-x", "--position"):
            xfile = arg
        elif opt in ("-v", "--velocity"):
            vfile = arg
        elif opt in ("-b", "--box"):
            bfile = arg
        elif opt in ("-m", "--masses"):
            massfile = arg
        elif opt in ("-t", "--threads"):
            nt = arg
        elif opt in ("-T", "--temperature"):
            temp = float(arg)
        
    masses = []
    with open(massfile, 'r') as file:
        for line in file:
            masses.append(float(line))

    return xfile,vfile,bfile,temp,masses,int(nt)



def cleanLogFile(file = logfile):
    try:
        os.remove(file)
    except OSError:
        pass


def Log(text, file = logfile, console = True):
    if console == True:
        print(text)
    with open(file, mode='a') as file:
        file.write(text)
        file.write('\n')


def loadxvg(filename):
    Log("- Loading {}...".format(filename))
    t1 = time.time()
    data = np.loadtxt(filename, comments=["#", "@"], dtype=np.float64)
    t2 = time.time()
    Log("- {} loaded (took {:.1f} s)!".format(filename, t2-t1))
    return data


def createMassArray(natoms, atomsPerMolecule, masses):
    m = np.empty(natoms)
    for j in range(atomsPerMolecule):
        m[j::atomsPerMolecule] = masses[j]
    Log('Atom masses = {}'.format(masses))
    # np.savetxt('m.txt', m) DEBUG
    return m


def powerspec(x):
    fourier = np.fft.rfft(x)
    power = np.abs(fourier)**2
    return power



def plotDOS(nu, DOS_tr, DOS_rot, DOS_vib, DOS_tot):
    plt.plot(nu/c, c*DOS_tr,  label = 'translational')
    plt.plot(nu/c, c*DOS_rot, label = 'rotational')
    plt.plot(nu/c, c*DOS_vib, label = 'vibrational')

    plt.gca().set(
        xlim = [0,4000],
        title = 'Density of States',
        xlabel = r'Wavenumber $(cm^{-1})$',
        ylabel = r'$S(\nu)$ (cm)',
    )
    plt.legend()
    plt.savefig('DoS_cm.png', dpi = 200, bbox_inches = 'tight')
    plt.clf()

    plt.plot(nu/c, c*DOS_tot, label = 'Total DoS')

    plt.gca().set(
        xlim = [0,4000],
        title = 'Density of States',
        xlabel = r'Wavenumber $(cm^{-1})$',
        ylabel = r'$S(\nu)$ (cm)',
    )
    plt.legend()
    plt.savefig('DoS_total.png', dpi = 200, bbox_inches = 'tight')
    plt.clf()