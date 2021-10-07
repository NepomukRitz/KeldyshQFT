import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
import time
from scipy.integrate import quad

'''
This script compiles and runs a C++ program to get an integrands.
It can only be used from the workstation.

Workflow:
1. run fRG flow, saving intermediate results
2. run this script to obtain integrands for an arbitrary set of external parameters
'''

def plot_integrand(ax, filename, kwargs={}):
    with h5py.File(filename, 'r') as f:
        freqs = np.array(f['v'])
        integrand_re = np.array(f['integrand_re'])
        integrand_im = np.array(f['integrand_im'])
    ax.plot(freqs, integrand_re, **kwargs)
    ax.plot(freqs, integrand_im, **kwargs)

def plot_integrand_alli2(ax, filename_pt1, filename_pt2, kwargs={}):
    nfreqs = len(freqs)
    ax.plot(freqs, integrand_re, **kwargs)
    ax.plot(freqs, integrand_im, **kwargs)

def compile_get_integrand_dGamma_1loop_Cpp():
    text = """
#!/bin/bash

module load hdf5/1.10.5-gcc

export LANG=C
export LC_ALL=C

MPI="-I/usr/include/mpi/"
HDF5="-I/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/include -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp"
FFTW="-lfftw3"

mpiCC --std=c++17  ../tests/get_integrand_dGamma_1loop.cpp -o ../tests/get_integrand_dGamma_1loop.o -fopenmp $MPI $HDF5 $FFTW

    """
    f = open("./compile_asc_saveIntegrand.sh", "w")
    f.write(text)
    f.close()
    print("Compiling get_integrand_dGamma_1loop")
    os.system("bash ./compile_asc_saveIntegrand.sh")
    os.system("rm ./compile_asc_saveIntegrand.sh")

def get_integrand_dGamma_1loop(exe_dir, data_dir, it_Lambda, rkStep, k_class_int, channel, i0, i2, w, v, vp, i_in):
    cmd = exe_dir + 'get_integrand_dGamma_1loop.o ' + data_dir + " " + str(it_Lambda) + " " + str(rkStep) + " " + str(k_class_int) + " " + str(channel) + " " + str(i0) + " " + str(i2) + " {}".format(w) + " {}".format(v) + " {}".format(vp) + " " + str(i_in)
    os.system(cmd)

def get_integrand_filename_dGamma_1loop(data_dir, it_Lambda, rkStep, k_class_int, channel, i0, i2, w, v, vp, i_in):
    if k_class_int == 0:
        k_class = '1'
    elif k_class_int == 1:
        k_class = '2'
    elif k_class_int == 3:
        k_class = '3'
    else:
        raise Exception("k_class_int " + str(k_class_int) + "not allowed!")
    filename = data_dir + "integrands/dGamma1Loop_iLambda" + str(1) + "_RKstep" + str(rkStep) + "integrand_K" + k_class + channel + "_i0=" + str(i0) + "_i2=" + str(i2) + "_w={:.6f}".format(w)
    if k_class == '2':
        filename = filename + "_v={:.6f}".format(v)
    if k_class == '3':
        filename = filename + "_vp={:.6f}".format(v)
    filename = filename + ".h5"
    return filename

if __name__ == '__main__':
    exe_dir = "/tmp/tmp.0CYjLbB99I/tests/"  # directory where the C++ executable is / will be compiled
    os.chdir(exe_dir)
    #compile_get_integrand_dGamma_1loop_Cpp()

    data_dir = "/tmp/tmp.0CYjLbB99I/Data_MF/"
    it_Lambda = 1 # iteration of ODE solver
    rkStep = 0 # Runge-Kutta step >= 0
    k_class_int = 3 # K1->0; K2->1; K3->3
    channel = 'a'
    i0 = 0
    i2 = 0
    w = 100
    v = 0
    vp = 0
    i_in = 0
    get_integrand_dGamma_1loop(exe_dir, data_dir, it_Lambda, rkStep, k_class_int, channel, i0, i2, w, v, vp, i_in)

    integrandFilename = get_integrand_filename_dGamma_1loop(data_dir, it_Lambda, rkStep, k_class_int, channel, i0, i2, w, v, vp, i_in)

    fig1, axs1 = plt.subplots(num=1, nrows=1, ncols=1)
    plot_integrand(axs1, integrandFilename, {"label": "testtest", "ls": "--"})
    fig1.savefig(data_dir + "integrands/test.png")



    #plt.show()