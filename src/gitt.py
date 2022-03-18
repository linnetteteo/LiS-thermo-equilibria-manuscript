import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import re
import pandas as pd

## GITT paper helper functions

def readFullRelaxFiles(folder_path):
    """
    Reads files for GITT relaxation experiment, simulated with 1D full model.
    Returns arrays for the run/filename, run number,
    dod at which cell is discharged too before current is turned off,
    and crate of discharge
    """

    run_arr = []
    Nrun_arr = []
    dod_arr = []
    crate_arr = []
    count=0

    # find number of files that starts with run
    # (this is the data file we want to read)
    for file in os.listdir(folder_path):
        if file.startswith("relaxrun"):
            count+=1

    # order the data files by run number, so we get descending crates
    Nrun=1
    for i in range(count+5):
        for file in os.listdir(folder_path):
            if file.startswith("relaxrun_"+str(Nrun)+"-"):
                run_arr.append(file)
                dod = re.search('dod=(.*).txt', file).group(1)
                crate = re.search('Crate=(.*)_',file).group(1)
                Nrun_arr.append(np.round(int(Nrun),decimals=0))
                dod_arr.append(float(dod))
                crate_arr.append(float(crate))
        Nrun+=1
    print(len(run_arr))

    return run_arr, Nrun_arr, dod_arr, crate_arr

def readNSRelaxFiles(folder_path):
    """
    Reads files for GITT relaxation experiment,
    simulated with 0D no spatial variation/thermo-kinetic model.
    Returns arrays for the run/filename, run number,
    dod at which cell is discharged too before current is turned off,
    and crate of discharge
    """

    run_arr = []
    Nrun_arr = []
    dod_arr = []
    crate_arr = []
    count=0

    # find number of files that starts with run
    # (this is the data file we want to read)
    for file in os.listdir(folder_path):
        if file.startswith("relaxrun"):
            count+=1

    # order the data files by run number, so we get descending crates
    Nrun=1
    for i in range(count+5):
        for file in os.listdir(folder_path):
            if file.startswith("relaxrun_"+str(Nrun)+"_"):
                run_arr.append(file)
                dod = re.search('dod=(.*).txt', file).group(1)
                crate = re.search('Crate=(.*)_',file).group(1)
                Nrun_arr.append(np.round(int(Nrun),decimals=0))
                dod_arr.append(float(dod))
                crate_arr.append(float(crate))
        Nrun+=1
    print(len(run_arr))

    return run_arr, Nrun_arr, dod_arr, crate_arr


class GittData:
    """Data type that has current on till given dod and
    current turned off during a following relaxation step.
    endt is time taken for full discharge (theoretically)"""

    def __init__(self, data, dod, crate, endt):
        self.data = data
        self.dod = dod
        self.crate = crate
        self.endt = endt
        # self.time = []
        # self.voltage = []

    def split_relax(self, zeroMode = False):
        """
        Split up voltage and time arrays into two arrays each.
        Indicated by dod when current is turned off.
        Dod is based on endt, total discharge time.
        zeroMode turned on means both arrays are shifted s.t. relax portion starts from 0.
        Assigns new attributes timeOn, timeOff, voltageOn, voltageOff to class
        """

        endt_dod = (self.dod/100*self.endt)
        arr = self.time - endt_dod
        idx = np.where(arr >= 0, arr, np.inf).argmin()

        self.timeOn = self.time[:idx]
        self.timeOff = self.time[idx:]
        self.voltageOn = self.voltage[:idx]
        self.voltageOff = self.voltage[idx:]

        if zeroMode:
            self.timeOn = arr[:idx]
            self.timeOff = arr[idx:]
            self.voltageOn = self.voltage[:idx] - self.voltage[idx]
            self.voltageOff = self.voltage[idx:] - self.voltage[idx]

        return


class GittDataNoSpatial(GittData):
    """This data is simulated using the 0D model with no spatial variation.
    data is a mxn array with m being time points and n=13.
    n columns:
    ["time", "C1[Li]", "C1[S8l]", "C1[S8m]",
    "C1[S6m]", "C1[S4m]", "C1[S2m]", "C1[Sm2]",
    "epc[S8solid]", "epc[Sm2solid]", "ec",
    "phiC2", "phi1"]
    tau0 is needed to to scale time since model uses dimensionless time.
    """

    def __init__(self, data, dod, crate, endt, tau0):
        GittData.__init__(self, data, dod, crate, endt)
        self.time = data[:,0]*tau0/3600
        self.voltage = data[:,12]
        self.tau0 = tau0


class GittDataFull1D(GittData):
    """This data is simulated using the 1D full model (no separator).
    data: n x (m+1) matrix, with n being number of time points,
    m being number of dependent variables
    (at variables points in x, so multiplied by Nx or Nx+2);
    col 0 is the time information
    Nx: integer, number of (internal) collocation points in x
    [see plotVariablesWithTimeAndX for more info about variables order]

    """

    def __init__(self, data, dod, crate, endt, Nx):
        GittData.__init__(self, data, dod, crate, endt)
        self.time = data[:,0]/3600
        self.voltage = data[:,-1]
        self.Nx = Nx

class GittDataFull1D_catsep(GittData):
    """This data is simulated using the 1D full model (with separator).
    data: n x (m+1) matrix, with n being number of time points,
    m being number of dependent variables
    (at variables points in x, so multiplied by Nx or Nx+2);
    col 0 is the time information
    Nx: integer, number of (internal) collocation points in x
    [see plotVariablesWithTimeAndX for more info about variables order]

    """

    def __init__(self, data, dod, crate, endt, Nx):
        GittData.__init__(self, data, dod, crate, endt)
        self.time = data[:,0]/3600
        ndepvars, neps = 23, 6
        self.voltage = data[:,ndepvars*Nx + (ndepvars-neps)*2 -1]
        self.Nx = Nx

class GittDataFull1D_catsep_AvgVars(GittData):
    """This data is simulated using the 1D full model (with separator).
    data: n x (m+1) matrix, with n being number of time points,
    m being number of dependent variables
    (at variables points in x, so multiplied by Nx or Nx+2);
    col 0 is the time information
    Nx: integer, number of (internal) collocation points in x
    [see plotVariablesWithTimeAndX for more info about variables order]

    """

    def __init__(self, data, dod, crate, endt, Nx):
        GittData.__init__(self, data, dod, crate, endt)
        self.time = data[:,0]/3600
        self.voltage = data[:,-2]
        self.Nx = Nx
