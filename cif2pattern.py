from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter
import os

##############################################################################
# INSTRUCTIONS: 
# 1) Open script in mantid workbench
# 2) Configure below (lines 25 to 42)
# 3) Execute script

# Expected output: 
# for each sample "sample", creates a workspace:
# "Peaks: sample" with a simulated powder diffraction spectrum
# "Ticks: sample" set of points at positions of expected peaks (and user specified height)
#
# "TotalSimulatedPattern" contains the total spectrum with weighted sum of 
# input sample spectra, simulated background and simulated noise.
# These are scaled to a reference measurement of nickel powder in a DAC  
##############################################################################

#all 3 lists below must have same number of entries...

sampleNames = ['Sapphire (ruby)','a-quartz GeO2'] # a list of sample names
sampleVolRatio = [1.0,1.0] # volume of each phase each ratio should be a float from 0 to 1.0, these needn't sum to 1
cifDir = 'cif/' #path to where cif files live
cifFilenames = ['ICSD_CollCode9770_ruby.cif','ICSD_CollCode59624_a-quartz_GeO2.cif']

#The values below are for the total sample volume, it will be occupied according to 
#the specified sampleVolRatio

totalSampleDiameter = 0.75 # illuminated diameter, enter in mm
totalSampleThickness = 0.05 #enter in mm
relPackingFraction = 1.67 #use this if you think simulated sample will be packed more than reference nickel (e.g. liquid loading)
measTime = 4 #count time in hours

#sets the height that tickmarks will be plotted at
tickYVal = 0.75

# controls additional output which may or may not be helpful 
verbose = False

##############################################################################
#Don't edit below
##############################################################################

#reference measurement dimensions
refDiameter = 0.350     # mm
refThickness = 0.080  # mm

#range for Bragg calculation (doesn't affect range of calculated spectrum)
dMin = 0.75
dMax = 3.5
#
nPts = 600
#bank scattering angle
ttheta = 90.0

#don't change numbers below, these were empirically fitted to the reference data
sigCoef = [0.00002,0.000001,0.00]
relCountTime = measTime/2.5 #reference nickel run was 2.5 hours
noiseAmp = 1e-2/np.sqrt(relCountTime)
scale = 2.0e-4#Chosen to fit reference nickel data)
offset = 0#0.25 #just shifts plot
outWS = 'TotalSimulatedPattern'

def calcSpecFromCif(dMin,dMax,nPts,sigCoef,cifPath,ttheta):
    '''calcSpecFromCif: takes input cif file, specified by cifPath and calculates a simulated spectrum as
    point data in d-space. x-values are specified by dMin (minimum d-spacing), dMax (maximum d-spacing) and nPt
    (number of points). sigCoef are gsas-esque peak width coefficients. Scattering angle ttheta is used in the 
    peak intensity calculation (Lorentz term) '''

    x = np.linspace(dMin,dMax,nPts)
    CreateSampleWorkspace(OutputWorkspace='Phase1')
    LoadCIF(Workspace='Phase1', InputFile=cifPath)
    ws = mtd['Phase1']
    crystal1 = ws.sample().getCrystalStructure()
    unitCell = crystal1.getUnitCell()
    V = unitCell.volume()
    #Generate reflections
    generator = ReflectionGenerator(crystal1)
    # Create list of unique reflections between 0.7 and 3.0 Angstrom
    hkls = generator.getUniqueHKLsUsingFilter(dMin, dMax, ReflectionConditionFilter.StructureFactor)
    # Calculate d and F^2
    dValues = generator.getDValues(hkls)
    fSquared = generator.getFsSquared(hkls)
    pg = crystal1.getSpaceGroup().getPointGroup()

    # Make list of tuples and sort by d-values, descending, include point group for multiplicity.
    reflections = sorted([(hkl, d, fsq, len(pg.getEquivalents(hkl))) for hkl, d, fsq in zip(hkls, dValues, fSquared)],
                                key=lambda x: x[1] - x[0][0]*1e-6, reverse=True)

    nRef = len(reflections)
    simSpec = np.zeros(nPts)

    d = []
    for i in range(nRef):

        Amp = reflections[i][2]*reflections[i][3] #Fsq times multiplicity
        dSpac = reflections[i][1]
        sig2 = sigCoef[0]+sigCoef[1]*dSpac**2+sigCoef[2]*dSpac**4
        sig = np.sqrt(sig2)
        
        if verbose:
            print(f"ref: {i} hkl: {reflections[i][0]}, d-spac = {reflections[i][1]:.3f} Ang, Fsq = {reflections[i][2]:.4f}")
        
        Height = dSpac**4*np.sin(np.radians(ttheta/2.0))*Amp/(np.sqrt(2*np.pi*sig))/V
        singlePeakY = Height*np.exp(-0.5*(x-dSpac)**2/(sig)**2)
        simSpec = simSpec + singlePeakY
        d.append(dSpac) 
    return x,simSpec,d,V


# reference data. This is a dataset that was used to empirically:
# * scale signal to background 
# * estimate noise
# * set peak width coefficients

# need to locate this file as paths are specified relative to it.
scriptLoc = os.path.dirname(__file__)

LoadNexus(Filename=f'{scriptLoc}/data/Nickel_DAC_1bar_48741.nxs', OutputWorkspace='nickel')
Rebin(InputWorkspace='nickel', OutputWorkspace='nickel', Params='0.75,-0.002,2.14')
ExtractSingleSpectrum(InputWorkspace='nickel',OutputWorkspace='nickel',WorkspaceIndex=1)
ConvertToPointData(InputWorkspace='nickel',OutputWorkspace='nickel')

LoadNexus(Filename=f'{scriptLoc}/data/Nickel_DAC_1bar_48741_fitted_bgnd.nxs', 
OutputWorkspace='nickel_bgnd')
ConvertToPointData(InputWorkspace='nickel_bgnd',OutputWorkspace='nickel_bgnd')

#phase 1: reference nickel data
cifPath = f'{scriptLoc}/cif/nickel_CollCode37502.cif'
x,niSpec,d1,V_ni = calcSpecFromCif(dMin,dMax,nPts,sigCoef,cifPath,ttheta)

#get background, rebin to match simulation   
ws = mtd['nickel_bgnd']
bgdx = ws.dataX(0)
bgdy = ws.dataY(0)
bgdprof = np.interp(x,bgdx,bgdy)
  
# calculate noise
noise = np.random.normal(scale=noiseAmp,size=nPts)

#very rough absorption calculation
mu1 = 0.7
gashAbs = np.exp(-mu1*x)

CreateWorkspace(OutputWorkspace='nickel_sim',
    DataX=x,
    DataY=scale*np.multiply(niSpec,gashAbs)+bgdprof+noise+offset)    

#ticks
yloc = np.ones(len(d1))

CreateWorkspace(OutputWorkspace='Ticks: Nickel',
    DataX=np.array(d1),
    DataY=1.4*yloc) 

#
# Simulate phase(s) and scale these to nickel
#
refVol = np.pi*(refDiameter/2)**2*refThickness
samVol = np.pi*(totalSampleDiameter/2)**2*totalSampleThickness
samVolumeFactor = relPackingFraction*(samVol/refVol)

backgroundFactor = (totalSampleDiameter/refDiameter)**2
noiseFactor = np.sqrt(backgroundFactor)

totalSimBragg = np.zeros_like(niSpec)
for i,cifFile in enumerate(cifFilenames):
    
    #calculate simulated Bragg intensities for phase_spectrum
    x,simSpec,d2,V_sim = calcSpecFromCif(dMin,dMax,nPts,sigCoef,f"{scriptLoc}/{cifDir}{cifFile}",ttheta)
    
    #apply scale factor and absorption
    simSampleSpec = scale*np.multiply(simSpec, gashAbs)

    #calculate factor to scale simultated sample to reference
    F = V_ni/V_sim

    #scale Bragg scattering for factor F and phase fractional_indexing and add to total spectrum
    simSampleSpec = samVolumeFactor*sampleVolRatio[i]*F*simSampleSpec
    totalSimBragg = totalSimBragg+simSampleSpec
    
    #create workspace with d-spacings of reflections (for plotting ticks)
    yloc = np.ones(len(d2))
    CreateWorkspace(OutputWorkspace=f"Ticks: {sampleNames[i]}",
    DataX=np.array(d2),
    DataY=tickYVal*yloc,
    UnitX='d-spacing')  

    #create workspace with calculated Bragg peaks for each phase
    CreateWorkspace(OutputWorkspace=f"Peaks: {sampleNames[i]}",
    DataX=x,
    DataY=simSampleSpec,
    UnitX='d-spacing')

#Create total pattern: adding background, noise and scaling to reference volume

totalSimSpec = totalSimBragg+backgroundFactor*bgdprof+noiseFactor*noise+offset #should correct for greater or larger number of unit cells

CreateWorkspace(OutputWorkspace=outWS,
    DataX=x,
    DataY=totalSimSpec)    

#clean up 
DeleteWorkspace(Workspace='Phase1')
DeleteWorkspace(Workspace='Nickel')

if verbose:
    print(f'ref volume: {refVol:.4f}')
    print(f'sim volume: {samVol:.4f}')
    print(f'sampleVolumeFactor: {samVolumeFactor:.4f}')
    print(f'backgroundFactor: {backgroundFactor:.4f}')
    print(f'noiseFactor: {noiseFactor:.4f}')
else:
    DeleteWorkspace(Workspace='nickel_bgnd')
    DeleteWorkspace(Workspace='nickel_sim')
