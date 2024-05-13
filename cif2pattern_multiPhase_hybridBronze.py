# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter

#Configure simluation
dMin = 0.75
dMax = 3.5#2.12
nPts = 600
#Gaussian peak width parameters
sigCoef = [0.00002,0.000001,0.00]
#bank scattering angle
ttheta = 90.0
noiseAmp = 1e-2
scale = 2.0e-4#Chosen to fit reference nickel data)
offset = 0#0.25 #just shifts plot
outWS = 'Multiphase sim'



def calcSpecFromCif(dMin,dMax,V,nPts,sigCoef,cifPath,ttheta):

    x = np.linspace(dMin,dMax,nPts)
    CreateSampleWorkspace(OutputWorkspace='Phase1')
    LoadCIF(Workspace='Phase1', InputFile=cifPath)
    ws = mtd['Phase1']
    crystal1 = ws.sample().getCrystalStructure()
    #Generate reflections
    generator = ReflectionGenerator(crystal1)
    # Create list of unique reflections between 0.7 and 3.0 Angstrom
    hkls = generator.getUniqueHKLsUsingFilter(dMin, dMax, ReflectionConditionFilter.StructureFactor)
    print(hkls)
    # Calculate d and F^2
    dValues = generator.getDValues(hkls)
    fSquared = generator.getFsSquared(hkls)
    pg = crystal1.getSpaceGroup().getPointGroup()

    # Make list of tuples and sort by d-values, descending, include point group for multiplicity.
    reflections = sorted([(hkl, d, fsq, len(pg.getEquivalents(hkl))) for hkl, d, fsq in zip(hkls, dValues, fSquared)],
                                key=lambda x: x[1] - x[0][0]*1e-6, reverse=True)

    # print('{0:<8}{1:>8}{2:>8}{3:>4}'.format('HKL', 'd', 'F^2', 'M'))
    # for reflection in reflections:
    #     print('{0!s:<8}{1:>8.5f}{2:>8.2f}{3:>4}'.format(*reflection))

    nRef = len(reflections)
    simSpec = np.zeros(nPts)

    d = []
    for i in range(nRef):

        Amp = reflections[i][2]*reflections[i][3] #Fsq times multiplicity
        dSpac = reflections[i][1]
        sig2 = sigCoef[0]+sigCoef[1]*dSpac**2+sigCoef[2]*dSpac**4
        sig = np.sqrt(sig2)
        
        # height set to equal model observed Bragg intensity using: 
        # L*Mult*Fsq/V following GSAS manual page 133.
        
        Height = dSpac**4*np.sin(np.radians(ttheta/2.0))*Amp/(np.sqrt(2*np.pi*sig))/V
        singlePeakY = Height*np.exp(-0.5*(x-dSpac)**2/(sig)**2)
        simSpec = simSpec + singlePeakY
        d.append(dSpac) 
    return x,simSpec,d

def calcSpecMan(dMin,dMax,V,nPts,sigCoef,cifPath,ttheta):

    x = np.linspace(dMin,dMax,nPts)
    CreateSampleWorkspace(OutputWorkspace='Phase1')
    LoadCIF(Workspace='Phase1', InputFile=cifPath)
    ws = mtd['Phase1']
    
    crystal1 = CrystalStructure("3.396 9.254 9.006","P 1 21/c 1",
    """
    N 0.714 0.129 0.413 1 0.05;
    D 0.821 0.042 0.380 1 0.05;
    D 0.792 0.133 0.528 1 0.05;
    D 0.406 0.146 0.401 1 0.05;
    O 0.039 0.867 0.238 1 0.05;
    D 0.135 0.783 0.185 1 0.05;
    D 0.150 0.956 0.195 1 0.05;
    N 0.665 0.875 0.008 1 0.05;
    D 0.751 0.876 0.008 1 0.05;
    D 0.715 0.772 0.855 1 0.05;
    D 0.380 0.896 0.890 1 0.05
    
    """)


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

    # print('{0:<8}{1:>8}{2:>8}{3:>4}'.format('HKL', 'd', 'F^2', 'M'))
    # for reflection in reflections:
    #     print('{0!s:<8}{1:>8.5f}{2:>8.2f}{3:>4}'.format(*reflection))

    nRef = len(reflections)
    simSpec = np.zeros(nPts)
    
    d = []
    for i in range(nRef):

        Amp = reflections[i][2]*reflections[i][3] #Fsq times multiplicity
        dSpac = reflections[i][1]
        sig2 = sigCoef[0]+sigCoef[1]*dSpac**2+sigCoef[2]*dSpac**4
        sig = np.sqrt(sig2)
        Height = dSpac**4*np.sin(np.radians(ttheta/2.0))*Amp/(np.sqrt(2*np.pi*sig))/V
        singlePeakY = Height*np.exp(-0.5*(x-dSpac)**2/(sig)**2)
        simSpec = simSpec + singlePeakY
        d.append(dSpac)    
    return x,simSpec,d

# reference data
LoadNexus(Filename='/Users/66j/Documents/ORNL/code/cif2pattern/Nickel_DAC_1bar_48741.nxs', OutputWorkspace='nickel')
Rebin(InputWorkspace='nickel', OutputWorkspace='nickel', Params='0.75,-0.002,2.14')
ExtractSingleSpectrum(InputWorkspace='nickel',OutputWorkspace='nickel',WorkspaceIndex=1)
ConvertToPointData(InputWorkspace='nickel',OutputWorkspace='nickel')
LoadNexus(Filename='/Users/66j/Documents/ORNL/code/cif2pattern/Nickel_DAC_1bar_48741_fitted_bgnd.nxs', 
OutputWorkspace='nickel_bgnd')
ConvertToPointData(InputWorkspace='nickel_bgnd',OutputWorkspace='nickel_bgnd')

#phase 1: reference nickel data
cifPath = '/Users/66j/Downloads/nickel_260169.cif'
V_ni = 43.9 #A^3
x,simSpec1,d1 = calcSpecFromCif(dMin,dMax,V_ni,nPts,sigCoef,cifPath,ttheta)

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

volFac = 4

#simulate nickel data to compare with measurement and to estimate scale factor
# simSpec1 = volFac*0.0*np.multiply(simSpec1,gashAbs)+bgdprof+noise
# simSPec1 bgdprof+noise

CreateWorkspace(OutputWorkspace='nickel_sim',
DataX=x,
DataY=scale*np.multiply(simSpec1,gashAbs)+bgdprof+noise+offset)    

#ticks
yloc = np.ones(len(d1))

CreateWorkspace(OutputWorkspace='Ticks: Nickel',
DataX=np.array(d1),
DataY=1.4*yloc) 

#
# Simulate other phase(s) and scale these to nickel
#

#phase 2: hybridBronze
cifPath = '/Users/66j/Downloads/Hx(4,4\'-bipyridine)0p5MoO3_290K_I4mmm.cif'
V_sim = 618.432 #A^3
x,simSpec2,d2 = calcSpecFromCif(dMin,dMax,V_sim,nPts,sigCoef,cifPath,ttheta)


simSpec2 = scale*np.multiply(simSpec2,gashAbs)

#IMPORTANT NOTE: calculation is based on number of unit cells in the beam.
#For Reference Nickel:
#
# have N_ni unit cells, with volume V_ni occupying a total volume of V_tot = N_ni*V_ni
#
# for the same illuminated volume, a simulated sample with unit cell volume V_sim we can say
#
# V_sim*N_sim = V_ni*N_ni
#
# Therefore:
#
# F = N_sim/N_ni = V_ni/V_sim
# 
# will be needed to scale simulation for variable number of unit cells in the beam and multiplying by this
# will represent the scattering from the same sample volume as the nickel reference.
#
#
F = V_ni/V_sim

#
# Here I simulate the situation for a doubling of sample area and 50% increase in sample thickness
# a factor of 6. Assume background increases with area = factor of 4. Assume noise is 100% dependent on backgrounf
# so will increase by a factor of sqrt(4) assuming equivalent measurement time.

simSpec2 = 6*F*simSpec2+4*bgdprof+2*noise+offset #should correct for greater or larger number of unit cells


#use this to scale to arbitrary phase   
CreateWorkspace(OutputWorkspace='sim_hybridBronze',
    DataX=x,
    DataY=simSpec2)    

#ticks
yloc = np.ones(len(d2))

CreateWorkspace(OutputWorkspace='Ticks: Hbronze',
DataX=np.array(d2),
DataY=1.0*yloc,
UnitX='d-spacing')    
