# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter

#basic parameters
dMin = 0.82
dMax = 2.12
nPts = 1000
sigCoef = [0.00004,0.000001,0.00]
bgdCoef = [1.038751,
-0.972576,
0.456113,
-0.194804,
0.069082,
-0.011982
]

#set up reference data: 

LoadNexus(Filename='/Users/66j/Documents/ORNL/code/cif2pattern/Nickel_DAC_1bar_48741.nxs', OutputWorkspace='nickel')
Rebin(InputWorkspace='nickel', OutputWorkspace='nickel_r', Params='0.75,-0.002,2.13')
SumSpectra(InputWorkspace='nickel_r', OutputWorkspace='nickel_r_sum', StartWorkspaceIndex=1, EndWorkspaceIndex=4, UseFractionalArea=False)


ttheta = 90.0
noiseAmp = 1e-1#1e-5
scale = 1.4e-5
outWS = 'sim3'
cifPath = '/Users/66j/Downloads/Hx(4,4\'-bipyridine)0p5MoO3_290K_I4mmm.cif'

x = np.linspace(dMin,dMax,nPts)
# y = np.zeros(nPts)

#Phase1
V1 = 618.432 #A^3
CreateSampleWorkspace(OutputWorkspace='Phase1')
LoadCIF(Workspace='Phase1', InputFile=cifPath)
ws = mtd['Phase1']
crystal1 = ws.sample().getCrystalStructure()
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

for i in range(nRef):

    Amp = reflections[i][2]*reflections[i][3] #Fsq times multiplicity
    dSpac = reflections[i][1]
    sig2 = sigCoef[0]+sigCoef[1]*dSpac**2+sigCoef[2]*dSpac**4
    sig = np.sqrt(sig2)
    Height = dSpac**4*np.sin(np.radians(ttheta/2.0))*Amp/(np.sqrt(2*np.pi*sig))
    singlePeakY = Height*np.exp(-0.5*(x-dSpac)**2/(sig)**2)
    simSpec = simSpec + singlePeakY
   
   
# chebychev background


# bgdprof = np.polynomial.chebyshev.chebval(x,bgdCoef)
 

# Fitted background 
  
ws = mtd['Nickel_r']

bgdx = ws.dataX(3)[0:-1]
bgdy = ws.dataY(3)

bgdprof = np.interp(x,bgdx,bgdy)
  
# noise
noise = np.random.normal(scale=noiseAmp,size=nPts)

mu1 = 0.7
gashAbs = simSpec*np.exp(-mu1*x)
   
CreateWorkspace(OutputWorkspace=outWS,
DataX=x,
# DataY=scale*simSpec+bgdprof+noise)    
DataY=scale*gashAbs+bgdprof+noise+1)    