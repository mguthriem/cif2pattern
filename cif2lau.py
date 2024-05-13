# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter
from mantid import config
config.setLogLevel(0, quiet=True)

#specify full path for input file
cifPath = '/Users/66j/Downloads/nickel_260169.cif'
#specify full path for output file
lauPath = '/Users/66j/Downloads/nickel_260169.lau'
#specify min and max d-spacing limits to process
dMin = 0.8
dMax = 2.2
##########################################################

#Phase1
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


lines = ['H K L Multiplicity Sin(Theta)/Lambda d_spacing |F|^2'] #header line    
for ref in reflections:
    stol = 1/(2*ref[1])
    h = int(ref[0][0])
    k = int(ref[0][1])
    l = int(ref[0][2])
    lines.append(f'{h:4d} {k:4d} {l:4d} {ref[3]:4d} {stol:.4f} {ref[1]:.4f} {ref[2]:.4f}')

allLines = '\n'.join(lines)
print(allLines)
with open(lauPath,'w+') as f:
    f.writelines(allLines)
print(f'Created file: {lauPath}')
    