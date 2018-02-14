#!/usr/bin/env python
import sys, re, os
from math import *

# Run as ./alkane N [RCC RCH]

nargin = len(sys.argv) - 1
if (nargin < 1):
    raise Exception('Need to specify N')

N = int(sys.argv[1])
RCC = 1.54
RCH = 1.10
TH = pi - acos(1.0/3.0)

# Bond lengths
if (nargin > 1):
    RCC = float(sys.argv[2])

if (nargin > 2):
    RCH = float(sys.argv[3])

# Geometry
Natom = 3 * (2 * N) + 2
mol = []

PHI = 0.5 * pi - 0.5 * TH
CX = (0.5 * RCC) * cos(PHI)
CZ = (0.5 * RCC) * sin(PHI)

HX = (RCH) * cos(PHI)
HY = (RCH) * cos(PHI)
HZ = (RCH) * sin(PHI) 

S = 1.0 # Phase
DX = CX

for k in range(N): 
    CR  = ['C',  DX,  0.0,  S*CZ]
    CL  = ['C', -DX,  0.0, -S*CZ]
    HR1 = ['H',  DX,   HY,  S*(CZ + HZ)]
    HR2 = ['H',  DX,  -HY,  S*(CZ + HZ)]
    HL1 = ['H', -DX,   HY,  -S*(CZ + HZ)]
    HL2 = ['H', -DX,  -HY,  -S*(CZ + HZ)]

    mol.append(CR)
    mol.append(CL)
    mol.append(HR1)
    mol.append(HR2)
    mol.append(HL1)
    mol.append(HL2)
    
    DX += 2.0 * CX
    S = -S

DX -= 2.0 * CX
DX += HX
HR = ['H',  DX, 0.0, -S*(CZ - HZ)]
HL = ['H', -DX, 0.0,  S*(CZ - HZ)]
mol.append(HR)
mol.append(HL)

# Write config
fh = open('./data/alkane_%d_config.dat' % (Natom), 'w')
fh.write('molecule {\n\n')
for atom in mol:
    fh.write('%1s %24.16f %24.16f %24.16f\n' %(atom[0], atom[1], atom[2], atom[3]))
fh.write('\n}\n')
fh.write('memory 4000 mb\n\n')
fh.write('plugin_load(\"./test_plugin.so\")\n\n')
fh.write('set {\n'); 
fh.write('    basis cc-pVDZ\n')
fh.write('    scf_type pk\n')
fh.write('}\n\n')
fh.write('set test_plugin {\n')
fh.write('    print 1\n')
fh.write('}\n\n')
fh.write('plugin(\"./test_plugin.so\")\n')
fh.close()

# Write NWChem
fh = open('./data/alkane_%d.nw' % (Natom), 'w')
fh.write('title alkane_%d\n' % (Natom))
fh.write('print medium \"screening statistics"\n')
fh.write('set scf:variable logical .true.\n')
fh.write('set fock:replicated logical .false.\n')
fh.write('geometry units an\n')
for atom in mol:
    fh.write('%1s %24.16f %24.16f %24.16f\n' %(atom[0], atom[1], atom[2], atom[3]))
fh.write('end\n')
fh.write('basis\n')
fh.write('    H library cc-pVDZ\n')
fh.write('    C library cc-pVDZ\n')
fh.write('end\n')
fh.write('scf\n')
fh.write('    direct\n')
fh.write('    sym off\n')
fh.write('    tol2e 1e-10\n')
fh.write('    profile\n')
fh.write('    maxiter 1\n')
fh.write('    thresh 100.0\n')
fh.write('end\n')
fh.write('task scf\n')
fh.close()