#!/usr/bin/env python3

import os
import sys

# Read the gromacs input file 
#groFile = "test.gro"
groFile = sys.argv[1]
f = open(groFile, "r")
header1 = f.readline()
nbAtoms = f.readline()
nbAtoms = int(nbAtoms)
print("Number of atoms in gromacs file= ", nbAtoms)

countAtom = 0
idInRes = list()
resId = list()
resName = list()
atomType = list()
atomNb = list()
x = list()
y = list()
z = list()
mapName = list()

res = -1
resOld = -1
idInResCount = 1

for line in f:
	if countAtom < nbAtoms:
		res = int(line[:5])
		if res != resOld:
			idInResCount = 1
		resId.append(res)
		resName.append(line[5:10].strip())
		atomType.append((line[10:15]))
		atomNb.append(int(line[15:20]))
		x.append(10.0*float(line[20:28]))
		y.append(10.0*float(line[28:36]))
		z.append(10.0*float(line[36:44]))
		idInRes.append(idInResCount)
		countAtom += 1
		idInResCount += 1
		resOld = res
		if not line[5:10].strip() in mapName:
			mapName.append(line[5:10].strip())

boxLine = line.split(' ')
boxLine = [c for c in boxLine if c != '']
print(boxLine)
xBox = 10.0*float(boxLine[0])
yBox = 10.0*float(boxLine[1])
zBox = 10.0*float(boxLine[2])		
print("Simulation box size= ", xBox, yBox, zBox)
print("Residues in simulation= ", mapName)	
f.close()


mapAll = {}
mapAllConj = {}
frag = ""
fragList = []
nbConj = 0


# Read the map file 
for map in mapName:
	mapFile = map + ".map"
	f = open(mapFile, "r")
	mapDico = {}
	mapDicoConj = {}
	#for line in f:
	while True:
		line = f.readline()
		if not line: break
		if len(line) > 1:
			if line[0] == "[":
				frag = line[1:-2]
				nbAtFrag = int(f.readline())
				for i in range(1,nbAtFrag+1):
					line = f.readline()
					line = line.split()
					key = int(line[0])
					mapDico[key] = [line[1], int(line[2]), frag]
				fragList.append(frag)
			if line[0] == "{":
				conj = line[1:-2]
				nbConj = nbConj + 1
				nbFrag = int(f.readline())
				for i in range(1,nbFrag+1):
					line = f.readline()
					mapDicoConj[line.strip()] = [conj, nbConj]
	f.close()
	mapAll[map] = mapDico
	mapAllConj[map] = mapDicoConj


# Create a pdb file with only the mapped fragments, and saturate bonds with vmd (create tcl script)
resOld = -1
### CUBIC mod
#pdbFile = groFile[:-4] + ".pdb"
pdbFile = "saturated.pdb"
satFile = "sat.tcl"
satCount = 0
atomList = list()

#print(mapAll)
#print(mapAllConj)
#print(fragList)


for i,res in enumerate(resName):		
	

	if i == 0:

		idxTot = 0
		countAtom = 1
		resTypeSaturated = list()
		#satFile = "sat" + str(satCount) + ".tcl"
		#pdbFile = groFile[:-4] + str(satCount) + ".pdb"
		tcl = open(satFile, "w")
		header = """
		# Load package
		package require molefacture

		# Load a molecule
		"""
		tcl.write(header)
		tcl.write('%s%s%s\n' % ("set mol [mol new ",pdbFile," type pdb waitfor all]"))
		header = """
		set sel0 [atomselect top all]

		# Open in molefacture
		::Molefacture::molefacture_gui $sel0
		"""
		tcl.write(header)

		f = open(pdbFile, "w")
		f.write('%s%9.3f%9.3f%9.3f%s\n' % ("CRYST1",xBox,yBox,zBox,"  90.00  90.00  90.00 P 1           1") )	



	
	if (resId[i] != resOld and i != 0) or (i==nbAtoms-1):

		if i==nbAtoms-1:
			if res in mapAll:
				if idInRes[i] in mapAll[res]:
					#if True:
					if mapAll[res][idInRes[i]][2] in mapDicoConj:
						#print(mapAll[res][idInRes[i]][2])
						f.write( '%4s%7d%5s%6s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%10s\n' % ("ATOM", countAtom, mapAll[res][idInRes[i]][0], mapAll[res][idInRes[i]][2][:3]+"0", resId[i], x[i], y[i], z[i] ,0.0,float(mapAll[res][idInRes[i]][2][3:]),mapAll[res][idInRes[i]][0]) )
						countAtom += 1
						idxTot += 1
			
						# Do I need to saturate bond?
						if mapAll[res][idInRes[i]][1] > 0:
							for r in range(mapAll[res][idInRes[i]][1]): 
								tcl.write( '%s %d \n %s \n' % ("Molefacture::select_atoms", countAtom-2, "Molefacture::add_hydrogen_gui") )
								idxTot += 1
								resTypeSaturated.append(float(mapAll[res][idInRes[i]][2][3:]))

		tcl.write('\n%s%d%s\n%s\n%s' % ("set sel1 [atomselect top index<", idxTot, "]", "$sel1 writepdb saturated.pdb", "exit") )
		tcl.close()			
		f.close()
		command = "vmd -e " + satFile
		print("la")
		print(command)
		print(satCount)
		### CUBIC mod		
		#os.system(command)	
		satCount += 1

		
		f = open("saturated.pdb", "r")
		f.readline()
		countSaturatedRes = 0
		for line in f:
			if len(line) > 4:
				# CUBIC mod
				#elem = line[13:14]
				elem = line[15:16]
				#resN = line[17:20]
				resN = line[18:21]
				resX = line[21:22]
				xt = float(line[30:38])
				yt = float(line[38:46])
				zt = float(line[46:54])
				rest = int(float(line[60:66]))
				chain = int(line[22:26])
				if resX == "X":
					rest = int(resTypeSaturated[countSaturatedRes])
					countSaturatedRes += 1
				crAt = resN+str(rest)
				#frAt = -1
				#if crAt in mapDicoConj:
				#      	frAt = list(mapDicoConj.keys()).index(crAt) + 1
				frAt = fragList.index(crAt) + 1
				atomList.append([rest, resN, elem, xt, yt, zt, chain,frAt])
		atomList.sort(key=lambda x: x[7])		
		f.close()
		
		idxTot = 0
		countAtom = 1
		resTypeSaturated = list()
		#satFile = "sat" + str(satCount) + ".tcl"
		#pdbFile = groFile[:-4] + str(satCount) + ".pdb"
		tcl = open(satFile, "w")
		header = """
		# Load package
		package require molefacture

		# Load a molecule
		"""
		tcl.write(header)
		tcl.write('%s%s%s\n' % ("set mol [mol new ",pdbFile," type pdb waitfor all]"))
		header = """
		set sel0 [atomselect top all]

		# Open in molefacture
		::Molefacture::molefacture_gui $sel0
		"""
		tcl.write(header)

		f = open(pdbFile, "w")
		f.write('%s%9.3f%9.3f%9.3f%s\n' % ("CRYST1",xBox,yBox,zBox,"  90.00  90.00  90.00 P 1           1") )		
				
	resOld = resId[i]





	if res in mapAll and i!= nbAtoms-1:
		if idInRes[i] in mapAll[res]:
			#if True:
			if mapAll[res][idInRes[i]][2] in mapDicoConj:
				#print(mapAll[res][idInRes[i]][2])
				f.write( '%4s%7d%5s%6s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%10s\n' % ("ATOM", countAtom, mapAll[res][idInRes[i]][0], mapAll[res][idInRes[i]][2][:3]+"0", resId[i], x[i], y[i], z[i] ,0.0,float(mapAll[res][idInRes[i]][2][3:]),mapAll[res][idInRes[i]][0]) )
				countAtom += 1
				idxTot += 1
			
				# Do I need to saturate bond?
				if mapAll[res][idInRes[i]][1] > 0:
					for r in range(mapAll[res][idInRes[i]][1]): 
						tcl.write( '%s %d \n %s \n' % ("Molefacture::select_atoms", countAtom-2, "Molefacture::add_hydrogen_gui") )
						idxTot += 1
						resTypeSaturated.append(float(mapAll[res][idInRes[i]][2][3:]))
				

			




#atomList.sort(key=lambda x: x[6])
### CUBIC mod
#finalPdb = pdbFile[:-4] + "_fragmented.pdb"
finalPdb = groFile[:-4] + "_fragmented.pdb"
#finalPdb = "fragmented.pdb"
f = open(finalPdb, "w")
f.write('%s%9.3f%9.3f%9.3f%s\n' % ("CRYST1",xBox,yBox,zBox,"  90.00  90.00  90.00 P 1           1") )

countAtom = 1
resId = 0
chainOld = -1
chainCurrent = 1
fragmentOld = -1
for atom in atomList:
	chain = atom[6]
	#fragment = atom[0]	
	if chain != chainOld and countAtom != 1:
		resId += nbConj
		chainCurrent += 1
	conjugatedSiteId = -1
	conjugatedName = "NUL"
	#if fragList[fragment-1] in mapDicoConj:
	currentFragment = atom[1]+str(atom[0])
	#fragment = -1
	fragment = fragList.index(currentFragment) + 1 
	if currentFragment in mapDicoConj:
		#fragment = list(mapDicoConj.keys()).index(currentFragment) + 1
		conjugatedSiteId = resId+mapDicoConj[currentFragment][1]
		conjugatedName = mapDicoConj[currentFragment][0]
	f.write( '%4s%7d%5s%6s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%10s\n' % ("ATOM", countAtom, atom[2], conjugatedName+"0", conjugatedSiteId % 10000, atom[3], atom[4], atom[5], float(fragment), chainCurrent % 10000, atom[2] ))
	countAtom += 1
	chainOld = chain
	fragmentOld = fragment

f.close()
