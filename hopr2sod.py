'''
Script to convert HOPR output format (https://hopr.readthedocs.io/en/latest/userguide/meshformat.html) to SOD2D.
'''

import numpy as np
import h5py
import matplotlib.pyplot as plt

from sides_mapping.hopr import *

##HOPR input file 
BASEDIR  = '/home/benet/Dropbox/UNIVERSITAT/PhD/windsor/examples/Test_hopr'
CASESTR  = 'CYLINDER_CURVED'
#BASEDIR  = '/home/benet/Dropbox/UNIVERSITAT/PhD/windsor/examples/Test_hopr/cylind_jordi'
#CASESTR  = 'CYLINDER'
hoprname = '%s/%s_mesh.h5' % (BASEDIR, CASESTR)

##Specify element mappings:
hex64_map = np.loadtxt('element_mapping/hopr_gmsh/hex64.txt', dtype=np.int32)-1
hex64_map = np.argsort(hex64_map)
qua16_map = np.loadtxt('element_mapping/hopr_gmsh/qua16.txt', dtype=np.int32)-1
qua16_map = np.argsort(qua16_map)

##Debug output to GMSH
visugmsh     = True
gmsh_voltype = 92
gmsh_boutype = 36

##Read file data
hoprfile = h5py.File(hoprname, 'r')
# Read attributes
group    = hoprfile["/"]
attrs    = group.attrs
nnod = np.array(attrs['nNodes'])[0]
print('Number of nodes: %i' % (nnod))
nnodU = np.array(attrs['nUniqueNodes'])[0]
print('Number of unique nodes: %i' % (nnodU))
nel = np.array(attrs['nElems'])[0]
print('Number of elements: %i' % (nel))
nbcs = np.array(attrs['nBCs'])[0]
print('Number of boundary zones: %i' %(nbcs))
nnodxel = 64    #All elements are HEX64 -> ElemInfo has to be read to particularize for element type
nsidxel = 6     #All elements are HEX64 -> ElemInfo has to be read to particularize for element type
nnodxelbou = 16 #All boundary elemets are QUA16
order = int(attrs['Ngeo'])
sidemap = hopr_hexa_sides(order)
nsid = nsidxel*nel
# Read data to construct boundaries and periodicity information
bcNames  = np.array(hoprfile['BCNames'])
bcNames  = [elem.decode('utf-8') for elem in bcNames]
bcCodes  = np.array(hoprfile['BCType'])[:,0]
perCodes = np.array(hoprfile['BCType'])[:,3]
bcIDs    = np.array(hoprfile['SideInfo'][:,4])
# Read data to construct connectivity
globalID = np.array(hoprfile['GlobalNodeIDs'])
myglobal = np.arange(nnod)
indexsID = np.unique(globalID,return_index=True)[1]
uniqueID = np.array([globalID[index] for index in sorted(indexsID)])
for nod in uniqueID:
    indices = np.where(globalID == nod)[0]
    if len(indices) > 1:
        for inod in indices:
            myglobal[inod] = indices[0]
# Read coordinates
coord = np.array(hoprfile['NodeCoords'][np.unique(myglobal),:])
nnodU  = coord.shape[0]
#Close HOPR file
hoprfile.close()

##Construct connectivity matrix
next_num = 0
dict_nums = {}
for i in range(myglobal.size):
    if myglobal[i] not in dict_nums:
        dict_nums[myglobal[i]] = next_num
        globalID[i] = next_num
        next_num += 1
    else:
        globalID[i] = dict_nums[myglobal[i]]
lnods = globalID.reshape(nel,nnodxel,order='C')

##Build boundary elements
bcmask  = bcIDs > 0
bcsides = np.array(np.where(bcmask), dtype=np.int32)[0,:]
codbo   = bcIDs[bcmask]
nelbou  = codbo.shape[0]
bcelems = np.array(bcsides/nsidxel, dtype=np.int32)
locsid  = np.array(np.mod(bcsides,nsidxel), dtype=np.int32)
lnodb   = np.zeros((nelbou,nnodxelbou), dtype=np.int32)
for iside, side in enumerate(locsid):
    volel    = bcelems[iside]
    lnodb[iside,:] = lnods[volel, sidemap[side]]

##Build periodic elements
'''
perIDs   = np.where(bcCodes == 1)
permask  = np.in1d(bcIDs,perIDs)
persides = np.array(np.where(permask), dtype=np.int32)[0,:]
codper   = bcIDs[permask]
nelper   = codper.shape[0]
print(nelper)
perelems = np.array(persides/nsidxel, dtype=np.int32)
locsid   = np.array(np.mod(persides,nsidxel), dtype=np.int32)
lnodp    = np.zeros((nelper,nnodxelbou), dtype=np.int32)
for iside, side in enumerate(locsid):
    volel    = perelems[iside]
    lnodp[iside,:] = lnods[volel, sidemap[side]]
'''

##Write the output
lnods = lnods[:,hex64_map]+1
lnodb = lnodb[:,qua16_map]+1
#lnodp = lnodp[:,qua16_map]+1
ntot  = nel + nelbou #+ nelper

#Write the output in SOD2D format
# Open HDF5 file
'''
sodname      = 'output_sod.h5'
sodfile      = h5py.File(sodname,'w')
dims_group   = sodfile.create_group('dims')
elems_dset   = dims_group.create_dataset('numElements',(1,),dtype='i8',data=nel)
bound_dset   = dims_group.create_dataset('numBoundaryFaces',(1,),dtype='i8',data=nelbou)
nodes_dset   = sodfile.create_dataset('coords',(nnodU,3),dtype='f8',data=coord,chunks=True,maxshape=(nnodU,3))
connec_dset  = sodfile.create_dataset('connec',(nel,lnods.shape[1]),dtype='i8',data=lnods,chunks=True,maxshape=(None,lnods.shape[1]))
bounds_dset  = sodfile.create_dataset('boundFaces',(nelbou,lnodb.shape[1]),dtype='i8',data=lnodb,chunks=True,maxshape=(None,lnodb.shape[1]))
boundId_dset = sodfile.create_dataset('boundFacesId',(nelbou,),dtype='i8',data=codbo,chunks=True,maxshape=(None))
#per_dset    = h5file.create_dataset('periodicFaces',(nel_periodic,lnodp_ndim),dtype='i8',data=lnodp,chunks=True,maxshape=(None,lnodp_ndim))
#per_dset    = dims_group.create_dataset('numPeriodicFaces',(1,),dtype='i8',data=nel_periodic)
sodfile.close()
'''

#Write to gmsh in case the user wants to visualize the mesh
if visugmsh:
    with open('mesh_output.msh', 'w') as f:
        # Write the header information
        f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')

        # Write the physical groups
        f.write('PhysicalNames\n')
        f.write(f'{nbcs+1}\n')
        for i in range(nbcs):
            f.write(f'{2} {i+2} "{bcNames[i]}"\n')
        f.write(f'{3} {1} "{"fluid"}"\n')
        f.write('$EndPhysicalNames\n')

        # Write the node coordinates
        f.write('$Nodes\n')
        f.write(f'{nnodU}\n')
        for i in range(nnodU):
            f.write(f'{i+1} {coord[i][0]} {coord[i][1]} {coord[i][2]}\n')
        f.write('$EndNodes\n')

        # Write the element connectivity
        f.write('$Elements\n')
        f.write(f'{ntot}\n')
        iel = 0
        for i in range(nelbou):
            iel = iel + 1
            f.write(f'{iel} {gmsh_boutype} 2 {codbo[i]} {codbo[i]} ')
            f.write(' '.join([str(x) for x in lnodb[i,:]]))
            f.write('\n')
        for i in range(nel):
            iel = iel + 1
            f.write(f'{iel} {gmsh_voltype} 2 1 1 ')
            f.write(' '.join([str(x) for x in lnods[i,:]]))
            f.write('\n')
        f.write('$EndElements\n')