import numpy as np
import h5py

def hopr_hexa_sides(p):
    '''
        Based on fig 2.10 on HOPR documentation (link at the code header)
        k = 0 -> side 1
        j = 0 -> side 2
        i = p -> side 3
        j = p -> side 4
        i = 0 -> side 5
        k = p -> side 6
    '''
    nsides = 6 #Always working with hexahedral meshes
    nodes = np.zeros((nsides,(p+1)**2),dtype=np.int32)
    for iside in range(nsides):
        inod = 0
        for j in range(p+1):
            for i in range(p+1):
                if iside == 0:
                    nodes[iside, inod] = 0*(p+1)**2 + j*(p+1) + i
                if iside == 1:
                    nodes[iside, inod] = j*(p+1)**2 + 0*(p+1) + i
                if iside == 2:
                    nodes[iside, inod] = j*(p+1)**2 + i*(p+1) + p
                if iside == 3:
                    nodes[iside, inod] = j*(p+1)**2 + p*(p+1) + i
                if iside == 4:
                    nodes[iside, inod] = j*(p+1)**2 + i*(p+1) + 0
                if iside == 5:
                    nodes[iside, inod] = p*(p+1)**2 + j*(p+1) + i
                inod = inod + 1
    return nodes