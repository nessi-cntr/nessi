# ----------------------------------------------------------------------

""" Helper functions for post processing for cntr hdf5-file based output
 data (Green's functions and observables).

Author: Hugo Strand, hugo.strand@gmail.com (2015)
Author of trunc routines: Christopher Stahl, christopher.stahl@fau.de (2022)

"""

# ----------------------------------------------------------------------

import h5py
import numpy as np

# ----------------------------------------------------------------------

from ReadCNTR import Dummy

# ----------------------------------------------------------------------
def read_group(group,key_lim=None):
    if key_lim==None:
        klist=list(group.items())
    else:
        klist=[]
        for key in key_lim:
            index=list(group.keys()).index(key)
            klist.append(list(group.items())[index])
    
    d = Dummy()
    for key, item in klist:

        #print key, type(item)
        
        if isgroup(item):

            if isgreensfunction(item):
                setattr(d, key, read_gf_group(item))
            else:
                setattr(d, key, read_group(item))

        elif isdataset(item):
            setattr(d, key, read_ndarray_from_data_set(item))

        else:
            print(key, type(item))
            raise NotImplementedError

    return d

def read_group_trunc(group,key_lim=None):
    if key_lim==None:
        klist=list(group.items())
    else:
        klist=[]
        for key in key_lim:
            index=list(group.keys()).index(key)
            klist.append(list(group.items())[index])
    
    d = Dummy()
    for key, item in klist:

        #print key, type(item)
        
        if isgroup(item):

            if isgreensfunction_trunc(item):
                setattr(d, key, read_gf_group_trunc(item))
            else:
                setattr(d, key, read_group(item))

        elif isdataset(item):
            setattr(d, key, read_ndarray_from_data_set(item))

        else:
            print(key, type(item))
            raise NotImplementedError

    return d


def read_group_slices(filename):

    fd = h5py.File(filename)

    imp = Dummy()
    for key, item in fd.items():
        if isgreensfunction_slices(item):
            setattr(imp, key, read_gf_slices_group(item))
        elif isgreensfunction(item):
            setattr(imp, key, read_gf_group(item))
        elif isgreensfunction_tavtrel(item):
            setattr(imp, key, read_gf_tavrel_group(item))
        else:
            setattr(imp, key, read_all_sets(item))
    fd.close()
    return imp

# ----------------------------------------------------------------------
def isdataset(group):
    return isinstance(group, h5py.Dataset)

# ----------------------------------------------------------------------
def isgroup(group):
    return isinstance(group, h5py.Group)

# ----------------------------------------------------------------------
def read_imp_h5file(filename,key_lim=None):

    fd = h5py.File(filename)
    imp = read_group(fd,key_lim)
    fd.close()

    return imp

def read_imp_trunc_h5file(filename,key_lim=None):

    fd = h5py.File(filename)
    imp = read_group_trunc(fd,key_lim)
    fd.close()

    return imp

# ----------------------------------------------------------------------
def read_imp_h5file_slices(filename,key_lim=None):

    fd = h5py.File(filename)
    if key_lim==None:
        klist=list(fd.items())
    else:
        klist=[]
        for key in key_lim:
            index=list(fd.keys()).index(key)
            klist.append(list(fd.items())[index])
    imp = Dummy()
    for key, item in klist:
        if isgreensfunction_slices(item):
            setattr(imp, key, read_gf_slices_group(item))
        elif isgreensfunction(item):
            setattr(imp, key, read_gf_group(item))
        elif isgreensfunction_tavtrel(item):
            continue
#setattr(imp, key, read_gf_tavrel_group(item))
        else:
            setattr(imp, key, read_all_sets(item))
    fd.close()
    return imp

# ----------------------------------------------------------------------
def read_imp_h5file_old(filename):

    fd = h5py.File(filename)

    imp = Dummy()
    for key, item in fd.items():
        
        if isgreensfunction(item):
            setattr(imp, key, read_gf_group(item))
        else:
            setattr(imp, key, read_all_sets(item))
    fd.close()
    return imp

# ----------------------------------------------------------------------
def read_gf_h5file(filename):

    fd = h5py.File(filename)

    if len(fd.keys()) > 1:
        raise NotImplementedError
    
    key = fd.keys()[0]
    group = fd[key]
    G = read_gf_group(group)
    fd.close() # done with hdf5 file
    
    return G
    
# ----------------------------------------------------------------------
def read_all_sets(group):

    # -- Read in data to dummy class
    d = Dummy()
    for key, item in group.items():        
        np_data = read_ndarray_from_data_set(group[key])
        setattr(d, key, np_data)

    return d

# ----------------------------------------------------------------------
def read_ndarray_from_data_set(hdf5_data):

    np_data = np.empty(hdf5_data.shape, dtype=hdf5_data.dtype)
    if( not hdf5_data.shape == (0,) ):
        hdf5_data.read_direct(np_data)
    return np_data
    
# ----------------------------------------------------------------------
def isgreensfunction_slices(group):

    attribs = set(('t-1', 't0'))
    return attribs.issubset(set(group.keys()))

def isgreensfunction_tavtrel(group):
    attribs = set(('gtr', 'les'))
    return attribs.issubset(set(group.keys()))

def isgreensfunction(group):

    # -- Sanity check for known Gf attributes
    attribs = set(('size1', 'size2', 'element_size', 'sig', 'les', 'mat'))
    return attribs.issubset(set(group.keys()))

def isgreensfunction_trunc(group):

    # -- Sanity check for known Gf attributes
    attribs = set(('size1', 'size2', 'element_size', 'sig', 'les'))
    return attribs.issubset(set(group.keys()))

# ----------------------------------------------------------------------
def isstrippedgreensfunction(group):

    if 'stripped' in group.keys() and group['stripped'].value == True:
        return True
    else:
        return False

# ----------------------------------------------------------------------
def read_gf_group(group):

    if not isgreensfunction(group):
        raise NotImplementedError
    
    #print group.keys()
    d = read_all_sets(group)

    # -- Don't reshape stripped Green's function
    if isstrippedgreensfunction(group): return d
    
    # -- Reshape the Green's function components
    
    # Recast to matrix form in time
    d.tv = d.tv.reshape(d.nt[0]+1, d.ntau[0]+1, d.size1[0], d.size2[0])

    # Recast vectors of trangular (upper/lower) to full matrices
    key_idx = [
        # NB! The np.triu_indices is row major, G_less is column major
        #('les', np.triu_indices(nt+1)),
        ('les', triu_indices_colmaj(d.nt+1)),
        ('ret', np.tril_indices(d.nt+1)),
        ]

    for key, idx in key_idx:
        data = getattr(d, key)

        new_data = np.zeros(
            (d.nt[0]+1, d.nt[0]+1, d.size1[0], d.size2[0]),
            dtype=complex)
        
        new_data[tuple(idx)] = data
        setattr(d, key, new_data)

    return d

def read_gf_group_trunc(group):

    if not isgreensfunction_trunc(group):
        raise NotImplementedError
    
    #print group.keys()
    d = read_all_sets(group)

    # -- Don't reshape stripped Green's function
    if isstrippedgreensfunction(group): return d
    
    # -- Reshape the Green's function components
    
    # Recast vectors of triangular (upper/lower) to full matrices
    key_idx = [
        ('les',trunc_idx(d.tc+1) ),
        ('ret',trunc_idx(d.tc+1) ),
        ]

    for key, idx in key_idx:
        data = getattr(d, key)

        new_data = np.zeros(
            (d.tc[0]+1, d.tc[0]+1, d.size1[0], d.size2[0]),
            dtype=np.complex)
        
        new_data[tuple(idx)] = data
        setattr(d, key, new_data)

    return d



# ----------------------------------------------------------------------
def read_gf_tavrel_group(group):

    if not isgreensfunction_tavtrel(group):
        raise NotImplementedError
    
    dout=Dummy()
    for key, item in group.items():
        length=len(item.items())
        maximum=max(max(p[1].shape for p in item.iteritems()))
        size1=item.items()[0][1].shape[1]
        size2=item.items()[0][1].shape[2]
        new_data = np.zeros((length,maximum,size1,size2), dtype=complex)
        t=0
        for skey, sitem in item.items():
            data = sitem[:][:][:]
            new_data[t][0:len(data)] = data
            t+=1
        setattr(dout, key, new_data)

    return dout

# ----------------------------------------------------------------------
def read_gf_slices_group(group):

    if not isgreensfunction_slices(group):
        raise NotImplementedError
    
    dall = []
    dmat = []
    for key, item in group.items():
        if key=="t-1":
            d=read_all_sets(item)
            dmat.append(d)
        else:
            d=read_all_sets(item)
            dall.append(d)
    # -- Reshape the Green's function components
    dout=Dummy()
    data = getattr(dmat[0], 'mat')
    setattr(dout, 'mat', data[:,0,0])

    key = ['les','ret']
    maximum=max(p.tstp for p in dall)
    
    
    for key in key:
        new_data = np.zeros((len(dall),maximum[0]+1, dall[len(dall)-1].size1[0], dall[len(dall)-1].size2[0]), dtype=complex)
        tstep = np.array([],dtype=int)
        for t in np.arange(len(dall)):
            data = getattr(dall[t], key)
            tstep=np.append(tstep,int(len(data)))
            new_data[t][0:len(data)] = data
        setattr(dout, key, new_data)
        setattr(dout, 'steps', tstep)

    return dout




# ----------------------------------------------------------------------
def triu_indices_colmaj(N):

    indices = [[], []]
    for idx1 in range(N[0]):
        for idx0 in range(0, idx1+1):
            indices[0].append(idx0)
            indices[1].append(idx1)

    return indices


def trunc_idx(N):
    indices =[[], []]
    for idx1 in range(N[0]):
        for idx0 in range(N[0]):
            indices[0].append(idx1)
            indices[1].append(idx0)

    return indices

# ----------------------------------------------------------------------
def strip_gf_group(group):

    g = read_gf_group(group)

    for key in group.keys():
        print(key)

    for a in ['tv', 'ret', 'les']:
        print(a, group[a].shape)

    if isstrippedgreensfunction(group):
        print('--> Trying to strip already stripped gf')
        return

    for key in ['tv', 'ret', 'les']:
        del group[key]

    group.create_dataset('tv', data=np.array([]))
    group.create_dataset('ret', data=g.ret[:, 0, :, :])
    group.create_dataset('les', data=g.les[0, :, :, :])

    group['stripped'] = True

# ----------------------------------------------------------------------
def strip_imp_h5file(filename):

    fd = h5py.File(filename, 'r+')
    
    print(dir(fd))
    
    for key in fd.keys():
        print(key)
        
    for key, item in fd.items():
        print(key, item)
        if isgreensfunction(item):
            print('--> got gf:', key)
            strip_gf_group(item)    

    fd.close()
