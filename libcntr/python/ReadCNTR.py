# ----------------------------------------------------------------------

""" Helper functions for post processing and input file generation
for cntr text-file based output data (Green's functions and observables).

Author: Hugo Strand, hugo.strand@gmail.com (2015) 

"""

# ----------------------------------------------------------------------

import os
import ast
import glob
import inspect
import numpy as np

# ----------------------------------------------------------------------
class Dummy(object):
    """ Dummy class used for storing data as members. """
    pass

# ----------------------------------------------------------------------
def write_input_file(d, filename):

    """ Generate an input file from dictionary or class

    Input: d -- dictionary with parameters (Nb! without underscores)
           filename -- input filename
    
    Returns: (writes the file to disk)

    """

    if type(d) == dict:
        pass
    elif hasattr(d, '__class__'):
        d = d.__dict__
    else:
        raise NotImplementedError
    
    f = open(filename, 'w')
    for key, value in d.items():
        f.write('__'+key+'=  ' + str(value) + '\n')
    f.close()

# ----------------------------------------------------------------------
def read_input_file(filename):

    """ Read cntr-type input file and return class with the parameters
    as members. """
    
    d = {}
    f = open(filename, 'r')
    for line in f.readlines():
        key, value = line.split('=')
        key = key.split('__')[1]

        if('txt' in value):
            d[key] = value
        elif('true' in value or 'false' in value):
            d[key] = bool(value)
        else:
            d[key] = ast.literal_eval(value.strip(' '))

    f.close()

    out = Dummy()
    out.__dict__.update(d)
    return out

# ----------------------------------------------------------------------
def read_obs_out(filename):

    """ General observables out file parser can take any format based on
    label1: 1.234 label2: 1245.0 2311.0 1231.3 label3: 1111.1 ...
    
    Returns: class instance with member ndarrays with the label names.
    """
    
    obs = np.genfromtxt(filename, dtype=None)
    d = Dummy()
    tmp = obs[0] # First row (is a tuple)

    # -- Get out string occurences
    names = []
    for el in tmp:
        if type(el) is np.string_:
            names.append(el.strip(':'))

    #print '--> Data:', names
    #print tmp

    # -- Locate the indices for the "names" labels in idx_list
    idx_list = []
    for name in names:
        mask = np.array([x == name+':' for x in tmp])
        idx = np.where(mask == True)[0][0]
        idx_list.append(idx)

    idx_list = np.array(idx_list)
    names = np.array(names)

    # -- Sort the "names" according to occurence in idx_list
    sidx = np.argsort(idx_list)
    idx_list, names = idx_list[sidx], names[sidx]

    # -- Pick out data for each "name" out of obs and put in dict "d"
    for idx, name in enumerate(names):
        #print '-'*72
        #print '--> Parsing: ', name

        s_idx = idx_list[idx]
        if idx+1 == len(idx_list):
            e_idx = len(tmp)
            ndata = len(tmp) - s_idx - 1
        else:
            e_idx = idx_list[idx+1]
            ndata = e_idx - s_idx - 1

        nrow = obs.shape[0]
        #print 'nrow =', nrow
        #print 'ndata =', ndata
        #print s_idx, e_idx
        
        data = np.zeros((nrow, ndata))
        for ridx in xrange(nrow):
            for cidx, didx in enumerate(xrange(s_idx+1, e_idx)):
                #print ridx, cidx, didx, obs[ridx], obs[ridx][didx]
                data[ridx, cidx] = obs[ridx][didx]

        setattr(d, name, data)

    #print d.__dict__

    # -- Cut out zero-time data
    for key, value in d.__dict__.items():
        setattr(d, key+'0', value[0])
        setattr(d, key, value[1:])

    #print d.__dict__.keys()
    return d

# ----------------------------------------------------------------------
def read_gf(filename):

    """ Read cntr-type text based Green's function file and return
    class instance with the components as members, i.e.

    > G = read_gf('G.out')
    > print type(G.mat)
    numpy.ndarray
    > print G.mat.shape
    (150, 2, 2)
    
    """

    d = Dummy()
    f = open(filename, 'r')

    # -- Dimensions in first row
    row = f.readline()
    row = row.strip('#')
    d.nt, d.ntau, d.dim, d.foo = np.fromstring(row, dtype=int, sep=' ')

    #for key, value in d.__dict__.items():
    #    print key, ' = ', value

    ntau, nt, dim = d.ntau, d.nt, d.dim

    G = { 'mat' : np.zeros((ntau+1, dim, dim), dtype=np.complex),
          'tv'  : np.zeros((nt+1, ntau+1, dim, dim), dtype=np.complex),
          'ret' : np.zeros((nt+1, nt+1, dim, dim), dtype=np.complex),
          'les' : np.zeros((nt+1, nt+1, dim, dim), dtype=np.complex), }
    
    for line in f.readlines():
        if not ':' in line:
            continue
        
        key, str_data = line.split(':')
        data = np.fromstring(str_data, dtype=float, sep=' ')

        if key == 'mat':
            idx = tuple(np.array([data[0]], dtype=np.int))
            data = data[1:]
        else:
            idx = tuple(np.array(data[:2], dtype=np.int))
            data = data[2:]

        data = data[::2] + 1.j*data[1::2]
        G[key][idx] = data.reshape(dim, dim)

    d.__dict__.update(G)
        
    return d

# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename1 = 'test_dict.txt'
    filename2 = 'test_class.txt'
    
    d = {'foo':1234.1}
    write_input_file(d, filename1)

    d = Dummy()
    d.foo = 1234.1
    write_input_file(d, filename2)
