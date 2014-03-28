'''
Created on 12/02/2014

@author: Steven

This is supposed to run the whole idea
'''

from hmf import sample_mf
import numpy as np
from place import place_halos
import os

#===============================================================================
# THE DRIVER ROUTINE
#===============================================================================
def fast(filename, ncells=1, max_mfrac=0.01,
         min_halo_mass=20, tol=3, verbose=True, exclusion=True,
         alg='stat', alpha=1.0, ** mf_kwargs):
    """
    Runs a FAST approximate synthetic galaxy simulation
    
    Parameters
    ----------
    filename : str
        Path to input simulation file
            
    ncells : int, default=1
        How many cells on a side for the placement.
        
    max_mfrac : float, optional, default None
        The log10 maximum mass of a halo to be sampled. If None, is set at 
        mass of entire simulation (or cell).
        
    min_halo_mass : float, optiona, default 20
        The minimum halo mass, in units of the particle mass.
    
    tol : int or float, default 3
        The tolerance with which to match the mass of the sim.
        If int, will be within tol*min_halo_mass of the full mass.
        If float, will be within tol*full_mass of the full mass
        
    exclusion : bool, default True
        Whether to use halo exclusion

    mf_kwargs : Any argument passed to MassFunction()
    
    Returns
    -------
    halopos : (N,3)-array 
        The x,y,z positions of the halos
        
    halomass : N-array
        The mass of each halo.
        
    """
    # Read in the dm particle positions from GADGET file
    if verbose:
        print "READING GADGET FILE...",
    dm_pos, header = _get_gadget(filename)
    mpart = header['massarr'][1] * 1e10
    boxsize = header['boxsize']
    npart = header['n_all'][1]

    # Make sure particles aren't on upper edge
    dm_pos[dm_pos == boxsize] = 0.0
    simvars = {'boxsize':boxsize,
               'npart': npart,
               'm_min':min_halo_mass,
               'mfrac':max_mfrac}

    # # Set mf_kwargs with things we know. This requires us to set omegab and omegac
    omegab_frac = 0.1666
    mf_kwargs["omegab"] = omegab_frac * header['omegam']
    mf_kwargs["omegac"] = (1 - omegab_frac) * header['omegam']
    mf_kwargs["omegav"] = header['omegav']

    halomasses, hmf, frac_in_bounds = sample_mf(simvars=simvars, tol=tol, match='mass',
                                                sort=True, **mf_kwargs)

    test = halomasses[1:] / halomasses[:-1]
    if test.max() > 1.0:
        print "halomasses weren't actually sorted!" + str(test[test > 1.0])
        halomasses = np.sort(halomasses)[::-1]

    test = halomasses[1:] / halomasses[:-1]
    print test.max()

    if verbose:
        print "  FRACTION OF MASS IN BOUNDS: ", frac_in_bounds

    if verbose:
        print "PLACING HALOS..."
    halopos, radii = place_halos(halomasses, dm_pos, boxsize, header['omegam'],
                                 ncells, mpart, frac_in_bounds, min_halo_mass,
                                 exclusion, verbose, alg, alpha)

    return halopos, halomasses, dm_pos, hmf, header, frac_in_bounds, radii

def _get_gadget(filename):
    """
    Reads a gadget binary specified by filename
    """
    positions_dtype = np.dtype([("x", 'f'),
                                ("y", 'f'),
                                ("z", 'f')
                                ])

    header_dtype = np.dtype([('npart', ('i', 6)),
                             ('massarr', ('d', 6)),
                             ('time', 'd'),
                             ('z', 'd'),
                             ('FlagSfr', 'i'),
                             ('FlagFeedback', 'i'),
                             ('n_all', ('i', 6)),
                             ('FlagCooling', 'i'),
                             ('num_files', 'i'),
                             ('boxsize', 'd'),
                             ('omegam', 'd'),
                             ('omegav', 'd'),
                             ('h', 'd'),
                             ('FlagAge', 'i'),
                             ('FlagMetals', 'i'),
                             ('n_all_HW', ('i', 6)),
                             ('flag_entr_ics', 'i')
                             ])

    # Check if its in multiple bits
    folder, actual_filename = os.path.split(filename)
    avail_files = os.listdir(folder)
    avail_files = [f for f in avail_files if f.startswith(actual_filename)]


    n_so_far = 0
    for i in range(len(avail_files)):
        if len(avail_files) > 1:
            extension = ".%s" % i
        else:
            extension = ""

        thisfile = filename + extension
        with open(thisfile, "rb") as f:
            f.read(4)
            x = str(f.read(4))
            if x == "HEAD":
                gtype = 2
            else:
                gtype = 1

        with open(thisfile, "rb") as f:
            if gtype == 2:
                f.read(16)

            # Header block
            f.read(4)
            header = np.fromfile(f, dtype=header_dtype, count=1)
            # So far only 196 out of 256 so read the rest
            f.read(256 - 196)
            f.read(4)


            # # Now we can allocate the final array of pos/id since we have n_all
            if i == 0:
                pos = np.empty(header['n_all'][0][1], dtype=positions_dtype)
                # ids = np.empty(header['n_all'][0][1])

            if gtype == 2:
                f.read(16)

            # Positions Block
            f.read(4)
            pos[n_so_far:n_so_far + header['npart'][0][1]] = np.fromfile(f, dtype=positions_dtype, count=header['npart'][0][1])
            f.read(4)

#             # SKIP Velocities Block (don't save it anywhere)
#             if gtype == 2:
#                 f.read(16)
#
#             f.seek(3 * 4 * header['npart'][0][1] + 2 * 4, 1)
#
#             # ID's block
#             if gtype == 2:
#                 f.read(16)
#             f.read(4)
#             ids[n_so_far:n_so_far + header['npart'][0][1]] = np.fromfile(f, dtype='i', count=header['npart'][0][1])

            n_so_far += header['npart'][0][1]

#     ids = np.array(ids - 1)
#
#     indices = np.argsort(ids)
#     pos = pos[indices]
#     ids = sorted(ids)

    header_dict = {}
    for name in header.dtype.names:
        header_dict[name] = header[name][0]

    return pos.view(np.float32).reshape(pos.shape + (-1,)) , header_dict


def _make_grid(dm_pos, boxsize, ncells):
    H = np.histogramdd(dm_pos, ncells)[0]
    return H

