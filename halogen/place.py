import os
from os.path import join
import numpy as np

from numpy.ctypeslib import ndpointer
import ctypes
LOCATION = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.cdll.LoadLibrary(join(LOCATION, 'place_halos.so'))

def place_halos(halomasses, dm_pos, boxsize, omegam, ncells,
                mp, alpha, mcuts, seed=-1, exclusion=True,
                verbose=True, alg="stat"):
    """
    A wrapper for the c-code which places halos spatially
    
    Parameters
    ----------
    halomasses : float array
        An array containing the mass of halos in the simulation. These are 
        ordered by cell number, and within this are ordered by mass (descending).
        
    dm_pos : (N_dm,3)-array, dtype=float32
        The positions of the dark matter particles. Used to get initial positions
        for the halos
        
    boxsize : float
        Size of the simulation box in Mpc/h
        
    Returns
    -------
    pos : (N_halos,3)-array, dtype=float32
        Final positions of the halos
        
    """
#     if nend is None:
#         nend = len(halomasses) - 1

    if alg == 'stat':
        return _place_halos_weighted(halomasses, dm_pos, ncells, boxsize,
                                     mp, alpha, mcuts, len(alpha), omegam,
                                     seed)

#         return halopos, r

    elif alg == "rank":
        return _place_halos_rank(halomasses, dm_pos, ncells, boxsize,
                                 mp, omegam)

def r200(m, omegam):
    return (3 * m / (4 * np.pi * 200 * 27.755 * 10 ** 10)) ** (1. / 3.)

def _place_halos_weighted(halomasses, dm_pos, ncells, boxsize, mp,
                          alpha, mcuts, seed, omegam):

    cplace_halos_mglobal = lib.place_halos  # _Mglobal
    cplace_halos_mglobal.restype = None
    cplace_halos_mglobal.argtypes = [ctypes.c_long, ndpointer(ctypes.c_float),
                                     ctypes.c_long, ctypes.c_long,
                                     ndpointer(ctypes.c_float),
                                     ndpointer(ctypes.c_float), ndpointer(ctypes.c_float),
                                     ctypes.c_float, ctypes.c_float,
                                     # ctypes.c_float, ctypes.c_float,
                                     ndpointer(ctypes.c_double), ndpointer(ctypes.c_double),
                                     ctypes.c_long, ctypes.c_long,
                                     ndpointer(ctypes.c_float), ndpointer(ctypes.c_float),
                                     ndpointer(ctypes.c_float)]
#     if halopos is None:
    x = np.zeros(len(halomasses)).astype('float32')
    y = np.zeros(len(halomasses)).astype('float32')
    z = np.zeros(len(halomasses)).astype('float32')
#     else:
#     x = halopos[:, 0].copy()
#     y = halopos[:, 1].copy()
#     z = halopos[:, 2].copy()

#     if grid is None:
#         grid = np.zeros(ncells ** 3).astype('float32')

#     alpha = alpha.astype('float32')
#     mcuts = mcuts.astype('float32')
    halomasses = halomasses.astype('float32')

    dmx = dm_pos[:, 0].copy()
    dmy = dm_pos[:, 1].copy()
    dmz = dm_pos[:, 2].copy()

    cplace_halos_mglobal(len(halomasses), halomasses, ncells,
                         dm_pos.shape[0], dmx, dmy, dmz, boxsize,
                         np.float32(mp),
                         alpha, mcuts, len(alpha), seed, x, y, z)

    return np.vstack((x, y, z)).T, r200(halomasses, omegam)  # , grid

def _place_halos_rank(halomasses, dm_pos, ncells, boxsize, frac_in_halos, mp, M_min,
                        omegam):

    cplace_halos_ranked = lib.place_halos_ranked
    cplace_halos_ranked.restype = None
    cplace_halos_ranked.argtypes = [ctypes.c_long, ndpointer(ctypes.c_float),
                                    ctypes.c_long, ctypes.c_long,
                                    ndpointer(ctypes.c_float),
                                    ndpointer(ctypes.c_float), ndpointer(ctypes.c_float),
                                    ctypes.c_float, ctypes.c_float, ctypes.c_float,
                                    ctypes.c_float, ndpointer(ctypes.c_float),
                                    ndpointer(ctypes.c_float), ndpointer(ctypes.c_float)]

    x = np.zeros(len(halomasses)).astype('float32')
    y = np.zeros(len(halomasses)).astype('float32')
    z = np.zeros(len(halomasses)).astype('float32')
    halomasses = halomasses.astype('float32')

    dmx = dm_pos[:, 0].copy()
    dmy = dm_pos[:, 1].copy()
    dmz = dm_pos[:, 2].copy()

    cplace_halos_ranked(len(halomasses), halomasses, ncells,
                         dm_pos.shape[0], dmx, dmy, dmz, boxsize,
                         np.float32(mp), np.float32(frac_in_halos),
                         np.float32(M_min), x, y, z)

    return np.vstack((x, y, z)).T, r200(halomasses, omegam)
