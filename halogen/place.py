import os
from os.path import join
import numpy as np

from numpy.ctypeslib import ndpointer
import ctypes
import sys
#===============================================================================
# SET UP INTERFACE
#===============================================================================
LOCATION = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.cdll.LoadLibrary(join(LOCATION, 'place_halos.so'))
cplace_halos = lib.place_halos_byparts  # _Mglobal
cplace_halos.restype = None
cplace_halos.argtypes = [ctypes.c_long, ctypes.c_long, ndpointer(ctypes.c_float),
                         ctypes.c_long, ctypes.c_long,
                         ndpointer(ctypes.c_float),
                         ndpointer(ctypes.c_float), ndpointer(ctypes.c_float),
                         ctypes.c_float, ctypes.c_float,
                         ctypes.c_long, ctypes.c_float,
                         ndpointer(ctypes.c_double), ndpointer(ctypes.c_double),
                         ctypes.c_long, ndpointer(ctypes.c_float),
                         ndpointer(ctypes.c_float), ndpointer(ctypes.c_float),
                         ndpointer(ctypes.c_float), ndpointer(ctypes.c_double)]

#===============================================================================
# PLACEMENT CLASS
#===============================================================================
class HaloPlacer(object):

    def __init__(self, halomasses, dm_pos, L, ncells,
                 mp, alpha, mcuts, seed=-1,
                 rho_ref=2.7755e11, keep_lastiter=False):

        self._halomasses = halomasses.astype('float32')
        self._L = np.float32(L)
        self._ncells = ncells
        self._mp = np.float32(mp)
        self._alpha = alpha
        self._mcuts = mcuts
        self._seed = seed
        self._rho_ref = np.float32(rho_ref)
        self._ndm = dm_pos.shape[0]
        self._keep_last_iter = keep_lastiter
        self._dmx = dm_pos[:, 0].copy()
        self._dmy = dm_pos[:, 1].copy()
        self._dmz = dm_pos[:, 2].copy()

        # # Save outputs
        self._x = np.empty(len(halomasses)).astype('float32')
        self._y = np.empty(len(halomasses)).astype('float32')
        self._z = np.empty(len(halomasses)).astype('float32')
        self._r = np.empty(len(halomasses)).astype('float32')

        self._massleft = np.empty(ncells ** 3)

        # The total number of halos placed so far.
        self.nhalos_placed = 0

        if keep_lastiter:
            self._old_massleft = np.empty(ncells ** 3)
            self._old_nhalos_placed = -1

    def place(self, nstart=0, nend=None):
        if nend is None:
            nend = len(self._halomasses)

        print "NSTART: ", nstart
        print "NEND: ", nend
        print "HALOS PLACED: ", self.nhalos_placed
        print "OLD HALOS PLACED: ", self._old_nhalos_placed
        print "XSIZE: ", self._x.nbytes
        print "YSIZE: ", self._y.nbytes
        print "ZSIZE: ", self._z.nbytes
        print "rSIZE: ", self._r.nbytes
        print "HALOMASSES size: ", self._halomasses.nbytes
        print "DM x,y,x size: ", self._dmx.nbytes, self._dmy.nbytes, self._dmz.nbytes
        print "massleft, old size: ", self._massleft.nbytes, self._old_massleft.nbytes
        print "alpha size: ", self._alpha.nbtyes
        print "mcuts size: ", self._mcuts.nbytes
        if nstart != self.nhalos_placed:
            if self._keep_last_iter:
                if nstart != self._old_nhalos_placed:
                    raise ValueError("nstart must be the same as the last halo placed.")
            else:
                raise ValueError("nstart must be the same as the last halo placed.")

        if not self._keep_last_iter:
            # Increment the total halos placed.
            self.nhalos_placed = nend
        else:
            if nstart == self._old_nhalos_placed:
                self._massleft = self._old_massleft
            else:
                self._old_nhalos_placed = self.nhalos_placed
                self.nhalos_placed = nend

        cplace_halos(nstart, nend, self._halomasses, self._ncells,
                     self._ndm, self._dmx, self._dmy, self._dmz, self._L,
                     self._rho_ref, self._seed,
                     self._mp, self._alpha, self._mcuts, len(self._alpha),
                     self._x, self._y, self._z, self._r, self._massleft)

    @property
    def halopos(self):
        return np.vstack((self._x[:self.nhalos_placed],
                          self._y[:self.nhalos_placed],
                          self._z[:self.nhalos_placed])).T

    @property
    def r(self):
        return self._r[np.logical_not(np.isnan(self._r))]

def place_halos(halomasses, dm_pos, boxsize, ncells,
                mp, alpha, mcuts, seed=-1,
                nstart=0, nend=None, rho_ref=2.7755e11,
                halopos=None, massleft=None, r=None):
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
    if nend is None:
        nend = len(halomasses) - 1


    if halopos is None:
        x = np.zeros(len(halomasses)).astype('float32')
        y = np.zeros(len(halomasses)).astype('float32')
        z = np.zeros(len(halomasses)).astype('float32')
    else:  # Some have been filled already
        x = halopos[:, 0].copy()
        y = halopos[:, 1].copy()
        z = halopos[:, 2].copy()

    if massleft is None:
        massleft = np.empty(ncells ** 3)
    if r is None:
        r = np.empty(len(halomasses)).astype('float32')

    halomasses = halomasses.astype('float32')



    cplace_halos(nstart, nend, halomasses, ncells,
                 dm_pos.shape[0], dmx, dmy, dmz, boxsize, rho_ref, seed,
                 np.float32(mp), alpha, mcuts, len(alpha), x, y, z, r, massleft)

    return np.vstack((x, y, z)).T, r, massleft
