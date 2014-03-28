'''
Created on 25/02/2014

@author: Steven
'''
import numpy as _np
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
import os
from os.path import join

def _nanmean(vec):
    x = vec[_np.logical_and(_np.logical_not(_np.isnan(vec)),
                            _np.logical_not(_np.isinf(vec)))]
    return _np.mean(x)

def _nanstd(vec):
    x = vec[_np.logical_and(_np.logical_not(_np.isnan(vec)),
                            _np.logical_not(_np.isinf(vec)))]
    return _np.std(x)

def test_overlap(rad, pos):
    k = 0
    for i in range(len(rad[1:])):
        p = _np.sqrt(_np.sum((pos[:i + 1] - pos[i + 1]) ** 2, axis=1))
        k += _np.sum((p - rad[:i + 1]) < rad[i + 1])

    return k

def test_mf_true(sample, true):
    """
    Test how close the sample mass function is to the predicted one
    
    TODO: Test closeness to n-body mf
    """
    spl = _spline(_np.log(true[:, 0]), _np.log(true[:, 1]))
    true_spl = _np.exp(spl(sample[:, 0]))
    err = _np.sqrt(_nanmean(_np.square((true_spl - sample[:, 1]) / true_spl)))
    return true_spl, err

def test_mf_nbody(sample, nbody):
    err = _np.sqrt(_nanmean(_np.square((sample - nbody) / nbody)))
    return err

def test_mass_cons(pos, masses, dm_pos, true_n, mpart, boxsize, frac_in_bounds,
                   nmax=100):
    """
    Test how mass is conserved within cells of differing size
    """

    ratio = []
    mean = []
    std = []
    for i in range(1, nmax + 1):
        edges = _np.linspace(0, boxsize, i + 1)
        edges_array = _np.vstack((edges, edges, edges))
        h = _np.histogramdd(dm_pos, bins=edges_array)[0] * mpart
        h *= frac_in_bounds
        hhalo = _np.histogramdd(pos, bins=edges_array, weights=masses)[0]

        ratio.append(_np.log(hhalo / h))

        mean.append(_nanmean(ratio[i - 1]))
        std.append(_nanstd(ratio[i - 1]))

    return ratio, mean, std


def tpcorr(data, bins, vol):
    """
    Really simple/dodge 2pcf for periodic boxes
    """
    hist = _np.zeros(len(bins) - 1)
    for i in range(data.shape[0]):
        rad = _np.sqrt(_np.sum((data - data[i]) ** 2, axis=1))
        hist += _np.histogram(rad, bins=bins)[0]

    binvol = 4 * _np.pi * (bins[1:] ** 3 - bins[:-1] ** 3) / 3
    print hist
    print binvol

    return hist * vol / (data.shape[0] * binvol)

def mk_input(infile, numlines, inform, output, boxsize):
    input_str = """
#input information
data_filename= %s
num_lines= %s
input_format= %s
output_filename= %s
box_size= %s

#tree stuff
use_tree= 0
max_tree_order= 7
max_tree_nparts= 100

#irrelevant yet
use_pm= 0
n_grid_side= 256
resolution_factor= 6
""" % (infile, numlines, inform, output, boxsize)

    return input_str

def run_2pcorr(cute_makefile, data_file, boxsize, nr=24, rmax=None,
               n_logint=8, numlines=1000, inform=0, remake=True,
               tmp_prefix=""):

    cute_dir = os.path.dirname(cute_makefile)
    os.chdir(cute_dir)

    if not rmax or rmax > boxsize:
        rmax = boxsize / 2

    if remake:
        with open(cute_makefile) as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if line.startswith("NB_R"):
                    lines[i] = "NB_R = %s\n" % nr
                elif line.startswith("R_MAX"):
                    lines[i] = "R_MAX = %s\n" % rmax
                elif line.startswith("N_LOGINT"):
                    lines[i] = "N_LOGINT = %s\n" % n_logint
        with open(cute_makefile, 'w') as f:
            for line in lines:
                f.write(line)

        os.system("make clean")
        os.system("make CUTE_box")

    # Create the data tables we need.
    with open(join(tmp_prefix, "tpcinput.in"), 'w') as f:
        f.write(mk_input(data_file, numlines, inform,
                         join(tmp_prefix, "tpc"), boxsize))

    os.system("./CUTE_box %s" % join(tmp_prefix, "tpcinput.in"))

    s2pcf = _np.genfromtxt(join(tmp_prefix, "tpc"))

    print s2pcf
    # os.system("rm %s" % join(tmp_prefix, "tpc*"))
    return s2pcf[:, :2]

def test_2pcorr(nbodyfile, samplefile, makefile, boxsize, nr=24, rmax=None,
                 n_logint=8, numlines=1000, nbody_corr=[], remake=True,
                 tmp_prefix=""):
    """
    Test closeness of the 2PCF
    """
    t2pcf = nbody_corr
    print "tmp_prefix: ", tmp_prefix
    if not os.path.exists(tmp_prefix): os.mkdir(tmp_prefix)

    s2pcf = run_2pcorr(makefile, samplefile, boxsize, nr, rmax,
                       n_logint, numlines, 0, remake=remake, tmp_prefix=tmp_prefix)

    print "t2pcf: ", t2pcf
    if len(t2pcf) == 0:
        t2pcf = run_2pcorr(makefile, nbodyfile, boxsize, nr, rmax,
                           n_logint, numlines, 0, remake=False)

    return s2pcf, t2pcf, _np.sqrt(_np.mean(_np.square(t2pcf[:, 1] - s2pcf[:, 1])))

