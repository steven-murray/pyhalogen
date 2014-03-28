import numpy as np
from scipy.integrate import cumtrapz
from hmf import sample_mf
import matplotlib.pyplot as plt

class GetSample(object):
    def __init__(self, masses, mnum, M):
        self.mass = masses
        self.num = mnum
        self.M = M

    def nltm(self, cons):
        """ Calculate sampled MF"""
        mass = getattr(self, cons)
        nltm = np.array([np.sum(mass < m) for m in self.M])
        return nltm

    def mltm(self, cons):
        mass = getattr(self, cons)
        mltm = np.zeros(len(self.M))
        for i, m in enumerate(self.M[1:]):
            if np.sum(mass < m) > 0:
                mltm[i + 1] = np.sum(mass[mass < m])
        return mltm

    def ngtm(self, cons):
        mass = getattr(self, cons)
        ngtm = np.array([np.sum(mass > m) for m in self.M])
        return ngtm

    def mgtm(self, cons):
        mass = getattr(self, cons)
        mgtm = np.zeros(len(self.M))
        for i, m in enumerate(self.M[1:]):
            if np.sum(mass > m) > 0:
                mgtm[i + 1] = np.sum(mass[mass > m])
        return mgtm

class GetTrue(object):
    def __init__(self, hmf, boxsize):
        self.hmf = hmf
        self.boxsize = boxsize

    def nltm(self):
        return cumtrapz(hmf.dndlnm, np.log10(hmf.M), initial=1e-20) * np.log(10) * self.boxsize ** 3

    def ngtm(self):
        nltm = cumtrapz(hmf.dndlnm[::-1], dx=np.log10(hmf.M[1]) - np.log10(hmf.M[0]), initial=1e-20) * np.log(10) * self.boxsize ** 3
        return nltm[::-1]

    def mltm(self):
        return cumtrapz(hmf.M * hmf.dndlnm, dx=np.log10(hmf.M[1]) - np.log10(hmf.M[0]), initial=1e-20) * np.log(10) * self.boxsize ** 3

    def mgtm(self):
        mgtm = cumtrapz((hmf.M * hmf.dndlnm)[::-1], dx=np.log10(hmf.M[1]) - np.log10(hmf.M[0]), initial=1e-20) * np.log(10) * boxsize ** 3
        return mgtm[::-1]

class GetAx(object):
    def __init__(self, figsize=(15, 10)):
        self.num_mltm_fig = plt.figure(figsize=figsize)
        self.num_mgtm_fig = plt.figure(figsize=figsize)
        self.num_nltm_fig = plt.figure(figsize=figsize)
        self.num_ngtm_fig = plt.figure(figsize=figsize)
        self.mass_mltm_fig = plt.figure(figsize=figsize)
        self.mass_mgtm_fig = plt.figure(figsize=figsize)
        self.mass_nltm_fig = plt.figure(figsize=figsize)
        self.mass_ngtm_fig = plt.figure(figsize=figsize)

    def getfig(self, cons, quant):
        return getattr(self, "%s_%s_fig" % (cons, quant))

def _plot_and_print(sample, true, fig, M, boxsize, mfrac, pl):
    if sample[0] < sample[-1]:
        sample = sample[1:]
        true = true[1:]
        M = M[1:]
    else:
        sample = sample[:-1]
        true = true[:-1]
        M = M[:-1]

    dev = np.log(sample / true)
    dev = dev[np.logical_not(np.isnan(dev))]
    dev = dev[np.logical_not(np.isinf(dev))]
    rms = np.exp(np.sqrt(np.mean(np.square(dev)))) - 1
    print "  RMS = %s" % rms
    print "  SYS ERR: %s" % (np.exp(np.mean(dev)) - 1)

    # col = np.random.random(3)
    col = ['r', 'b', 'g']
    ax1 = fig.add_subplot(211)
    ax1.plot(M, true, color=col[pl])
    ax1.scatter(M, sample, color=col[pl])
    ax2 = fig.add_subplot(212)
    ax2.plot(M, sample / true, color=col[pl],
            label="L=%s, mfrac=%s, RMS=%g" % (boxsize, mfrac, rms))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log', basey=2)
    ax1.legend(loc=0)

#     ax.plot(M, true / true.max(), color=col)
#     ax.plot(M, 1 + np.sqrt(true) / true, color=col, linestyle='--')
#     ax.plot(M, 1 - np.sqrt(true) / true, color=col, linestyle='--')

if __name__ == "__main__":
    import logging
    hmflog = logging.getLogger('hmf')
    hmflog.addHandler(logging.NullHandler())
    omegab = 0.05
    omegac = 0.25

    boxsizes = [500, 1000, 2000]
    nparts = [256, 512, 1024]
    mfracs = [0.01]
    outdir = "/Users/Steven/Documents/PhD/sample_mf_plots"
    from os.path import join

    mf_kwargs = {"omegab":omegab,
                 "omegac":omegac,
                 "lAccuracyBoost":2,
                 "transfer__kmax":20
                 }
    M_max = 100

    ax = GetAx()

    pl = 0
    for boxsize, npart in zip(boxsizes, nparts):
        for mfrac in mfracs:

            print "BOXSIZE: ", boxsize, " MFRAC: ", mfrac

            simvars = {'boxsize':boxsize,
                       'npart': npart ** 3,
                       'm_min':20,
                       'mfrac':mfrac}
            mtot = boxsize ** 3 * (omegab + omegac) * 2.7755e11
            nvars = {'m':mtot,
                     'm_min':np.log10(20 * mtot / npart ** 3),
                     'm_max':np.log10(mtot * mfrac)}

            m, hmf, frac = sample_mf(simvars, tol=3, **mf_kwargs)
            mnum = sample_mf(nvars=nvars, tol=3, **mf_kwargs)[0]
            print len(m), len(mnum)
            s = GetSample(m, mnum, hmf.M)
            t = GetTrue(hmf, boxsize)
            for cons in ["mass", "num"]:
                for quantity in ["nltm", "mltm", "ngtm", "mgtm"]:
                    print "  " + "="*15 + " %s-C %s " % (cons.upper(), quantity.upper()) + "="*15
                    _plot_and_print(getattr(s, quantity)(cons), getattr(t, quantity)(),
                                    ax.getfig(cons, quantity), hmf.M, boxsize, mfrac, pl)
            pl += 1
    for cons in ["mass", "num"]:
        for quantity in ["nltm", "mltm", "ngtm", "mgtm"]:
#             x = ax.getax(cons, quantity)
#             x.set_xscale('log')
#             x.set_yscale('log', basey=2)
#             # box = x.get_position()
#             # x.set_position([box.x0, box.y0, box.width * 0.6, box.height])
#             x.legend(loc=0)  # , bbox_to_anchor=(1, 0.5))
            ax.getfig(cons, quantity).savefig(join(outdir, "%s_%s.pdf" % (cons, quantity)))





