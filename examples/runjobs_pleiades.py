#!/usr/bin/python

import subprocess
import sys
import os

import json
from ConfigParser import SafeConfigParser as cfg
from os.path import join
from argparse import ArgumentParser

#===============================================================================
# Configuration
#===============================================================================
def main():
    parser = ArgumentParser()

    # parser.add_argument("nodes", type=int, help="the number of nodes requested")
#     parser.add_argument("-q", "--queue", default='routequeue', choices=['debugq', "routequeue"], help="the queue to put it on")
    parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)")
#     parser.add_argument("-W", '--walltime', default="01:00:00")

    # Arguments for the actual run
    parser.add_argument("filename", help="Gadget snapshot to be populated")
    parser.add_argument("configfile", help="Configuration file for the run")
    args = parser.parse_args()

    #===========================================================================
    # Run the configuration file
    #===========================================================================
    config = cfg()
    config.read(args.configfile)

    queue = config.get("PBS", 'queue')
    walltime = config.get("PBS", 'walltime')

    outdir_base = config.get("IO", "outdir_base")
    outdir_anl = config.get("IO", "outdir_anl")

    tol = json.loads(config.get("RunOptions", 'tol'))
    ncells = json.loads(config.get("RunOptions", 'ncells'))
    min_halo_mass = json.loads(config.get("RunOptions", 'min_halo_mass'))
    max_mfrac = json.loads(config.get("RunOptions", 'max_mfrac'))
    alphas = json.loads(config.get("RunOptions", 'alpha'))
    nthreads = config.getint("RunOptions", "nthreads")

    density = config.getboolean("RunFlags", "density")
    exclusion = config.getboolean("RunFlags", "exclusion")
    analyse = config.getboolean("RunFlags", "analyse")
    plot = config.getboolean("RunFlags", "plot")
    post_analyse = config.getboolean("RunFlags", "post_analyse")
    alg = config.get("RunFlags", "alg")


    bins = config.getint("AnalyseOptions", "bins")
    nmax = config.getint("AnalyseOptions", 'nmax')
    rmax = config.getfloat("AnalyseOptions", "rmax")
    n_logint = config.getint("AnalyseOptions", "n_logint")
    makefile = config.get("AnalyseOptions", "makefile")
    numhalos = config.getint("AnalyseOptions", "numhalos")
    nr = config.getint("AnalyseOptions", "nr")

    mf_kwargs = config.get("MFOptions", "mf_kwargs")

    JOBNAME = "fast"

    subprocess.call(["python", "setup.py", "install"])
    subprocess.call(["rm", "-rf", "build"])
    subprocess.call(["rm", "-rf", "dist"])
    #===============================================================================
    # SEND THE JOB
    #===============================================================================
    i = 0
    jobs = []
    fulloutdir = []
    for mfrac in max_mfrac:
        for m in min_halo_mass:
            for N in ncells:
                for a in alg:
                    for alpha in alphas:

                        job_name = JOBNAME + "_" + str(i)

                        folder = os.path.split(args.filename)[0]
                        fname = os.path.split(args.filename)[1]
                        if "NB" in fname:
                            simtype = "NB"
                        elif "2LPT" in fname:
                            simtype = "2LPT"

                        folder = os.path.split(folder)[1]
#                         parts = folder.split("_")
#                         # Get boxsize, npart and redshift
#                         boxsize = '1000'
#                         npart = '512'
#                         z = '0.0'
#                         for part in parts:
#                             if "L" in part:
#                                 boxsize = part[1:]
#                             elif "N" in part:
#                                 npart = part[1:]
#                             elif 'z' in part:
#                                 z = part[1:]
                        # Get a string rep of the original input which goes through
                        outdir_part = "N" + str(N) + "_m" + str(m) + "_alpha" + str(alpha)\
                                       + "_mfrac" + str(mfrac) + "_" + a + "_" + \
                                       folder + "_" + simtype
                        if exclusion:
                            outdir_part += "_" + "ex"

                        fulloutdir.append(join(outdir_base, outdir_part))

                        if not os.path.exists(fulloutdir[i]):
                            os.system("mkdir " + fulloutdir[i])

                        # Make the config file for the runfast script
                        with open(join(fulloutdir[i], "config"), 'w') as f:
                            CFG = mkconfig(N, m, mfrac, density, exclusion,
                                           analyse, plot, bins, nmax,
                                           nr, rmax, n_logint, makefile, numhalos,
                                           a, tol, alpha, mf_kwargs)

                            CFG.write(f)

                        # # ALSO COPY IN THE WHOLE SOURCE CODE
                        os.system("cp -r pyfast %s" % join(fulloutdir[i]))
                        job_string = """
        #!/bin/bash
        
        #PBS -N %s
        #PBS -q %s
        #PBS -l select=1:ncpus=12:mem=23gb
        #PBS -l walltime=%s
        #PBS -o %s
        #PBS -e %s
        
        cd $PBS_O_WORKDIR
        source ~/.bashrc
        export OMP_NUM_THREADS=%s
        time runfast %s %s %s 
        """ % (job_name, queue, walltime, join(fulloutdir[i], "output.txt"),
               join(fulloutdir[i], "error.txt"), nthreads, args.filename, fulloutdir[i],
               join(fulloutdir[i], "config"))

                        with open(join(fulloutdir[i], "job"), "w") as f:
                            f.write(job_string)
                        print "="*20
                        print job_string
                        jobs.append(subprocess.Popen(["qsub", join(fulloutdir[i], "job")], stdout=subprocess.PIPE))
                        i += 1

    #===========================================================================
    # ANALYSIS PART
    #===========================================================================
    # Get the number for each created job
    if post_analyse:
        jobs = [job.stdout.read().replace("\n", "").replace(".epic", "") for job in jobs]
        jobs = ":".join(jobs)
        outdir = join(outdir_base, outdir_anl)
        analysis_jobstring = """
    #!/bin/bash
    
    #PBS -N analysis
    #PBS -W group_list=partner712
    #PBS -q %s
    #PBS -l select=1:ncpus=12:mem=23gb
    #PBS -l walltime=%s
    #PBS -o %s
    #PBS -e %s
    #PBS -W depend=afterok:%s
    
    cd $PBS_O_WORKDIR
    source ~/.bashrc
    time analyse %s %s %s %s --nmax %s --nr %s --rmax %s --n_logint %s --numhalos %s
    """ % (queue, walltime, join(outdir, "output.txt"),
           join(outdir, "error.txt"), jobs, " ".join(fulloutdir),
           args.filename, outdir, makefile, nmax, nr, rmax, n_logint, numhalos)

        if not os.path.exists(outdir):
            os.system("mkdir " + outdir)

        with open(join(outdir, "job_an"), "w") as f:
            f.write(analysis_jobstring)
            print "="*20
            print analysis_jobstring
        subprocess.Popen(["qsub", join(outdir, "job_an")])


def mkconfig(ncells, min_halo_mass, max_mfrac, density,
             exclusion, analyse, plot, bins, nmax,
             nr, rmax, n_logint, makefile, numhalos, alg, tol,
             alpha, mf_kwargs):
    config = cfg()

    config.add_section("RunOptions")
    config.set("RunOptions", "tol", str(tol))
    config.set("RunOptions", "ncells", str(ncells))
    config.set("RunOptions", "min_halo_mass", str(min_halo_mass))
    config.set("RunOptions", "max_mfrac", str(max_mfrac))
    config.set("RunOptions", "alpha", str(alpha))

    config.add_section("RunFlags")
    config.set("RunFlags", "density", str(density))
    config.set("RunFlags", "exclusion", str(exclusion))
    config.set("RunFlags", "analyse", str(analyse))
    config.set("RunFlags", "plot", str(plot))
    config.set("RunFlags", "alg", str(alg))

    config.add_section("AnalyseOptions")
    config.set("AnalyseOptions", "bins", str(bins))
    config.set("AnalyseOptions", "nmax", str(nmax))
    config.set("AnalyseOptions", "nr", str(nr))
    config.set("AnalyseOptions", "rmax", str(rmax))
    config.set("AnalyseOptions", "makefile", str(makefile))
    config.set("AnalyseOptions", "n_logint", str(n_logint))
    config.set("AnalyseOptions", "numhalos", str(numhalos))

    config.add_section("MFOptions")
    config.set("MFOptions", "mf_kwargs", str(mf_kwargs))

    return config
if __name__ == "__main__":
    sys.exit(main())

