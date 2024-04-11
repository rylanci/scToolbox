from itertools import repeat
from multiprocessing import Pool
from datetime import datetime
import subprocess
import os 
import logging


def run_hmmratac(ipath, opath):
    file=os.path.basename(ipath)
    name=file[:-4] 

    logging.info(f"Run hmmratac: {name}")    
    cmd=f"macs3 hmmratac -i {ipath} -f BEDPE -n {name} --outdir {opath}"
    with open(opath + name + '_log.txt', 'w') as log:
        log.write(name)
        log.write("command: " + cmd)
        subprocess.run(cmd, shell=True, stdout=log, stderr=log)



def runp_hmmratac(p, idir, odir):
    files=os.listdir(idir)
    nfile=len(files)
    file_paths=[idir + file for file in files]
    file_outs=list(repeat(odir, nfile))
    # format args as 2d list for starmap
    args=list(zip(file_paths, file_outs))

    with Pool(processes=p) as pool:
        pool.starmap(func=run_hmmratac, iterable=args)


def run_callpeak(ipath, opath):
    file=os.path.basename(ipath)
    name=file[:-4] 

    logging.info(f"Run callpeak: {name}")    
    cmd=f"macs3 callpeak -t {ipath} -f BED -n {name} -g hs --nomodel --shift -100 --extsize 200 --call-summits --outdir {opath}"
    with open(opath + name + '_log.txt', 'w') as log:
        subprocess.run(cmd, shell=True, stdout=log, stderr=log)


def runp_callpeak(p, idir, odir):
    files=os.listdir(idir)
    nfile=len(files)
    file_paths=[idir + file for file in files]
    file_outs=list(repeat(odir, nfile))
    # format args as 2d list for starmap
    args=list(zip(file_paths, file_outs))

    with Pool(processes=p) as pool:
        pool.starmap(func=run_callpeak, iterable=args)



if __name__ == '__main__':
    
    logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.INFO)
    logging.info('Start.')
    startTime = datetime.now()


    idir="/tscc/projects/ps-epigen/users/rlan/sandbox/CallPeaks/splitfrags_out/"
    ipath="/tscc/projects/ps-epigen/users/rlan/sandbox/CallPeaks/splitfrags_out/Cones.bed"
    opath="/tscc/projects/ps-epigen/users/rlan/sandbox/CallPeaks/macs3_outv3/"
    copath="/tscc/projects/ps-epigen/users/rlan/sandbox/CallPeaks/callpeak_out/"

    #idir="/home/rlan/SSHFS/users/rlan/sandbox/CallPeaks/splitfrags_out/"
    #opath="/home/rlan/SSHFS/users/rlan/sandbox/CallPeaks/macs3_outv2/"

    #### Test 
    #subprocess.run(f"macs3 hmmratac -i {ipath} -f BEDPE -n Cones --outdir {opath}", shell=True)

    #run_hmmratac(ipath=ipath, opath=opath)
    #runp_hmmratac(p=8, idir=idir, odir=opath)
    #run_callpeak(ipath=ipath, opath=copath)
    runp_callpeak(p=8, idir=idir, odir=copath)


    logging.info((datetime.now() - startTime)) 
    logging.info('Finish.')


