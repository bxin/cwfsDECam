import os
import sys
import argparse
import multiprocessing

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from cwfsAlgo import cwfsAlgo
from cwfsInstru import cwfsInstru
from cwfsImage import cwfsImage

intraname=['FN1','FN3','FS2','FS4']
extraname=['FN2','FN4','FS1','FS3']
instruFile = 'decam15'
algoFile = 'exp'
model = 'paraxial'
znmax = 22

def main():
    parser = argparse.ArgumentParser(
        description='----- runNstampWFS.py ---------')
    parser.add_argument('snrcut', type=float, help='threshold on SNR')
    parser.add_argument('-snroff', help='w/o recalculating SNR of images',
                        action='store_true')
    parser.add_argument('-cwfsoff', help='w/o running cwfs',
                        action='store_true')
    parser.add_argument('-plotsoff', help='w/o making plots',
                        action='store_true')
    parser.add_argument('-p', dest='numproc', default=1, type=int,
                        help='Number of Processors Phosim uses')
    parser.add_argument('-d', dest='debugLevel', type=int,
                        default=0, choices=(-1, 0, 1),
                        help='debug level, -1=quiet, 0=normal, \
                        1=verbose, default=0')
    args = parser.parse_args()
    if args.debugLevel >= 1:
        print(args)

    rvddate = '151209'
    dataset = '20140613s1'

    expidList0 = os.listdir(os.path.join(rvddate, dataset))
    expidList = expidList0.copy()
    for expid in expidList0:
        if not expid.isdigit():
            expidList.remove(expid)

    # read in the DECam parameters
    inst = cwfsInstru(instruFile, 128) #128 is tentative, doesn't matter
    outerR = inst.donutR/inst.pixelSize
    if args.debugLevel >= 0:
        print(outerR)
            
    for expid in expidList:
        imgdir = os.path.join(rvddate, dataset, expid)
        if args.debugLevel >= 0:
            print(expid, imgdir)

        snrfile = os.path.join(imgdir, 'snr.txt')
        if not args.snroff:
            getSNRandType(snrfile,
                imgdir, outerR, inst.obscuration, args.debugLevel)
        fileSNR, fileType, fileList = readSNRfile(snrfile)

        if not args.plotsoff:
            plotSNR(fileSNR, fileType, args.snrcut,
                    os.path.join(imgdir, 'snr.png'))
            plotExampleDonuts(fileList, fileType, fileSNR, args.snrcut,
                              os.path.join(imgdir, 'donuts.png'), imgdir)

        for isenGrp in range(4):
            zcfile = os.path.join(imgdir, 'cwfs_grp%d.txt'%isenGrp)
            if args.cwfsoff:
                # read in cwfs results
                a=1
            else:
                # run cwfs
                stamp = cwfsImage(os.path.join(imgdir, fileList[0]), [0, 0], '')
                stampSize = stamp.sizeinPix
                parallelCwfs(imgdir, fileList, isenGrp, fileSNR, fileType, args.snrcut, outerR,
                                stampSize, args.numproc, args.debugLevel, zcfile)
                sys.exit()
            ax = plt.subplot(2, 2 , isenGrp)
            # plot zc
            
        sys.exit()

def readSNRfile(snrfile):
    fileSNR = []
    fileType =[]
    fileList=[]
    fidr = open(snrfile, 'r')
    for line in fidr:
        fileList.append(line.split()[0])
        fileType.append(line.split()[1])
        fileSNR.append(float(line.split()[2]))
    fidr.close()
    fileSNR = np.array(fileSNR)
    return fileSNR, fileType, fileList
            
def parallelCwfs(imgdir, fileList, isenGrp, fileSNR, fileType, snrcut, outerR,
                 stampSize, numproc, debugLevel, zcfile):
    inst = cwfsInstru(instruFile, stampSize)
    algo = cwfsAlgo(algoFile, inst, debugLevel)
    jobs = []
    counter = 0

    intraIdx = np.array([x=='intra%d'%isenGrp for x in fileType])
    extraIdx = np.array([x=='extra%d'%isenGrp for x in fileType])
    npair = min(sum(fileSNR[intraIdx]>snrcut),
                sum(fileSNR[extraIdx]>snrcut))
    aa = np.array(range(len(fileSNR)))
    intraIdxS = sorted(aa[intraIdx],key=lambda x:fileSNR[x], reverse=True)
    extraIdxS = sorted(aa[extraIdx],key=lambda x:fileSNR[x], reverse=True)

    zcarray = np.zeros((npair, znmax-3))
    for ipair in range(npair):
        stamp = cwfsImage(os.path.join(imgdir, fileList[intraIdxS[ipair]]), [0, 0], '')
        stamp.rmBkgd(outerR, debugLevel)
        stamp.normalizeI(outerR, inst.obscuration)
        I1Array = stamp.image
        stamp = cwfsImage(os.path.join(imgdir, fileList[extraIdxS[ipair]]), [0, 0], '')
        stamp.rmBkgd(outerR, debugLevel)
        stamp.normalizeI(outerR, inst.obscuration)
        I2Array = stamp.image
        # test, pdb cannot go into the subprocess
        # runcwfs(ipair, I1Array, I2Array, inst, algo, zcarray)
        p = multiprocessing.Process(
            target=runcwfs, name='cwfs%d' % ipair, args=(ipair,
                I1Array, I2Array, inst, algo, zcarray))
        jobs.append(p)
        p.start()
        counter += 1
        if (counter == numproc) or (ipair == npair - 1):
            for p in jobs:
                p.join()
            counter = 0
            jobs = []
    np.savetxt(zcfile, zcarray)

def runcwfs(ipair, I1File, I2File, inst, algo, zcarray):
    I1Field = [0, 0]
    I2Field = [0, 0]
    
    I1 = cwfsImage(I1File, I1Field, 'intra')
    I2 = cwfsImage(I2File, I2Field, 'extra')
    algo.reset(I1, I2)
    algo.runIt(inst, I1, I2, model)

    zcarray[ipair, :] = algo.zer4UpNm
    
def plotExampleDonuts(fileList, fileType, fileSNR, snrcut, pngfile, imgdir):
    # take a look at some example images 
    # plot 3 pairs from each sensor group    
    nplot = 3

    for isenGrp in range(4):
        intraIdx = np.array([x=='intra%d'%isenGrp for x in fileType])
        extraIdx = np.array([x=='extra%d'%isenGrp for x in fileType])
        npair = min(sum(fileSNR[intraIdx]>snrcut),
                     sum(fileSNR[extraIdx]>snrcut))
        aa = np.array(range(len(fileSNR)))
        intraIdxS = sorted(aa[intraIdx],key=lambda x:fileSNR[x], reverse=True)
        extraIdxS = sorted(aa[extraIdx],key=lambda x:fileSNR[x], reverse=True)
        for i in range(min(nplot, npair)):
            I1name = os.path.join(imgdir, fileList[intraIdxS[i]])
            IHDU = fits.open(I1name)
            I1 = IHDU[0].data
            IHDU.close()
            I1namesplt = I1name.split('.')
            I1title = '%s %s'%(I1namesplt[1], I1namesplt[2])
            
            I2name = os.path.join(imgdir, fileList[extraIdxS[i]])
            IHDU = fits.open(I2name)
            I2 = IHDU[0].data
            IHDU.close()
            I2namesplt = I2name.split('.')
            I2title = '%s %s'%(I2namesplt[1], I2namesplt[2])

            if i==0:
                myvmin = min(np.min(I2), np.min(I1))
                myvmax = max(np.max(I2), np.max(I1))
            isub=i+1+nplot*isenGrp
            if isenGrp>=2:
                isub=isub+2*nplot
            ax = plt.subplot(4,2*nplot,isub)
            plt.imshow(I1, origin='lower', vmin=myvmin, vmax=myvmax)
            plt.axis('off')
            plt.title(I1title)
            ax = plt.subplot(4,2*nplot,isub+2*nplot)
            plt.imshow(I2, origin='lower', vmin=myvmin, vmax=myvmax)
            plt.axis('off')
            plt.title(I2title)
    plt.savefig(pngfile)
    
def getSNRandType(snrfile, imgdir, outerR, obsR, debugLevel):
    fileList0 = os.listdir(imgdir)
    fileList = fileList0.copy()
    for filename in fileList0:
        if not filename.endswith('fits'):
            fileList.remove(filename)

    fidw = open(snrfile, 'w')
    for filename in fileList:
        if debugLevel >= 1:
            print(filename)
        stamp = cwfsImage(os.path.join(imgdir, filename), [0, 0], '')
        stamp.rmBkgd(outerR, debugLevel)
        stamp.getSNR(outerR, obsR)
        for isenGrp in range(len(intraname)):
            if intraname[isenGrp] in filename:
                fileType = ('intra%d'%isenGrp)
                break
            elif extraname[isenGrp] in filename:
                fileType = ('extra%d'%isenGrp)
                break

        fidw.write('%s\t %s \t%7.2f \t %10.1f\t %10.1f\n'%
                   (filename, fileType, stamp.SNR, stamp.SNRsig, stamp.SNRbg))
    fidw.close()
            
def plotSNR(fileSNR, fileType, snrcut, pngfile):
    for isenGrp in range(4):
        ax = plt.subplot(2, 2, isenGrp+1)
        intraIdx = np.array([x=='intra%d'%isenGrp for x in fileType])
        plt.plot(sorted(fileSNR[intraIdx], reverse=True), marker='o',
                 color='r', label='intra%d'%isenGrp)
        extraIdx = np.array([x=='extra%d'%isenGrp for x in fileType])
        plt.plot(sorted(fileSNR[extraIdx], reverse=True), marker='o',
                 color='b', label='extra%d'%isenGrp)
        plt.plot([0, max(sum(intraIdx), sum(extraIdx))],[snrcut, snrcut],
                 color='k', linestyle='-')
        plt.xlabel('donut number')
        plt.ylabel('SNR')
        plt.grid()
        plt.legend()
        npair = min(sum(fileSNR[intraIdx]>snrcut),
                     sum(fileSNR[extraIdx]>snrcut))
        plt.title('%s - %s, %d pairs>20'%(intraname[isenGrp],
                                          extraname[isenGrp], npair))
    plt.tight_layout()
    plt.savefig(pngfile)
    
if __name__ == "__main__":
    main()
    
