import os
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
    parser.add_argument('-startid', dest='startid', default=-1, type=int,
                        help='exposure ID; default = -1, run over everything')
    parser.add_argument('-endid', dest='endid', default=1e10, type=int,
                        help='exposure ID; default = 1e10, run over everything')
    parser.add_argument('-snroff', help='w/o recalculating SNR of images',
                        action='store_true')
    parser.add_argument('-cwfsoff', help='w/o running cwfs',
                        action='store_true')
    parser.add_argument('-plotsoff', help='w/o making plots',
                        action='store_true')
    parser.add_argument('-deszcfitsoff', help='w/o reading DES zc fits files',
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

    iexp = 0
    for expid in expidList:
        iexp = iexp + 1
        if (args.startid>0 and int(expid) <args.startid):
            continue
        if (int(expid) > args.endid):
            continue
        imgdir = os.path.join(rvddate, dataset, expid)
        outdir = os.path.join('output', imgdir)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        if args.debugLevel >= 0:
            print(iexp, ': -----', expid, imgdir)

        snrfile = os.path.join(outdir, 'snr.txt')
        if not args.snroff:
            getSNRandType(snrfile,
                imgdir, outerR, inst.obscuration, args.debugLevel)
        continue
        fileSNR, fileType, fileList = readSNRfile(snrfile)

        if not args.plotsoff:
            plotSNR(fileSNR, fileType, fileList, args.snrcut,
                    os.path.join(outdir, 'snr.png'))
            plotExampleDonuts(fileSNR, fileType, fileList, args.snrcut,
                              os.path.join(outdir, 'donuts.png'), imgdir)

        for isenGrp in range(4):
            pairListFile = os.path.join(outdir, 'pairs_grp%d.txt'%isenGrp)
            if not args.cwfsoff:
                zcfile = os.path.join(outdir, 'cwfs_grp%d.txt'%isenGrp)
                # run cwfs
                stamp = cwfsImage(os.path.join(imgdir, fileList[0]), [0, 0], '')
                stampSize = stamp.sizeinPix
                parallelCwfs(imgdir, fileList, isenGrp, fileSNR, fileType, args.snrcut, outerR,
                                stampSize, args.numproc, args.debugLevel, zcfile, pairListFile)
            for isolu in range(2):
                if not args.deszcfitsoff:
                    zcfile = os.path.join(outdir, 'deszc%d_grp%d.txt'%(isolu, isenGrp))
                    readDESzcfits(pairListFile, zcfile)
                        
        # plot zc
        failStatFile = os.path.join(outdir, 'failStat.txt')
        plotComp(outdir, failStatFile)
        
            
def plotComp(outdir, failStatFile):
    
    failStat = np.zeros((12, 3))
    # caustic, 0, total
    # zc0fail, zc0missing, zc0total
    # zc1fail, zc1missing, zc1total
    # repeat above 4 times for 4 sensor groups
    
    x = range(4, znmax + 1)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,10))
    for isenGrp in range(4):
            
        #plt.subplot(2, 2 , isenGrp+1)
        icol = isenGrp%2
        irow = np.int8(np.floor(isenGrp/2))
        # read in cwfs results
        zcfile = os.path.join(outdir, 'cwfs_grp%d.txt'%isenGrp)
        zcarray = np.loadtxt(zcfile)
        if (zcarray.shape[0]>0):
            if zcarray.ndim == 1:
                temp = np.zeros((1, zcarray.shape[0]))
                temp[0, :] = zcarray
                zcarray = temp
            failStat[isenGrp*3, 0] = sum(zcarray[:,-1]>0.5) #caustic flag
            failStat[isenGrp*3, 2] = zcarray.shape[0]
        
            #filter out those with caustic warning
            zcarray = zcarray[zcarray[:,-1]==0, 0:-1]
        else:
            failStat[isenGrp*3:isenGrp*3+3, 0] = 0
            failStat[isenGrp*3:isenGrp*3+3, 2] = 1e-3 #so that later in Matlab we don't divide by 0

            npair = 0
            if isenGrp == 0:
                axes[irow, icol].set_title('%s/%s, %s - %s, npair=%d'%(
                    outdir.split('/')[2],outdir.split('/')[3],
                    intraname[isenGrp], extraname[isenGrp], npair))
            else:
                axes[irow, icol].set_title('%s - %s, npair=%d'%(
                    intraname[isenGrp], extraname[isenGrp], npair))
            
            continue
        npair = zcarray.shape[0]
        # read in DES results
        deszc0file = os.path.join(outdir, 'deszc0_grp%d.txt'%(isenGrp))
        deszc1file = os.path.join(outdir, 'deszc1_grp%d.txt'%(isenGrp))
        deszc0 = np.loadtxt(deszc0file)
        deszc1 = np.loadtxt(deszc1file)
        failStat[isenGrp*3+1, 0] = sum((deszc0[:,0]>1e8) & (deszc0[:,0]<1.5e10))
        failStat[isenGrp*3+1, 1] = sum(deszc0[:,0]>1.5e10) #missing files
        failStat[isenGrp*3+1, 2] = deszc0.shape[0]
        failStat[isenGrp*3+2, 0] = sum((deszc1[:,0]>1e8) & (deszc1[:,0]<1.5e10))
        failStat[isenGrp*3+2, 1] = sum(deszc1[:,0]>1.5e10) #missing files
        failStat[isenGrp*3+2, 2] = deszc1.shape[0]
        
        deszc0 = deszc0*1e3 #convert to nm
        deszc1 = deszc1*1e3 #convert to nm
        #filter out the xe10
        deszc0 = deszc0[deszc0[:,0]<1e8,:]
        deszc1 = deszc1[deszc1[:,0]<1e8,:]
        nzc0 = deszc0.shape[0]
        nzc1 = deszc1.shape[0]
        
        for ipair in range(npair):
            if ipair == 0:
                axes[irow, icol].plot(x, zcarray[ipair, :], label = 'LSST All Pairs',
                        marker='.', color='r', markersize=5, linestyle='--')
            else:
                axes[irow, icol].plot(x, zcarray[ipair, :], # label = '',
                        marker='.', color='r', markersize=5, linestyle='--')
        for izc0 in range(nzc0):
            if izc0 == 0:
                axes[irow, icol].plot(x, deszc0[izc0, :], label = 'DES All z4-11',
                        marker='.', color='b', markersize=5, linestyle='--')
            else:
                axes[irow, icol].plot(x, deszc0[izc0, :], # label = '',
                        marker='.', color='b', markersize=5, linestyle='--')
        for izc1 in range(nzc1):
            if izc1 == 0:
                axes[irow, icol].plot(x, deszc1[izc1, :], label = 'DES All z4-11,14,15',
                        marker='.', color='g', markersize=5, linestyle='--')
            else:
                axes[irow, icol].plot(x, deszc1[izc1, :], # label = '',
                        marker='.', color='g', markersize=5, linestyle='--')
        axes[irow, icol].plot(x, np.mean(zcarray, axis=0), label = 'LSST Ave.',
                marker='v', color='r', markersize=8, linewidth=4)
        axes[irow, icol].plot(x, np.mean(deszc0, axis=0), label = 'DES Ave. z4-11',
                marker='*', color='b', markersize=15, linewidth=4)
        axes[irow, icol].plot(x, np.mean(deszc1, axis=0), label = 'DES Ave. z4-11,14,15',
                marker='o', color='g', markersize=8, linewidth=4)
        
        np.savetxt(os.path.join(outdir, 'ave_grp%d.txt'%isenGrp),
                   np.vstack((np.mean(zcarray, axis=0),
                              np.mean(deszc0, axis=0),
                              np.mean(deszc1, axis=0))))

        axes[irow, icol].grid()
        axes[irow, icol].set_xlabel('Zernike index')
        axes[irow, icol].set_ylabel('Coefficients (nm)')
        axes[irow, icol].legend(loc='best', framealpha = 0.5)
        
        if isenGrp == 0:
            axes[irow, icol].set_title('%s/%s, %s - %s, npair=%d'%(
                outdir.split('/')[2],outdir.split('/')[3],
                intraname[isenGrp], extraname[isenGrp], npair))
        else:
            axes[irow, icol].set_title('%s - %s, npair=%d'%(
                intraname[isenGrp], extraname[isenGrp], npair))
            
    plt.savefig(os.path.join(outdir, 'solution.png'))
    np.savetxt(failStatFile, failStat)
    
def readDESzcfits(pairListFile, zcfile):
    fidr = open(pairListFile, 'r')
    npair = sum(1 for line in fidr) #number of lines/image pairs
    fidr.close()
    fidr = open(pairListFile, 'r')
    zcI1I2 = np.zeros((npair, znmax-3, 2))*np.nan #I1,I2 separate; 
    ipair = 0
    for line in fidr:
        ipair += 1
        for ifits in range(2):
            filename = line.split()[ifits]
            filename = filename.replace('151209','forwardSolutions')
            filename = filename.replace('DECam','v20/DECam')
            if 'zc0' in zcfile:
                isolu = 0
                filename = filename.replace('stamp.fits','first.donut.fits.fz')
            elif 'zc1' in zcfile:
                isolu = 1
                filename = filename.replace('stamp.fits','second.donut.fits.fz')
            try:
                IHDU = fits.open(filename)
                for zi in range(4, znmax+1):
                    try:
                        zcI1I2[ipair-1, zi-4, ifits] = IHDU[0].header['ZERN%d'%zi]
                    except KeyError:
                        zcI1I2[ipair-1, zi-4, ifits] = 0

                if not (np.sqrt(sum(zcI1I2[ipair-1,5-4:10-4,ifits]**2))<3 and
                        abs(zcI1I2[ipair-1,4-4,ifits])>5 and
                        (abs(zcI1I2[ipair-1,4-4,ifits])<15)):
                    zcI1I2[ipair-1, :, ifits] = 1e10
                    print('bad DES solution for -----------  %s'%filename)
                    plt.figure(figsize=(3, 6))
                    plt.subplot(3, 1, 1)
                    plt.imshow(IHDU[1].data, origin='lower')
                    plt.axis('off')
                    plt.title('fit')
                    plt.colorbar()
                    plt.subplot(3, 1, 2)
                    plt.imshow(IHDU[2].data, origin='lower')
                    plt.axis('off')
                    plt.title('data')
                    plt.colorbar()
                    plt.subplot(3, 1, 3)
                    plt.imshow(IHDU[3].data, origin='lower')
                    plt.axis('off')
                    plt.title('data-fit')
                    plt.colorbar()
                    plt.savefig(line.split()[ifits].replace('stamp.fits','DESzc%d.png'%isolu))
                
                IHDU.close()
            except FileNotFoundError:
                print('missing DES file -----------  %s'%filename)
                for zi in range(4, znmax+1):
                    zcI1I2[ipair-1, zi-4, ifits] = 2e10          

    fidr.close()
    #get rid of the large z4 terms in Roodman's results                
    avez4 = np.zeros(2)*np.nan
    for ifits in range(2):
        #average all the z4 for I1, then for I2, across npair
        idx = (zcI1I2[:, 4-4, ifits]<1e8)
        avez4[ifits]=np.mean(zcI1I2[idx, 4-4, ifits])
    for ifits in range(2):
        idx = (zcI1I2[:, 4-4, ifits]<1e8)
        zcI1I2[idx, 4-4, ifits] = np.mean(avez4)

    zcI1I2DES = np.zeros((npair*2, znmax-3))
    for ipair in range(npair):
        for ifits in range(2):
            zcI1I2DES[ipair*2+ifits, :] = zcI1I2[ipair, :, ifits]

    np.savetxt(zcfile, zcI1I2DES)

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
                 stampSize, numproc, debugLevel, zcfile, pairListFile):
    inst = cwfsInstru(instruFile, stampSize)
    algo = cwfsAlgo(algoFile, inst, debugLevel)
    
    intraIdx = np.array([x=='intra%d'%isenGrp for x in fileType])
    extraIdx = np.array([x=='extra%d'%isenGrp for x in fileType])
    npair = min(sum(fileSNR[intraIdx]>snrcut),
                sum(fileSNR[extraIdx]>snrcut))
    aa = np.array(range(len(fileSNR)))
    intraIdxS = sorted(aa[intraIdx],key=lambda x:fileSNR[x], reverse=True)
    extraIdxS = sorted(aa[extraIdx],key=lambda x:fileSNR[x], reverse=True)

    argList = []
    fidw = open(pairListFile, 'w')
    for ipair in range(npair):
        I1File = os.path.join(imgdir, fileList[intraIdxS[ipair]])
        stamp = cwfsImage(I1File, [0, 0], '')
        stamp.rmBkgd(outerR, debugLevel)
        stamp.normalizeI(outerR, inst.obscuration)
        I1Array = stamp.image
        
        I2File = os.path.join(imgdir, fileList[extraIdxS[ipair]])
        stamp = cwfsImage(I2File, [0, 0], '')
        stamp.rmBkgd(outerR, debugLevel)
        stamp.normalizeI(outerR, inst.obscuration)
        I2Array = stamp.image
        
        # test, pdb cannot go into the subprocess
        # runcwfs(I1Array, I2Array, inst, algo)
        argList.append((I1Array, I2Array, inst, algo))
        fidw.write('%s\t%s\n'%(I1File, I2File))
    fidw.close()
    pool = multiprocessing.Pool(numproc)
    zcarray = pool.map(runcwfs, argList)
    pool.close()
    pool.join()
        
    np.savetxt(zcfile, zcarray)

def runcwfs(argList):
    I1File = argList[0]
    I2File = argList[1]
    inst = argList[2]
    algo = argList[3]
    
    I1Field = [0, 0]
    I2Field = [0, 0]
    
    I1 = cwfsImage(I1File, I1Field, 'intra')
    I2 = cwfsImage(I2File, I2Field, 'extra')
    algo.reset(I1, I2)
    algo.runIt(inst, I1, I2, model)

    return np.append(algo.zer4UpNm, algo.caustic)
    
def plotExampleDonuts(fileSNR, fileType, fileList, snrcut, pngfile, imgdir):
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
            plt.subplot(4,2*nplot,isub)
            plt.imshow(I1, origin='lower', vmin=myvmin, vmax=myvmax)
            plt.axis('off')
            plt.title(I1title)
            plt.subplot(4,2*nplot,isub+2*nplot)
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
            
def plotSNR(fileSNR, fileType, fileList, snrcut, pngfile):

    # set the saturated images to SNR=1e-5
    idx = fileSNR<0
    nsat = sum(idx)
    for i in range(nsat):
        print(fileList[i])
        
    for isenGrp in range(4):
        plt.subplot(2, 2, isenGrp+1)
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
    
