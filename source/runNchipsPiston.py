import os
import argparse
import multiprocessing

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from cwfsAlgo import cwfsAlgo
from cwfsInstru import cwfsInstru
from cwfsImage import cwfsImage

traname = ['intra', 'extra']
senGrpName = ['N', 'S']
senGrpMarker = ['*','o']
nsensor = 31
algoFile = 'exp'
model = 'paraxial'
znmax = 22

def main():
    parser = argparse.ArgumentParser(
        description='----- runNchipsPiston.py ---------')
    parser.add_argument('snrcut', type=float, help='threshold on SNR')
    parser.add_argument('dl', type=float, help='image space offset')
    parser.add_argument('-startid0', dest='startid0', default=-1, type=int,
                        help='exposure ID; default = -1, run over everything')
    parser.add_argument('-endid0', dest='endid0', default=1e10, type=int,
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
    dataset = ['', '']
    if abs(args.dl - 3.0)< 1e-4:
        dataset[0] = '20140807s2'
        dataset[1] = '20140807s3'
    elif abs(args.dl - 1.5)< 1e-4:
        dataset[0] = '20140807s5'
        dataset[1] = '20140807s4'

    instruFile = 'decamHex%d'%(args.dl*10)

    expidList = [[], []]
    for itra in range(2):
        expidList0 = os.listdir(os.path.join(rvddate, dataset[itra]))
        expidList[itra] = expidList0.copy()
        for expid in expidList0:
            if not expid.isdigit():
                expidList[itra].remove(expid)

    # read in the DECam parameters
    inst = cwfsInstru(instruFile, 128) #128 is tentative, doesn't matter
    outerR = inst.donutR/inst.pixelSize
    if args.debugLevel >= 0:
        print(outerR)

    expid = [0, 0]
    imgdir = ['', '']
    outdir = ['','']
    snrfile = ['', '']
    for iexp in range(len(expidList[0])):

        for itra in range(2):
            expid[itra] = expidList[itra][iexp]
            imgdir[itra] = os.path.join(rvddate, dataset[itra], expid[itra])
            outdir[itra] = os.path.join('output', imgdir[itra])
            if not os.path.isdir(outdir[itra]):
                os.makedirs(outdir[itra])
            snrfile[itra] = os.path.join(outdir[itra], 'snr.txt')
        if (args.startid0>0 and int(expid[0]) <args.startid0):
            continue
        if (int(expid[0]) > args.endid0):
            continue
        if args.debugLevel >= 0:
            print(iexp+1, ': ----Intra', imgdir[0], '  -----Extra', imgdir[1])

        for itra in range(2):
            if not args.snroff:
                getSNRandType(snrfile[itra], traname[itra],
                    imgdir[itra], outerR, inst.obscuration, args.debugLevel)
        fileSNR, fileType, fileList = readSNRfile(snrfile)

        if not args.plotsoff:
            plotSNR(fileSNR, fileType, args.snrcut,
                    os.path.join(outdir[0], 'snr.png')) #this is only in the intra dir
            plotExampleDonuts(fileList, fileType, fileSNR, args.snrcut,
                              os.path.join(outdir[0], 'donuts.png'), imgdir)

        pairListFile = os.path.join(outdir[0], 'pairs.txt')
        if not args.cwfsoff:
            zcfile = os.path.join(outdir[0], 'cwfs.txt')
            # run cwfs
            stamp = cwfsImage(os.path.join(imgdir[0], fileList[0]), [0, 0], '')
            stampSize = stamp.sizeinPix
            parallelCwfs(instruFile, imgdir, fileList, fileSNR, fileType, args.snrcut, outerR,
                            stampSize, args.numproc, args.debugLevel, zcfile, pairListFile)

        for isolu in range(2):
            if not args.deszcfitsoff:
                zcfile = os.path.join(outdir[0], 'deszc%d.txt'%(isolu))
                readDESzcfits(pairListFile, zcfile, args.dl)
                        
        # plot zc
        failStatFile = os.path.join(outdir[0], 'failStat.txt')
        plotComp(outdir, pairListFile, failStatFile, 31)
        # plotComp(outdir, pairListFile, failStatFile, 4)        

            
def plotComp(outdir, pairListFile, failStatFile, nplot):

    failStat = np.zeros((3, 3))
    # caustic, 0, total
    # zc0fail, zc0missing, zc0total
    # zc1fail, zc1missing, zc1total
    
    # read in cwfs results
    zcfile = os.path.join(outdir[0], 'cwfs.txt')
    zcarray = np.loadtxt(zcfile)
    failStat[0, 0] = sum(zcarray[:,-1]>0.5) #caustic flag
    failStat[0, 2] = zcarray.shape[0]
    # filter out those with caustic warning
    zcarray[zcarray[:,-1]>0.5, 0:-1] = 1e10
    zcarray = zcarray[:, 0:-1]
    npairTot = zcarray.shape[0]
    
    # read in DES results
    deszc0file = os.path.join(outdir[0], 'deszc0.txt')
    deszc1file = os.path.join(outdir[0], 'deszc1.txt')
    deszc0 = np.loadtxt(deszc0file)
    deszc1 = np.loadtxt(deszc1file)
    failStat[1, 0] = sum((deszc0[:,0]>1e8) & (deszc0[:,0]<1.5e10))
    failStat[1, 1] = sum(deszc0[:,0]>1.5e10) #missing files
    failStat[1, 2] = deszc0.shape[0]
    failStat[2, 0] = sum((deszc1[:,0]>1e8) & (deszc1[:,0]<1.5e10))
    failStat[2, 1] = sum(deszc1[:,0]>1.5e10) #missing files
    failStat[2, 2] = deszc1.shape[0]
    np.savetxt(failStatFile, failStat)
    # return

    intraFileList = []
    extraFileList = []
    fidr = open(pairListFile, 'r')
    for line in fidr:
        intraFileList.append(line.split()[0])
        extraFileList.append(line.split()[1])        
    fidr.close()
            
    x = range(4, znmax + 1)
    nrow = 7
    ncol =5
    if nplot == 4:
        nrow = 2
        ncol = 2
        
    for isenGrp in range(2):
        plt.figure(figsize=(15,15))
        for isensor in range(nplot):

            idx = np.array(['%s%d.'%(
                senGrpName[isenGrp], isensor+1) in intraFileList[i] for i in range(npairTot)])
            idx2 = np.repeat(idx, 2)
            
            plt.subplot(nrow, ncol, isensor+1)
            npair = sum(idx)
            zcsub = zcarray[idx, :]
            for ipair in range(npair):
                if ipair == 0:
                    plt.plot(x, zcsub[ipair, :], label = 'LSST All Pairs',
                            marker='.', color='r', markersize=5, linestyle='--')
                else:
                    plt.plot(x, zcsub[ipair, :], # label = '',
                            marker='.', color='r', markersize=5, linestyle='--')
            #filter out the xe10
            zc0sub = deszc0[idx2, :]
            zc1sub = deszc1[idx2, :]
            zc0sub = zc0sub[zc0sub[:,0]<1e8,:]
            zc1sub = zc1sub[zc1sub[:,0]<1e8,:]
            nzc0 = zc0sub.shape[0]
            nzc1 = zc1sub.shape[0]
            zc0sub = zc0sub*1e3
            zc1sub = zc1sub*1e3
            for izc0 in range(nzc0):
                if izc0 == 0:
                    plt.plot(x, zc0sub[izc0, :], label = 'DES All z4-11',
                            marker='.', color='b', markersize=5, linestyle='--')
                else:
                    plt.plot(x, zc0sub[izc0, :], # label = '',
                            marker='.', color='b', markersize=5, linestyle='--')
            for izc1 in range(nzc1):
                if izc1 == 0:
                    plt.plot(x, zc1sub[izc1, :], label = 'DES All z4-11,14,15',
                            marker='.', color='g', markersize=5, linestyle='--')
                else:
                    plt.plot(x, zc1sub[izc1, :], # label = '',
                            marker='.', color='g', markersize=5, linestyle='--')
            plt.plot(x, np.mean(zcsub, axis=0), label = 'LSST Ave.',
                    marker='v', color='r', markersize=5, linewidth=2)
            plt.plot(x, np.mean(zc0sub, axis=0), label = 'DES Ave. z4-11',
                    marker='*', color='b', markersize=7, linewidth=2)
            plt.plot(x, np.mean(zc1sub, axis=0), label = 'DES Ave. z4-11,14,15',
                    marker='o', color='g', markersize=5, linewidth=2)
            
            np.savetxt(os.path.join(outdir[0], 'ave_%s%d.txt'%(senGrpName[isenGrp], isensor)),
                    np.vstack((np.mean(zcsub, axis=0),
                                np.mean(zc0sub, axis=0),
                                np.mean(zc1sub, axis=0))))
    
            plt.grid()
            #plt.xlabel('Zernike index')
            #plt.ylabel('Coefficients (nm)')
            l = plt.legend(loc='best', framealpha = 0.5, fontsize=5,
                       title='%s%d, npair=%d'%(
                    senGrpName[isenGrp],isensor+1, npair))
            l.get_title().set_fontsize(5)
            if isensor == 0:
                plt.title('%s - %s'%(
                    outdir[0].split('/')[3], outdir[1].split('/')[3]))

        #plt.tight_layout()
        pngname = 'solution%d_%s.png'%(nplot, senGrpName[isenGrp])
        plt.savefig(os.path.join(outdir[0], pngname), dpi=300)
        plt.close('all')
    
def readDESzcfits(pairListFile, zcfile, dl):
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
    
                z4low = 5
                z4high = 15
                if abs(dl - 3.0)< 1e-4:
                    z4low = 10
                    z4high = 30
                if not (np.sqrt(sum(zcI1I2[ipair-1,5-4:10-4,ifits]**2))<3 and
                        abs(zcI1I2[ipair-1,4-4,ifits])>z4low and
                        (abs(zcI1I2[ipair-1,4-4,ifits])<z4high)):
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
                    plt.close('all')
                    
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
    for itra in range(2):
        fidr = open(snrfile[itra], 'r')
        for line in fidr:
            fileList.append(line.split()[0])
            fileType.append(line.split()[1])
            fileSNR.append(float(line.split()[2]))
        fidr.close()
    fileSNR = np.array(fileSNR)
    return fileSNR, fileType, fileList
            
def parallelCwfs(instruFile, imgdir, fileList, fileSNR, fileType, snrcut, outerR,
                 stampSize, numproc, debugLevel, zcfile, pairListFile):
    inst = cwfsInstru(instruFile, stampSize)
    algo = cwfsAlgo(algoFile, inst, debugLevel)

    argList = []
    fidw = open(pairListFile, 'w')
    for isensor in range(nsensor):
        for isenGrp in range(2):
    
            intraIdx = np.array([x=='intra%s%d'%(
                senGrpName[isenGrp], isensor+1) for x in fileType])
            extraIdx = np.array([x=='extra%s%d'%(
                senGrpName[isenGrp], isensor+1) for x in fileType])

            npair = min(sum(fileSNR[intraIdx]>snrcut),
                sum(fileSNR[extraIdx]>snrcut))
            aa = np.array(range(len(fileSNR)))
            intraIdxS = sorted(aa[intraIdx],key=lambda x:fileSNR[x], reverse=True)
            extraIdxS = sorted(aa[extraIdx],key=lambda x:fileSNR[x], reverse=True)

            for ipair in range(npair):
                I1File = os.path.join(imgdir[0], fileList[intraIdxS[ipair]])
                stamp = cwfsImage(I1File, [0, 0], '')
                stamp.rmBkgd(outerR, debugLevel)
                stamp.normalizeI(outerR, inst.obscuration)
                I1Array = stamp.image
        
                I2File = os.path.join(imgdir[1], fileList[extraIdxS[ipair]])
                stamp = cwfsImage(I2File, [0, 0], '')
                stamp.rmBkgd(outerR, debugLevel)
                stamp.normalizeI(outerR, inst.obscuration)
                I2Array = stamp.image
                
                argList.append((I1Array, I2Array, inst, algo))
                fidw.write('%s\t%s\n'%(I1File, I2File))
        
                # test, pdb cannot go into the subprocess
                # aa = runcwfs(argList[0])
                
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
    
def plotExampleDonuts(fileList, fileType, fileSNR, snrcut, pngfile, imgdir):
    # take a look at some example images 
    # plot 1 pairs from each sensor group    

    nplot = 6
    plt.figure(figsize=(10, 8))
    for isenGrp in range(2):
        for isensor  in range(nplot):
        
            intraIdx = np.array([x=='intra%s%d'%(
                senGrpName[isenGrp], isensor+1) for x in fileType])
            extraIdx = np.array([x=='extra%s%d'%(
                senGrpName[isenGrp], isensor+1) for x in fileType])

            aa = np.array(range(len(fileSNR)))
            intraIdxS = sorted(aa[intraIdx],key=lambda x:fileSNR[x], reverse=True)
            extraIdxS = sorted(aa[extraIdx],key=lambda x:fileSNR[x], reverse=True)

            I1name = os.path.join(imgdir[0], fileList[intraIdxS[0]])
            IHDU = fits.open(I1name)
            I1 = IHDU[0].data
            IHDU.close()
            I1namesplt = I1name.split('.')
            I1title = '%s %s'%(I1namesplt[1], I1namesplt[2])
            
            I2name = os.path.join(imgdir[1], fileList[extraIdxS[0]])
            IHDU = fits.open(I2name)
            I2 = IHDU[0].data
            IHDU.close()
            I2namesplt = I2name.split('.')
            I2title = '%s %s'%(I2namesplt[1], I2namesplt[2])
            
            isub=isensor+1+2*nplot*isenGrp
            plt.subplot(4,nplot,isub)
            plt.imshow(I1, origin='lower')
            plt.axis('off')
            plt.title(I1title)
            plt.subplot(4,nplot,isub+nplot)
            plt.imshow(I2, origin='lower')
            plt.axis('off')
            plt.title(I2title)
    plt.savefig(pngfile)
    
def getSNRandType(snrfile, fileType, imgdir, outerR, obsR, debugLevel):
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
        
        fidw.write('%s\t %s \t%7.2f \t %10.1f\t %10.1f\n'%
                   (filename, '%s%s'%(fileType, filename.split('.')[1]),
                    stamp.SNR, stamp.SNRsig, stamp.SNRbg))
    fidw.close()
            
def plotSNR(fileSNR, fileType, snrcut, pngfile):
    
    plt.figure(figsize=(15, 15))
    for isensor in range(nsensor):
        plt.subplot(7, 5, isensor+1)
        npair = [0, 0]
        for isenGrp in range(2):
            intraIdx = np.array([x=='intra%s%d'%(
                senGrpName[isenGrp], isensor+1) for x in fileType])
            plt.plot(sorted(fileSNR[intraIdx], reverse=True), marker=senGrpMarker[isenGrp],
                    color='r', label='intra%s%d'%(senGrpName[isenGrp], isensor+1))
            extraIdx = np.array([x=='extra%s%d'%(
                senGrpName[isenGrp], isensor+1) for x in fileType])
            plt.plot(sorted(fileSNR[extraIdx], reverse=True), marker=senGrpMarker[isenGrp],
                    color='b', label='extra%s%d'%(senGrpName[isenGrp], isensor+1))
            plt.plot([0, max(sum(intraIdx), sum(extraIdx))],[snrcut, snrcut],
                    color='k', linestyle='-')
            npair[isenGrp] = min(sum(fileSNR[intraIdx]>snrcut),
                        sum(fileSNR[extraIdx]>snrcut))
        plt.xlabel('donut number')
        plt.ylabel('SNR')
        plt.grid()
        plt.legend(fontsize=8)
        plt.title('%s: %d, %s: %d pairs>20'%(
            senGrpName[0], npair[0], senGrpName[1], npair[1]))

    plt.tight_layout()
    plt.savefig(pngfile)
    
if __name__ == "__main__":
    main()
    
