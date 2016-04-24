function [] = diff_zc(plotStyle)

% dataType can only be 's1'. 
% for 's2' and 's5', use diff_zc_NS()

% plotStyle = 'diff', 'ratio', 'histogram', 'scatter', 'history'

% dataset = 'output/skymap/20140613s1';
dataset = 'output/nights/20140315s1';

sensorName = {'FN1/N2','FN3/N4','FS2/S1','FS4/S3'};
expIdList = dir(dataset);
nexp = 0;
for i = 1:size(expIdList,1)
    filename = sprintf('%s/%s/snr.txt',dataset,expIdList(i).name);
    if exist(filename, 'file')
        nexp = nexp + 1;
    end
end

znmax = 15;
cwfs = zeros(4, nexp, znmax-3);
zc0 = zeros(4, nexp, znmax-3);
zc1 = zeros(4, nexp, znmax-3);
iexp = 0;
for i = 1:size(expIdList,1)
    filename = sprintf('%s/%s/snr.txt',dataset,expIdList(i).name);
    if exist(filename, 'file')
        iexp = iexp + 1;
        for isenGrp=0:3
            filename=sprintf('%s/%s/ave_grp%d.txt',dataset,expIdList(i).name,isenGrp);
            if exist(filename, 'file')
                data = load(filename);
                cwfs(isenGrp+1, iexp, :) =  data(1,1:znmax-3);
                zc0(isenGrp+1, iexp, :) =  data(2,1:znmax-3);
                zc1(isenGrp+1, iexp, :) =  data(3,1:znmax-3);
            else
                cwfs(isenGrp+1, iexp, :) =  nan;
                zc0(isenGrp+1, iexp, :) =  nan;
                zc1(isenGrp+1, iexp, :) =  nan;
            end
        end
    end
end

if strcmp(plotStyle, 'diff')
    % first - cwfs
    figure(1);clf;
    for isenGrp=0:3
        subplot(2,2,isenGrp+1);
        plot(4:znmax, squeeze(zc0(isenGrp+1, :, :) - cwfs(isenGrp+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        xlabel('Zernike Number');ylabel('coefficient (in nm)');
        title('FM (z4-11) - cwfs');
    end
    % second - cwfs
    figure(2);clf;
    for isenGrp=0:3
        subplot(2,2,isenGrp+1);
        plot(4:znmax, squeeze(zc1(isenGrp+1, :, :) - cwfs(isenGrp+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        xlabel('Zernike Number');ylabel('coefficient (in nm)');
        title('FM (z4-11, 14,15) - cwfs');
    end
elseif strcmp(plotStyle, 'ratio')
    % first/cwfs
    figure(1);clf;
    for isenGrp=0:3
        subplot(2,2,isenGrp+1);
        plot(4:znmax, squeeze(zc0(isenGrp+1, :, :)./cwfs(isenGrp+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title('FM (z4-11) / cwfs');
    end
    % second/cwfs
    figure(2);clf;
    for isenGrp=0:3
        subplot(2,2,isenGrp+1);
        plot(4:znmax, squeeze(zc1(isenGrp+1, :, :)./cwfs(isenGrp+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title('FM (z4-11, 14,15) / cwfs');
    end
elseif strcmp(plotStyle, 'histogram')
    % first - cwfs, by Zernike
    xbin = -500:50:500;
    figure(1);clf; %zc0
    for iz=4:znmax
        subplot(5,4,iz-3);
        aa=hist(reshape(squeeze(zc0(:, :, iz-3) - cwfs(:, :, iz-3)),[], ...
            1),xbin);
        aa=aa/(sum(aa))*100;
        bar(xbin,aa);
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in %%)',iz));
    end
    figure(2);clf; %zc1
    for iz=4:znmax
        subplot(5,4,iz-3);
        aa=hist(reshape(squeeze(zc1(:, :, iz-3) - cwfs(:, :, iz-3)),[], ...
            1),xbin);
        aa=aa/(sum(aa))*100;
        bar(xbin,aa);
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in %%)',iz));
    end
elseif strcmp(plotStyle, 'scatter')
    % scatter(first, cwfs), by Zernike
    xlow = -1000;
    xhigh = 1000;
    figure(1);clf; %zc0
    for iz=4:znmax
        subplot(5,4,iz-3);
        scatter(reshape(squeeze(zc0(:, :, iz-3)),[],1), reshape(cwfs(:, :, iz-3),[],1),180,'.');
        line([xlow, xhigh],[xlow, xhigh],'color','r');
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in nm)',iz));
        xlabel('FM'); ylabel('LSST');
    end
    figure(2);clf; %zc1
    for iz=4:znmax
        subplot(5,4,iz-3);
        scatter(reshape(squeeze(zc1(:, :, iz-3)),[],1), reshape(cwfs(:, :, iz-3),[],1),180,'.');
        line([xlow, xhigh],[xlow, xhigh],'color','r');
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in nm)',iz));
        xlabel('FM'); ylabel('LSST');
    end
elseif strcmp(plotStyle, 'history')
    x=1:nexp;
    for iz=4:znmax
        figure(1);clf; %zc0
        ymax=max(max(max(squeeze(zc0(:, :, iz-3)))),max(max(squeeze(cwfs(:,:,iz-3)))));
        ymin=min(min(min(squeeze(zc0(:, :, iz-3)))),min(min(squeeze(cwfs(:,:,iz-3)))));
        for isenGrp=1:4
            subplot(4,1,isenGrp);
            fm=(squeeze(zc0(isenGrp, :, iz-3)));
            cw=(squeeze(cwfs(isenGrp,:, iz-3)));
            idxfm = ~isnan(fm);
            idxcw = ~isnan(cw);
            idx = (idxfm & idxcw);
            myrms = rms(fm(idx)-cw(idx));
            plot(x, fm, '-r.', x, cw, '-b.','markersize',10);
            grid on;
            if isenGrp ==1
                legend({'FM','LSST'},'location','best');
            end
            text(0.81,0.9, sprintf('z%d (rms diff=%3.0fnm)',iz, myrms), 'units','Normalized');
            text(0.01,0.9, sensorName{isenGrp}, 'units','Normalized');
            ylabel(sprintf('Z%d (in nm)', iz));
            ylim([ymin ymax]);
        end
        samexaxis('xmt','on','ytac','join','yld',1);
        xlabel('Exposure Number');
    end
    figure(2);clf; %zc1
    for iz=4:znmax
        subplot(5,4,iz-3);
        fm=mean(squeeze(zc1(:, :, iz-3)));
        cw=mean(squeeze(cwfs(:,:, iz-3)));
        idxfm = ~isnan(fm);
        idxcw = ~isnan(cw);
        idx = (idxfm & idxcw);
        myrms = rms(fm(idx)-cw(idx));
        xfm = x(idxfm);
        fm=fm(idxfm);
        xcw=x(idxcw);
        cw=cw(idxcw);
        plot(xfm, fm, '-ro', xcw, cw, '-b*');
        grid on;
        title(sprintf('z%d (rms diff=%3.0fnm)',iz, myrms));
        legend({'FM','LSST'},'location','best');
    end
    
elseif strcmp(plotStyle, 'history1')
    x=1:nexp;
    figure(1);clf; %zc0
    for iz=4:znmax
        subplot(3,4,iz-3);
        fm=mean(squeeze(zc0(:, :, iz-3)));
        cw=mean(squeeze(cwfs(:,:, iz-3)));
        idxfm = ~isnan(fm);
        idxcw = ~isnan(cw);
        idx = (idxfm & idxcw);
        myrms = rms(fm(idx)-cw(idx));
        xfm = x(idxfm);
        fm=fm(idxfm);
        xcw=x(idxcw);
        cw=cw(idxcw);
        plot(xfm, fm, '-ro', xcw, cw, '-b*');
        grid on;
        title(sprintf('z%d (rms diff=%3.0fnm)',iz, myrms));
        legend({'FM','LSST'},'location','best');
    end
    figure(2);clf; %zc1
    for iz=4:znmax
        subplot(3,4,iz-3);
        fm=mean(squeeze(zc1(:, :, iz-3)));
        cw=mean(squeeze(cwfs(:,:, iz-3)));
        idxfm = ~isnan(fm);
        idxcw = ~isnan(cw);
        idx = (idxfm & idxcw);
        myrms = rms(fm(idx)-cw(idx));
        xfm = x(idxfm);
        fm=fm(idxfm);
        xcw=x(idxcw);
        cw=cw(idxcw);
        plot(xfm, fm, '-ro', xcw, cw, '-b*');
        grid on;
        title(sprintf('z%d (rms diff=%3.0fnm)',iz, myrms));
        legend({'FM','LSST'},'location','best');
    end

end

end
