function [] = diff_zc(plotStyle,znmax)

% dataType can only be 's1'. 
% for 's2' and 's5', use diff_zc_NS()

% plotStyle = 'diff', 'ratio', 'histogram', 'scatter', 'history'
% quadrafoil = 0: use z4-z11; =1: use z4-z11,14,15
% znmax = 11, 15, or 22

if znmax==11
    fmlabel = 'FM (z4-11)';
else
    fmlabel = 'FM (z4-11,14,15)';
end
% dataset = 'output/skymap/20140613s1';
dataset = 'output/nights/20140315s1';
% dataset = 'output/nights/20150801s1';

expIdList = dir(dataset);
%% count the number of exposures
nexp = 0;
for i = 1:size(expIdList,1)
    filename = sprintf('%s/%s/snr.txt',dataset,expIdList(i).name);
    if exist(filename, 'file')
        nexp = nexp + 1;
    end
end

cwfs = zeros(4, nexp, znmax-3);
fmzc = zeros(4, nexp, znmax-3);
iexp = 0;
for i = 1:size(expIdList,1)
    filename = sprintf('%s/%s/snr.txt',dataset,expIdList(i).name);
    if exist(filename, 'file')
        iexp = iexp + 1;
        for isenGrp=0:3
            filename=sprintf('%s/%s/ave_grp%d.txt',dataset,expIdList(i).name,isenGrp);
            if exist(filename, 'file')
                data = load(filename);
                zread = [0 0 0 data(1,:)];
                if znmax == 11
                    znew = z2z(zread,0.3396, 0 , znmax);
                    cwfs(isenGrp+1, iexp, :) = znew(4:end);
                    fmzc(isenGrp+1, iexp, :) =  data(2,1:znmax-3);
                elseif znmax == 15
                    znew = z2z(zread,0.3396, 0 , [1:11 14:15]);
                    cwfs(isenGrp+1, iexp, :) = znew(4:end);
                    fmzc(isenGrp+1, iexp, :) =  data(3,1:znmax-3);
                elseif znmax == 22
                    cwfs(isenGrp+1, iexp, :) = zread(4:end);
                    fmzc(isenGrp+1, iexp, :) =  data(3,1:znmax-3);
                end
            else
                cwfs(isenGrp+1, iexp, :) =  nan;
                fmzc(isenGrp+1, iexp, :) =  nan;
            end
        end
    end
end

if strcmp(plotStyle, 'diff')
    figure(1);clf;
    for isenGrp=0:3
        subplot(2,2,isenGrp+1);
        plot(4:znmax, squeeze(fmzc(isenGrp+1, :, :) - cwfs(isenGrp+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        xlabel('Zernike Number');ylabel('coefficient (in nm)');
        title(sprintf('%s - cwfs',fmlabel));
    end
elseif strcmp(plotStyle, 'ratio')
    figure(1);clf;
    for isenGrp=0:3
        subplot(2,2,isenGrp+1);
        plot(4:znmax, squeeze(fmzc(isenGrp+1, :, :)./cwfs(isenGrp+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title(sprintf('%s / cwfs',fmlabel));
    end
elseif strcmp(plotStyle, 'histogram')
    % fm - cwfs, by Zernike
    xbin = -500:50:500;
    figure(1);clf; %zc0
    for iz=4:znmax
        subplot(5,4,iz-3);
        aa=hist(reshape(squeeze(fmzc(:, :, iz-3) - cwfs(:, :, iz-3)),[], ...
            1),xbin);
        aa=aa/(sum(aa))*100;
        bar(xbin,aa);
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in %%)',iz));
    end
elseif strcmp(plotStyle, 'scatter')
    % scatter(fm, cwfs), by Zernike
    xlow = -1000;
    xhigh = 1000;
    figure(1);clf; %zc0
    for iz=4:znmax
        subplot(5,4,iz-3);
        scatter(reshape(squeeze(fmzc(:, :, iz-3)),[],1), reshape(cwfs(:, :, iz-3),[],1),50,'.');
        line([xlow, xhigh],[xlow, xhigh],'color','r');
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in nm)',iz));
        xlabel('FM'); ylabel('LSST');
    end
elseif strcmp(plotStyle, 'history')
    for iz=4:znmax
        plot_zn_history(1,iz, fmzc, cwfs);
        if iz==5
            fmzc = subtract_misalign(fmzc);
            cwfs = subtract_misalign(cwfs);
        end
        if (iz==5 || iz==6)
            plot_zn_history(2,iz, fmzc, cwfs);
        end
    end
    
elseif strcmp(plotStyle, 'history1')
    x=1:nexp;
    figure(1);clf; %zc0
    for iz=4:znmax
        subplot(3,4,iz-3);
        fm=mean(squeeze(fmzc(:, :, iz-3)));
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

function zc0new = subtract_misalign(zc0)
theta=[pi-atan(0.8); pi; atan(0.8); 0];
nexp = size(zc0,2);
zc0new = zc0;
options = optimset('Display', 'off');
for iexp=1:nexp
%     zc0(1,iexp,5-3) = 300;     zc0(1, iexp, 6-3) = -600;
%     zc0(2,iexp,5-3) = 1;       zc0(2, iexp, 6-3) = -700;
%     zc0(3,iexp,5-3) = 300;       zc0(3, iexp, 6-3) = 500;
%     zc0(4,iexp,5-3) = 1;       zc0(4, iexp, 6-3) = 700;
    
    ytan = zc0(:,iexp,5-3)./zc0(:, iexp, 6-3);
    idx = ~isnan(ytan);
    fun = @(phil, theta)tan(theta+phil);
    phil = lsqcurvefit(fun, 0, theta(idx),ytan(idx),[],[],options);
    yG2 = (zc0(:,iexp,6-3).^2+zc0(:, iexp, 5-3).^2);
    idx = ~isnan(yG2);
    fun = @(G, theta)(G^2+theta*0);
    G = lsqcurvefit(fun, 300, theta(idx), yG2(idx),[],[],options);
    zc0new(:, iexp, 5-3) = zc0(:, iexp, 5-3)-G*sin(theta+phil);
    zc0new(:, iexp, 6-3) = zc0(:, iexp, 6-3)-G*cos(theta+phil);
    
%     clf;
%     centerx = 1*cos(theta);
%     centery = 1*sin(theta);
%     plot(centerx, centery, '.','markersize',50);
%     ylim([-1.5 1.5]);axis square;grid;
%     
%     angle = atan(ytan);
%     angle(ytan<0)=angle(ytan<0)+pi;
%     halfx = sqrt(yG2)*2e-4.*cos(angle/2);
%     halfy = sqrt(yG2)*2e-4.*sin(angle/2);
%     for i=1:4
%         line([centerx(i)-halfx(i) centerx(i)+halfx(i)],[centery(i)-halfy(i) centery(i)+halfy(i)],'linewidth',5);
%     end
%     % phil=0;
%     halfx = (G)*2e-4.*cos((theta+phil)/2);
%     halfy = (G)*2e-4.*sin((theta+phil)/2);
%     for i=1:4
%         line([centerx(i)-halfx(i) centerx(i)+halfx(i)],[centery(i)-halfy(i) centery(i)+halfy(i)],'color','r','linewidth',2);
%     end
end
end

function [] = plot_zn_history(figureN, iz, zc0, cwfs)
sensorName = {'FN1/N2','FN3/N4','FS2/S1','FS4/S3'};
nexp = size(zc0,2);
x=1:nexp;
figure(figureN);clf; %zc0
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