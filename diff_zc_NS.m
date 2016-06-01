function [] = diff_zc_NS(dataType, plotStyle, fieldcenter, znmax, writez2z)

% NS means the N sensors and S sensors
% the difference between this and diff_zc() is that 
% diff_zc() only deals with 4 pairs of WF sensors
% diff_zc_NS() deals with 31x2 pairs of sensors (each pair is made of two exposures)

% dataType = 's2/s5'
% plotStyle = 'diff', 'ratio', 'histogram', 'scatter', 'history'
% znmax = 11, or 15

if nargin<3
    fieldcenter=0;
end
centerNames = [2, 3, 4, 5, 6, 9, 10, 11, 12, 15, 16, 17, 18, 21, 22, 23];
if fieldcenter 
    nsensor=length(centerNames);
else
    nsensor = 31;
end

dataset = sprintf('output/camPiston/20140807%s', dataType);
nrow = 3;
ncol = 4;
if znmax==11
    fmlabel = 'FM (z4-11)';
else
    fmlabel = 'FM (z4-11,14,15)';
    nrow = 5;
    ncol = 2;
end

expIdList = dir(dataset);
nexp = 0;
for i = 1:size(expIdList,1)
    filename = sprintf('%s/%s/snr.txt',dataset,expIdList(i).name);
    if exist(filename, 'file')
        nexp = nexp + 1;
    end
end

cwfs = zeros(2, nsensor, nexp, znmax-3); %first 2 is for north and south sensors
fmzc = zeros(2,nsensor, nexp, znmax-3);
iexp = 0;
if fieldcenter
    z2zfile = sprintf('%s/z2z_%d_center.mat',dataset, znmax);    
else
    z2zfile = sprintf('%s/z2z_%d.mat',dataset, znmax);
end
if writez2z
    for i = 1:size(expIdList,1)
        filename = sprintf('%s/%s/snr.txt',dataset,expIdList(i).name);
        if exist(filename, 'file')
            iexp = iexp + 1;
            for isensor=0:nsensor-1
                if fieldcenter
                    isensorN = centerNames(isensor+1)-1;
                else
                    isensorN = isensor;
                end
                filename={sprintf(['%s/%s/' ...
                    'ave_N%d.txt'],dataset,expIdList(i).name,isensorN), ...
                    sprintf(['%s/%s/' ...
                    'ave_S%d.txt'],dataset,expIdList(i).name,isensorN)};
                for ins=1:2 %north or south sensor
                    if exist(filename{ins}, 'file')
                        data = load(filename{ins});
                        zread = [0 0 0 data(1,:)];
                        if znmax == 11
                            znew = z2z(zread,0.3396, 0 , znmax);
                            cwfs(ins, isensor+1, iexp, :)= znew(4:end);
                            fmzc(ins, isensor+1, iexp, :) =  data(2,1:znmax-3);
                        elseif znmax == 15
                            znew = z2z(zread,0.3396, 0 , [1:11 14:15]);
                            cwfs(ins, isensor+1, iexp, :) = znew(4:end);
                            fmzc(ins, isensor+1, iexp, :) =  data(3,1:znmax-3);
                        elseif znmax == 22
                            cwfs(ins, isensor+1, iexp, :) = zread(4:end);
                            fmzc(ins, isensor+1, iexp, :) =  data(3,1:znmax-3);
                        end
                    else
                        cwfs(ins, isensor+1, iexp, :) =  nan;
                        fmzc(ins, isensor+1, iexp, :) =  nan;
                    end
                end
            end
        end
    end
    save(z2zfile,'cwfs','fmzc');
else
    load(z2zfile,'cwfs','fmzc');
end

% close all;
if strcmp(plotStyle, 'diff')
    % fm - cwfs (N)
    figure(1);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(fmzc(1,isensor+1, :, :) - cwfs(1,isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        title(sprintf('N: %s - cwfs',fmlabel));
    end

    % fm - cwfs (S)
    figure(3);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(fmzc(2, isensor+1, :, :) - cwfs(2, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        title(sprintf('S: %s - cwfs',fmlabel));
    end

elseif strcmp(plotStyle, 'ratio')
    % fm/cwfs (N)
    figure(1);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(fmzc(1, isensor+1, :, :)./cwfs(1, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title(sprintf('N: %s / cwfs',fmlabel));
    end
 
    % fm/cwfs (S)
    figure(3);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(fmzc(2, isensor+1, :, :)./cwfs(2, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title(sprintf('S: %s / cwfs',fmlabel));
    end

elseif strcmp(plotStyle, 'histogram')
    % fm - cwfs, by Zernike
    xbin = -500:50:500;
    figure(1);clf;
    for iz=4:znmax
        subplot(5,4,iz-3);
        aa=hist(reshape(squeeze(fmzc(:, :, :, iz-3) - cwfs(:, :, :, iz-3)),[], ...
            1),xbin);
        aa=aa/(sum(aa))*100;
        bar(xbin,aa);
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in %%) %s',iz, fmlabel));
    end

elseif strcmp(plotStyle, 'scatter')
    % scatter(fm, cwfs), by Zernike
    xlow = -1000;
    xhigh = 1000;
    figure(1);clf; 
    for iz=4:znmax
        isub = iz-3;
        if iz==12 || iz==13
           continue;
        elseif iz==14 || iz==15
            isub = isub-2;
        end
        subplot(nrow,ncol,isub);
        scatter(reshape(squeeze(fmzc(:, :, :, iz-3)),[],1), reshape(cwfs(:, :, :, iz-3),[],1),180,'.');
        line([xlow, xhigh],[xlow, xhigh],'color','r');
        xlim([-1000 1000]);
        ylim([-1000 1000]);
        grid on;
        text(0.1,0.9,sprintf('z%d',iz),'unit','normalized','fontweight','bold','fontsize',12);
        if mod(isub,ncol)==1
            ylabel('LSST (in nm)');
        end
        if iz>=znmax-ncol+1
            xlabel('DECam (in nm)');
        end
    end
    samexyaxis('xmt','on','ytac','join','yld',1);

end




end
