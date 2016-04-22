function [] = diff_zc_NS(dataType, plotStyle, fieldcenter)

% NS means the N sensors and S sensors
% the difference between this and diff_zc() is that 
% diff_zc() only deals with 4 pairs of WF sensors
% diff_zc_NS() deals with 31x2 pairs of sensors (each pair is made of two exposures)

% dataType = 's2/s5'
% plotStyle = 'diff', 'ratio', 'histogram', 'scatter', 'history'

if nargin<3
    fieldcenter=0;
end
centerNames = [3, 4, 5, 10, 11, 16, 17];
if fieldcenter 
    nsensor=7;
else
    nsensor = 31;
end

dataset = sprintf('output/camPiston/20140807%s', dataType);

expIdList = dir(dataset);
nexp = 0;
for i = 1:size(expIdList,1)
    filename = sprintf('%s/%s/snr.txt',dataset,expIdList(i).name);
    if exist(filename, 'file')
        nexp = nexp + 1;
    end
end

znmax = 22;
cwfs = zeros(2, nsensor, nexp, znmax-3);
zc0 = zeros(2, nsensor, nexp, znmax-3);
zc1 = zeros(2, nsensor, nexp, znmax-3);
iexp = 0;
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
            for ins=1:2
                if exist(filename{ins}, 'file')
                    data = load(filename{ins});
                    cwfs(ins, isensor+1, iexp, :) =  data(1,:);
                    zc0(ins, isensor+1, iexp, :) =  data(2,:);
                    zc1(ins, isensor+1, iexp, :) =  data(3,:);
                else
                    cwfs(ins, isensor+1, iexp, :) =  nan;
                    zc0(ins, isensor+1, iexp, :) =  nan;
                    zc1(ins, isensor+1, iexp, :) =  nan;                    
                end
            end
        end
    end
end

% close all;
if strcmp(plotStyle, 'diff')
    % first - cwfs (N)
    figure(1);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(zc0(1,isensor+1, :, :) - cwfs(1,isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        title('N: FM (z4-11) - cwfs');
    end
    % second - cwfs (N)
    figure(2);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(zc1(1, isensor+1, :, :) - cwfs(1, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        title('N: FM (z4-11, 14,15) - cwfs');
    end
    % first - cwfs (S)
    figure(3);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(zc0(2, isensor+1, :, :) - cwfs(2, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        title('S: FM (z4-11) - cwfs');
    end
    % second - cwfs (N)
    figure(4);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(zc1(2, isensor+1, :, :) - cwfs(2, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        grid on;
        title('S: FM (z4-11, 14,15) - cwfs');
    end
elseif strcmp(plotStyle, 'ratio')
    % first/cwfs (N)
    figure(1);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(zc0(1, isensor+1, :, :)./cwfs(1, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title('N: FM (z4-11) / cwfs');
    end
    
    % second/cwfs (N)
    figure(2);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(zc1(1, isensor+1, :, :)./cwfs(1, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title('N: FM (z4-11, 14,15) / cwfs');
    end
    % first/cwfs (S)
    figure(3);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(zc0(2, isensor+1, :, :)./cwfs(2, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title('S: FM (z4-11) / cwfs');
    end
    
    % second/cwfs (S)
    figure(4);clf;
    for isensor=0:nsensor-1
        subplot(6,6,isensor+1);
        plot(4:znmax, squeeze(zc1(2, isensor+1, :, :)./cwfs(2, isensor+1, :, :))','o');
        xlim([4-0.5 znmax+0.5]);
        ylim([-10 10]);
        grid on;
        title('S: FM (z4-11, 14,15) / cwfs');
    end
elseif strcmp(plotStyle, 'histogram')
    % first - cwfs, by Zernike
    xbin = -500:50:500;
    figure(1);clf;
    for iz=4:znmax
        subplot(5,4,iz-3);
        aa=hist(reshape(squeeze(zc0(:, :, :, iz-3) - cwfs(:, :, :, iz-3)),[], ...
            1),xbin);
        aa=aa/(sum(aa))*100;
        bar(xbin,aa);
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in %%)',iz));
    end
    % second - cwfs, by Zernike
    xbin = -500:50:500;
    figure(2);clf;
    for iz=4:znmax
        subplot(5,4,iz-3);
        aa=hist(reshape(squeeze(zc1(:, :, :, iz-3) - cwfs(:, :, :, iz-3)),[], ...
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
        scatter(reshape(squeeze(zc0(:, :, :, iz-3)),[],1), reshape(cwfs(:, :, :, iz-3),[],1),180,'.');
        line([xlow, xhigh],[xlow, xhigh],'color','r');
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in nm)',iz));
        xlabel('FM'); ylabel('LSST');
    end
    figure(2);clf; %zc1
    for iz=4:znmax
        subplot(5,4,iz-3);
        scatter(reshape(squeeze(zc1(:, :, :, iz-3)),[],1), reshape(cwfs(:, :, :, iz-3),[],1),180,'.');
        line([xlow, xhigh],[xlow, xhigh],'color','r');
        xlim([-1000 1000]);
        grid on;
        title(sprintf('z%d (in nm)',iz));
        xlabel('FM'); ylabel('LSST');
    end
end




end
