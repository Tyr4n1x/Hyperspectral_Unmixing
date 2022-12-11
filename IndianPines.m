clear all, close all, clc
%% Load the data

addpath(genpath('./Data'))
addpath(genpath('./Images'))
addpath(genpath('./Functions'))

load('IndianPines.mat'); % default from Hyperspectral Toolbox
load('IndianPines_corrected.mat');
load('Indian_pines_gt.mat');

data = indian_pines_corrected;

%% False Color and RGB plot

wavelengths = linspace(0.4,2.5,size(indian_pines,3))*10^3; % [nm]
if data == indian_pines_corrected
    wavelengths([104:108,150:163,220]) = []; % [nm]
    endmembers([104:108,150:163,220],:) = [];
end

hcube = hypercube(data, wavelengths);
rgb = colorize(hcube,'Method','RGB');
[falsecolor,idx] = colorize(hcube,'Method','falsecolor');

figure('WindowState','maximized');
tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
ax1 = nexttile;
imagesc(falsecolor)
axis off
title('False Color','FontSize',20)
ax2 = nexttile;
imagesc(rgb)
axis off
title('RGB','FontSize',20)

exportgraphics(gcf,'./Images/RGB-False_Color.png')

%% Extract endmembers

    %% Average knowing Ground Truth

m = size(data,1);
n = size(data,2);
L = size(data,3); % number of bands
p = length(unique(indian_pines_gt(indian_pines_gt~=0))); % number of endmembers

M = zeros(L,p);
freq = zeros(1,p);
for i = 1:m
    for j = 1:n
        k = indian_pines_gt(i,j); % which endmember is in the pixel
        if k>0 % if k = 0 --> no endmember detected
            r = squeeze( data(i,j,:) ); % measurement at one pixel
            M(:,k) = M(:,k) + r; % add the endmember
            freq(k) = freq(k) + 1; % count the number of pixels per endmember
        end
    end
end
M = M./freq; % to take average

    %% Count Endmembers HFC

figure(); set(gca,'XScale','log','XDir','reverse'); hold on
ylim([15 34])
for PFA = 10.^(-10:1:-1)
    numEndmembers = countEndmembersHFC(hcube,'PFA',PFA,'NoiseWhiten',true);
    plot(PFA,numEndmembers,'r*')
end
xlabel('PFA')
xticks(10.^(-10:1:-1))
xlim([10^-11 10^0])
title('Number of endmembers','FontSize',16)
exportgraphics(gcf,'./Images/Convergence_Endmembers.png')

numEndmembers = countEndmembersHFC(hcube,'PFA',10^-5);

    %% Comparison of the different methods to the reference
    
endmembers1 = ppi(hcube.DataCube,numEndmembers,'NumVectors',10^4,'ReductionMethod','MNF');
endmembers2 = nfindr(hcube.DataCube,numEndmembers,'NumIterations',3*numEndmembers,'ReductionMethod','MNF');
endmembers3 = fippi(hcube.DataCube,numEndmembers,'ReductionMethod','MNF');

figure(); imagesc( corr(endmembers) );
ax = gca; ax.YDir = 'normal'; ax.FontSize = 12;
colorbar;
title('Autocorrelation of reference','FontSize',16)

exportgraphics(gcf,'./Images/Correlation_Endmembers.png')

figure(); imagesc( corr(endmembers,M) );
ax = gca; ax.YDir = 'normal'; ax.FontSize = 12;
colorbar;
title('Average knowing ground truth','FontSize',16)

exportgraphics(gcf,'./Images/Correlation_Endmembers1.png')

figure(); imagesc( corr(endmembers,endmembers1) );
ax = gca; ax.YDir = 'normal'; ax.FontSize = 12;
colorbar;
title('PPI','FontSize',16)

exportgraphics(gcf,'./Images/Correlation_Endmembers2.png')

figure(); imagesc( corr(endmembers,endmembers2) );
ax = gca; ax.YDir = 'normal'; ax.FontSize = 12;
colorbar;
title('N-FINDR','FontSize',16)

exportgraphics(gcf,'./Images/Correlation_Endmembers3.png')

figure(); imagesc( corr(endmembers,endmembers3) );
ax = gca; ax.YDir = 'normal'; ax.FontSize = 12;
colorbar;
title('FPPI','FontSize',16)

exportgraphics(gcf,'./Images/Correlation_Endmembers4.png')

    %%

figure('WindowState','maximized');
t = tiledlayout(4,4,'TileSpacing','Compact','Padding','Compact');

for tile = 1:numEndmembers
    nexttile; hold on
    plot(wavelengths, endmembers(:,tile),'LineWidth',2.5)
    plot(wavelengths, M(:,tile))
    plot(wavelengths, endmembers1(:,tile))
    plot(wavelengths, endmembers2(:,tile))
    plot(wavelengths, endmembers3(:,tile))
    xlim([min(wavelengths)-50, max(wavelengths)+50])
    xlabel('Wavelength')
    title(sprintf('Endmember = %d',tile),'FontSize',14)
end
l = legend('Reference','Ground Truth','PPI','N-FINDR','FPPI','Fontsize',16);
l.Layout.Tile = 'East';
linkaxes(t.Children,'xy')

exportgraphics(gcf,'./Images/Determination_Endmembers.png')

%% Calculate covariance matrix

K = m*n; % total number of pixels per band
R = zeros(L, L);
for i = 1:m
    for j = 1:n
        r = squeeze( data(i,j,:) ); % measurement at one pixel
        R = R + r*r';
    end
end
R = R./K;

%% Calculate the abundance for every endmember

    %% Pseudo-inverse (assume M is known)
    
[alphas,classification] = pseudo_inverse(data,M,indian_pines_gt,false);

C = confusionmat( reshape(indian_pines_gt,1,[]) , reshape(classification,1,[]) );

figure('WindowState','maximized');
confusionchart(C,0:16,'RowSummary','row-normalized','Title','Confusion Chart','FontSize',20);
exportgraphics(gcf,'./Images/Pseudo_Inverse/ConfusionChart.png')

plot_classes(classification);
title('Classification','FontSize',20)
exportgraphics(gcf,'./Images/Pseudo_Inverse/Classification.png')

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
exportgraphics(gcf,'./Images/Pseudo_Inverse/AbundanceMap.png')

fprintf("Determinant of M'*M = ")
disp( det(M'*M) )

% not good, endmembers not spectrally independent (<-> assumption pseudo-inverse) + rounding !

    %% Optimum detection (assume U is known)

[alphas,classification] = optimum_detection(data,M,indian_pines_gt,false);

C = confusionmat( reshape(indian_pines_gt,1,[]) , reshape(classification,1,[]) );

figure('WindowState','maximized');
confusionchart(C,0:16,'RowSummary','row-normalized','Title','Confusion Chart','FontSize',20);
exportgraphics(gcf,'./Images/Optimum_Detection/ConfusionChart.png')

plot_classes(classification);
title('Classification','FontSize',20)
exportgraphics(gcf,'./Images/Optimum_Detection/Classification.png')

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
exportgraphics(gcf,'./Images/Optimum_Detection/AbundanceMap.png')

% should be same as previous method

    %% Detection with unknown U

[alphas,classification] = unknownU(data,M,R,indian_pines_gt,false);

C = confusionmat( reshape(indian_pines_gt,1,[]) , reshape(classification,1,[]) );

figure('WindowState','maximized'); tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile; confusionchart(C,0:16,'RowSummary','row-normalized','Title','Confusion Chart','FontSize',20);
exportgraphics(gcf,'./Images/UnknownU/ConfusionChart.png')

plot_classes(classification);
title('Classification','FontSize',20)
exportgraphics(gcf,'./Images/UnknownU/Classification.png')

figure(); tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
exportgraphics(gcf,'./Images/UnknownU/AbundanceMap.png')
