clear all, close all, clc
%% Load the data

addpath(genpath('./Data'))
addpath(genpath('./Images'))
addpath(genpath('./Functions'))

load('Indian_pines.mat');
load('Indian_pines_corrected.mat');
load('Indian_pines_gt.mat');

data = indian_pines_corrected;

%% False Color and RGB plot

wavelengths = linspace(0.4,2.5,size(indian_pines,3))*10^3; % [nm]
if data == indian_pines_corrected
    wavelengths([104:108,150:163,220]) = []; % [nm]
end

hcube = hypercube(data, wavelengths);
% hcube = hypercube( denoiseNGMeet(hcube.DataCube), hcube.Wavelength);
rgb = colorize(hcube,'Method','RGB');
[falsecolor,idx] = colorize(hcube,'Method','falsecolor');

figure('WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','Compact');
ax1 = nexttile;
imagesc(falsecolor)
axis off
title('False Color','FontSize',16)
ax2 = nexttile;
imagesc(rgb)
axis off
title('RGB','FontSize',16)
%% Extract endmembers

    %% Average knowing Ground Truth

% need to normalize indian_pines ?
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
    numEndmembers = countEndmembersHFC(hcube,'PFA',PFA);
    plot(PFA,numEndmembers,'r*')
end
xlabel('PFA')
xticks(10.^(-10:1:-1))
title('Number of endmembers','FontSize',14)
saveas(gcf,'./Images/Convergence_Endmembers.png')

numEndmembers = countEndmembersHFC(hcube,'PFA',10^-5);
%%
endmembers = ppi(hcube.DataCube,numEndmembers,'NumVectors',10^5,'ReductionMethod','MNF');
endmembers2 = nfindr(hcube.DataCube,numEndmembers,'NumIterations',3*numEndmembers,'ReductionMethod','MNF');
endmembers3 = fippi(hcube.DataCube,numEndmembers,'ReductionMethod','MNF');

% figure('WindowState','maximized');
% 
% subplot(2,2,1)
% plot(wavelengths, M)
% xlabel('Wavelength')
% xlim([min(wavelengths)-50, max(wavelengths)+50])
% title('Average knowing Ground Truth')
% 
% subplot(2,2,2)
% plot(wavelengths, endmembers)
% xlabel('Wavelength')
% xlim([min(wavelengths)-50, max(wavelengths)+50])
% title('PPI')
% 
% subplot(2,2,3)
% plot(wavelengths, endmembers2)
% xlabel('Wavelength')
% xlim([min(wavelengths)-50, max(wavelengths)+50])
% title('N-FINDR')
% 
% subplot(2,2,4)
% plot(wavelengths, endmembers3)
% xlabel('Wavelength')
% xlim([min(wavelengths)-50, max(wavelengths)+50])
% title('FIPPI')
% 
% sgtitle(['Number of Endmembers: ' num2str(p)])
% saveas(gcf,'./Images/Determination_Endmembers.png')

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

% abundanceMap = estimateAbundanceLS(hcube_c,endmembers);

    %% Pseudo-inverse (assume M is known)
    
[alphas,classification] = pseudo_inverse(data,M,indian_pines_gt,false);

C = confusionmat( reshape(indian_pines_gt,1,[]) , reshape(classification,1,[]) );

figure('WindowState','maximized');
confusionchart(C,0:16,'RowSummary','row-normalized','Title','Confusion Chart','FontSize',14);
saveas(gcf,'./Images/Pseudo_Inverse/ConfusionChart.png')

plot_classes(classification);
title('Classification','FontSize',14)
saveas(gcf,'./Images/Pseudo_Inverse/Classification.png')

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
saveas(gcf,'./Images/Pseudo_Inverse/AbundanceMap.png')

% not good, endmembers not spectrally independent (<-> assumption pseudo-inverse) + rounding !

    %% Optimum detection (assume U is known)

[alphas,classification] = optimum_detection(data,M,indian_pines_gt,false);

C = confusionmat( reshape(indian_pines_gt,1,[]) , reshape(classification,1,[]) );

figure('WindowState','maximized');
confusionchart(C,0:16,'RowSummary','row-normalized','Title','Confusion Chart','FontSize',14);
saveas(gcf,'./Images/Optimum_Detection/ConfusionChart.png')

plot_classes(classification);
title('Classification','FontSize',14)
saveas(gcf,'./Images/Optimum_Detection/Classification.png')

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
saveas(gcf,'./Images/Optimum_Detection/AbundanceMap.png')

% should be same as previous method

    %% Detection with unknown U

[alphas,classification] = unknownU(data,M,R,indian_pines_gt,false);

C = confusionmat( reshape(indian_pines_gt,1,[]) , reshape(classification,1,[]) );

figure('WindowState','maximized');
confusionchart(C,0:16,'RowSummary','row-normalized','Title','Confusion Chart','FontSize',14);
saveas(gcf,'./Images/UnknownU/ConfusionChart.png')

plot_classes(classification);
title('Classification','FontSize',14)
saveas(gcf,'./Images/UnknownU/Classification.png')

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
saveas(gcf,'./Images/UnknownU/AbundanceMap.png')

%%

counter = 0;
for i = 1:m
    for j = 1:n
        r = squeeze( alphas(i,j,:) ); % measurement at one pixel
        if sum(r==1) > 1
            counter = counter + 1;
            fprintf('%d: %d,%d \n',counter,i,j)
        end
    end
end