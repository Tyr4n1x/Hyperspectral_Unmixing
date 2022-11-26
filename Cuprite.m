clear all, close all, clc
%% Load the data

addpath(genpath('./Data'))
addpath(genpath('./Images'))
addpath(genpath('./Functions'))

load('Cuprite.mat')

data = double(X);
data(:,:,[1:2,221:224,104:113,148:167]) = []; % [nm]

%% False Color and RGB plot

wavelengths = linspace(0.4,2.5,size(X,3))*10^3; % [nm]
wavelengths([1:2,221:224,104:113,148:167]) = []; % [nm]

hcube = hypercube(data, wavelengths);
% hcube = hypercube( denoiseNGMeet(hcube.DataCube), hcube.Wavelength);
rgb = colorize(hcube,'Method','RGB');
[falsecolor,idx] = colorize(hcube,'Method','falsecolor');

figure('WindowState','maximized');
tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
ax1 = nexttile;
imagesc(falsecolor)
axis off
title('False Color','FontSize',16)
ax2 = nexttile;
imagesc(rgb)
axis off
title('RGB','FontSize',16)

%% Extract endmembers

m = size(data,1);
n = size(data,2);
L = size(data,3); % number of bands
p = 12;

%%
endmembers = ppi(hcube.DataCube,p,'NumVectors',10^5,'ReductionMethod','MNF');
endmembers2 = nfindr(hcube.DataCube,p,'NumIterations',3*numEndmembers,'ReductionMethod','MNF');
endmembers3 = fippi(hcube.DataCube,p,'ReductionMethod','MNF');

figure('WindowState','maximized');
t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');

nexttile;
plot(wavelengths, endmembers)
xlim([min(wavelengths)-50, max(wavelengths)+50])
xlabel('Wavelength')
title('PPI','FontSize',14)

nexttile;
plot(wavelengths, endmembers2)
xlim([min(wavelengths)-50, max(wavelengths)+50])
xlabel('Wavelength')
title('N-FINDR','FontSize',14)

nexttile;
plot(wavelengths, endmembers3)
xlim([min(wavelengths)-50, max(wavelengths)+50])
xlabel('Wavelength')
title('FIPPI','FontSize',14)

linkaxes(t.Children,'xy')

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

M = endmembers;

    %% Pseudo-inverse (assume M is known)
    
alphas = zeros(m,n,p);
M_pseudo_inverse = inv(M'*M)*M';

for i = 1:m
    for j = 1:n
        r = squeeze( data(i,j,:) ); % measurement at one pixel
        alpha = M_pseudo_inverse*r;
        alphas(i,j,:) = alpha;
    end
end
    
[~,classification] = max(alphas,[],3);
    
alphas = max(0,min(1,alphas)); % 0 <= alpha <= 1
for endmember = 1:p
    figure();
    t = tiledlayout(1,1,'TileSpacing','Compact');
    title(t,'Endmember ' + string(endmember),'VerticalAlignment', 'bottom', 'FontSize',16)

    ax1 = nexttile;
    imagesc(alphas(:,:,endmember))
    colormap(ax1,'default')
    axis off
    title('Abundance Map')

    c = colorbar(ax1,'FontSize',12);
    c.Layout.Tile = 'west';
    ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)-c.Position(4)*4;
end

plot_classes(classification);
title('Classification','FontSize',14)

% not good, endmembers not spectrally independent (<-> assumption pseudo-inverse) + rounding !

    %% Optimum detection (assume U is known)

[alphas,classification] = optimum_detection(data,M,indian_pines_gt,false);

C = confusionmat( reshape(indian_pines_gt,1,[]) , reshape(classification,1,[]) );

figure('WindowState','maximized');
confusionchart(C,0:16,'RowSummary','row-normalized','Title','Confusion Chart','FontSize',14);

plot_classes(classification);
title('Classification','FontSize',14)

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)

% should be same as previous method

    %% Detection with unknown U

[alphas,classification] = unknownU(data,M,R,indian_pines_gt,false);

C = confusionmat( reshape(indian_pines_gt,1,[]) , reshape(classification,1,[]) );

figure('WindowState','maximized');
confusionchart(C,0:16,'RowSummary','row-normalized','Title','Confusion Chart','FontSize',14);

plot_classes(classification);
title('Classification','FontSize',14)

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)

    %% Matlab method
    
abundanceMap = estimateAbundanceLS(hcube,M,'Method','FCLS');
figure()
montage(abundanceMap,'Size',[4 4],'BorderSize',[10 10]);
colormap default

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