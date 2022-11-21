clear all, close all, clc
%% Load the data

addpath(genpath('./Data'))
addpath(genpath('./Images'))

load('Indian_pines.mat');
load('Indian_pines_corrected.mat');
load('Indian_pines_gt.mat');

%% Extract endmembers

% need to normalize indian_pines ?
m = size(indian_pines,1);
n = size(indian_pines,2);
L = size(indian_pines,3); % number of bands
p = length(unique(indian_pines_gt(indian_pines_gt~=0))); % number of endmembers

M = zeros(L,p);
freq = zeros(1,p);
for i = 1:m
    for j = 1:n
        k = indian_pines_gt(i,j); % which endmember is in the pixel
        if k>0 % if k = 0 --> no endmember detected
            r = squeeze( indian_pines(i,j,:) ); % measurement at one pixel
            M(:,k) = M(:,k) + r; % add the endmember
            freq(k) = freq(k) + 1; % count the number of pixels per endmember
        end
    end
end
M = M./freq; % to take average

wavelengths = linspace(0.4,2.5,L)*10^3; % [nm]

hcube = hypercube(indian_pines, wavelengths);
numEndmembers = countEndmembersHFC(hcube,'PFA',10^-10); % asymptotic behavior
endmembers = ppi(hcube.DataCube,numEndmembers);
endmembers2 = nfindr(hcube.DataCube,numEndmembers);
endmembers3 = fippi(hcube.DataCube,numEndmembers);

figure('WindowState','maximized');

subplot(2,2,1)
plot(wavelengths, M)
xlabel('Wavelength')
xlim([min(wavelengths)-50, max(wavelengths)+50])
title('Average knowing Ground Truth')

subplot(2,2,2)
plot(wavelengths, endmembers)
xlabel('Wavelength')
xlim([min(wavelengths)-50, max(wavelengths)+50])
title('PPI')

subplot(2,2,3)
plot(wavelengths, endmembers2)
xlabel('Wavelength')
xlim([min(wavelengths)-50, max(wavelengths)+50])
title('N-FINDR')

subplot(2,2,4)
plot(wavelengths, endmembers3)
xlabel('Wavelength')
xlim([min(wavelengths)-50, max(wavelengths)+50])
title('FIPPI')

sgtitle(['Number of Endmembers: ' num2str(p)])
saveas(gcf,'./Images/Determination_Endmembers.png')

%% Calculate covariance matrix

K = m*n; % total number of pixels per band
R = zeros(L, L);
for i = 1:m
    for j = 1:n
        r = squeeze( indian_pines(i,j,:) ); % measurement at one pixel
        R = R + r*r';
    end
end
R = R./K;

%% Calculate the abundance for every endmember

    %% Pseudo-inverse (assume M is known)
    
alphas = zeros(m,n,p);
M_pseudo_inverse = inv(M'*M)*M';

for i = 1:m
    for j = 1:n
        r = squeeze( indian_pines(i,j,:) ); % measurement at one pixel
        alpha = M_pseudo_inverse*r;
        alpha = max(0,min(1,alpha)); % 0 <= alpha <= 1
        
        alphas(i,j,:) = alpha;
    end
end

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
saveas(gcf,'./Images/Pseudo_Inverse/AbundanceMap.png')


for endmember = 1:p
    figure();
    t = tiledlayout(1,2,'TileSpacing','Compact');
    title(t,'Endmember ' + string(endmember),'VerticalAlignment', 'bottom', 'FontSize',16)
    
    ax1 = nexttile;
    imagesc(alphas(:,:,endmember))
    colormap(ax1,'default')
    axis off
    title('Abundance Map')
    
    ax2 = nexttile;
    imagesc(indian_pines_gt==endmember)
    colormap(ax2,[1,1,1; 0.9769,0.9839,0.0805]) % to get two colors
    axis off
    title('Ground Truth')
    
    c = colorbar(ax1,'FontSize',12);
    c.Layout.Tile = 'west';
    ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)-c.Position(4)*4;
    
    saveas(gcf,'./Images/Pseudo_Inverse/Abundance_endmember_' + string(endmember)  + '.png')
end

% not good, endmembers not spectrally independent (<-> assumption pseudo-inverse) + rounding !

    %% Optimum detection (assume U is known)

alphas = zeros(m,n,p);

for endmember = 1:p
    d = M(:,endmember);
    U = M; U(:,endmember)=[];
    
    P_u = eye(L) - U*inv(U'*U)*U';

    for i = 1:m
        for j = 1:n
            r = squeeze( indian_pines(i,j,:) ); % measurement at one pixel
            alpha = (d'*P_u*r)/(d'*P_u*d);
            alpha = max(0,min(1,alpha)); % 0 <= alpha <= 1

            alphas(i,j,endmember) = alpha;
        end
    end
    
    figure();
    t = tiledlayout(1,2,'TileSpacing','Compact');
    title(t,'Endmember ' + string(endmember),'VerticalAlignment', 'bottom', 'FontSize',16)
    
    ax1 = nexttile;
    imagesc(alphas(:,:,endmember))
    colormap(ax1,'default')
    axis off
    title('Abundance Map')
    
    ax2 = nexttile;
    imagesc(indian_pines_gt==endmember)
    colormap(ax2,[1,1,1; 0.9769,0.9839,0.0805]) % to get two colors
    axis off
    title('Ground Truth')
    
    c = colorbar(ax1,'FontSize',12);
    c.Layout.Tile = 'west';
    ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)-c.Position(4)*4;
    
    saveas(gcf,'./Images/Optimum_Detection/Abundance_endmember_' + string(endmember)  + '.png')
end

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
saveas(gcf,'./Images/Optimum_Detection/AbundanceMap.png')

    %% Detection with unknown U
    
alphas = zeros(m,n,p);

for endmember = 1:p
    d = M(:,endmember);
    w = ( inv(R)*d ) / (d'*inv(R)*d); % calculate the appropriate filter

    for i = 1:m
        for j = 1:n
            r = squeeze( indian_pines(i,j,:) ); % measurement at one pixel
            alpha = w'*r;
            alpha = max(0,min(1,alpha)); % 0 <= alpha <= 1

            alphas(i,j,endmember) = alpha;
        end
    end
    
    figure();
    t = tiledlayout(1,2,'TileSpacing','Compact');
    title(t,'Endmember ' + string(endmember),'VerticalAlignment', 'bottom', 'FontSize',16)
    
    ax1 = nexttile;
    imagesc(alphas(:,:,endmember))
    colormap(ax1,'default')
    axis off
    title('Abundance Map')
    
    ax2 = nexttile;
    imagesc(indian_pines_gt==endmember)
    colormap(ax2,[1,1,1; 0.9769,0.9839,0.0805]) % to get two colors
    axis off
    title('Ground Truth')
    
    c = colorbar(ax1,'FontSize',12);
    c.Layout.Tile = 'west';
    ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)-c.Position(4)*4;
    
    saveas(gcf,'./Images/UnknownU/Abundance_endmember_' + string(endmember)  + '.png')
end

figure();
montage(alphas,'Size',[4 4],'BorderSize',[20 20])
colormap default
c = colorbar('FontSize',14);
ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)+c.Position(4)*5;
title('Abundance Map','FontSize',14)
saveas(gcf,'./Images/UnknownU/AbundanceMap.png')
