clear all, close all, clc
%% Extract the data

% 0.4-2.5 µm
load('Indian_pines.mat');
load('Indian_pines_corrected.mat');

wavelengths = linspace(0.4,2.5,size(indian_pines,3))*10^3; % [nm]
wavelengths_corrected = wavelengths; wavelengths_corrected([104:108,150:163,220]) = []; % [nm]

%%

figure(); set(gcf, 'Position',  [400, 100, 1000, 800])
counter = 0;
for i = linspace(1,size(indian_pines,3),4)
    counter = counter + 1;
    subplot(2, 2, counter)
    imagesc(indian_pines(:,:,i)./max(indian_pines(:,:,i))*255)
    hold on
    plot(40,100,'ro', 'LineWidth',2)
    title( sprintf('$\\lambda = %d nm$', round(wavelengths(i)) ),'Interpreter', 'Latex')
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
end

figure();
plot(wavelengths, squeeze(indian_pines(40,100,:)./max(indian_pines(40,100,:))))
xlabel('$\lambda [nm]$','Interpreter','Latex')
ylabel('Normalized Intensity')
title('Spectrum of pixel (40,100)')
xlim([min(wavelengths)-50, max(wavelengths)+50])

%% 
hcube = hypercube(indian_pines, wavelengths);
hcube_denoise = hypercube( denoiseNGMeet(hcube.DataCube), hcube.Wavelength);

rgb = colorize(hcube,'Method','RGB');
rgb_denoise = colorize(hcube_denoise,'Method','RGB');

figure();
subplot(1,2,1); imagesc(rgb)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
subplot(1,2,2); imagesc(rgb_denoise)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])

%%
% hyperspectralViewer(hcube);

%% 
hcube_c = hypercube(indian_pines_corrected, wavelengths_corrected);

rgb_c = colorize(hcube_c,'Method','RGB');

figure();
subplot(1,2,1); imagesc(rgb)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
subplot(1,2,2); imagesc(rgb_c)
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])

%%

numEndmembers = countEndmembersHFC(hcube_c,'PFA',10^-7);
endmembers = ppi(hcube_c.DataCube,numEndmembers);

figure();
plot(endmembers)
xlabel('Band Number (Corrected)')
ylabel('Pixel Values')
title(['Number of Endmembers: ' num2str(size(endmembers,2))])

%%
abundanceMap = estimateAbundanceLS(hcube_c,endmembers);

figure()
montage(abundanceMap,'Size',[4 4],'BorderSize',[20 20])
colormap default
colorbar
title('Abundance Map','FontSize',14)

%%
load('Indian_pines_gt.mat');

figure();
imagesc(indian_pines_gt)
title('Ground Thruth')
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])