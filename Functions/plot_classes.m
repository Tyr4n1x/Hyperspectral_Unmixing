function plot_classes(data)

    n_colors = length(unique(data(data~=0)));

%     my_colors = [255 255 255;
%                  100 200 255;
%                  255 140 0;
%                  255 200 0;
%                  255 255 0;
%                  50  200 50;
%                  0   128 0;
%                  120 255 0;
%                  160 40 40;
%                  128 0   128;
%                  255 100 70;
%                  255 160 120;
%                  255 0   0;
%                  255 20  150;
%                  0   0   128;
%                  0   0   0;
%                  170 170 170]/255;

    
    if any( find(data==0) ) % if contains at least one zero
        my_colors = distinguishable_colors(n_colors+1,'white');
        my_colors(1,:) = [255 255 255]/255; % white background
    else
        my_colors = distinguishable_colors(n_colors,'white');
    end
    
    figure('WindowState','maximized');
    colormap(my_colors)
    imagesc(data)
    c = colorbar('TickLabels',{'Alfalfa',...
                           'Corn-notill',...
                           'Corn-mintill',...
                           'Corn',...
                           'Grass-pasture',...
                           'Grass-trees',...
                           'Grass-pasture-mowed',...
                           'Hay-windrowed',...
                           'Oats',...
                           'Soybean-notill',...
                           'Soybean-mintill',...
                           'Soybean-clean',...
                           'Wheat',...
                           'Woods',...
                           'Buildings-Grass-Trees-Drives',...
                           'Stone-Steel-Towers'},...
              'FontSize',14);
    axis equal; axis off
    
    d = n_colors/( 2*(n_colors+1) );
    if any( find(data==0) ) % if contains at least one zero
        set(c,'Ticks',c.Limits(1)+3*d:2*d:c.Limits(2)-d)
    else
        set(c,'Ticks',c.Limits(1)+d:2*d:c.Limits(2)+d)
    end
    
end

