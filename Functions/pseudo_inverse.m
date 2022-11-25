function [alphas,classification] = pseudo_inverse(data,M,ground_truth,do_plot)
%PSEUDO_INVERSE Calculates the abundance of the different endmembers
%using the method of the pseudo-inverse and create a confusion matrix based on
%the maximal abundance per pixel.
    

    m = size(data,1);
    n = size(data,2);
    L = size(data,3);
    p = size(M,2);
    
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
    
    if do_plot
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
            imagesc(ground_truth==endmember)
            colormap(ax2,[1,1,1; 0.9769,0.9839,0.0805]) % to get two colors
            axis off
            title('Ground Truth')

            c = colorbar(ax1,'FontSize',12);
            c.Layout.Tile = 'west';
            ylabel(c,'\alpha', 'FontSize',20, 'Rotation',0); c.Label.Position(1) = c.Position(2)-c.Position(4)*4;

            saveas(gcf,'./Images/Pseudo_Inverse/Abundance_endmember_' + string(endmember)  + '.png')
        end
    end
end
