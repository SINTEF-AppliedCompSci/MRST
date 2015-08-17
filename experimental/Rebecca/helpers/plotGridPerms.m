function [ hfig ] = plotGridPerms( G, rock )

% G and rock could be either 3D or 2D (top-surface) grids and rock
% properties. However, if 2D, only permX = permY is plotted.


    [~,n] = size(rock.perm);


    figure; 

    subplot(1,n,1)
    plotCellData(G, rock.perm(:,1), 'EdgeColor', 'none')
    view(3); colorbar; title('Perm X')
    axis equal
    set(gca,'DataAspect',[1 1 0.02]); grid
    
    if n==3

        subplot(1,n,2)
        plotCellData(G, rock.perm(:,2), 'EdgeColor', 'none')
        view(3); colorbar; title('Perm Y')
        axis equal
        set(gca,'DataAspect',[1 1 0.02]); grid

        subplot(1,n,3)
        plotCellData(G, rock.perm(:,3), 'EdgeColor', 'none')
        view(3); colorbar; title('Perm Z')
        axis equal
        set(gca,'DataAspect',[1 1 0.02]); grid
        
        set(gcf,'Position',[1 1 3000 500])
        
    else
        title('Perm X and Y')
        set(gcf,'Position',[1 1 800 500])

    end


    hfig = gcf;

end

