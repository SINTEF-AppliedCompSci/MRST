function [ hfig, hax ] = plot3DandTopGrids( G, Gt )

% Inspect orientation of grid (compare with orientation shown in Singh
    % et al 2010, Cavanagh 2013, etc.)
    figure; set(gcf,'Position',[1 1 1000 1000])
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(3)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x'); ylabel('y'); zlabel('z');
    % Create textarrow
    view(44,22)
    annotation(gcf,'textarrow',[0.4734 0.5391],[0.7825 0.81],'String',{'North'});
    set(gca,'FontSize',14)
    
    
    % Inspect orientation of grid: x-sectional views
    figure; set(gcf,'Position',[1 1 1100 500])
    % facing 
    subplot(2,2,1)
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))                                 % ,'EdgeColor','none')
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(270,0)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Facing East')
    
    % facing 
    subplot(2,2,3)
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(90,0)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Facing West')
    
    % facing 
    subplot(2,2,2)
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(180,0)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Facing South')
    
    % facing 
    subplot(2,2,4)
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(0,0)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Facing North')
    
    
    % Inspect Top Surface Grid and H
    figure; set(gcf,'Position',[1 1 1000 1000])
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r');
    plotCellData(Gt, Gt.cells.z)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title({'Top Surface elevation';['average thickness ',num2str(mean(Gt.cells.H)),' meters']})
    % Create textarrow
    view(44,22)
    annotation(gcf,'textarrow',[0.4734 0.5391],[0.7825 0.81],'String',{'North'});
    set(gca,'FontSize',14)
    
    hfig = gcf;
    hax  = gca;

end

