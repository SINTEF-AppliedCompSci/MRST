function [plumes,topsurface, topfit, hCO2] = makeSurfaceDataAndPlots(plumes,Gt)
if(isfield(Gt,'cartDims') && Gt.cells.num==prod(Gt.cartDims))
X=reshape(Gt.cells.centroids(:,1),Gt.cartDims(1),Gt.cartDims(2));
Y=reshape(Gt.cells.centroids(:,2),Gt.cartDims(1),Gt.cartDims(2));
Z=reshape(Gt.cells.z,Gt.cartDims(1),Gt.cartDims(2));
topsurface=@(coord) interp2(X',Y',Z',coord(:,1),coord(:,2));
else
    F = scatteredInterpolant(Gt.cells.centroids(:,1),Gt.cells.centroids(:,2),Gt.cells.z);
    topsurface =@(coord) F(coord(:,1),coord(:,2));
end
for i=1:numel(plumes)
    line_coord=plumes{i}.outline;
    %%{
    nd=size(line_coord,1);
    A=[ones(nd,1),line_coord];
    rhs=topsurface(line_coord);
    vec=A\rhs;
    topfit=@(coord) [ones(size(coord,1),1),coord]*vec;
    
    figure; set(gcf,'Position',[1 1 1300 500])
    subplot(1,2,1)
    
        plumeInside_coord_tmp = [Gt.cells.centroids(:,1).*insidePolygon(line_coord, Gt.cells.centroids), ...
        Gt.cells.centroids(:,2).*insidePolygon(line_coord, Gt.cells.centroids)];
        tmp1 = plumeInside_coord_tmp(:,1);
        tmp2 = plumeInside_coord_tmp(:,2);
        plumeInside_coord = [tmp1(plumeInside_coord_tmp(:,1)~=0), tmp2(plumeInside_coord_tmp(:,2)~=0)];
        
        
        hold on
        plot3(plumeInside_coord(:,1), plumeInside_coord(:,2), topsurface(plumeInside_coord), 'xr')
        plot3(plumeInside_coord(:,1), plumeInside_coord(:,2), topfit(plumeInside_coord), 'sqb')
        plot3(line_coord(:,1), line_coord(:,2), topsurface(line_coord), 'r', 'LineWidth',2)
        plot3(line_coord(:,1), line_coord(:,2), topfit(line_coord), 'b', 'LineWidth',2)
        view(3)
        
        title({['Elevations of (x,y) inside plume as determined by'];['functions (wrt ',num2str(plumes{i}.year),' plume).']})
        legend('topsurface(x,y)','topfit(x,y)','topsurface(x_p,y_p)','topfit(x_p,y_p)','Location','NorthEast')
        axis equal tight;
        set(gca,'DataAspect',[1 1 1/50])
        grid
        zlabel('depth (meters)')
    
      % whole grid
%     hold on
%     plot3(Gt.cells.centroids(:,1), Gt.cells.centroids(:,2), topsurface(Gt.cells.centroids), 'xr')
%     plot3(Gt.cells.centroids(:,1), Gt.cells.centroids(:,2), topfit(Gt.cells.centroids), 'sqb')
%     plot3(line_coord(:,1), line_coord(:,2), topsurface(line_coord), 'r', 'LineWidth',2)
%     plot3(line_coord(:,1), line_coord(:,2), topfit(line_coord), 'b', 'LineWidth',2)
%     view(3)
%     
%     title({['Elevations at (x,y) and (x_p,y_p) obtained using'];['functions topsurface and topfit (wrt ',num2str(plumes{i}.year),' plume)']})
%     legend('topsurface(x,y)','topfit(x,y)','topsurface(x_p,y_p)','topfit(x_p,y_p)','Location','NorthEast')
%     axis equal;
%     set(gca,'DataAspect',[1 1 1/50])
%     grid
%     zlabel('depth (meters)')
    %}
%{
    [ydata,ia,ic]=unique(line_coord(:,2));
    zt =@(y) interp1(ydata,topsurface(line_coord(ia,:)),y);
    topfit =@(coord)  zt(coord(:,2));
%}
    hCO2_tmp =@(coord)(topfit(coord)-topsurface(coord));
    hCO2 =@(coord) hCO2_tmp(coord).*(hCO2_tmp(coord)>0).*insidePolygon(line_coord,coord);
    plumes{i}.hCO2=hCO2;
    plumes{i}.h = hCO2(Gt.cells.centroids);
    
    subplot(1,2,2)
    plotCellData(Gt, hCO2(Gt.cells.centroids)); colorbar
    line(line_coord(:,1), line_coord(:,2), topsurface(line_coord)-1, 'LineWidth',2, 'Color','r'); % -1 to ensure visibility of plume outline against other plots
    view(2)
    
    title('hCO2: (topfit - topsurface) > 0, inside plume')
    legend('hCO2(x,y)','topsurface(x_p,y_p)', 'Location','NorthEast')
    axis equal tight;
    set(gca,'DataAspect',[1 1 1/50])
    grid
    zlabel('depth (meters)')
    %{
    figure(44),clf
    plotCellData(Gt,hCO2(Gt.cells.centroids));colorbar
    line(line_coord(:,1), line_coord(:,2),topsurface(line_coord), 'LineWidth',3, 'Color','r')
    %}
    
    %if i==numel(plumes)
%         plumeInside_coord_tmp = [Gt.cells.centroids(:,1).*insidePolygon(line_coord, Gt.cells.centroids), ...
%             Gt.cells.centroids(:,2).*insidePolygon(line_coord, Gt.cells.centroids)];
%         tmp1 = plumeInside_coord_tmp(:,1);
%         tmp2 = plumeInside_coord_tmp(:,2);
%         plumeInside_coord = [tmp1(plumeInside_coord_tmp(:,1)~=0), tmp2(plumeInside_coord_tmp(:,2)~=0)];
%         
%         figure;
%         hold on
%         plot3(plumeInside_coord(:,1), plumeInside_coord(:,2), topsurface(plumeInside_coord), 'xr')
%         plot3(plumeInside_coord(:,1), plumeInside_coord(:,2), topfit(plumeInside_coord), 'sqb')
%         plot3(line_coord(:,1), line_coord(:,2), topsurface(line_coord), 'r', 'LineWidth',2)
%         plot3(line_coord(:,1), line_coord(:,2), topfit(line_coord), 'b', 'LineWidth',2)
%         view(3)
%         
%         title('Elevations of (x,y) inside plume as determined by functions.')
%         legend('topsurface(x,y)','topfit(x,y)','topsurface(x_p,y_p)','topfit(x_p,y_p)','Location','NorthEast')
%         axis equal tight;
%         set(gca,'DataAspect',[1 1 1/50])
%         grid
%         zlabel('depth (meters)')
        % note: in the above plot, the differences are taken, and the
        % positive values are then assigned as hCO2. hCO2 is thus an
        % indication of how much the top surface differs from what the top
        % surface should be such that the plume migrates under a planar
        % caprock
    %end
    
    
end
end