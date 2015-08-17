function plumes = makeSurfaceData(plumes,Gt)
    X=reshape(Gt.cells.centroids(:,1),Gt.cartDims(1),Gt.cartDims(2));
    Y=reshape(Gt.cells.centroids(:,2),Gt.cartDims(1),Gt.cartDims(2));
    Z=reshape(Gt.cells.z,Gt.cartDims(1),Gt.cartDims(2));
    topsurface=@(coord) interp2(X',Y',Z',coord(:,1),coord(:,2));
    for i=1:numel(plumes)
        line_coord=plumes{i}.outline;
    nd=size(line_coord,1);
    A=[ones(nd,1),line_coord];
    rhs=topsurface(line_coord);
    vec=A\rhs;
    topfit=@(coord) [ones(size(coord,1),1),coord]*vec;
    hCO2_tmp =@(coord)(topfit(coord)-topsurface(coord));
    hCO2 =@(coord) hCO2_tmp(coord).*(hCO2_tmp(coord)>0).*insidePolygon(line_coord,coord);
    plumes{i}.hCO2=hCO2;
    plumes{i}.h = hCO2(Gt.cells.centroids);
    %{
    figure(44),clf
    plotCellData(Gt,hCO2(Gt.cells.centroids));colorbar
    line(line_coord(:,1), line_coord(:,2),topsurface(line_coord), 'LineWidth',3, 'Color','r')
    %}    
    end   
end