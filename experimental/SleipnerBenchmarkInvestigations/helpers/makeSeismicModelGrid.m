function [G, Gt, rock, rock2D] = makeSeismicModelGrid( refineLevel )
% Construct SEISMIC grid; a grid constructed from seismic data

    % get seismic data
    tdata = readSeismicAsc('data/1994-TopLayer9-kriging_Depth.asc');
    
    % use original grid to help construct seismic grid
    [ ~, Gt_tmp, ~, rock2D_tmp ] = makeSleipnerModelGrid('modelName', 'ORIGINALmodel', 'refineLevel',refineLevel);
    
    X = reshape(Gt_tmp.cells.centroids(:,1), Gt_tmp.cartDims);
    Y = reshape(Gt_tmp.cells.centroids(:,2), Gt_tmp.cartDims);
    H = reshape(Gt_tmp.cells.H, Gt_tmp.cartDims);
    FH = scatteredInterpolant(X(:), Y(:), H(:));
    
    % making 3D and 2D grids
    G      = ascData2Grid(tdata, 'cartdims',Gt_tmp.cartDims, 'FH',FH);  % num cells doesn't match rock, rather matches rock2D
    Gt     = topSurfaceGrid(computeGeometry(G));                        % num cells matches rock2D
    rock   = rock2D_tmp;
    rock2D = rock2D_tmp;


end

