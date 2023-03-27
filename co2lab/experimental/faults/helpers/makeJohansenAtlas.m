function [ grdecl, G, Gt, rock2D ] = makeJohansenAtlas()
% Construct NPD's Atlas model of Johansen

    % We could call 'getFormationTopGrid', however we want to get grdecl
    % and we also want G which may have changed during call to
    % 'topSurfaceGrid' yet wasn't updated in Gt.parent.
    
    %[Gt_atls, rock2D_atls] = getFormationTopGrid('Johansenfm',1);
    
    [grdecl, dataset, petroinfo] = getAtlasGrid('Johansenfm');
    grdecl = grdecl{1};
    G = processGRDECL(grdecl);
    G = mcomputeGeometry(G);
    [Gt, G] = topSurfaceGrid(G); % @@
    
    %% rock properties
    petinfo = petroinfo{1};
    rock2D.perm = petinfo.avgperm*ones(Gt.cells.num,1);
    rock2D.poro = petinfo.avgporo*ones(Gt.cells.num,1);
    
end

