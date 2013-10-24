function T = computeCoupledTrans(G, rock)
% Computes transmissibilities for coupled systems. Wrapper for computeTrans

    tags = G.cells.faces(:,2);
    % Find top/bottom half faces
    vf  = (tags == 5 | tags == 6) | ~ismember(G.cells.faces(:,1), G.facesBnd.index);
    
    % Modify centroids by setting z-component to zero. This is because
    % z-components are not meaningful in the transition between 2D and 3D
    % grids.
    cc_h = G.cells.centroids;
    cc_h(:,3) = 0;
    
    cfc_h = G.faces.centroids(G.cells.faces(:,1), :);
    cfc_h(:,3) = 0;
    
    % Calculate horizontal and vertical transmissibilities. 
    T_h = computeTrans(G, rock, 'cellCenters', cc_h, 'cellFaceCenters', cfc_h);
    T_v = computeTrans(G, rock);
    n = size(G.cells.faces,1);
    
    ind     = false(n,1);
    ind(vf) = true;

    % Use the vertical transmissibilities where suitable
    T = zeros(numel(ind),1);
    T(ind)  = T_v(ind);
    T(~ind) = T_h(~ind);
end
