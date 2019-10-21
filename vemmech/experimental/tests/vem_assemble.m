G = cartGrid([2, 3, 4]);

T_dim_ind = SmartTensor(ones(G.griddim, 1), {'dim'});
T_node_coords = SmartTensor(G.nodes.coords, {'node', 'dim'});

T_face_node_ind = ...
    SmartTensor([], ...
            [rldecode((1:G.faces.num)', diff(G.faces.nodePos)), G.faces.nodes], ...
                {'face', 'node'});

T_cell_face_ind = ...
    SmartTensor([], [rldecode((1:G.cells.num)', diff(G.cells.facePos)), ...
                     G.cells.faces(:,1)], ...
                {'cell', 'face'});

T_cell_node_ind = SmartTensor.ind(T_face_node_ind * T_cell_face_ind);


T_cell_node_coords = T_node_coords.semi_contract_with(T_cell_node_ind);


T_cell_numnodes = T_cell_node_ind.contract_one('node');

% ----------------------------- compute centroids -----------------------------
T_cell_centroids = (T_node_coords * T_cell_node_ind) ./ (T_cell_numnodes * T_dim_ind);



% ---------------------- compute node-centroid distances ----------------------

T_cell_node_distance = ...
    T_cell_node_coords - T_cell_centroids.semi_contract_with(T_cell_node_ind);

