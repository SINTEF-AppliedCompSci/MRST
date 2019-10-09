G = cartGrid([2,3,4]);

T_dim_ind = Tensor(ones(G.griddim, 1), 'd');
T_node_coords = Tensor(G.nodes.coords(:), {'node', G.nodes.num, 'd', G.griddim});

T_face_node_ind = ...
    Tensor(struct('face', rldecode((1:G.faces.num)', diff(G.faces.nodePos)),...
                  'node', G.faces.nodes));
T_cell_face_ind = ...
    Tensor(struct('cell', rldecode((1:G.cells.num)', diff(G.cells.facePos)),...
                  'face', G.cells.faces(:,1)));

T_cell_node_ind = Tensor.toInd(T_face_node_ind * T_cell_face_ind);

T_cell_numnodes = T_cell_node_ind.contract('node');

T_cell_node_coords = T_node_coords.project(T_cell_node_ind);

T_cell_centroids = T_cell_node_coords.contract('node') ./ ...
                   T_cell_numnodes.project(T_dim_ind);

T_cell_node_distance = ...
    T_cell_node_coords - T_cell_centroids.project(T_cell_node_ind);

Nc = [1 0 0 2 0 3;
      0 2 0 1 3 0;
      0 0 3 0 2 1];
T_Nc = Tensor(Nc, {'d', 'dlform'});

T_cell_node_Nc = T_Nc * T_cell_node_ind;