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


T_cell_node_coords = T_node_coords ^ T_cell_node_ind;


T_cell_numnodes = T_cell_node_ind.contract('node');

% ----------------------------- compute centroids -----------------------------
fprintf('Compute centroids\n');

T_cell_centroids = (T_node_coords * T_cell_node_ind) ./ (T_cell_numnodes * T_dim_ind);



% ---------------------- compute node-centroid distances ----------------------
fprintf('Compute node-centroid distances\n');
T_cell_node_distance = ...
    T_cell_node_coords - (T_cell_centroids ^ T_cell_node_ind);

% change to homogeneous coordinates (to allow translation)
T_hom = SmartTensor([0,0,0,1], {'dim'}); % fourth component of homogeneous coord.
T_cell_node_distance_hom = T_cell_node_distance + (T_hom ^ T_cell_node_ind);

% -------------------------------- compute Nc --------------------------------
fprintf('Compute Nc\n');
% Nc = [1 0 0 2 0 3;
%       0 2 0 1 3 0;
%       0 0 3 0 2 1];

NcInd.d     = [1, 1, 1, 2, 2, 2, 3, 3, 3];
NcInd.lform = [1, 4, 6, 2, 4, 5, 3, 5, 6];
NcInd.k     = [1, 2, 3, 2, 1, 3, 3, 2, 1];

T_Nc = SmartTensor([], [NcInd.d', NcInd.lform', NcInd.k'], {'dim', 'lform', 'kc'});

% -------------------------------- compute q_i --------------------------------
fprintf('Compute q_i\n');

G = createAugmentedGrid(computeGeometry(G));
[qc, qf, qcvol] = calculateQC(G);

T_cell_node_d_qc = T_cell_node_coords.sortIndices({'dim', 'cell', 'node'});
T_cell_node_d_qc.components{1}.coefs = qc(:);

T_cell_node_d_qc_hom = T_cell_node_d_qc + T_hom ^ T_cell_node_ind;

% -------------------------------- compute Wc --------------------------------
fprintf('Compute Wc.\n');

T_Wc = SmartTensor([2, 1, 1, 2, 1, 1, 2, 1, 1], ...
                   [NcInd.d', NcInd.lform', NcInd.k'], {'dim', 'lform', 'kc'});

% -------------------------------- compute Nr --------------------------------
fprintf('Compute Nr.\n');
% Nr = [_1  0  0 2 0 -3;
%        0 _1  0 -1 3 0;
%        0  0 _1 0 -2 1]

NrInd = NcInd;
NrInd.k(1:3:8) = 4; % translation represented by 4. component of hom. coord

T_Nr = SmartTensor([1, 1, -1, 1, -1, 1, 1, -1, 1], ...
                   [NrInd.d', NrInd.lform', NrInd.k'], {'dim', 'lform', 'kr'});

% -------------------------------- Compute Wr --------------------------------
fprintf('Compute Wr.\n');


