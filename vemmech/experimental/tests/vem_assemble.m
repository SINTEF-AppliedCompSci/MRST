%G = cartGrid([1, 1, 1]);
%G = cartGrid([2,2,2]);
%G = cartGrid([5, 10, 10]);
%G = cartGrid([10, 10, 10]);
%G = cartGrid([15, 15, 15]);
%G = cartGrid([20, 20, 20]);
%G = cartGrid([40, 40, 40]);
%G = cartGrid([22, 22, 22]);
G = cartGrid([30 30 30]);
%G = cartGrid([50, 50, 50]);

T_dim_ind = SparseTensor(ones(G.griddim, 1), {'dim'});
T_node_coords = SparseTensor(G.nodes.coords, {'node', 'dim'});

T_face_node_ind = ...
    SparseTensor([], ...
            [rldecode((1:G.faces.num)', diff(G.faces.nodePos)), G.faces.nodes], ...
                {'face', 'node'});

T_cell_face_ind = ...
    SparseTensor([], [rldecode((1:G.cells.num)', diff(G.cells.facePos)), ...
                     G.cells.faces(:,1)], ...
                {'cell', 'face'});

T_cell_node_ind = SparseTensor.ind(T_face_node_ind * T_cell_face_ind);

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
T_hom = SparseTensor([0,0,0,1], {'dim'}); % fourth component of homogeneous coord.
T_cell_node_distance_hom = T_cell_node_distance + (T_hom ^ T_cell_node_ind);

% -------------------------------- compute Nc --------------------------------
fprintf('Compute Nc\n');
% Nc = [1 0 0 2 0 3;
%       0 2 0 1 3 0;
%       0 0 3 0 2 1];

NcInd.d     = [1, 1, 1, 2, 2, 2, 3, 3, 3];
NcInd.lform = [1, 4, 6, 2, 4, 5, 3, 5, 6];
NcInd.k     = [1, 2, 3, 2, 1, 3, 3, 2, 1];

T_Nc = SparseTensor([], [NcInd.d', NcInd.lform', NcInd.k'], {'dim', 'lform', 'k'});

T_Nc = T_Nc * T_cell_node_distance.changeIndexName('dim', 'k');

% -------------------------------- compute q_i --------------------------------
fprintf('Compute q_i\n');

G = createAugmentedGrid(computeGeometry(G));
[qc, qf, qcvol] = calculateQC(G);


% use 'toInd()' to ensure coordinates with physical value '0' aren't eliminated.
T_cell_node_d_qc = (T_node_coords.toInd() ^ T_cell_node_ind);
T_cell_node_d_qc = T_cell_node_d_qc.sortIndices({'dim', 'cell', 'node'});
T_cell_node_d_qc.components{1}.coefs = qc(:);

T_cell_node_d_qc_hom = T_cell_node_d_qc + T_hom ^ T_cell_node_ind;

% -------------------------------- compute Wc --------------------------------
fprintf('Compute Wc.\n');

T_Wc = SparseTensor([2, 1, 1, 2, 1, 1, 2, 1, 1], ...
                   [NcInd.d', NcInd.lform', NcInd.k'], {'dim', 'lform', 'k'});

T_Wc = T_Wc * T_cell_node_d_qc.changeIndexName('dim', 'k');

% -------------------------------- compute Nr --------------------------------
fprintf('Compute Nr.\n');
% Nr = [_1  0  0 2 0 -3;
%        0 _1  0 -1 3 0;
%        0  0 _1 0 -2 1]

NrInd = NcInd;
NrInd.k(1:3:8) = 4; % translation represented by 4. component of hom. coord

T_Nr = SparseTensor([1, 1, -1, 1, -1, 1, 1, -1, 1], ...
                   [NrInd.d', NrInd.lform', NrInd.k'], {'dim', 'lform', 'k'});

T_Nr = T_Nr * T_cell_node_distance_hom.changeIndexName('dim', 'k');

% -------------------------------- Compute Wr --------------------------------
fprintf('Compute Wr.\n');

% Wr = 1/n   0   0   q2   0 -q3
%        0 1/n   0  -q1  q3   0 
%        0   0 1/n    0 -q2  q1

T_Wr = SparseTensor([1, 1, -1, 1, -1, 1, 1, -1, 1], ...
                   [NrInd.d', NrInd.lform', NrInd.k'], ...
                   {'dim', 'lform', 'k'});

T_cell_numnodes_inv = T_cell_numnodes.toInd() ./ T_cell_numnodes;
T_cell_node_numnodes_inv_hom = (T_hom * T_cell_numnodes_inv) ^ T_cell_node_ind;

T_cell_node_d_qc_numnodes_inv_hom = (T_cell_node_d_qc + T_cell_node_numnodes_inv_hom);

T_Wr = T_Wr * T_cell_node_d_qc_numnodes_inv_hom.changeIndexName('dim', 'k');

% --------------------------------- compute P ---------------------------------

T_Pr = T_Wr * T_Nr.changeIndexName({'dim', 'node'}, {'dim2', 'node2'});
T_Pc = T_Wc * T_Nc.changeIndexName({'dim', 'node'}, {'dim2', 'node2'});
tic; T_Pp = T_Pr + T_Pc; toc