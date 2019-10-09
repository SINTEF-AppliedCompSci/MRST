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

% ========================== Compute cell centroids ==========================

T_cell_centroids = T_cell_node_coords.contract('node') ./ ...
                   T_cell_numnodes.project(T_dim_ind);

% ============ Computing cell-node distances (x_node - x_centroid) ============

T_cell_node_distance = ...
    T_cell_node_coords - T_cell_centroids.project(T_cell_node_ind);

% =============================== Computing Nc ===============================

% Nc = [1 0 0 2 0 3;
%       0 2 0 1 3 0;
%       0 0 3 0 2 1];

NcInd.d     = [1, 1, 1, 2, 2, 2, 3, 3, 3];
NcInd.lform = [1, 4, 6, 2, 4, 5, 3, 5, 6];
NcInd.k     = [1, 2, 3, 2, 1, 3, 3, 2, 1];

T_cell_node_Nc = Tensor(NcInd) * T_cell_node_distance.changeIndexName('d', 'k');

% ================================ Compute q_i ================================

G = createAugmentedGrid(computeGeometry(G));
[qc, qf, qcvol] = calculateQC(G);

% make the qc-tensor (which has same structure as T_cell_node_coords)
T_cell_node_d_qc = T_cell_node_coords.sortIndices({'d', 'cell', 'node'});
T_cell_node_d_qc.nzvals = qc(:);

% - we need the face normals in cell_face_dim space (beware of the sign!)
% - we also need cell volumes
% - q in cell_node_dim space  (a 3D vector associated with each cell corner)
% - best option: use Xavier's function =calculateQC= to get the values, then ...
%     make a copy of the T_cell_node_coords tensor and write the values into the ...
%     nzvals (careful to use the correct strides)  (d has longest stride, cell ...
%                                                   medium and node shortest)

% ================================ Compute Wc ================================

WcInd = NcInd; % they have exactly the same structure
T_cell_node_Wc = ...
    Tensor([2, 1, 1, 2, 1, 1, 2, 1, 1], WcInd) * ...
    T_cell_node_d_qc.changeIndexName({'d', 'k'});
