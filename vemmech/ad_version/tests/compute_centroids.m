%G = cartGrid([20,30,10]);
G = cartGrid([10, 10, 10]);
%G = cartGrid([7, 7, 7]);
%G = cartGrid([2,3,4]);
%G = cartGrid([2,1, 1]);

T_dim_ind = Tensor(ones(G.griddim, 1), 'd');
T_node_coords = Tensor(G.nodes.coords(:), {'node', G.nodes.num, 'd', G.griddim});

T_face_node_ind = ...
    Tensor(struct('face', rldecode((1:G.faces.num)', diff(G.faces.nodePos)),...
                  'node', G.faces.nodes));
T_cell_face_ind = ...
    Tensor(struct('cell', rldecode((1:G.cells.num)', diff(G.cells.facePos)),...
                  'face', G.cells.faces(:,1)));

T_cell_node_ind = Tensor.ind(T_face_node_ind * T_cell_face_ind);

T_cell_numnodes = T_cell_node_ind.contract('node');

T_cell_node_coords = T_node_coords.project(T_cell_node_ind);

% ========================== Compute cell centroids ==========================

fprintf('Computing cell centroids\n');

T_cell_centroids = T_cell_node_coords.contract('node') ./ ...
                   T_cell_numnodes.project(T_dim_ind);

% ============ Computing cell-node distances (x_node - x_centroid) ============

T_cell_node_distance = ...
    T_cell_node_coords - T_cell_centroids.project(T_cell_node_ind);

% =============================== Computing Nc ===============================

fprintf('Computing Nc.\n');
% Nc = [1 0 0 2 0 3;
%       0 2 0 1 3 0;
%       0 0 3 0 2 1];

NcInd.d     = [1, 1, 1, 2, 2, 2, 3, 3, 3];
NcInd.lform = [1, 4, 6, 2, 4, 5, 3, 5, 6];
NcInd.k     = [1, 2, 3, 2, 1, 3, 3, 2, 1];

T_cell_node_Nc = Tensor(NcInd) * T_cell_node_distance.changeIndexName('d', 'k');

% ================================ Compute q_i ================================

fprintf('Computing q_i.\n');

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

fprintf('Computing Nc.\n');

WcInd = NcInd; % they have exactly the same structure
T_cell_node_Wc = ...
    Tensor([2, 1, 1, 2, 1, 1, 2, 1, 1], WcInd) * ...
    T_cell_node_d_qc.changeIndexName('d', 'k');

% ================================ Compute Nr ================================


fprintf('Compute Nr.\n');
% Nr = [_1  0  0 2 0 -3;
%        0 _1  0 -1 3 0;
%        0  0 _1 0 -2 1]

NrIndRot = Tensor([0 1 -1 0 -1 1 0 -1 1], NcInd); % zero coefs for translation part
NrIndTrans = Tensor([1 0 0 1 0 0 1 0 0], rmfield(NcInd, 'k'));

T_cell_node_Nr = ...
    NrIndRot   * T_cell_node_distance.changeIndexName('d', 'k') + ...
    NrIndTrans * T_cell_node_ind;

% ================================ Compute Wr ================================

fprintf('Compute Wr.\n');

% Wr = 1/n   0   0   q2   0 -q3
%        0 1/n   0  -q1  q3   0 
%        0   0 1/n    0 -q2  q1

T_cell_node_d_lform_ind = T_cell_node_Nr.toInd(); % create indicator function

WrIndRot = NrIndRot; % identical to Nr
WrIndTrans = NrIndTrans; % identical to Nr

T_cell_node_Wr = ...
    WrIndRot * T_cell_node_d_qc.changeIndexName('d', 'k') + ...
    ((WrIndTrans * T_cell_node_ind) ./ ...
     (T_cell_numnodes.project(T_cell_node_d_lform_ind)));


% ================================= Compute P =================================

fprintf('Compute P.\n');

T_Pr = T_cell_node_Wr * T_cell_node_Nr.changeIndexName({'d', 'node'}, {'d2', 'node2'});
fprintf('.\n');
T_Pc = T_cell_node_Wc * T_cell_node_Nc.changeIndexName({'d', 'node'}, {'d2', 'node2'});
fprintf('.\n');
Pp = T_Pr + T_Pc;
fprintf('.\n');
% ============================== Compute (I - P) ==============================

fprintf('Compute I-P.\n');

% identity tensor in dxd
Id = Tensor(struct('d', [1,2,3], 'd2', [1,2,3]));

% identity tensor in d x d x cell x cell x node x node
Id_node_d = Id.project(Pp.toInd());

%spy(Id_node_d.asMatrix({{'d', 'node'}, {'d2', 'node2'}}))

IminusP = Pp - Id_node_d;

%spy(IminusP.pruneZeros().asMatrix({{'d', 'node'}, {'d2', 'node2'}}))


% ================================= Compute D =================================

fprintf('Compute D.\n');

E = 0.51e9;
nu = 0.15;

Ev  = repmat(E, G.cells.num, 1);
nuv = repmat(nu, G.cells.num, 1);
C   = Enu2C(Ev, nuv, G);
D   = C2D(C, G);

T_D = Tensor(D(:), {'cell', G.cells.num, 'lform1', 6, 'lform2', 6});

% ============================== Compute Wc D Wc ==============================

fprintf('Compute Wc D Wc.\n');

tmp = ...
    T_cell_node_Wc.project(T_cell_node_Wc.changeIndexName({'d', 'lform', 'node'}, ...
                                                          {'d2', 'lform2', 'node2'}));

fprintf('Compute WDW\n');
WDW = tmp * T_D.changeIndexName('lform1', 'lform');
fprintf('.\n');

spy(WDW.asMatrix({{'d', 'node'}, {'d2', 'node2'}}))


% ================================= Compute S =================================

% compute alpha
%T_alpha = ;


% ======================== Compute full system matrix ========================


% Compute system matrix
% M = |E| WDW + IminusP S IminusP


% compute alpha