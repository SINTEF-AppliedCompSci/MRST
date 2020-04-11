function [S, A, C, D] = vem_system(G, E_ad, nu_ad)

assert(G.griddim == 3);

%% preparing grid if necessary
if isempty(intersect('computeGeometry', G.type))
   G = computeGeometry(G);
end
if isempty(intersect('createAugmentedGrid', G.type))
   G = createAugmentedGrid(G);
end

%% construct C-tensor
ixnames = {'i', 'j', 'cell'}; % 'c' indicates cell subscript
C = Enu2C_AD(E_ad, nu_ad, G, ixnames);

%% construct D-tensor
D = C ^ ...
    SparseTensor([1,1,1,2,2,2]', ixnames{1}) ^ ...
    SparseTensor([1,1,1,2,2,2]', ixnames{2});

%% compute cell centroids and node-centroid distances

% dimension tensor
DIM = SparseTensor(ones(G.griddim, 1), {'dim'});

% node coordinates
NCOORD = SparseTensor(G.nodes.coords, {'node', 'dim'}); 

% node-cell indicators
NCIND = SparseTensor([], ...
                     [rldecode((1:G.cells.num)', diff(G.cells.nodePos)), ...
                      G.cells.nodes], ...
                     {'cell', 'node'});

% number of incident nodes per cell 
CELL_NUMNODES = NCIND.contract('node'); 

% cell centroids
CTROID = NCOORD * NCIND ./ (CELL_NUMNODES * DIM); 

% node-cell centroid distances
NCDIST = (NCOORD ^ NCIND) - (CTROID ^ NCIND);

%% compute volumes and areas
[qc, qf, qcvol] = calculateQC(G);

ED = D ^ SparseTensor(qcvol, {'cell'}); % matrix |E| D

% setup QC tensor (3D QC value for each cell-node pair)
QC = DIM ^ NCIND;
QC = QC.sortIndices({'dim', 'cell', 'node'}).expandall();
QC.components{1}.coefs = qc(:);

%% compute Wc

NcInd.d     = [1, 1, 1, 2, 2, 2, 3, 3, 3];
NcInd.lform = [1, 4, 6, 2, 4, 5, 3, 5, 6];
NcInd.k     = [1, 2, 3, 2, 1, 3, 3, 2, 1];

WC = SparseTensor([2, 1, 1, 2, 1, 1, 2, 1, 1], ...
                  [NcInd.d', NcInd.lform', NcInd.k'], ...
                  {'dim', 'lform', 'k'});
                  
WC = WC * QC.changeIndexName('dim' ,'k');



%% Assembly

WC2 = WC.changeIndexName({'dim', 'node', 'lform'}, {'dim2', 'node2', 'lform2'});

%WC; WC2; ED;

% S = WC.sub('cell', 1) * ED.sub('cell', 1) * WC2.sub('cell', 1);
% S = S.contract('i', 'lform').contract('j', 'lform2');
% for i = 2:G.cells.num
%    tmp = WC.sub('cell', i) * ED.sub('cell', i) * WC2.sub('cell', i);
%    S = S + tmp.contract('i', 'lform').contract('j', 'lform2').expandall();
%    i
% end
   

% S = WC ^ ED;
% S = S.contract('i', 'lform');
% S.expandall();
% S = S * WC2.changeIndexName('lform2', 'j');
% S.expandall();

S = WC ^ ED ^ WC2;
S = S.contract('i', 'lform').contract('j', 'lform2').contract('cell').expandall();

%% prepare return values
%S = ED;
A = [];
end