%% testing update to trapAnalysis for closed boundaries

[Gt, rock2D] = getFormationTopGrid('Stofm',1);


%% Get face indexes that are closed.
% To do this, we first grab all the boundary face indexes. Then, we
% determine which ones are east or south facing, and are located within a
% certain part of the grid (i.e., the south-east boundary of Sto).
bfinx_all = boundaryFaces(Gt);

% v = [1, -0.7];
% tmp = ones(size(Gt.faces.neighbors,1),1);
% tmp(Gt.faces.neighbors(:,1) == 0) = -1;
% dotProd = bsxfun(@times, tmp, Gt.faces.normals) * v';
% dotProd_of_bdryFaces = zeros(Gt.faces.num,1);
% dotProd_of_bdryFaces(bfinx_all) = dotProd(bfinx_all);
% bfinx_se = find(dotProd_of_bdryFaces > 0);

% The closed boundaries of Sto lies within two rectangular areas of size W
% x H, where W is a line given by points (Wx1,Wy1) and (Wx2,Wy2), and where
% H is a line given by points (Hx1,Hy1) and (Hx2,Hy2).

% The first rectangular area is given by:
W1 = getCellIndex(Gt, 8.852e5, 7.911e6);
W2 = getCellIndex(Gt, 8.977e5, 7.911e6);
H1 = getCellIndex(Gt, 8.902e5, 7.906e6);
H2 = getCellIndex(Gt, 8.967e5, 7.911e6);
% And the boundary faces are facing south, east, and west.
W1_ij = [Gt.cells.ij(W1,1); Gt.cells.ij(W1,2)];
W2_ij = [Gt.cells.ij(W2,1); Gt.cells.ij(W2,2)];
H1_ij = [Gt.cells.ij(H1,1); Gt.cells.ij(H1,2)];
H2_ij = [Gt.cells.ij(H2,1); Gt.cells.ij(H2,2)];
bfinx_s = boundaryFaceIndices(Gt, 'South', W1_ij(1):W2_ij(1), H1_ij(2):H2_ij(2));
bfinx_e = boundaryFaceIndices(Gt, 'East',  W1_ij(1):W2_ij(1), H1_ij(2):H2_ij(2));
bfinx_w = boundaryFaceIndices(Gt, 'West',  W1_ij(1):W2_ij(1), H1_ij(2):H2_ij(2));

% The second rectangular area is given by:
cinx1 = getCellIndex(Gt, 8.967e5, 7.911e6);
cinx2 = getCellIndex(Gt, 9.642e5, 7.959e6);

% And the boundary faces are facing south and east.
ij_inx_1 = [Gt.cells.ij(cinx1,1); Gt.cells.ij(cinx1,2)];
ij_inx_2 = [Gt.cells.ij(cinx2,1); Gt.cells.ij(cinx2,2)];
bfinx_s = [bfinx_s; boundaryFaceIndices(Gt, 'South', ij_inx_1(1):ij_inx_2(1), ij_inx_1(2):ij_inx_2(2))];
bfinx_e = [bfinx_e; boundaryFaceIndices(Gt, 'East',  ij_inx_1(1):ij_inx_2(1), ij_inx_1(2):ij_inx_2(2))];

%% Compare trapAnalysis with and without closed boundaries
% We obtain the traps in the Sto formation, using open boundaries and using
% closed boundaries along the south-east boundary.
ta    = trapAnalysis(Gt,false);
ta_cb = trapAnalysis(Gt,false,'closed_boundary_edges',[bfinx_s; bfinx_e; bfinx_w]);

figure;
mapPlot(gcf, Gt, 'traps',ta.traps, 'rivers',ta.cell_lines)
h1 = plotFaces(Gt, bfinx_all, 'EdgeColor','g', 'FaceColor','g','LineWidth',2);
legend([h1],{'open boundary'})
axis equal tight off

figure;
mapPlot(gcf, Gt, 'traps',ta_cb.traps, 'rivers',ta_cb.cell_lines)
h1 = plotFaces(Gt, bfinx_all, 'EdgeColor','g', 'FaceColor','g', 'LineWidth',2);
h2 = plotFaces(Gt, [bfinx_s; bfinx_e; bfinx_w], 'EdgeColor','r', 'FaceColor','r', 'LineWidth',2);
legend([h1, h2],{'open boundary'; 'closed boundary'})
axis equal tight off

%% De-bug mapPlot for trap along edge
figure;
trap1 = zeros(Gt.cells.num,1);
trap1(ta_cb.traps == 2) = ta_cb.traps(ta_cb.traps == 2);
mapPlot(gcf, Gt, 'traps',trap1, 'rivers',ta_cb.cell_lines)





