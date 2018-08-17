function fault_faces = faults_on_boundaries_Sto(Gt)
% Returns (face) index of fault faces on boundaries in Stø formation
%
% These indices should be combined with the other internal fault faces
% found within the formation.

    % ------------------------------------
    %   Fault faces along south boundary
    % ------------------------------------
    % The southern faulted boundaries of Sto lies within two rectangular
    % areas of size W x H, where W is a line given by points (Wx1,Wy1) and
    % (Wx2,Wy2), and where H is a line given by points (Hx1,Hy1) and
    % (Hx2,Hy2).

    % The first rectangular area is given by:
    cinxAA = getCellIndex(Gt, 8.902e5, 7.906e6);
    cinxBB = getCellIndex(Gt, 8.967e5, 7.911e6);

    % And the boundary faces are facing south, east, and west.
    cinxAA_ij = [Gt.cells.ij(cinxAA,1); Gt.cells.ij(cinxAA,2)];
    cinxBB_ij = [Gt.cells.ij(cinxBB,1); Gt.cells.ij(cinxBB,2)];
    bfinx_s = boundaryFaceIndices(Gt, 'South', cinxAA_ij(1):cinxBB_ij(1), cinxAA_ij(2):cinxBB_ij(2));
    bfinx_e = boundaryFaceIndices(Gt, 'East',  cinxAA_ij(1):cinxBB_ij(1), cinxAA_ij(2):cinxBB_ij(2));


    % The second rectangular area is given by:
    cinx1 = getCellIndex(Gt, 8.96e5, 7.91e6);
    cinx2 = getCellIndex(Gt, 9.642e5, 7.959e6);

    % And the boundary faces are facing south and east.
    ij_inx_1 = [Gt.cells.ij(cinx1,1); Gt.cells.ij(cinx1,2)];
    ij_inx_2 = [Gt.cells.ij(cinx2,1); Gt.cells.ij(cinx2,2)];
    bfinx_s = [bfinx_s; boundaryFaceIndices(Gt, 'South', ij_inx_1(1):ij_inx_2(1), ij_inx_1(2):ij_inx_2(2))];
    bfinx_e = [bfinx_e; boundaryFaceIndices(Gt, 'East',  ij_inx_1(1):ij_inx_2(1), ij_inx_1(2):ij_inx_2(2))];


    % ------------------------------------
    %   Fault faces along north boundary
    % ------------------------------------
    % Here, boundary faces are either north-facing or west-facing.
    cinxN1 = getCellIndex(Gt, 9.308e5, 8.016e6);
    cinxN2 = getCellIndex(Gt, 9.672e5, 8.043e6);
    ij_inx_N1 = [Gt.cells.ij(cinxN1,1); Gt.cells.ij(cinxN1,2)];
    ij_inx_N2 = [Gt.cells.ij(cinxN2,1); Gt.cells.ij(cinxN2,2)];
    bfinx_n = boundaryFaceIndices(Gt, 'North', ij_inx_N1(1):ij_inx_N2(1), ij_inx_N1(2):ij_inx_N2(2));
    bfinx_w = boundaryFaceIndices(Gt, 'West',  ij_inx_N1(1):ij_inx_N2(1), ij_inx_N1(2):ij_inx_N2(2));

    
    % ------------------------------------
    %   Putting face indices together:
    % ------------------------------------
    fault_faces = [bfinx_s; bfinx_e; bfinx_n; bfinx_w];

end