function bfinx = getFormationsClosedBdryFaces(fmName, Gt)
% A fault line may exist along the boundary of a CO2 Atlas formation.
% By specifying the formation name, the indexes of the boundary faces of Gt
% will be returned

switch lower(fmName),
   case {'stofm'},
    % The closed boundaries of Sto lies within two rectangular areas of
    % size W x H, where W is a line given by points (Wx1,Wy1) and
    % (Wx2,Wy2), and where H is a line given by points (Hx1,Hy1) and
    % (Hx2,Hy2).

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
    
    % Return boundary face indexes that are considered to be 'closed'
    bfinx = [bfinx_s; bfinx_e; bfinx_w];
    
    
    
   otherwise
      error(1,['Either incorrect formation name, ', ...
          'or no closed boundary faces reported for that formation yet']);
end

    

end