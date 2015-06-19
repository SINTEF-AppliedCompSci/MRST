function e = boundedges ( p, t )

%*****************************************************************************80
%
%% BOUNDEDGES finds the boundary edges in a triangular mesh.
%
%  Discussion:
%
%    You may need this routine if you need to enforce boundary
%    conditions in a PDE.
%
%    The 3D version of this code is called SURFTRI.
%
%  Licensing:
%
%    (C) 2004 Per-Olof Persson. 
%    See COPYRIGHT.TXT for details.
%
%  Reference:
%
%    Per-Olof Persson and Gilbert Strang,
%    A Simple Mesh Generator in MATLAB,
%    SIAM Review,
%    Volume 46, Number 2, June 2004, pages 329-345.
%
%  Parameters:
%
%    Input, real P(NP,2), the coordinates of a set of nodes.
%
%    Input, integer T(NT,1:3), a list of the nodes which make up each triangle
%    of a triangulation of the nodes in P.
%
%    Output, integer E(*,*), ?
%

%
%  Form all edges, non-duplicates are boundary edges
%
  edges = [t(:,[1,2]);
           t(:,[1,3]);
           t(:,[2,3])];

  node3 = [t(:,3);t(:,2);t(:,1)];
  edges = sort(edges,2);
  [foo,ix,jx] = unique(edges,'rows');
  vec = histc(jx,1:max(jx));
  qx = find(vec==1);
  e = edges(ix(qx),:);
  node3 = node3(ix(qx));
%
%  Orientation
%
  v1 = p(e(:,2),:)-p(e(:,1),:);
  v2 = p(node3,:)-p(e(:,1),:);
  ix = find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
  e(ix,[1,2]) = e(ix,[2,1]);

  return
end
