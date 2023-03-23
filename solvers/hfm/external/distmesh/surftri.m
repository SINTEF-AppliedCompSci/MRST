function tri = surftri ( p, t )

%*****************************************************************************80
%
%% SURFTRI finds the surface triangles in a tetrahedral mesh.
%
%  Discussion:
%
%    The 2D version of this code is called BOUNDEDGES.
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
%    Input, real P(NP,3), the coordinates of a set of nodes.
%
%    Input, integer T(NT,1:4), a list of the nodes which make up each 
%    tetrahedron in the mesh.
%

%
% Form all faces, non-duplicates are surface triangles
%
  faces = [t(:,[1,2,3]);
           t(:,[1,2,4]);
           t(:,[1,3,4]);
           t(:,[2,3,4])];
  node4 = [t(:,4);t(:,3);t(:,2);t(:,1)];
  faces = sort(faces,2);
  [foo,ix,jx] = unique(faces,'rows');
  vec = histc(jx,1:max(jx));
  qx = find(vec==1);
  tri = faces(ix(qx),:);
  node4 = node4(ix(qx));
%
% Orientation
%
  v1 = p(tri(:,2),:)-p(tri(:,1),:);
  v2 = p(tri(:,3),:)-p(tri(:,1),:);
  v3 = p(node4,:)-p(tri(:,1),:);
  ix = find(dot(cross(v1,v2,2),v3,2)>0);
  tri(ix,[2,3]) = tri(ix,[3,2]);

  return
end
