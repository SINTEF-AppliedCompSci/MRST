function u = uniformity ( p, t, fh, varargin )

%*****************************************************************************80
%
%% UNIFORMITY determines the uniformity of a mesh.
%
%  Discussion:
%
%    The measure of quality is the comparison of the radius of the
%    circumcenter of each triangle to the mesh size function evaluated
%    at the center of the circumscribed circle.
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
%    Input, external FH, the mesh size function.
%
%    Input, VARARGIN, room for optional arguments to be passed to FH.
%
%    Output, real U(NT), the uniformity measure.
%
  [ pc, r ] = circumcenter ( p, t );
  hc = feval ( fh, pc, varargin{:} );

  sz = r ./ hc;
  u = std ( sz ) / mean ( sz );

  return
end
