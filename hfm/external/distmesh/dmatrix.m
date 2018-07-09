function d = dmatrix ( p, xx, yy, dd, varargin )

%*****************************************************************************80
%
%% DMATRIX returns the signed distance by interpolation from known values on a Cartesian grid.
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
%    Input, real P(NP,2), the point coordinates.
%
%    Input, real XX(NDATA), YY(NDATA), the points at which the distance was evaluated.
%
%    Input, real DD(NDATA), the signed distance function at the data points.
%
%    Input, VARARGIN, optional arguments.
%
%    Output, real D(NP), the signed distance to the region, estimated by interpolation.
%
  d = interp2 ( xx, yy, dd, p(:,1), p(:,2), '*linear' );

  return
end
