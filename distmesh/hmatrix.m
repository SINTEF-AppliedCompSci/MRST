function h = hmatrix ( p, xx, yy, dd, hh, varargin )

%*****************************************************************************80
%
%% HMATRIX computes the mesh size function from values specified on a Cartesian grid.
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
%    Input, real XX(NDATA), YY(NDATA), the coordinates of points at which the
%    mesh size function has been specified.
%
%    Input, real DD, a dummy parameter included for consistency with the DMATRIX routine.
%
%    Input, real HH(NDATA), the value of the mesh size function at the data points.
%
%    Input, VARARGIN, room for extra arguments.
%
%    Output, real H(NP), the interpolated value of the mesh size function at each
%    of the input points.
%
  h = interp2 ( xx, yy, hh, p(:,1), p(:,2), '*linear' );

  return
end
