function p = pshift ( p, x0, y0 )

%*****************************************************************************80
%
%% PSHIFT shifts a set of points by a given increment.
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
%    Input/output, real P(NP,2), a set of points to be shifted.
%
%    Input, real X0, Y0, the X and Y increments to be added to each point.
%
  p(:,1) = p(:,1) + x0;
  p(:,2) = p(:,2) + y0;

  return
end
