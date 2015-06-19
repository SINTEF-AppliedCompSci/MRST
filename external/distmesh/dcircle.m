function d = dcircle ( p, xc, yc, r )

%*****************************************************************************80
%
%% DCIRCLE returns the signed distance of one or more points to a circle.
%
%  Discussion:
%
%    The corresponding routine in 3D is DSPHERE.
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
%    Input, real XC, YC, the coordinates of the center of the circle.
%
%    Input, real R, the radius of the circle.
%
%    Output, real D(NP), the signed distance of each point to the
%    circle.  The point is inside, on, or outside the circle depending
%    on whether D is negative, zero, or positive.
%
  d = sqrt ( ( p(:,1) - xc ).^2 + ( p(:,2) - yc ).^2 ) - r;

  return
end
