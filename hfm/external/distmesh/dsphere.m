function d = dsphere ( p, xc, yc, zc, r )

%*****************************************************************************80
%
%% DSPHERE returns the signed distance of one or more points to a sphere.
%
%  Discussion:
%
%    The corresponding routine in 2D is called DCIRCLE.
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
%    Input, real P(NP,3), the point coordinates.
%
%    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
%
%    Input, real R, the radius of the sphere.
%
%    Output, real D(NP), the signed distance of each point to the
%    sphere.  The point is inside, on, or outside the sphere depending
%    on whether D is negative, zero, or positive.
%
  d = sqrt ( ( p(:,1) - xc ).^2 ...
           + ( p(:,2) - yc ).^2 ...
           + ( p(:,3) - zc ).^2 ) - r;

  return
end
