function d = drectangle ( p, x1, x2, y1, y2 )

%*****************************************************************************80
%
%% DRECTANGLE returns the signed distance of one or more points to a rectangle.
%
%  Discussion:
%
%    The formula used here is not quite correct.  In particular, it is wrong
%    at the corners.  For a true distance function, compare DRECTANGLE0.
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
%    Input, real X1, X2, Y1, Y2, the left, right, bottom and top
%    coordinates of the rectangle.
%
%    Output, real D, the signed distance of the points to the rectangle,
%    which is negative, 0 or positive depending on whether each point
%    is inside, or, or outside the rectangle.
%
  d = - min ( min ( min ( -y1+p(:,2), y2-p(:,2) ), -x1+p(:,1) ), x2-p(:,1) );

  return
end
