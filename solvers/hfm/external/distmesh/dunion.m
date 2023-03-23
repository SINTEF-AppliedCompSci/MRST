function d = dunion ( d1, d2 )

%*****************************************************************************80
%
%% DUNION returns the signed distance to a union of two regions.
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
%    Input, real D1, D2, the signed distance of one or more points
%    to two regions.
%
%    Output, real D, the signed distance of one or more points to
%    the region formed by the union of the two regions.
%
  d = min ( d1, d2 );

  return
end
