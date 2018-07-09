function d = ddiff ( d1, d2 )

%*****************************************************************************80
%
%% DDIFF returns the signed distance to a region that is the difference of two regions.
%
%  Discussion:
%
%    The author comments that this is not the true signed distance function for
%    the difference set, in particular, around corners.
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
%    Input, real D1, D2, the signed distances to region 1 and 2.
%
%    Output, real D, the signed distance to the region formed by
%    removing from region 1 its intersection with region 2.
%
  d = max ( d1, -d2 );

  return
end

