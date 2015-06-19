function d = dintersect ( d1, d2 )

%*****************************************************************************80
%
%% DINTERSECT sets the signed distance to the intersection of two regions.
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
%    Input, real D1, D2, the signed distance to each of the regions.
%
%    Output, real D, the signed distance to the region formed by the
%    intersection of the two regions.
%
  d = max ( d1, d2 );

  return
end
