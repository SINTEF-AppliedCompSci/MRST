function [ pc, r ] = circumcenter ( p, t )

%*****************************************************************************80
%
%% CIRCUMCENTER computes the circumcenters of a set of triangles.
%
%  Discussion:
%
%    The circumcenter of a triangle is the circle which passes through
%    all three vertices of the triangle.
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
%  Modified:
%
%    11 June 2004
%
%  Parameters:
%
%    Input, real P(NP,2), the coordinates of a set of nodes.
%
%    Input, integer T(NT,1:3), a list of the nodes which make up each triangle
%    of a triangulation of the nodes in P.
%
%    Output, real PC(NT,2), the centers of the circumcircles.
%
%    Output, real R(NT,1), the radii of the circumcircles.
%
  nt = size ( t, 1 );

  pc = zeros ( nt, 2 );
  r  = zeros ( nt, 1 );

  for it = 1 : nt

    ct = t(it,:);
    dp1 = p(ct(2),:) - p(ct(1),:);
    dp2 = p(ct(3),:) - p(ct(1),:);
  
    mid1 = ( p(ct(2),:) + p(ct(1),:) ) / 2;
    mid2 = ( p(ct(3),:) + p(ct(1),:) ) / 2;
  
    s = [ -dp1(2), dp2(2); dp1(1), -dp2(1) ] \ [ -mid1 + mid2 ]';
  
    cpc = mid1 + s(1) * [ -dp1(2), dp1(1) ];
    cr = norm ( p(ct(1),:) - cpc );
  
    pc(it,:) = cpc;
    r(it,1) = cr;
    
  end

  return 
end
