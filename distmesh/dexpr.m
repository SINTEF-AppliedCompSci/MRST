function d = dexpr ( p, fin, nit, alpha )

%*****************************************************************************80
%
%% DEXPR computes the signed distance for a general implicit expression FIN.
%
%  Discussion:
%
%    The assumption is that the boundary is defined as a level set function,
%    f(x,y) = 0, with f(x,y) < 0 for points inside the region, and 0 < f(x,y)
%    for points outside the region.
%
%    To find the point (X,Y) on the level set that is nearest to a given point 
%    (X0,Y0), we require that f(X,Y) be zero and that the vector from (X,Y) to
%    (X0,Y0) be parallel to the gradient (fx,fy) at (X,Y), which defines a pair
%    of nonlinear equations L(X,Y):
%
%      f(X,Y) = 0
%      ( X - X0 ) * fy - ( Y - Y0 ) * fx = 0
%
%    These equations are solved by Newton's method.
%
%    The use of this function requires access to the Symbolic Toolbox, although
%    it would be possible to rewrite it to accept derivatives of FIN as input
%    as well.  This routine is quite slow relative to a C version.
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
%    Input, string FIN, the expression that defines F(X,Y).
%
%    Input, integer NIT, the number of iterations.  The default value of NIT is 20.
%
%    Input, real ALPHA, the damping parameter for the Newton method.  
%    The default value of ALPHA is 0.1.
%
%    Output, real D(NP), the signed distance of each point to the boundary.
%
  if ( nargin < 3 ) 
    nit = 20; 
  end

  if ( nargin < 4 )
    alpha = 0.1;
  end

  fx = inline ( vectorize(maple(['diff(',fin,',x)'])), 'x', 'y' );
  fy = inline ( vectorize(maple(['diff(',fin,',y)'])), 'x', 'y' );
  fxx = inline ( vectorize(maple(['diff(',fin,',x$2)'])), 'x', 'y' );
  fyy = inline ( vectorize(maple(['diff(',fin,',y$2)'])), 'x', 'y' );
  fxy = inline ( vectorize(maple(['diff(diff(',fin,',x),y)'])), 'x', 'y' );
  f = inline ( vectorize(fin), 'x', 'y' );

  x0 = p(:,1);
  y0 = p(:,2);
  x = x0;
  y = y0;

  for it = 1 : nit

    cf = feval(f,x,y);
    cfx = feval(fx,x,y);
    cfy = feval(fy,x,y);
    cfxx = feval(fxx,x,y);
    cfxy = feval(fxy,x,y);
    cfyy = feval(fyy,x,y);
%
%  Evaluate L:
%
    F1 = cf;
    F2 = (x-x0) .* cfy - (y-y0) .* cfx;
%
%  Compute the jacobian of L:
%
%    fx                                            fy
%    fy + ( x - x0 ) * fxy - ( y - y0 ) * fxx     -fx - ( y - y0 ) * fxy + ( x - x0 ) * fyy
%    
    J11 = cfx;
    J12 = cfy;
    J21 =  cfy + (x-x0) .* cfxy - (y-y0) .* cfxx;
    J22 = -cfx - (y-y0) .* cfxy + (x-x0) .* cfyy;
  
    detJ = J11 .* J22 - J12 .* J21;
    detJ(detJ==0) = inf;
%
%  Apply the Newton correction to (X,Y).
%  
    x = x - alpha * (  J22 .* F1 - J21 .* F2 ) ./ detJ;
    y = y - alpha * ( -J12 .* F1 + J11 .* F2 ) ./ detJ;

  end
%
%  Having found the nearest point (X,Y) on the level set, we can compute the signed distance.
%
  d = sqrt((x-x0).^2+(y-y0).^2) .* sign(feval(f,x0,y0));

  return
end
