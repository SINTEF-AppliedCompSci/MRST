function [ p, t ] = distmesh_2d ( fd, fh, h0, box, iteration_max, pfix, ...
  varargin )

%*****************************************************************************80
%
%% DISTMESH_2D is a 2D mesh generator using distance functions.
%
%  Example: 
%
%    Uniform Mesh on Unit Circle:
%
%      fd = inline('sqrt(sum(p.^2,2))-1','p');
%      [p,t] = distmesh_2d(fd,@huniform,0.2,[-1,-1;1,1],100,[]);
%
%    Rectangle with circular hole, refined at circle boundary:
%
%      fd = inline('ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5))','p');
%      fh = inline('min(4*sqrt(sum(p.^2,2))-1,2)','p');
%      [p,t] = distmesh_2d(fd,fh,0.05,[-1,-1;1,1],500,[-1,-1;-1,1;1,-1;1,1]);
%
%  Licensing:
%
%    (C) 2004 Per-Olof Persson. 
%    See COPYRIGHT.TXT for details.
%
%  Modified:
%
%    09 June 2012
%
%  Author:
%
%    Per-Olof Persson
%    Modifications by John Burkardt
%
%  Reference:
%
%    Per-Olof Persson, Gilbert Strang,
%    A Simple Mesh Generator in MATLAB,
%    SIAM Review,
%    Volume 46, Number 2, June 2004, pages 329-345.
%
%  Parameters:
%
%    Input, function FD, signed distance function d(x,y).
%
%    Input, function FH, scaled edge length function h(x,y).
%
%    Input, real H0, the initial edge length.
%
%    Input, real BOX(2,2), the bounding box [xmin,ymin; xmax,ymax].
%
%    Input, integer ITERATION_MAX, the maximum number of iterations.
%    The iteration might terminate sooner than this limit, if the program decides
%    that the mesh has converged.
%
%    Input, real PFIX(NFIX,2), the fixed node positions.
%
%    Input, VARARGIN, aditional parameters passed to FD and FH.
%
%    Output, real P(N,2), the node positions.
%
%    Output, integer T(NT,3), the triangle indices.
%
%  Local parameters:
%
%    Local, real GEPS, a tolerance for determining whether a point is "almost" inside
%    the region.  Setting GEPS = 0 makes this an exact test.  The program currently 
%    sets it to 0.001 * H0, that is, a very small multiple of the desired side length
%    of a triangle.  GEPS is also used to determine whether a triangle falls inside
%    or outside the region.  In this case, the test is a little tighter.  The centroid
%    PMID is required to satisfy FD ( PMID ) <= -GEPS.
%
  dptol = 0.001;
  ttol = 0.1;
  Fscale = 1.2;
  deltat = 0.2;
  geps = 0.001 * h0;
  deps = sqrt ( eps ) * h0;
  iteration = 0;
  triangulation_count = 0;
%
%  1. Create the initial point distribution by generating a rectangular mesh
%  in the bounding box.
%
  [ x, y ] = meshgrid ( box(1,1) : h0           : box(2,1), ...
                        box(1,2) : h0*sqrt(3)/2 : box(2,2) );
%
%  Shift the even rows of the mesh to create a "perfect" mesh of equilateral triangles.
%  Store the X and Y coordinates together as our first estimate of "P", the mesh points
%  we want.
%
  x(2:2:end,:) = x(2:2:end,:) + h0 / 2;
  p = [ x(:), y(:) ];
%
%  Instead of a regular mesh, you can initialize P with random values here.
%
%  2. Remove mesh points that are outside the region, 
%  then satisfy the density constraint.
%
%  Keep only points inside (or almost inside) the region, that is, FD(P) < GEPS.
%
  p = p( feval ( fd, p, varargin{:} ) < geps, : );
%
%  Set R0, the relative probability to keep a point, based on the mesh density function.
%
  r0 = 1 ./ feval ( fh, p, varargin{:} ).^2;
%
%  Apply the rejection method to thin out points according to the density.
%
  p = [ pfix; p(rand(size(p,1),1) < r0 ./ max ( r0 ),: ) ];
%
%  The following obscure commands are written so that:
%  * We ALWAYS keep the fixed points at the beginning of the P array.
%  * We remove any points which are duplicates of these points.
%  That way, we do not later allow the fixed points to move,
%  because we know they are at the beginning of the array,
%  and DELAUNAYN stops complaining about duplicate points.
%  JVB, 09 June 2012.
%
  [ q, i, j ] = unique ( p, 'rows', 'first' );
  k = unique ( i );
  p = p(k,:);

  N = size ( p, 1 );
%
%  If ITERATION_MAX is 0, we're almost done.
%  For just this case, do the triangulation, then exit.
%  Setting ITERATION_MAX to 0 means that we can see the initial mesh
%  before any of the improvements have been made.
%
  if ( iteration_max <= 0 )
    t = delaunayn ( p );
    triangulation_count = triangulation_count + 1;
    return
  end

  pold = inf;

  while ( iteration < iteration_max )
      
    iteration = iteration + 1;

    if ( mod ( iteration, 10 ) == 0 )
      fprintf ( 1, '  %d iterations, %d triangulations\n', ...
        iteration, triangulation_count );
    end
%
%  3. Retriangulation by the Delaunay algorithm.
%
%  Was there large enough movement to retriangulate?
%
%  If so, save the current positions, get the list of
%  Delaunay triangles, compute the centroids, and keep
%  the interior triangles (whose centroids are within the region).
%
    if ( ttol < max ( sqrt ( sum ( ( p - pold ).^2, 2 ) ) / h0 ) )
      N = size ( p, 1 );
      pold = p;  
      t = delaunayn ( p );
      triangulation_count = triangulation_count + 1;
      pmid = ( p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:) ) / 3; 
      t = t( feval ( fd, pmid, varargin{:} ) < -geps, : );
%
%  4. Describe each bar by a unique pair of nodes.
%
      bars = [ t(:,[1,2]); t(:,[1,3]); t(:,[2,3]) ];
      bars = unique ( sort ( bars, 2 ), 'rows' );
%
%  5. Graphical output of the current mesh.
%
      trimesh ( t, p(:,1), p(:,2), zeros(N,1), 'EdgeColor', 'b', 'Linewidth', 2 )
      title ( sprintf ( 'Iteration %d', iteration ) );
      view(2), axis equal, axis off, drawnow
%
%  Put a "pause" command here if you'd like to see each new mesh.
%
    end
    
%
%  6. Move mesh points based on bar lengths L and forces F.
%
%  Make a list of the bar vectors and lengths.
%  Set L0 to the desired lengths, F to the scalar bar forces,
%  and FVEC to the x, y components of the bar forces.
%
%  At the fixed positions, reset the force to 0.
%
    barvec = p(bars(:,1),:) - p(bars(:,2),:); 
    L = sqrt ( sum ( barvec.^2, 2 ) ); 
    hbars = feval ( fh, (p(bars(:,1),:)+p(bars(:,2),:))/2, varargin{:} );
    L0 = hbars * Fscale * sqrt ( sum(L.^2) / sum(hbars.^2) );
    F = max ( L0 - L, 0 );                               
    Fvec = F ./ L * [1,1] .* barvec;                     
    Ftot = full ( sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2) );
    Ftot(1:size(pfix,1),:) = 0;                         
    p = p + deltat * Ftot;      
%
%  7. Bring outside points back to the boundary.
%
%  Use the numerical gradient of FD to project points back to the boundary.
%
    d = feval ( fd, p, varargin{:} );
    ix = d > 0;
    dgradx = ( feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:}) - d(ix) ) / deps;
    dgrady = ( feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:}) - d(ix) ) / deps;
    p(ix,:) = p(ix,:) - [ d(ix) .* dgradx, d(ix) .* dgrady ];
%
%  8. Termination criterion: All interior nodes move less than dptol (scaled).
%
    if ( max ( sqrt ( sum ( deltat * Ftot ( d < -geps,:).^2, 2 ) ) / h0 ) < dptol )
      break; 
    end
%   if iteration == iteration_max
%       hold on
%        plot(pfix(5:24,1),pfix(5:24,2),'mo-',...
%                 'LineWidth',3,...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor',[.49 1 .63],...
%                 'MarkerSize',3);
%             hold on
%             plot(pfix(25:44,1),pfix(25:44,2),'mo-',...
%                 'LineWidth',3,...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor',[.49 1 .63],...
%                 'MarkerSize',3);
%             hold off
%   end
  end
  

  return
end
