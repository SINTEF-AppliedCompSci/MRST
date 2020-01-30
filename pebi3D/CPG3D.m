function [G,optPts,f,g] = CPG3D(pts, bnd, varargin)
% Construct a 3D centroidal Pebi Grid(CPG). The CPG is found by
% minimizing the CPG energy function using the lbfgs algorithm.
%
% SYNOPSIS:
%   [G, optPts, f, g] = CPG3D(pts, bnd)
%   [...] = CPG3D(..., 'Name1', Value1,'Name2', Value2,...)
%
% PARAMETERS:
%   p        - A n X 3 array of coordinates. Initial guess for Voronoi
%              sites.
%   bnd      - A k X 3 array of coordinates. Each coordinate coresponds to 
%              a vertex of the polyhedron boundary. The boundary is assumed
%              to be the convex hull of bnd
%   fixedPts - OPTIONAL.
%              Default value []. A m X 3 array of fixed points. These
%              points are not moved during the optimization algorithm.
%   rho      - OPTIONAL.
%              Default value @(pts) ones(size(pts,1),1). Function that
%              gives the mass density function for the Voronoi cells
%   tol      - OPTIONAL.
%              Default value 1e-4. The convergence tolerance passed on to
%              the lbfgs algorithm. The convergence requirement is 
%              grad F_k <= tol * grad F_0
%   storedVec- OPTIONAL
%              Default value 7. Number of iterations remembered by the
%              lbfgs algorithm. Normal values are between 3 and 20. A
%              higher value will in general result in fewer iterations to
%              convergence, but it is a trade off with higher computational
%              cost at each iteration.
%   maxIt    - OPTIONAL.
%              Default value 200. Maximum number of iteration in the 
%              lbfgs algorithm.
%   minStep  - OPTIONAL
%              DefaultValue 10*eps. if step length is smaler than minStep,
%              the lbfgs returns
% RETURNS:
%   G        - Valid MRST grid definition.
%   optPts   - A n + m X 3 array with the optimal points. 
%   f        - The CVD energy function at each iteration of the lbfgs
%              algortihm. 
%   g        - The two norm of the gradient of the CVD energy function at
%              each step of the lbfgs algorithm
%
% EXAMPLE:
%   p = rand(30,3);
%   bnd = [0,0,0;0,1,0;1,1,0;1,0,0;...
%          0,0,1;0,1,1;1,1,1;1,0,1];
%   G = CPG3D(p,bnd);
%   figure();
%   plotGrid(G);
%   axis equal tight
%
% SEE ALSO
%   compositePebiGrid3D, pebiGrid2D, clippedPebi3D, mirroredPebi3D, CPG2D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


assert (size(pts, 2) == 3, ...
      ['Function ''%s'' is only supported in three ', ...
       'space dimensions.'], mfilename);
assert (size(bnd, 2) == 3, ...
      ['Function ''%s'' is only supported in three ', ...
       'space dimensions.'], mfilename);

opt = struct('density',   @(x) ones(size(x,1),1),...
             'storedVec', 5,                     ...
             'tol',       1e-6,                  ...         
             'maxIt',     40,                    ...
             'minStep',   10*eps,                ...
             'fixedPts',  []);
opt = merge_options(opt, varargin{:});

nf = size(opt.fixedPts,1);
dt = delaunayTriangulation(bnd);
pts = [opt.fixedPts;pts];
F = @(pts) objectiveFunc(pts, bnd,nf, opt.density);
inDomain = @(p) repmat(~isnan(pointLocation(dt,reshape(p,3,[])')),3,1); % The returned logical is shuffled, but it does not matter as lbfgs checks for any point outside domain.
pts = reshape(pts',[],1);

[optPts, f, g] = lbfgs(pts, F, 'storedVec',opt.storedVec, ...
                               'maxIt',    opt.maxIt,...
                               'minStep',  opt.minStep,  ...
                               'tol',      opt.tol, ...
                               'inDomain', inDomain);

optPts = reshape(optPts,3,[])';
G = mirroredPebi3D(optPts, bnd);
end

function [f, g] = objectiveFunc(pts, bndr, nf, rho)
    pts = reshape(pts,3,[])';

    G = mirroredPebi3D(pts, bndr);
    G = computeGeometry(G);
    G = mrstGridWithFullMappings(G);

    intFun = @(x,i) sum(repmat(rho(x),1,3).*(x-repmat(pts(i,:),size(x,1),1)).^2,2);
    f = sum(polyhedronInt_upr(G,1:G.cells.num, intFun,3));

    massFun = @(x,i) rho(x);
    masses = polyhedronInt_upr(G,1:G.cells.num,massFun,3);
    g = reshape((2*repmat(masses,1,3).*(pts - G.cells.centroids))',[],1);
    g(1:3*nf) = 0;
end
