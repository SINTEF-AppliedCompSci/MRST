function [G,optPts,f,g] = CPG2D(pts,bnd,varargin)
% Construct a 2D centroidal Voronoi Diagram (CVD). The CVD is found by
% minimizing the CVD energy function using the lbfgs algorithm.
%
% SYNOPSIS:
%   [G, optPts, f, g] = createCVD(pts, bnd)
%   [...] = createCVD(..., 'Name1', Value1,'Name2', Value2,...)
%
% PARAMETERS:
%   p        - A n X 2 array of coordinates. Each coordinate coresponds to 
%              one Voronoi site.
%   bnd      - A k X 2 array of coordinates. Each coordinate coresponds to a
%              vertex in the polygon boundary. The coordinates must be
%              ordered clockwise or counter clockwise. 
%   fixedPts - OPTIONAL.
%              Default value []. A m X 2 array of fixed points. These
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
%              DefaultValue 10*eps. If step length is smaler than minStep,
%              the lbfgs returns
% RETURNS:
%   G        - Valid MRST grid definition.
%   optPts   - A n + m X 2 array with the optimal points. 
%   f        - The CVD energy function at each iteration of the lbfgs
%              algortihm. 
%   g        - The two norm of the gradient of the CVD energy function at
%              each step of the lbfgs algorithm
%
% EXAMPLE:
%   p = rand(30,2);
%   bnd = [0,0;0,1;1,1;1,0];
%   G = CPG2D(p,bnd);
%   plotGrid(G);
%   axis equal tight
%
% SEE ALSO
%   compositePebiGrid2D, pebi, pebiGrid2D, clippedPebi2D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

assert (size(pts, 2) == 2, ...
      ['Function ''%s'' is only supported in two ', ...
       'space dimensions.'], mfilename);
assert (size(bnd, 2) == 2, ...
      ['Function ''%s'' is only supported in two ', ...
       'space dimensions.'], mfilename);
   
opt = struct('fixedPts', [],                       ...
             'rho',     @(pts) ones(size(pts,1),1),...
             'tol',     1e-4,                      ...
             'maxIt',     40,                      ...
             'minStep',   10*eps,                  ...
             'storedVec',7);
opt = merge_options(opt,varargin{:});
fixedPts = opt.fixedPts;
rho = opt.rho;

inDomain = @(p) inpolygon(p(1:2:end),p(2:2:end),bnd(:,1), bnd(:,2));
if any(~inDomain(reshape(pts', [], 1)))
    error('Located sites outside domain')
end

nf = size(fixedPts,1);

p = [fixedPts; pts];
F = @(pts) objFunc(pts, bnd,nf,rho);

p = reshape(p',[],1);
[optPts,f,g] = lbfgs(p, F, 'tol',      opt.tol,      ...
                           'storedVec',opt.storedVec,...
                           'maxIt',    opt.maxIt,    ...
                           'minStep',  opt.minStep,  ...
                           'inDomain', inDomain);
optPts = reshape(optPts,2,[])';

G = clippedPebi2D(optPts, bnd);

end


function [f, g] = objFunc(p, bnd,nf,rho)
    pts = reshape(p,2,[])';

    G = clippedPebi2D(pts, bnd);
    G = computeGeometry(G);
    G = mrstGridWithFullMappings(G);

    intFun = @(x,i) sum(repmat(rho(x),1,2).*(x-repmat(pts(i,:),size(x,1),1)).^2,2);
    f = sum(polygonInt_v2(G,1:G.cells.num, intFun,7));

    massFun = @(x,i) rho(x);
    masses = polygonInt_v2(G,1:G.cells.num,massFun,7);
    g = reshape((2*repmat(masses,1,2).*(pts - G.cells.centroids))',[],1);
    g(1:2*nf) = 0;
end
