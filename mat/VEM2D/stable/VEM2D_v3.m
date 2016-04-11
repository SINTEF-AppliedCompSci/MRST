function [sol, varargout] = VEM2D_v3(G, f, k, bc, varargin)

nN = G.nodes.num;
nE = G.faces.num;
nK = G.cells.num;

opt = struct('alpha'       , ones(nK,1), ...
             'src'         , []        , ...
             'projectors'  , false     , ...
             'cellAverages', false     );
opt = merge_options(opt, varargin{:});
alpha = opt.alpha;
src = opt.src;
projectors = opt.projectors;
cellAverages = opt.cellAverages;

assert(size(alpha,1) == nK & size(alpha,2) == 1, ...
       'Dimensions of paramter matrix alpha must be G.cells.num x 1')

if cellAverages
    projectors = true;
end

[A,b,PNstarT] = VEM2D_glob_v3(G, f, k, bc, alpha, src, projectors);

fprintf('Solving linear system ...\n')
tic;

U = A\b;

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

nodeValues  = full( U( 1:nN)                            );
edgeValues  = full( U((1:nE*(k-1)) + nN)                );
cellMoments = full( U((1:nK*k*(k-1)/2) + nN + nE*(k-1)) );

sol = struct(...
             'nodeValues' , {nodeValues} , ...
             'edgeValues' , {edgeValues} , ...
             'cellMoments', {cellMoments}     );
if projectors
    G.('PNstarT') = PNstarT;
    if k == 1
        PNstarPos = [1, cumsum( diff(G.cells.nodePos')) + 1];
    elseif k == 2
        PNstarPos = [1, cumsum( diff(G.cells.nodePos') + ...
                                diff(G.cells.facePos') + 1) + 1];
    end
    G.PNstarPos = PNstarPos;
    varargout(1) = {G};
end

if cellAverages && k == 1
    cellMoments = calculateCellAverages(G, nodeValues);
    sol.cellMoments = cellMoments;
end
         
end