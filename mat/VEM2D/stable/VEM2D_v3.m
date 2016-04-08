function [sol, varargout] = VEM2D_v3(G, f, k, bc, varargin)

nN = G.nodes.num;
nE = G.faces.num;
nK = G.cells.num;

opt = struct('alpha',  ones(nK,1), ...
             'src', []            );         
opt = merge_options(opt, varargin{:});
alpha = opt.alpha;
src = opt.src;

assert(size(alpha,1) == nK & size(alpha,2) == 1, ...
       'Dimensions of paramter matrix alpha must be G.cells.num x 1')

projectors = false;
if nargout == 2
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
         
if nargout == 2
    varargout(1) = {PNstarT};
end
         
end