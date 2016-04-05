function sol = VEM2D_v3(G, f, k, bc, varargin)

nN = G.nodes.num;
nE = G.faces.num;
nK = G.cells.num;

nk = (k+1)*(k+2)/2;

if nargin < 5
    alpha = ones(nK,1);
else
    alpha = varargin{1};
    assert(size(alpha,1) == nK & size(alpha,2) == 1, ...
           ['Diemnsions of paramter matrix alpha must be ' ...
           'G.cells.num x 1'])
    G.cells.('alpha') = varargin{1};
end

[A,b] = VEM2D_glob_v3(G, f, k, bc, alpha);

% fprintf('Solving linear system ...\n')
% tic;

U = A\b;

% stop = toc;
% fprintf('Done in %f seconds.\n\n', stop);

nodeValues  = full( U( 1:nN)                            );
edgeValues  = full( U((1:nE*(k-1))       + nN)          );
cellMoments = full( U((1:nK*k*(k-1)/2) + nN + nE*(k-1)) );

sol = struct(...
             'nodeValues' , {nodeValues} , ...
             'edgeValues' , {edgeValues} , ...
             'cellMoments', {cellMoments}     );
         
end