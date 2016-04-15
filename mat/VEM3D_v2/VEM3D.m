function [sol, A, b,G] = VEM3D(G, f, bc, k, varargin)

addpath('../');

nN = G.nodes.num;
nE = G.edges.num;
nF = G.faces.num;
nK = G.cells.num;

nk = (k+1)*(k+2)*(k+3)/6;

% if nargin < 5
%     G.cells.('alpha') = repmat(G.cells.diameters, 1, nAlpha);
% else
%     alpha = varargin{1};
%     assert(size(alpha,1) == nK & size(alpha,2) == nk, ...
%            ['Diemnsions of paramter matrix alpha must be ' ...
%            'G.cells.num x k*(k+1)*(k+2)/6'])
%     G.cells.('alpha') = varargin{1};
% end

[A,b,G] = VEM3D_glob(G, f, bc, k);

% fprintf('Solving linear system ...\n')
% tic;

U = A\b;

% stop = toc;
% fprintf('Done in %f seconds.\n\n', stop);

nodeValues  = full( U( 1:nN)                                             );
edgeValues  = full( U((1:nE*(k-1))       + nN)                           );
faceMoments = full( U((1:nF*k*(k-1)/2)   + nN + nE*(k-1))                );
cellMoments = full( U((1:nK*k*(k^2-1)/6) + nN + nE*(k-1) + nF*k*(k-1)/2) );

sol = struct(...
             'nodeValues' , {nodeValues} , ...
             'edgeValues' , {edgeValues} , ...
             'faceMoments', {faceMoments}, ...
             'cellMoments', {cellMoments}     );
         
end