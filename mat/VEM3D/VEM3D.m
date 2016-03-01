function [sol, G] = VEM3D(G, f, bc, k)


[A,b] = VEM3D_glob(G,f,bc,k);

fprintf('Solving linear system ...\n')
tic;

U = A\b;

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

nN = G.nodes.num;
nE = G.edges.num;
nF = G.faces.num;
nK = G.cells.num;

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