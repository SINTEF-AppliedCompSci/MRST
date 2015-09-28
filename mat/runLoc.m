run('../../matlab/project-mechanics-fractures/mystartup.m')

% G = cartGrid([1,1]);
% G = mrstGridWithFullMappings(G);
% G = computeGeometry(G);
% 
% iNodes = [G.cells.nodePos(1):G.cells.nodePos(2)-1];
% nodes = G.cells.nodes(iNodes);
% X = G.nodes.coords(nodes,:);
% 
% Sl = locS(X);

d = sqrt(2);
X = [0.0, 0.0
     1.0, 0.0;
     1.0, 1.0;
     0.0, 1.0];

locS(X)