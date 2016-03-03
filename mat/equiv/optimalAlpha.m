clc; clear; close all;
addpath('../')
addpath('/home/strene/Documents/master/coop/pebiGridding/voronoi3D')

n = 3;
gridLim = [1,1,1];

G = cartGrid([n,n,n],gridLim);
% G = voronoiCube(50, gridLim);


%--------------------------------------------------------------------------
%   -\delta u = 0,
%           u = 1/(2\pi||x-C||)
%--------------------------------------------------------------------------
f = @(X) zeros(size(X,1),1);
C = -[.2,.2,.2];
gD = @(X) -1./(2*pi*sqrt(sum((X-repmat(C,size(X,1),1)).^2,2)));
k = 1;


G = computeVEMGeometry(G,f,k);

boundaryFaces = (1:G.faces.num)';
                           
boundaryFaces = boundaryFaces( G.faces.neighbors(:,1) == 0 | ...
                               G.faces.neighbors(:,2) == 0 );

bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryFaces}}, 'bcType', {{'dir'}});

u = gD(G.nodes.coords);

beta = 4*1.8395265*10e-5;

alphaOpt = gridLim(1)/n*(1/20*beta +1/5)*3*ones(G.cells.num,1);

alphaVec = (.01:0.01:1)'*sqrt(3);
nAlph = size(alphaVec,1);
errVec = zeros(nAlph,1);

for i = 1:nAlph

fprintf('Iteration %d\n\n', i)
    
alpha = alphaVec(i)*G.cells.diameters;
G = computeVEMGeometry(G,f,k);

sol = VEM3D(G,f,bc,k, alpha);
U = [sol.nodeValues; sol.edgeValues; sol.faceMoments; sol.cellMoments];

errVec(i) = norm(U - u,2);

end

plot(alphaVec, errVec)
alphaMin = alphaVec(errVec == min(errVec));

%   RESULTS
%   n = 5, xMax = 1: alpha = 0.132, alphaOpt = 0.12, beta = 4.