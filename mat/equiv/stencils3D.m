function [A, epsx, epsy, epsz] = stencils3D(G,w1,w2,w3)

assert(strcmp(G.type(2), 'cartGrid') & G.griddim == 3, '3D cartesian grid required.')

hx = (max(G.nodes.coords(:,1))-min(G.nodes.coords(:,1)))/G.cartDims(1);
hy = (max(G.nodes.coords(:,2))-min(G.nodes.coords(:,2)))/G.cartDims(2);
hz = (max(G.nodes.coords(:,3))-min(G.nodes.coords(:,3)))/G.cartDims(3);

assert(w1 + w2 + w3 == 1, 'Wheights must sum to 1');

assert((w2 ~= 0 & hx == hy & hx == hz) | w2 == 0, ...
  'Distances betwwen nodes must be the same in all coordinate directions.')

if w1 == 1
    epsx = hy*hz/hx;
    epsy = hx*hz/hy;
    epsz = hx*hy/hz;
else
    epsx = 1;
    epsy = 1;
    epsz = 1;
end

nx = G.cartDims(1); ny = G.cartDims(2); nz = G.cartDims(3);
nN = G.nodes.num;

sub = ones(nN,1);
main = -2*ones(nN,1);

Ax = -spdiags([sub,main,sub], [-1,0,1], nN, nN);
Ay = -spdiags([sub,main,sub], [-(nx+1),0,nx+1], nN, nN);
Az = -spdiags([sub,main,sub], [-(nx+1)*(ny+1),0,(nx+1)*(ny+1)], nN, nN);

A = epsx*Ax + epsy*Ay + epsz*Az;