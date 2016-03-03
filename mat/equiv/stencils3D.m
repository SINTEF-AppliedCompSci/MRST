function [A_FD, epsx, epsy, epsz] = stencils3D(G,w1,w2,w3)

assert(strcmp(G.type(2), 'cartGrid') & G.griddim == 3, '3D cartesian grid required.')

hx = (max(G.nodes.coords(:,1))-min(G.nodes.coords(:,1)))/G.cartDims(1);
hy = (max(G.nodes.coords(:,2))-min(G.nodes.coords(:,2)))/G.cartDims(2);
hz = (max(G.nodes.coords(:,3))-min(G.nodes.coords(:,3)))/G.cartDims(3);

tol = 10e-10;
assert((w2 ~= 0 & abs(hx - hy)/hx < tol & abs(hx - hz)/hx < tol) | w2 == 0, ...
  'Distances betwwen nodes must be the same in all coordinate directions.')

epsx = (hy*hz)/hx;
epsy = (hx*hz)/hy;
epsz = (hx*hy)/hz;

nx = G.cartDims(1); ny = G.cartDims(2); nz = G.cartDims(3);
nN = G.nodes.num;

sub = ones(nN,1);
main = -2*ones(nN,1);

x = 1;
y = (nx+1);
z = (nx+1)*(ny+1);

Ax = -spdiags([sub,main,sub], [-x,0,x], nN, nN);
Ay = -spdiags([sub,main,sub], [-y,0,y], nN, nN);
Az = -spdiags([sub,main,sub], [-z,0,z], nN, nN);

A = epsx*Ax + epsy*Ay + epsz*Az;

yx =  (nx+1) + (nx+1)*(ny+1);
zx = -(nx+1) + (nx+1)*(ny+1);
xy = 1 + (nx+1)*(ny+1);
zy = -1 + (nx+1)*(ny+1);
xz = 1 + (nx+1);
yz = - 1 + (nx+1);


dyxdyx = spdiags([sub,main,sub], [-yx, 0, yx], nN, nN);
dzxdzx = spdiags([sub,main,sub], [-zx, 0, zx], nN, nN);
dxydxy = spdiags([sub,main,sub], [-xy, 0, xy], nN, nN);
dzydzy = spdiags([sub,main,sub], [-zy, 0, zy], nN, nN);
dxzdxz = spdiags([sub,main,sub], [-xz, 0, xz], nN, nN);
dyzdyz = spdiags([sub,main,sub], [-yz, 0, yz], nN, nN);

Bx = -epsx/2*(dyxdyx + dzxdzx) + epsx*Ax;
By = -epsx/2*(dxydxy + dzydzy) + epsx*Ay;
Bz = -epsx/2*(dxzdxz + dyzdyz) + epsz*Az;

d1 = 1 + (nx+1) + (nx+1)*(ny+1);
d2 = 1 + (nx+1) - (nx+1)*(ny+1);
d3 = 1 - (nx+1) + (nx+1)*(ny+1);
d4 = 1 - (nx+1) - (nx+1)*(ny+1);

d1d1 = spdiags([sub,main,sub], [-d1, 0, d1], nN, nN);
d2d2 = spdiags([sub,main,sub], [-d2, 0, d2], nN, nN);
d3d3 = spdiags([sub,main,sub], [-d3, 0, d3], nN, nN);
d4d4 = spdiags([sub,main,sub], [-d4, 0, d4], nN, nN);

d1d2 = spdiags([sub,-sub,-sub,sub], [-xz, -z, z, xz], nN, nN);
d1d3 = spdiags([sub,-sub,-sub,sub], [-zy, -y, y, zy], nN, nN);
d1d4 = spdiags([sub,-sub,-sub,sub], [-x, -yx, yx, x], nN, nN);
d2d3 = spdiags([sub,-sub,-sub,sub], [-x, -zx, zx, x], nN, nN);
d2d4 = spdiags([sub,-sub,-sub,sub], [-xy, -y, y, xy], nN, nN);
d3d4 = spdiags([sub,-sub,-sub,sub], [-yz, -z, z, yz], nN, nN);

B1 = -epsx/2*(d1d1+d2d2+d3d3+d2d3-d1d3-d1d2);
B2 = -epsx/2*(d1d1+d2d2+d4d4+d1d4-d2d4-d1d2);
B3 = -epsx/2*(d1d1+d3d3+d4d4+d1d4-d1d3-d3d4);
B4 = -epsx/2*(d2d2+d3d3+d4d4+d2d3-d2d4-d3d4);
            
A_FD = w1*A + w2/3*(Bx + By + Bz) + w3/4*(B1 + B2 + B3 + B4);

% d2d1 = spdiags([sub,-sub,-sub,sub], ...
%                [-xz, -(nx+1)*(ny+1), (nx+1)*(ny+1), xz], nN, nN);
% d3d1 = spdiags([sub,-sub,-sub,sub], ...
%                [-zy, -(nx+1), (nx+1), zy], nN, nN);
% d4d1 = spdiags([sub,-sub,-sub,sub], ...
%                [-1, -yx, yx, 1], nN, nN);
% d3d2 = spdiags([sub,-sub,-sub,sub], ...
%                [-1, zx, -zx, 1], nN, nN);
% d4d2 = spdiags([sub,-sub,-sub,sub], ...
%                [-xy, -(nx+1), (nx+1), xy], nN, nN);
% d4d3 = spdiags([sub,-sub,-sub,sub], ...
%                [yz, -(nx+1)*(ny+1), (nx+1)*(ny+1), -yz], nN, nN);
% 
% B1 = -1/2*epsx*(spdiags([sub,sub,sub,3*main,sub,sub,sub], ...
%                  [-d1, -d2, -d3, 0, d1, d2, d3], ...
%                    nN, nN) ...
%                 - d2d1 - d3d1 + d3d2);
% 
% B2 = -1/2*epsx*(spdiags([sub,sub,sub,3*main,sub,sub,sub], ...
%                  [-d1, -d2, -d4, 0, d1, d2, d4], ...
%                    nN, nN) ...
%                 + d4d1 - d4d2 - d2d1);
% 
% B3 = -1/2*epsx*(spdiags([sub,sub,sub,3*main,sub,sub,sub], ...
%                  [-d1, -d3, -d4, 0, d1, d3, d4], ...
%                    nN, nN) ...
%                 + d4d1 - d3d1 - d4d3);
% 
% B4 = -1/2*epsx*(spdiags([sub,sub,sub,3*main,sub,sub,sub], ...
%                  [-d2, -d3, -d4, 0, d2, d3, d4], ...
%                    nN, nN) ...
%                 + d3d2 - d4d2 -d4d3);

% B11 = -1/3*epsx*(spdiags([sub,sub,sub,3*main,sub,sub,sub], ...
%                  [-d1, -d2, -d3, 0, d1, d2, d3], ...
%                    nN, nN));
%                
% B22 = -1/3*epsx*(spdiags([sub,sub,sub,3*main,sub,sub,sub], ...
%                  [-d1, -d2, -d4, 0, d1, d2, d4], ...
%                    nN, nN));
%                
% B33 = -1/3*epsx*(spdiags([sub,sub,sub,3*main,sub,sub,sub], ...
%                  [-d1, -d3, -d4, 0, d1, d3, d4], ...
%                    nN, nN));
%                
% B44 = -1/3*epsx*(spdiags([sub,sub,sub,3*main,sub,sub,sub], ...
%                  [-d2, -d3, -d4, 0, d2, d3, d4], ...
%                    nN, nN));
% B1+B2+B3+B4-(B11+B22+B33+B44)


% Bx = -.5*epsx*spdiags([sub,sub,2*main,sub,sub], ...
%                  [-yx, -zx, 0, yx, zx], nN, nN) + epsx*Ax;
%  
% By = -.5*epsx*spdiags([sub,sub,2*main,sub,sub], ...
%                  [-xy, -zy, 0, xy, zy], nN, nN) + epsx*Ay;
% 
% Bz = -.5*epsx*spdiags([sub,sub,2*main,sub,sub], ...
%                  [-xz, -yz, 0, xz, yz], nN, nN) + epsx*Az;
