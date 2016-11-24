%% Simple Model of an Anticline
mrstModule add spe10 coarsegrid

%% Cartesian grid and rock parameters
[nx,ny,nz] = deal(30,30,15);
G = cartGrid([nx ny nz], [nx ny nz].*[20 10 2]*ft);

rock = getSPE10rock(1:nx,1:ny,nz:-1:1);
rock.poro(rock.poro==0) = 1e-5;
%rock.perm = 10*milli*darcy*ones(G.cells.num,1);
%rock.poro = 0.2*ones(G.cells.num,1);

%% Make anticline structure
x = G.nodes.coords(:,1);
y = G.nodes.coords(:,2);

normalize = @(x) (x-min(x))/(max(x)-min(x));
x = 2*normalize(x)-1;
y = 2*normalize(y)-1;
r = min(sqrt(x.^2 + y.^2),1);
dz = r.^2 ./ (r.^2 + (1-r).^4);

G.nodes.coords(:,3) = G.nodes.coords(:,3) + dz*(nz+1)*2*ft;
G = computeGeometry(G);

%% Set initial saturation
s = zeros(G.nodes.num,1);
s(G.nodes.coords(:,3)>=nz*2*ft) = 1;
[nodes,pos] = gridCellNodes(G, (1:G.cells.num).');
c = rldecode(1:G.cells.num, diff(pos), 2).';
A = sparse(c,nodes,1)*[s, ones([G.nodes.num,1])];
s = bsxfun(@rdivide, A(:,1:(end-1)), A(:,end));

%% Set wells
% producers
args = {'Type', 'bhp', 'Val', 100*barsa, 'Comp_i', [0 1]};
W = verticalWell([], G, rock, 14, 15, [], args{:}, 'name', 'P1');
W = verticalWell(W,  G, rock, 17, 16, [], args{:}, 'name', 'P2');

% injectors
args = {'Type', 'bhp', 'Val', 150*barsa, 'Comp_i', [1 0]};
W = verticalWell(W,  G, rock,  9,  9, [], args{:}, 'name', 'I1');
W = verticalWell(W,  G, rock,  9, 22, [], args{:}, 'name', 'I2');
W = verticalWell(W,  G, rock, 22,  9, [], args{:}, 'name', 'I3');
W = verticalWell(W,  G, rock, 22, 22, [], args{:}, 'name', 'I4');

