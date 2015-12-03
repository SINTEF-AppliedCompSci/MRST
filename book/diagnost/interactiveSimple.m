%% Set up the geomodel and specify wells
[nx,ny] = deal(64);
G = cartGrid([nx,ny,1],[500,250,10]);
G = computeGeometry(G);
p = gaussianField(G.cartDims(1:2), [0.2 0.4], [11 3], 2.5);
K = p.^3.*(1.5e-5)^2./(0.81*72*(1-p).^2);
rock.poro = p(:);
rock.perm = K(:);

%% Set up and solve flow problem, compute diagnostics
rate  = sum(poreVolume(G,rock))/(30*year);
n = 12;
W = addWell([],  G, rock, nx*n+n+1, ...
    'Type', 'rate', 'Comp_i', 1, 'name', 'I1', 'Val', rate/2);
W = addWell(W, G, rock, nx*n+n+1+nx-2*n, ...
    'Type','rate',  'Comp_i', 1, 'name', 'I2', 'Val', rate/2);
W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx-.3*nx), ...
    'Type','rate',  'Comp_i', 0, 'name', 'P1', 'Val', -rate/4);
W = addWell(W, G, rock, G.cells.num-(n-.5)*nx, ...
    'Type','rate',  'Comp_i', 0, 'name', 'P2', 'Val', -rate/2);
W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx+.3*nx), ...
    'Type','rate',  'Comp_i', 0, 'name', 'P3', 'Val', -rate/4);

%%
close all;
interactiveDiagnostics(G, rock, W);
axis normal tight; set(gca,'dataaspect',[40 25 20]); view(5,40); zoom(1.5);
set(gca,'Projection','Perspective');