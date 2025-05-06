function showMonotonicityExample(G, K, seed, flag)
% Compare monotonicity and grid orientation for TPFA, mimetic, and MPFA
%
% SYNOPSIS: 
%   showMonotonicityExample(G, K, seed)
%   showMonotonicityExample(G, K, seed, quiver)
%
% PARAMETERS:
%   G     - grid structure
%   K     - 2x2 matrix with anisotropic permeability [mD]
%   seed  - start point for tracing streamlines
%   flag  - show quiver plot of cell velocities (default: false)
%
% DESCRIPTION:
%   Simulate single-phase flow in a square domain given by the grid G for a
%   homogeneous, anisotropic permeability K. The flow is driven by a
%   pressure drop from left to right (Dirichlet boundary conditions pL = 1
%   bar and pR = 0) and no-flow boundary conditions along the top and
%   bottom. Solutions are computed for the TPFA method, the mimetic method,
%   and the MPFA-O method and visualized by tracing streamlines backward
%   and forward from the positions in seed and by plotting pressure values
%   at cell centers.

MODS = mrstModule;

if nargin<4, flag = false; end

G = computeGeometry(G);
rock.perm = repmat([K(1,1), K(1,2), K(2,2)]*milli*darcy,[G.cells.num, 1]);

% Boundary conditions
bc  = pside([], G, 'left',  2);
bc  = pside(bc, G, 'right', 1);

% Fluid object
fluid = initSingleFluid('mu' ,    1*centi*poise     , ...
                        'rho', 1014*kilogram/meter^3);

% TPFA method
mrstModule add incomp
xr0     = initResSol(G, 0, 0);
hT      = computeTrans(G, rock);
xr{1}   = incompTPFA(xr0, G, hT, fluid, 'bc', bc);
name{1} = 'TPFA';

% Mimetic method
mrstModule add mimetic
S       = computeMimeticIP(G, rock,'InnerProduct', 'ip_quasitpf');
xr{2}   = incompMimetic(xr0, G, S, fluid, 'bc', bc);
name{2} = 'Mimetic';

% MPFA-O method
mrstModule add mpfa
mT      = computeMultiPointTrans(G, rock);
xr{3}   = incompMPFA(xr0, G, mT, fluid, 'bc', bc);
name{3} = 'MPFA-O';

%% Plot pressure solutions
clf, set(gcf,'Position',[300 250 1000 500]);
c = G.cells.centroids;
for i=1:3, 
   subplot(2,3,i)
   plotCellData(G, xr{i}.pressure, 'EdgeColor', 'none');
   title(name{i});
   set(gca,'XTick',[],'YTick',[]);
   axis equal tight
   if i==1, cax = caxis; ax = axis; end
   caxis(cax), axis(ax);
   if flag,
      v = faceFlux2cellVelocity(G, xr{i}.flux);
      hold on, 
      quiver(c(seed,1),c(seed,2),v(seed,1),v(seed,2),'Color','k','LineWidth',1);
      hold off
   end
end
colorbar('Position',[.92 .58 .02 .34])

% Plot streamlines
mrstModule add streamlines
for i=1:3
   subplot(2,3,3+i)
   Sf = pollock(G, xr{i}, seed, 'substeps', 1);
   Sb = pollock(G, xr{i}, seed, 'substeps', 1, 'reverse', true);
   hf=streamline(Sf);
   hb=streamline(Sb);
   set([hf; hb],'Color','k');
   p=get(gca,'Position'); p(2)=p(2)+.15; set(gca,'Position', p);
   box on;
   set(gca,'XTick',[],'YTick',[]); axis equal tight; axis(ax);
end

% Clean up module list
mrstModule clear
mrstModule('add', MODS{:})
