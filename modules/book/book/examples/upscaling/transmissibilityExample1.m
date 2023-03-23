%% Upscale a problem with a diagonal trend
% In this example we will upscale a problem with diagonal trend and compute
% the transmissibility associated with the coarse interface between two
% blocks in the x-direction

mrstModule add coarsegrid incomp;

%% Set up model
[Lx,Ly] = deal(200,100);
[nx,ny] = deal(30,15);
Dp      = [1 0]*barsa;
G = computeGeometry(cartGrid([nx ny], [Lx Ly]));
rock = makeRock(G, 30*milli*darcy, 0.3);
rock.perm( sin(4*sum(G.cells.centroids(:,:),2)*pi/Lx)>0 ) = 10*milli*darcy;
rock.perm = rock.perm .* (1+.2*rand(size(rock.perm)));
hT    = computeTrans(G, rock);
fluid = initSingleFluid('mu' ,1, 'rho', 1);

clf, plotCellData(G,rock.perm/(milli*darcy),'EdgeColor','none'); colorbar

%% Solve fine-scale problem
bc   = pside([], G, 'XMin', Dp(1));
bc   = pside(bc, G, 'XMax', Dp(2));
xr   = initResSol(G, 100*barsa, 1);
xr   = incompTPFA(xr, G, hT, fluid, 'bc', bc);

%% Loop over various partititions
pow=4:-1:-2;
x = G.cells.centroids(:,1); 
y = G.cells.centroids(:,2);
sp = [1 2 3 4 5 6 11];
clf, set(gcf,'Position',[360 500 930 300]);
for avg=[false, true]
   T = zeros(size(pow));
   for n=1:numel(pow);

      % Find faces on coarse interface
      q =  (y> 50 + 2^pow(n)*(100-x))+1; q = q(:);
      CG = generateCoarseGrid(G, q);
      CG = coarsenGeometry(CG);
      i = find(~any(CG.faces.neighbors==0,2));
      faces = CG.faces.fconn(CG.faces.connPos(i):CG.faces.connPos(i+1)-1);
      sgn   = 2*(CG.faces.neighbors(i, 1) == 1) - 1;

      % Upscale transmissibility
      flux = sgn*sum(xr.flux(faces));
      mu   = fluid.properties();
      if avg
         P = accumarray(q,xr.pressure)./accumarray(q,1);
      else
         cells = findEnclosingCell(G,CG.cells.centroids);
         P = xr.pressure(cells);
      end
      T(n) = mu*flux/(P(1) - P(2));

      if avg, continue; end
      
      % Plot grid
      subplot(3,5,sp(n));
      plotCellData(G,rock.perm(:,1),'EdgeColor','none');
      plotGrid(CG,'FaceColor','none','LineWidth',1);
      plotGrid(G, cells, 'FaceColor','k', 'EdgeColor','none');
      axis equal off; axis([-1 201 -1 101]); title(num2str(n));
   end
   flag = true(15,1); flag(sp) = false;
   subplot(3,5,find(flag)); hold on
   if avg
      plot(1:numel(pow),T,'--sb','LineWidth',1.5, 'MarkerSize',8, ...
         'MarkerEdgeColor','k','MarkerFaceColor','r');
   else
      plot(1:numel(pow),T,'--or','LineWidth',1.5, 'MarkerSize',8, ...
         'MarkerEdgeColor','k','MarkerFaceColor','b');
   end
   hold off
end
%axis([.5 7.5 [1.6 2.8]*1e-14]); legend('Centroid','Average',2);
colormap(.5*parula+.5*ones(size(parula)));
