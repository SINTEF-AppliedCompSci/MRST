%% Assessment of upscaling accuracy on 8x8x8 grid
% We upscale a small 8x8x8 cube with three different permeability
% realizations and measure the ratio between the outflows computed on the
% upscaled and the original fine-scale model.

mrstModule add coarsegrid incomp spe10 upscaling

%% Set up grid
warning('off','mrst:periodic_bc');
G = computeGeometry(cartGrid([8 8 8]));
coarse = [1 1 1];
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
fl = {'East', 'North', 'Top'};
fr = {'West', 'South', 'Bottom'};


%% Structures with boundary conditions
d = G.griddim;
[bcl,bcr,Dp]=deal(cell(d,1));
bcsides = {'XMin', 'XMax'; 'YMin', 'YMax'; 'ZMin', 'ZMax'};
for j = 1:d
   bcl{j} = pside([], G, bcsides{j, 1}, 0);
   bcr{j} = pside([], G, bcsides{j, 2}, 0);
   Dp{j}  = 0;
end
Dp{1}     = 4*barsa;
L         = max(G.faces.centroids)-min(G.faces.centroids);
[Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, Dp);
ofaces    = cell(d,1);
for j=1:d, ofaces{j} = bcp.face(bcp.tags==j); end

clc;
figure('Position',[440 460 800 300]);
disp('                      Harm/arith      Axial drop        Periodic');
for n=1:3
   disp('--------------------------------------------------------------------');
   switch n
      case 1
         rock.perm = repmat(300*milli*darcy, [G.cells.num, 1]);
         rock.perm( sin(4*sum(G.cells.centroids(:,[1 3]),2)*pi/norm(L([1 3])))>0 ) = 50*milli*darcy;
         disp('Layered:');
      case 2
         disp('Tarbert');
         rock = getSPE10rock(49:56,9:16,1:8);
      case 3
         disp('Upper Ness');
         rock = getSPE10rock(9:16,9:16,43:50);
   end
   hT = computeTrans(G, rock);

   subplot(1,3,n);
   plotCellData(G,log10(rock.perm(:,1)/(milli*darcy))); view(3); axis off; zoom(1.3);
   set(gca,'Clipping','off')
   
   %% Upscale values
   crock = cell(3,1);

   % Harmonic-arithmetic
   q   = ones(G.cells.num,1);
   vol = G.cells.volumes;
   for i=1:size(rock.perm,2)
      dims = G.cartDims; dims(i)=coarse(i);
      qq = partitionUI(G, dims);
      K = accumarray(qq,vol)./accumarray(qq,vol./rock.perm(:,i));
      crock{1}.perm(:,i) = accumarray(q,K(qq).*vol)./accumarray(q,vol);
   end
  
   % Pressure drop
   [v,dp] = deal(zeros(d, 1));
   for i=1:d
      bc = addBC([], bcl{i}.face, 'pressure', Dp{1});
      bc = addBC(bc, bcr{i}.face, 'pressure', Dp{2});

      xr = initResSol(G, 100*barsa, 1);
      xr = incompTPFA(xr, G, hT, fluid, 'bc', bc);

      v(i)  = sum(xr.flux(bcr{i}.face)) / ...
          sum(G.faces.areas(bcr{i}.face));
       dp(i) = Dp{1}/L(i);
   end
   crock{2}.perm = fluid.properties()*(v./dp)';
   
   % Periodic
   v  = nan(d);
   dp = Dp{1}*eye(d);
   nbcp = bcp;
   for i=1:d
      for j=1:d, nbcp.value(bcp.tags==j) = dp(j,i); end
      xr = initResSol(Gp, 100*barsa, 1);
      xr = incompTPFA(xr, Gp, hT, fluid, 'bcp', nbcp);
      for j=1:d
         v(j,i) = sum(xr.flux(ofaces{j})) / ...
            sum(Gp.faces.areas(ofaces{j}));
      end
   end
   dp = bsxfun(@rdivide, dp, L);
   K = fluid.properties()*v/dp;
   crock{3}.perm = abs(K([1 2 3 5 6 9]));

   %% Compare upscaling
   for j=1:numel(fl)
      % Fine-scale problem
      bc    = pside([], G, fl{j}, barsa);
      faces = bc.face;
      if j<3
         bc    = pside(bc, G, fr{j}, 0);
      else
         cf = gridCellFaces(G,1);
         cf = cf(any(G.faces.neighbors(cf,:)==0,2));
         bc = addBC([], cf, 'pressure', barsa);

         cf = gridCellFaces(G,G.cells.num);
         cf = cf(any(G.faces.neighbors(cf,:)==0,2));
         bc = addBC(bc, cf, 'pressure', 0);
      end
      x0    = initResSol(G,0);
      xr    = incompTPFA(x0, G, hT, fluid, 'bc', bc);
      flux  = sum(xr.flux(faces));

      fprintf('  %5s -> %6s: ', fl{j}, fr{j});
      for i=1:3
         nrock.perm = crock{i}.perm(ones(G.cells.num,1),:);
         chT    = computeTrans(G, nrock);
         x      = incompTPFA(x0, G, chT, fluid, 'bc', bc);
         fprintf('\t%f', sum(x.flux(faces)) / flux);
      end
      fprintf('\n');
   end
end
disp('--------------------------------------------------------------------');
cax = caxis;
for i=1:3, subplot(1,3,i), caxis(cax); end
h = colorbar;
set(h,'Position',[0.94 0.1100 0.02 0.8150], ...
   'XTick', .5, 'XTickLabel','[mD]', ...
   'YTickLabel',num2str(10.^(-2:3)'));