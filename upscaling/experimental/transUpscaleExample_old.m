require coarsegrid agglom

%% Set up grid and create partition
g=cartGrid([4,1,1]);
g=computeGeometry(g);
rock.perm=ones(g.cells.num,1);
% Calculate transmisibility on original grid
T=computeTrans(g,rock);
% Make partition for coarse grid
%p = partitionUI(g, floor(g.cartDims/2));
p = partitionUI(g, [g.cartDims(1), 1 1]);

p = processPartition(g, p, 'Verbose', true);
%% Generate coarse grid and define wells
cg  = generateCoarseGrid(g, p);
W=[];

% Add four wells which will be used for upscaling in linearly independent
% combinations.
if(false)
W = verticalWell(W, g, rock,  1,   1, (1:g.cartDims(3)),     ...
                 'Type', 'bhp', 'Val', 300*barsa, ...
                  'Radius', 0.125, 'Name', 'P1');
W = verticalWell(W, g, rock,  1,   g.cartDims(2), (1:g.cartDims(3)),     ...
                 'Type', 'bhp', 'Val', 300*barsa, ...
                  'Radius', 0.125, 'Name', 'P2');

W = verticalWell(W, g, rock, 1, g.cartDims(2), (1:g.cartDims(3)),     ...
                 'Type', 'bhp', 'Val', 300*barsa, ...
                  'Radius', 0.125, 'Name', 'P3');
W = verticalWell(W, g, rock, g.cartDims(1), g.cartDims(2), (1:g.cartDims(3)),     ...
                 'Type', 'bhp', 'Val', 300*barsa, ...
                  'Radius', 0.125, 'Name', 'P4');
W = verticalWell(W, g, rock, floor(g.cartDims(1)/2)+1, floor(g.cartDims(2)/2)+1, (1:floor(g.cartDims(3)/2)+1),     ...
                 'Type', 'bhp', 'Val', 500*barsa, ...
                  'Radius', 0.125, 'Name', 'I4');
              else
W = verticalWell(W, g, rock,  1,   1, (1:g.cartDims(3)),     ...
                 'Type', 'bhp', 'Val', 300*barsa, ...
                 'Radius', 0.125, 'Name', 'P1');
W = verticalWell(W, g, rock, g.cartDims(1), g.cartDims(2), (1:g.cartDims(3)),     ...
                 'Type', 'bhp', 'Val', 500*barsa, ...
                  'Radius', 0.125, 'Name', 'P4');
              end
clf;
plotGrid(g, 'FaceAlpha', 0)
plotWell(g, W)
%% Upscale transmissibility and plot the result
require agglom
[HT_cg, T_cg, Wcg, upscaled]=upscaleTrans(cg, T, 'match_method','lsq_flux',...
                                                 'bc_method','wells_simple',...
                                                 'wells',W);
figure(1),clf
for i=1:numel(upscaled)
   subplot(numel(upscaled),1,i)
   plotCellData(cg.parent,upscaled(i).state.pressure(cg.partition)./barsa);colorbar,view(3)
   outlineCoarseGrid(cg.parent, cg.partition,'LineWidth',3);
end

%%
% test options
% get an other well configuration

g=cg.parent;
W2 = verticalWell([], g, rock, g.cartDims(1), g.cartDims(2), (1:g.cartDims(3)),     ...
                 'Type', 'bhp', 'Val', 300*barsa, ...
                  'Radius', 0.125, 'Name', 'P4');
W2 = verticalWell(W2, g, rock, floor(g.cartDims(1)/2)+1, floor(g.cartDims(2)/2)+1, (1:floor(g.cartDims(3)/2)+1),     ...
                 'Type', 'bhp', 'Val', 500*barsa, ...
                  'Radius', 0.125, 'Name', 'I4');
% bc for setting outside
bc = cell(3,1);
bc{1} = pside([], cg.parent, 'West', 100*barsa);
bc{1} = pside(bc{1}, cg.parent, 'East', 200*barsa);
bc{2} = pside([], cg.parent, 'North', 100*barsa);
bc{2} = pside(bc{2}, cg.parent, 'South', 200*barsa);
bc{3} = pside([], cg.parent, 'Top', 100*barsa);
bc{3} = pside(bc{3}, cg.parent, 'Bottom', 200*barsa);
match_methods={'lsq_flux','max_flux'};
bc_methods={'bc_simple','bc','wells_simple','wells','bc'};
bcs={[],bc,[],[],{bc{1:2}}};
wells={[],[],W,{W,W2},[]};
%%
for i=1:numel(match_methods)
   for j=1:numel(bc_methods)
      disp(['Testing ', match_methods{i},' with ', bc_methods{j}]);
      [HT_cg,T_cg,Wcg,upscaled]=upscaleTrans(cg,T,'match_method',match_methods{i},...
         'bc_method',bc_methods{j},'wells',wells{j},'bc',bcs{j});
     disp(['Number of trans ',num2str(numel(T_cg))]);
     disp(['Number of negative trans ',num2str(sum(T_cg<0))]);
   end
end
