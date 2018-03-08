mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil blackoil-sequential gasinjection

%%

[state0, model, schedule]  = setupSPE10_AD('layers', 10);
G = model.G;
rock = model.rock;

xmax = max(G.nodes.coords);
G = cartGrid(G.cartDims(1:2), xmax(1:2));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);
rock.perm = model.rock.perm(:, 1:2);
fluid = model.fluid;
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1000, 800]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1]                 );

model = TwoPhaseOilWaterModel(G, rock, fluid);

modelFV = getSequentialModelFromFI(model);
modelDG = modelFV;
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid);
% modelDG.transportModel.AutoDiffBackend = DiagonalAutoDiffBackend();
disc = DGDiscretization(modelDG.transportModel, 2, 'degree', 1, 'basis', 'legendre');

modelDG.transportModel.disc = disc;

state0.s = repmat([0,1], G.cells.num,1);
state0 = disc.assignDofFromState(state0);

%%

[wsDG, statesDG, rep] = simulateScheduleAD(state0, modelDG, schedule);

%%

[wsFV, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

figure('position', [-2000, 0, 500, 1000]);
plotToolbar(G, statesDG);

figure('position', [-2000, 0, 500, 1000]);
plotToolbar(G, statesFV);


%%

close all

fig = figure('Position', [-2000, 0, 1500, 1000]);

states = {statesDG, statesFV};
titles = {'dG(1)', 'FV'};

for mNo = 1:numel(states)
    subplot(1,3,mNo)
    h(mNo) = plotCellData(G, states{mNo}{1}.s(:,1), 'edgec', 'none');
    colormap jet
    axis equal off
    text(290, 20, titles{mNo}, 'fontsize', 25, 'color', 'w'); 
    caxis([0,1])
end

% set(fig, 'Units', 'pixels');
% pos = get(fig, 'Position');
% set(fig, 'Units', 'normalized');
M = struct('cdata',[],'colormap',[]);

for sNo = 1:numel(schedule.step.val)
    
    for mNo = 1:numel(states)
        h(mNo).CData = states{mNo}{sNo}.s(:,1);
        pause(0.01);
    end
    
%     dx = 10;
%     dy = 10;
    rect = [0, 0, fig.Position(3:4)];
%     rect = pos;
    M(sNo) = getframe(fig, rect);
    
end

%%

pth = mrstPath('dg');
name = 'spe10';
duration = 10;
vo = VideoWriter(fullfile(pth, name));
vo.FrameRate = numel(states{1})/duration;
open(vo);

writeVideo(vo, M);

close(vo)

