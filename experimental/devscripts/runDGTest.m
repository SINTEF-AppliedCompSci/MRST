mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista upr compositional spe10
mrstVerbose on

%%

setup  = getDGTestCase('simple1d', 'n', 20, 'nkr', 2); %#ok

%%

setup  = getDGTestCase('qfs_wo_2d', 'useMomentFitting', false, 'nkr', 2, 'k', {[0,0], [0,0; 1,0; 0,1; 1,1]}, 'degree', 1, 'rotate', true, 'n', 20); %#ok

%%

setup = getDGTestCase('qfs_wog_3d', 'useMomentFitting', false, 'degree', [0,0,0; 1,1,0; 1,1,1]); %#ok

%%

setup = getDGTestCase('spe10_wo', 'useMomentFitting', true, 'I', 1:30, 'J', 1:110, 'k', {[0,0], [0,0; 1,0; 0,1; 1,1]}); %#ok

%%

setup = getDGTestCase('spe1', 'ijk', [5, 5, 1], 'useMomentFitting', false, 'degree', {0,1}); %#ok

%%

setup = getDGTestCase('spe9', 'ijk', [Inf, Inf, Inf], 'useMomentFitting', false, 'k', {[0,0,0], [0,0,0; 1,0,0; 0,1,0]}); %#ok

%%

setup = getDGTestCase('qfs_co2_2d', 'n', 3, 'useOverall', true); %#ok

%%

setup = getDGTestCase('viscous_fingers', 'rotate', true, 'k', {[0,0], [0,0; 1,0; 0,1], [0,0; 1,0; 0,1; 1,1]}, 'n', 20);

%%

sim = @(model,inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);
[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

[wsFI, stFI, repFI] = sim('modelFI', 1);

%%

ix = 2;
for i = ix
    [wsDG{i}, stDG{i}, repDG{i}] = sim('modelDG', i);
end

%%

six = 1:40;
schedule = setup.schedule;
schedule.step.val = schedule.step.val(six);
schedule.step.control = schedule.step.control(six);
[wsDG0, stDG0, repDG0] = simulateScheduleAD(setup.state0, setup.modelDG{1}, schedule);
[wsDG1, stDG1, repDG1] = simulateScheduleAD(setup.state0, setup.modelDG{2}, schedule);
[wsFV, stFV, repFV] = simulateScheduleAD(setup.state0, setup.modelFV{1}, schedule);

%%

inx = 2;
po = {'edgecolor', 'k'};
fn = plotLimiter(setup.modelDG{inx}.transportModel, 'plot1d', setup.plot1d, po{:}, 'n', 500, 'zlim', [-0.2, 1.2]);
setup.modelDG{inx}.transportModel.storeUnlimited = true;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.modelDG{inx}, setup.schedule, 'afterStepFn', fn);

%%

inx = 2;
po = {'edgecolor', 'k'};
fn = plotLimiter(setup.modelDG{inx}.transportModel, 'plot1d', false, po{:}, 'n', 100, 'zlim', [-0.2, 1.2]);
setup.modelDG{inx}.transportModel.storeUnlimited = true;
st0 = setup.state0;
st0.flux = stFV{2}.flux;
st0.wellSol = stFV{2}.wellSol;
[ws, st, rep] = simulateScheduleAD(st0, setup.modelDG{inx}.transportModel, setup.schedule, 'afterStepFn', fn);

%%

inx = 2;
po = {'edgecolor', 'none'};
fn = plotLimiter(setup.modelDG{inx}.transportModel, 'plot1d', true, po{:}, 'n', 100);
setup.modelDG{inx}.transportModel.storeUnlimited = true;
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.modelDG{inx}, setup.schedule, 'afterStepFn', fn);

%%

plotToolbar(setup.modelFV{1}.G, stDG{2});

%%

plotWellSols({wsDG0, wsDG1, wsFV}, schedule.step.val);

%%

coords = getPlotCoordinates(setup.modelFV{1}.G, 'n', 200);

%%

close all
hf = figure('Position', [0,0,2000,500]);

%%

set(gca,'nextplot','replacechildren');
for t = 1:numel(stDG{ix(1)})
    for i = 1:numel(ix)
        subplot(1, numel(ix), i)
        disc = setup.modelDG{ix(i)}.transportModel.discretization;
        h = plotSaturationDG(disc, stDG{ix(i)}{t}, 'edgecolor', 'none', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 1);
        axis tight
        zlim([0,1])
        pbaspect([1,1,0.25]);
        view([100,50]);
        camlight
        lighting gouraud
    end
    pause(0.2);
end


%%

ix = 2;
close all
figure('Position', [0,0,1000,500]);
azel = [64,-10];
view(azel);
disc = setup.modelDG{ix(1)}.transportModel.discretization;
h = plotSaturationDG(disc, stDG{ix(1)}{1}, 'edgecolor', 'k', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 3);
for t = 1:numel(stDG{ix(1)})
    for i = 1:numel(ix)
        delete(h);
%         h.Position = [-h.Position(1:2), h.Position(3)];
%         h.Position = -h.Position;
        subplot(1, numel(ix), i)
        h = plotSaturationDG(disc, stDG{ix(i)}{t}, 'edgecolor', 'k', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 3);
        ax = gca;
        plotWell(setup.modelFV{1}.G, setup.schedule.control(1).W);
%         plotGrid(setup.modelFV{1}.G, 'facec', 'none', 'facealpha', 0.2);
        ax.ZDir = 'reverse';
        box on
        view(azel);
        if t == 1
            hl = camlight;
        end
%         pbaspect([1,1,0.5]);
    end
    pause(0.1);
end

%%

vo = VideoWriter(fullfile(mrstPath('dg'), 'animations', 'spe10-subset.avi'));
vo.FrameRate = 5;
vo.open();

set(gca,'nextplot','replacechildren');
for t = 1:numel(stDG{ix(1)})
%     clf
    for i = 1:numel(ix)
        subplot(1, numel(ix), i)
        disc = setup.modelDG{ix(i)}.transportModel.discretization;
        h = plotSaturationDG(disc, stDG{ix(i)}{t}, 'edgecolor', 'none', 'edgealpha', 0.2, 'coords', coords, 'phaseNo', 1);
        axis tight
        zlim([0,1])
        pbaspect([1,2,0.25]);
        view([100,50]);
        camlight
        lighting gouraud
    end
    pause(0.2);
    M = getframe(hf);
    vo.writeVideo(M);
end
vo.close();

%% Test the triangle cubature
disc = setup.modelDG{1}.transportModel.discretization;
cubature = disc.cellCubature;
basis = disc.basis;
c = 5;
[W, x, w, cells] = cubature.getCubature(c, 'cell'); % Get points/weights
xr = cubature.transformCoords(x, cells); % Transform to reference coords
% The integral is easily evalauted using the matrix W, which has the
% cubature weights w placed so that int(psi)dx = W*psi(x). The weights sum
% to one, so that by construction, the integral of the first basis function
% (which is constant) should equal one, whereas the linear basis functions
% (1 and 2) should be zero
I = cellfun(@(psi) full(W*psi(xr)), basis.psi);
disp(I);

I = zeros(basis.nDof);
for i = 1:basis.nDof
    for j = 1:basis.nDof
        I(i,j) = W*(basis.psi{i}(xr).*basis.psi{j}(xr));
    end
end

%% Plot basis

% close all
disc = setup.modelDG{1}.transportModel.discretization;
plotDGBasis(disc.G, disc.basis, 5);

%%

plotCubature(disc.G, disc.faceCubature, 20, 'plotBoundingBox', false, 'facecolor', 'none')

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
