%% Conformance improvement for areal SPE10 subset
% This example considers an areal subset of the SPE10 bechmark and runs two
% simulations, one with waterflooding and one with polymer flooding. The
% two solutions are then compared to show the improved conformance and how
% polymer distributes in the reservoir.
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat

%% Set up and run polymer flooding
gravity reset on;
bookdir = getDatasetPath('EOR_Book_II');
fn      = fullfile(bookdir,'conformanceSPE10','SPE10_MODEL2.DATA');
[state0, model, schedule, nlsolver] = initEclipseProblemAD(fn);
model.rock.poro = max(model.rock.poro,1e-4);
[wellSolsSP, statesSP, reportsSP] =  ...
    simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);

%% Set up and run waterflooding
scheduleW = schedule;
scheduleW.control.W(1).cp = 0;
[wellSolsW, statesW, reportsW] = ...
    simulateScheduleAD(state0, model, scheduleW, 'NonLinearSolver', nlsolver);

%% Compare distribution of displacing fluid(s) at different times
figure('Position',[670 380 860 400]);
cells   = vertcat(scheduleW.control.W.cells);
[I,J,K] = gridLogicalIndices(model.G);
stepNo  = [8 16 24 30];
kx      = reshape(log10(model.rock.perm(:,1)),model.G.cartDims)';
for i=1:length(stepNo)
    bx1 = subplot(2,4,i);
    bx1.Position = bx1.Position+[0 -.03 .03 .05];
    s  = reshape(statesW{stepNo(i)}.s(:,1), model.G.cartDims)';
    c  = reshape(statesW{stepNo(i)}.cp, model.G.cartDims)';
    
    image(kx, 'CDataMapping','scaled'); set(gca,'YDir','normal'), axis equal off
    colormap(bx1,.85*flipud(pink)+.15);
    title(['Time: ' num2str(reportsW.ReservoirTime(stepNo(i))/day), ' days'], ...
        'FontWeight','normal');
    %
    axlim = axis(bx1);
    ax1 = axes('Position',bx1.Position);
    hs = image(ax1, s, 'CDataMapping','scaled'); set(gca,'YDir','normal')
    hs.AlphaData = s>.35;
    colormap(ax1,flipud(winter));
    axis(ax1,axlim); axis off
    
    %
    cx1 = axes('Position',bx1.Position);
    hc = image(cx1, c, 'CDataMapping','scaled'); set(gca,'YDir','normal')
    hc.AlphaData = .5*c;
    colormap(cx1,flipud(autumn));
    axis(cx1,axlim); axis off
    hold(cx1,'on');
    plot(cx1,I(cells)-.5,J(cells)-.5,'.k','MarkerSize',18)

    %
    bx2 = subplot(2,4,i+4); bx2.Position = bx2.Position+[0 0 .03 .05];
    s  = reshape(statesSP{stepNo(i)}.s(:,1), model.G.cartDims)';
    c  = reshape(statesSP{stepNo(i)}.cp, model.G.cartDims)';
    
    image(kx, 'CDataMapping','scaled'); set(gca,'YDir','normal'), axis equal off
    colormap(bx2,.85*flipud(pink)+.15);
    
    %
    axlim = axis(bx2);
    ax2 = axes('Position',bx2.Position);
    hs = image(ax2, s, 'CDataMapping','scaled'); set(gca,'YDir','normal')
    hs.AlphaData = s>.35;
    colormap(ax2,flipud(winter));
    axis(ax2,axlim); axis off
    
    %
    cx2 = axes('Position',bx2.Position);
    hc = image(cx2, c, 'CDataMapping','scaled'); set(gca,'YDir','normal')
    hc.AlphaData = .5*c;
    colormap(cx2,flipud(autumn));
    axis(cx2,axlim); axis off
    hold(cx2,'on');
    plot(cx2,I(cells)-.5,J(cells)-.5,'.k','MarkerSize',18)
end

%% Plot well solutions
plotWellSols({wellSolsW, wellSolsSP}, ...
    reportsW.ReservoirTime, 'datasetnames',{'Water','Polymer'});

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
