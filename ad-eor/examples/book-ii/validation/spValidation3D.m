%% Validation Against ECLIPSE on a 3D Sector Model
% This example is an extension of a polymer-flooding case reported
% previously by Bao et al. (Comput. Geosci. doi:10.1007/s10596-017-9624-5).
% Here, we consider four different strategies: waterflooding, polymer
% flooding, surfactant flooding, and surfactant-polymer flooding. These
% scenarios are simulated with the full three-phase, five-component model
% in MRST and compared with corresponding results from ECLIPSE.

mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui coarsegrid
mrstVerbose off

%% Surfactant-polymer flooding
% We load the data file and use MRST's feature for automatically selecting
% the correct model type and suitable solver defaults.
gravity reset on;
bookdir  = getDatasetPath('eor_book_ii');
fn       = fullfile(bookdir,'validation', '3D', '3D_SP.DATA');
[state0, model, schedule, nlsolver] = ...
    initEclipseProblemAD(fn, 'useCPR', false, ...
                         'AutoDiffBackend', 'diagonal-mex');
nlsolver.useRelaxation = 1;
[wellSolsSP, statesSP, reportsSP] =  ...
    simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);

%% Waterflooding
% To simplify setup, we use the full three-phase, five-component SP
% simulator class to simulate waterflooding.
scheduleW = schedule;
for i=1:3
    scheduleW.control(i).W(1).cs = 0;
    scheduleW.control(i).W(1).cp = 0;
end
[wellSolsW, statesW, reportsW] = ...
    simulateScheduleAD(state0, model, scheduleW, 'NonLinearSolver', nlsolver);

%% Polymer flooding
scheduleP = schedule;
for i=1:3
    scheduleP.control(i).W(1).cs = 0;
end
[wellSolsP, statesP, reportsP] = ...
    simulateScheduleAD(state0, model, scheduleP, 'NonLinearSolver', nlsolver);

%% Surfactant flooding
scheduleS = schedule;
for i=1:3
    scheduleS.control(i).W(1).cp = 0;
end
[wellSolsS, statesS, reportsS] = ...
    simulateScheduleAD(state0, model, scheduleS, 'NonLinearSolver', nlsolver);

%% Plot well solutions
plotWellSols({wellSolsW, wellSolsP, wellSolsS, wellSolsSP}, ...
    cumsum(schedule.step.val), 'datasetnames',{'Water', ...
    'Polymer', 'Surfactant', 'Surfactant+polymer'})

% -------------------------------------------------------------------------
% POST-PROCESSING AND VISUALIZATON OF RESULTS
%
% The following code is quite verbose and contains a lot of nitty gritty to
% tweak plots of setup, comparison of saturation/concentration profiles at
% different times, and comparing bottom-hole pressures and oil rates
% predicted by MRST and ECLIPSE
% -------------------------------------------------------------------------

%% Show initial setup
G     = model.G;
rock  = model.rock;
W     = schedule.control.W;

figure('Position',[350 450 1080 290])

% Horizontal permeability
ax(1) = subplot(1,2,1);
K = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(G, K, 'EdgeAlpha',.1); 

% Initial saturation
ax(2) = subplot(1,2,2);
plotCellData(G,state0.s(:,[3 2 1]),'EdgeAlpha',.1);

% Plot wells 
for i=1:2
    axes(ax(i)); %#ok<LAXES>
    camlight, material dull
    [~,ht] = plotWell(G,W,'FontSize',8,'height',50,'Color','b','Color2','r');
    set(ht,'BackgroundColor',[.9 .9 .9],'EdgeColor','k','HorizontalAlignment','right');
    view(135,23), set(gca,'DataAspectRatio',[1.5 3 1]);
    axis tight off
end

% Add colorbar with histogram to the first plot
colormap(ax(1),.9*flipud(pink)+.1)
[hc,hh]=mrstColorbar(ax(1), K, 'west');
hh.Children.NumBins   = 20;
hh.Children.FaceColor = [.92 .92 .74];
hc.TickLength = .05;
for i=1:2
    ax(i).Position(3) = .5;
    ax(i).Position(1) = ax(i).Position(1)-.055*(i+1);
end

% Add ternary triangle to the second plot
axes('Position',[.73 .19 .05 .15]);
patch('Vertices', [0 0; 2 0; 1 2*sin(pi/3)], 'Faces',1:3, ...
    'FaceVertexCData', [0 0 1; 0 1 0; 1 0 0],'FaceColor','interp','EdgeColor','none');
text(-0.05,0,'S_w','HorizontalAlignment','right'); 
text(2.05,0,'S_o','HorizontalAlignment','left');
text(1,2*sin(pi/3)+.25,'S_g','HorizontalAlignment','center');
set(gca,'FontSize',10); axis tight off

%% Setup well results from MRST and ECLIPSE
% To make the plotting setup as compact as possible, we store well data
% from MRST and ECLIPSE in cell arrays.
scenarioNames = {'Water','Polymer','Surfactant','Surfactant+polymer'};

% Compact setup of data from MRST       
time       = cumsum(schedule.step.val)/day;
reportMRST = {reportsW, reportsP, reportsS, reportsSP};
wsMRST     = {wellSolsW, wellSolsP, wellSolsS, wellSolsSP};
timeMRST   = cellfun(@(x)x.ReservoirTime/day, reportMRST, 'UniformOutput', false);

% Load precomputed results from ECLIPSE
data    = load(fullfile(bookdir, 'validation', '3D', 'results_ECL.mat'));
timeECL = data.results_ECL.timeECL;
wsECL = cellfun(@(x) data.results_ECL.wsECL.(x), ...
           {'WaterFlooding', 'PolymerFlooding', 'SurfactantFlooding', 'SurfactantPolymerFlooding'}, 'UniformOutput', false);
% wsECL   = data.results_ECL.wsECL;
clear data;

%% Plot the injector bhp
wellNum = 1;
field   = 'bhp';
unit    =  barsa;
bhpMRST = cellfun(@(y) cellfun(@(x)x(wellNum).(field)/ unit, y), ...
              wsMRST, 'UniformOutput', false);
bhpECL  = cellfun(@(y) cellfun(@(x)x(wellNum).(field)/ unit, y),  ...
              wsECL, 'UniformOutput', false);

figure, hold on
cols = {'[0,0.5,1]', '[0.7,0,0.6]', '[1,0.8,0]', '[1,0.25,0]'};
linewidthMRST = [1.5;1.5;1.5;1.5];
linewidthECL = [4;4;4;4];
database = [1,2,3,4];
patch(time([47 47 88 88]), [250.1 294.9 294.9 250.1],[1 1 1],...
    'FaceColor',[.95 .95 .95],'EdgeColor','none')
arrayfun(@(x)plot(timeMRST{x}, bhpMRST{x}, '-', 'color', cols{x}, ...
    'LineWidth', linewidthMRST(x)), database);
arrayfun(@(x)plot(timeECL, bhpECL{x}, '--', 'color', cols{x}, ...
    'LineWidth', linewidthECL(x)), database);
plot([0 timeMRST{1,1}(end)], [290 290],'--','Color',[.5 .5 .5]);
set(gca,'FontSize',12)
xlabel('Time [days]');
ylabel('Bottom-hole pressure [bar]');
axis([0 timeMRST{1,1}(end) 250 295]);
h=get(gca,'Children');
legend(h(end-1:-1:end-4),scenarioNames{:}, ...
    'Location','SouthOutside','Orientation','Horizontal');
box on
ax = gca; ax.Position([2 4]) = ax.Position([2 4]) - [.07 -.07];

%% Plot field oil production
wn1 = 2;
wn2 = 3;
unit    = -(meter^3/day);
figure, hold on
[coM,coE,h] = deal(zeros(4,1));
for i=1:4
    cumOilMRST = cumtrapz(timeMRST{i}, ...
        cellfun(@(x) x(wn1).qOs/unit + x(wn2).qOs/unit, wsMRST{i}));
    h(i)=plot(timeMRST{i},cumOilMRST,'-','Color',cols{i}, ...
        'LineWidth',linewidthMRST(i));
    coM(i)=cumOilMRST(end);
    if i>1
        text(timeMRST{i}(end), coM(i), ...
            [num2str((coM(i)-coM(1))/coM(1)*100,'+%.1f') '%'],'Color',[.4 .4 .4]);
    end
    
    cumOilECL = cumtrapz(timeECL, ...
        cellfun(@(x) x(wn1).qOs/unit + x(wn2).qOs/unit, wsECL{i}));
    plot(timeECL,cumOilECL,'--','Color',cols{i}, ...
        'LineWidth',linewidthECL(i));
    coE(i) = cumOilECL(end);
    text(.7*timeMRST{i}(end), (.5-i*.08)*coM(1), [scenarioNames{i} ': ' ...
        num2str((coM(i)-coE(i))/coM(i)*100,'%.1f') '%'],'Color',[.4 .4 .4]);
end
text(.65*timeMRST{1}(end), .5*coM(1),'Simulator discrepancy',...
    'FontWeight','bold','Color',[.4 .4 .4]);
plot([.63 .63 1.15 1.15 .63]*timeMRST{1}(end),...
    [.13 .54 .54 .13 .13]*coM(1),'Color',[.4 .4 .4]);
set(gca,'FontSize',12)
xlabel('Time [days]');
ylabel('Oil production [m^3]');
set(gca,'XLim',[0 1.2*timeMRST{1,1}(end)]);
legend(h,scenarioNames{:}, ...
    'Location','SouthOutside','Orientation','Horizontal');
box on
ax = gca; ax.Position([2 4]) = ax.Position([2 4]) - [.07 -.07];

%% Plot water rate
wn = 2;
unit = -meter^3/day;
figure, hold on
patch(time([47 47 88 88]), [0 1510 1510 0],[1 1 1],...
    'FaceColor',[.95 .95 .95],'EdgeColor','none')
for i=1:4
    qWsM = cellfun(@(x) x(wn).qWs/unit, wsMRST{i});
    qWsE = cellfun(@(x) x(wn).qWs/unit, wsMRST{i});
    h(i)=plot(timeMRST{i},qWsM,'-','Color',cols{i}, 'LineWidth',linewidthMRST(i));
    plot(timeECL,qWsE,'--','Color',cols{i},'LineWidth',linewidthECL(i));
end
set(gca,'FontSize',12)
xlabel('Time [days]');
ylabel('Water production rate [m^3/day]');
set(gca,'XLim',[0 timeMRST{1,1}(end)]);
legend(h,scenarioNames{:}, ...
    'Location','SouthOutside','Orientation','Horizontal');
box on
ax = gca; ax.Position([2 4]) = ax.Position([2 4]) - [.07 -.07];

%%
figure('Position',[775 20 620 700]);
state(1) = statesP{131};
state(2) = statesS{131};
for i=1:2
    ax(i)=subplot(2,1,i);
    plotCellData(G, state(i).s(:,[3 2 1]),state(i).s(:,1)>0.55, 'EdgeAlpha',.1);
    camlight, material dull
    [~,ht] = plotWell(G,W,'FontSize',10,'height',50,'Color','b','Color2','r');
    set(ht,'BackgroundColor',[.9 .9 .9],'EdgeColor','k','HorizontalAlignment','right');
    view(156,47), set(gca,'DataAspectRatio',[1.5 3 1]);
    axis tight off
    cg = generateCoarseGrid(G,ones(G.cells.num,1));
    plotFaces(cg,1:cg.faces.num,'EdgeColor',[.5 .5 .5],'FaceColor','none','LineWidth',1.5);
    camlight headlight
end
for i=1:2
    ax(i).Position = [-.5 -.325+i*.38 2 .55];
end

% Add ternary triangle to the second plot
axes('Position',[.78 .15 .1 .08]);
patch('Vertices', [0 0; 2 0; 1 2*sin(pi/3)], 'Faces',1:3, ...
    'FaceVertexCData', [0 0 1; 0 1 0; 1 0 0],'FaceColor','interp','EdgeColor','none');
text(-0.05,0,'S_w','HorizontalAlignment','right'); 
text(2.05,0,'S_o','HorizontalAlignment','left');
text(1,2*sin(pi/3)+.25,'S_g','HorizontalAlignment','center');
set(gca,'FontSize',10); axis tight off

%% Copyright Notice
%
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
