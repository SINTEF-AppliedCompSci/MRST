% 2D Validation Against ECLIPSE
%
% We consider a vertial 4000 m-by-200 m-by-125 m reservoir cross-section
% described on a 20-by-1-by-5 uniform Cartesian grid. The reservoir is
% produced by a single injector-producer pair, located to the lower right
% and the upper right, respectively, in the formation. The injection well
% is under rate control with target rate 1000 m3/day and upper limit of 800
% bar on the bottom-hole pressure (bhp), whereas the production well is
% under pressure control with target bottom-home pressure 260 bar.
%
% We consider a standard EOR setup consisting of a water preflush period of
% 1260 days, followed by a chemical slug injected over 1700 days, which is
% washed out by 8000 days of chase water injection. We consider three
% alternative chemical slugs: polymer, surfactant, and a combination of the
% two. For reference, we also simulate a pure waterflooding scenario. All
% simulations are performed with the full three-phase, five-component SP
% simulator class.
%
% We end the example by comparing bottom-hole pressures and oil rates
% predicted by MRST with similar results computed by ECLIPSE, that have
% been precomputed and stored on file.
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui coarsegrid

mrstVerbose off

%% Surfactant-polymer flooding
% We load the data file and use MRST's feature for automatically selecting
% the correct model type and suitable solver defaults.
gravity reset on;
bookdir  = getDatasetPath('eor_book_ii');
fn       = fullfile(bookdir,'validation', '2D', 'SURFACTANTPOLYMER2D.DATA');
[state0, model, schedule, nlsolver] = initEclipseProblemAD(fn);
[wellSolsSP, statesSP, reportSP] =  ...
    simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);

%% Waterflooding
% To simplify setup, we use the full three-phase, five-component SP
% simulator class to simulate waterflooding.
scheduleW = schedule;
for i=1:3
    scheduleW.control(i).W(1).cs = 0;
    scheduleW.control(i).W(1).cp = 0;
end
[wellSolsW, statesW, reportW] = simulateScheduleAD(state0, model, scheduleW);

%% Polymer flooding
scheduleP = schedule;
for i=1:3
    scheduleP.control(i).W(1).cs = 0;
end
[wellSolsP, statesP, reportP] = simulateScheduleAD(state0, model, scheduleP);

%% Surfactant flooding
scheduleS = schedule;
for i=1:3
    scheduleS.control(i).W(1).cp = 0;
end
[wellSolsS, statesS, reportS] = simulateScheduleAD(state0, model, scheduleS);

%% Plot well solutions
plotWellSols({wellSolsW, wellSolsP, wellSolsS, wellSolsSP}, ...
    cumsum(schedule.step.val), 'datasetnames',{'Water flooding', ...
    'Polymer flooding', 'Surfactant flooding', 'Surfactant-polymer flooding'})


% -------------------------------------------------------------------------
% POST PROCESSING AND VISUALIZATON OF RESULTS
%
% The following code is quite verbose and contains a lot of nitty gritty to
% tweak plots of setup, comparison of saturation/concentration profiles at
% different times, and comparing bottom-hole pressures and oil rates
% predicted by MRST and ECLIPSE
% -------------------------------------------------------------------------

%% Permeability and intial saturation
G    = model.G;
rock = model.rock;
W    = schedule.control(1).W;
figure('Position',[100 500 1000 190])

ax1 = subplot(1,2,1);
plotCellData(G, rock.perm(:,1)/(milli*darcy),'EdgeAlpha',.1); 
colorbar(gca,'EastOutside','FontSize',10);

ax2 = subplot(1,2,2);
plotCellData(G,state0.s(:,[3 2 1]),'EdgeAlpha',.1);

% Plot wells, fine-tune axes, and set colormap
wcoord = vertcat(G.cells.centroids(vertcat(W.cells),:));
for i=1:2
    subplot(1,2,i),
    view(0,0), axis tight
    set(gca,'FontSize',10);
    hold on
    plot3(wcoord(1:2,1),-[1 1],wcoord(1:2,3)+[-12;12],'-','LineWidth',3,'Color',[.4 .4 .5]);
    plot3(wcoord(3:4,1),-[1 1],wcoord(3:4,3)+[-12;12],'-','LineWidth',3,'Color',[.4 .4 .5]);
    text(wcoord(1,1),-1, wcoord(1,3)-15,'I','FontSize',14,'Color', [.4 .4 .5], ...
        'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Century Schoolbook')
    text(wcoord(4,1),-1, wcoord(4,3)+15,'P','FontSize',14,'Color', [.4 .4 .5], ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontName','Century Schoolbook')
    hold off
end
colormap(ax1,.75*flipud(pink)+.25)
drawnow;
ax2.Position(3) = ax1.Position(3);

% Add ternary triangle
axes('Position',[.85 .14 .06 .25]);
patch('Vertices', [0 0; 2 0; 1 2*sin(pi/3)], 'Faces',1:3, ...
    'FaceVertexCData', [0 0 1; 0 1 0; 1 0 0],'FaceColor','interp','EdgeColor','none');
text(-0.05,0,'S_w','HorizontalAlignment','right'); 
text(2.05,0,'S_o','HorizontalAlignment','left');
text(1,2*sin(pi/3)+.25,'S_g','HorizontalAlignment','center');
set(gca,'FontSize',10); axis tight off

%% Evolution of the three-phase water saturation and chemical concentrations
% We plot this as a 5x4 matrix of color plots in which the three-phase
% saturation is shown in RGB and outlines of surfactant and/or polymer are
% shown as lines. For this, we use two different axes overlayed on each
% other, one for the RGB plot and one for the contour plots.
figure('Position',[480 420 840 340])
stepNo = [40 80 148 188 228];
nrow = length(stepNo);
minC = min(G.cells.centroids);
maxC = max(G.cells.centroids);
xl = [minC(1) maxC(1)]; xc = unique(G.cells.centroids(:,1));
yl = [minC(2) maxC(2)]; yc = unique(G.cells.centroids(:,2));
zl = [minC(3) maxC(3)]; zc = unique(G.cells.centroids(:,3));
mapCS = min(bsxfun(@times, flipud(gray(128).^.25), [0.494 0.184 0.556])+.4,1);
mapCP = min(bsxfun(@times, flipud(gray(128).^.25), [0.929 0.694 0.125])+.4,1);
names = {'W flooding','P flooding','S flooding','SP flooding'};
for i=1:nrow
    spno = i;
    for j=1:4
        switch j
            case 1, state = statesW{stepNo(i)};
            case 2, state = statesP{stepNo(i)};
            case 3, state = statesS{stepNo(i)};
            case 4, state = statesSP{stepNo(i)};
        end
        
        % Plot 3-phase saturation (and add title if upper row)
        ax1 = subplot(4,nrow,spno); spno = spno + nrow;
        ax1.Position = (ax1.Position + 2*ax1.OuterPosition)/3-[.05 0 0 0];
        s   = permute(squeeze(reshape(state.s, [G.cartDims 3])),[2 1 3]);       
        image(ax1,xl,zl,s(:,:,[3 2 1])+.1); set(ax1,'XTick',[],'YTick',[]);
        if j==1
            title(['Time: ' num2str(reportW.ReservoirTime(stepNo(i))/year,3), ' yrs'], ...
                'FontWeight','normal','FontSize',10);
        end
        if i==1, ylabel(names{j},'FontSize',10), end
        
        % Add outline of polymer
        axlim = axis(ax1);
        bx1 = axes('Position',ax1.Position);
        if any(state.cp>.05)
            cg = generateCoarseGrid(G,(state.cp>.05)+1);
            plotFaces(cg,find(~any(cg.faces.neighbors==1,2)),...
                'EdgeColor',[1 0 0.8],'FaceColor','none','LineWidth',1);
        end
 
        % Add outline of surfactant
        if any(state.cs>.05)
            cg = generateCoarseGrid(G,(state.cs>.05)+1);
            plotFaces(cg, find(~any(cg.faces.neighbors==1,2)),...
                'EdgeColor',[1 0.8 0],'FaceColor','none','LineWidth',1);
        end
        set(bx1,'XLim',axlim(1:2),'ZLim',axlim(3:4),'ZDir','reverse'); view(0,0)
        axis(bx1,'off');
        
        % Add wells
        hold(bx1,'on');
        plot3(wcoord(1:2,1),-[1 1],wcoord(1:2,3)+[-12;12],'-w','LineWidth',1);
        plot3(wcoord(3:4,1),-[1 1],wcoord(3:4,3)+[-12;12],'-w','LineWidth',1);
        hold(bx1,'off')
    end
end
% Add ternary triangle
axes('Position',[.89 .45 .06 .15]);
patch('Vertices', [0 0; 2 0; 1 2*sin(pi/3)], 'Faces',1:3, ...
    'FaceVertexCData', [0 0 1; 0 1 0; 1 0 0]+.1,'FaceColor','interp','EdgeColor','none');
text(-0.05,0,'S_w','HorizontalAlignment','right'); 
text(2.05,0,'S_o','HorizontalAlignment','left');
text(1,2*sin(pi/3)+.25,'S_g','HorizontalAlignment','center');
axis tight off

% Add legend for polymer and surfactant
axes('Position',[.88 .325 .06 .1]);
patch([0 1 1 0],[0 0 1 1],[1 1 1 1],'FaceColor',[1 .8 0],'EdgeColor','none');
text(1.2,.5,'Surfactant','FontSize',10);
axis([0 3 -3 4]); axis off

axes('Position',[.88 .27 .06 .1]);
patch([0 1 1 0],[0 0 1 1],[1 1 1 1],'FaceColor',[1 0 .8],'EdgeColor','none');
text(1.2,.5,'Polymer','FontSize',10);
axis([0 3 -3 4]); axis off

%% Setup well results from MRST and ECLIPSE
% To make the plotting setup as compact as possible, we store well data
% from MRST and ECLIPSE in cell arrays.
scenarioNames = {'Water','Polymer','Surfactant','Surfactant+polymer'};

% Compact setup of data from MRST       
time       = cumsum(schedule.step.val)/year;
reportMRST = {reportW, reportP, reportS, reportSP};
wsMRST     = {wellSolsW, wellSolsP, wellSolsS, wellSolsSP};
timeMRST   = cellfun(@(x)x.ReservoirTime /year, reportMRST, 'UniformOutput', false);

% Load precomputed results from ECLIPSE
load(fullfile(bookdir, 'validation', '2D', 'ws2D_ecl.mat'));
% get parameters of injection well
wellNum = 1;
field   = 'bhp';
unit    = barsa;
bhpMRST = cellfun(@(y) cellfun(@(x)x(wellNum).(field)/ unit, y), ...
               wsMRST, 'UniformOutput', false);
bhpECL  = cellfun(@(y) cellfun(@(x)x(wellNum).(field)/ unit, y),  ...
               wsECL, 'UniformOutput', false);
% get parameters of production well
wellNum = 2;
field   = 'qOs';
unit    = -(meter^3/day);
qOsMRST  = cellfun(@(y) cellfun(@(x)x(wellNum).(field)/ unit, y),  ...
               wsMRST, 'UniformOutput', false);
qOsECL   = cellfun(@(y) cellfun(@(x)x(wellNum).(field)/ unit, y),  ...
               wsECL, 'UniformOutput', false);

%% Plot the injector bhp
figure, hold on
cols = {'[0,0.5,1]', '[0.7,0,0.6]', '[1,0.8,0]', '[1,0.25,0]'};
linewidth_MRST = [1.5;1.5;1.5;1.5];
linewidth_Ecl = [4;4;4;4];
database = [1,2,3,4];
patch(time([80 80 148 148]), [253 847 847 253],[1 1 1],...
    'FaceColor',[.95 .95 .95],'EdgeColor','none')
arrayfun(@(x)plot(timeMRST{x}, bhpMRST{x}, '-', 'color', cols{x}, ...
    'LineWidth', linewidth_MRST(x)), database);
arrayfun(@(x)plot(timeECL{x}, bhpECL{x}, '--', 'color', cols{x}, ...
    'LineWidth', linewidth_Ecl(x)), database);
plot([0 time(end)], [800 800], '--','Color',[.5 .5 .5]);
set(gca,'FontSize',12)
xlabel('Time [years]');
ylabel('Bottom-hole pressure [bar]');
title('Bottom-hole pressure in injector');
axis([0,timeMRST{1,1}(end),250,850]);
h=get(gca,'Children');
legend(h(end-1:-1:end-4),scenarioNames{:},'Location','SouthEast');
box on

%% Plot the oil-production rate
figure, hold on
patch(time([80 80 148 148]), [5 795 795 5],[1 1 1],...
    'FaceColor',[.95 .95 .95],'EdgeColor','none')
arrayfun(@(x)plot(timeMRST{x}, qOsMRST{x}, '-', 'color', cols{x}, ...
    'LineWidth', linewidth_MRST(x)), database);
arrayfun(@(x)plot(timeECL{x}, qOsECL{x}, '--', 'color', cols{x}, ...
    'LineWidth', linewidth_Ecl(x)), database);
set(gca,'FontSize',12)
xlabel('Time [years]');
ylabel('Oil rate [m^3/day]');
title('Oil rate in producer');
axis([0,timeMRST{1,1}(end),0,800]);
h=get(gca,'Children');
legend(h(end-1:-1:end-4),scenarioNames{:},'Location','NorthEast')
box on


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
