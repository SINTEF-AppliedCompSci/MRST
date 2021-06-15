function X = SPE9_2(eclout)
%Undocumented Helper Function for Examples

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

mrstVerbose true
mrstModule add deckformat
casenm  = 'SPE9B.DATA';

%

%eclout = '~/eclout/SPE9/more_output';
%eclout = 'c:\Documents and Settings\steink\My Documents\data\SPE9\more_output\';

%if compare&&isempty(who('rstrt'))
    [rstrt, rsspec] = readRestartLocal(fullfile(eclout, 'SPE9B'));
    [smry, smspec]  = readSummaryLocal(fullfile(eclout, 'SPE9B'));

%end

deck = readEclipseDeck(fullfile(eclout,casenm));

deck = resetSchedule(deck, smry);
deck = convertDeckUnits(deck);

%compare = ~isempty(dir(eclout));



%if compare
%end
init = readEclipseOutputFileUnFmt(fullfile(eclout, 'SPE9B.INIT'));
grid = readEclipseOutputFileUnFmt(fullfile(eclout, 'SPE9B.EGRID'));

[G, rock, N, T] = eclOut2mrst(init, grid);

s = setupSimComp(G, rock, 'neighbors', N, 'trans', T);

% f = initDeckADIFluid(deck);

gravity on
%gravity off

% init (hardcode conversion from field-units)
%ii = readEclipseOutputFileUnFmt('SPE9.X0000');
state.pressure = convertFrom(rstrt.PRESSURE{1}, psia);
sw = rstrt.SWAT{1};
sg = rstrt.SGAS{1};
state.s  = [sw, 1-sw-sg, sg];
state.rs = convertFrom(rstrt.RS{1}, (1000*ft^3)/stb);

%
schedule = deck.SCHEDULE;
%schedule = setControlToBHP(schedule);
%schedule = scheduleFromSummary(schedule, smry);


% schedule.control(1).WCONINJE{7} = 4500*psia;

% wellSols = runScheduleBlackOil3(state, G, rock, s, f, schedule);

% stepFunction = @(state0, state, meta, dt, W) stepBlackOil(state0, state, meta, dt, G, W, s, f, ...
%                'maxIts', 100, 'tol_relaxation', .1);
fluid = initDeckADIFluid(deck);

%%
system = initADISystem(deck, G, rock, fluid,'simComponents',s);
system.nonlinear.maxIterations = 30;
system.nonlinear.relaxRelTol = .1;
system.nonlinear.cpr = 1;
system.nonlinear.linesearch = false;
system.nonlinear.cprEllipticSolver = @(A,b) agmg(A,b, 5e-2);
%system.nonlinear.cprEllipticSolver = @mldivide;
system.nonlinear.cprRelTol = 5e-3;

% system.nonlinear.
system.nonlinear.bhpcontrols = false;
system.nonlinear.changeWells = false;
system.nonlinear.cprBlockInvert = false;
% system stepOptions
system.stepOptions.dpMax  = inf;
system.stepOptions.dsMax  = .02;
system.stepOptions.drsMax = inf;
system.stepOptions.solveWellEqs = true;

system.stepFunction = @(state0,state,meta,dt,W,G,system)stepBlackOil(state0,state,meta,dt,G,W,system,fluid);
[wellSols, states, iterations] = runScheduleADI(state, G, rock, system, schedule);

%system.stepFunction = @(state0,state,meta,dt,W,G,system)stepBlackOil(state0,state,meta,dt,G,W,system,fluid);
%[wellSols, states, iterations] = runScheduleADI(state, G, rock, system, schedule);


wwww = whos;
for v = { wwww.name },
   X.(v{1}) = eval(v{1});
end

return
%%
figure(1);
W = processWells(G, rock, schedule.control(1));
for i = 1:numel(states)
    clf;
%     subplot(2,1,2)
%     plotCellData(G, states{i}.s(:,3))
%     subplot(2,1,1)
    plotGrid(G, 'facec', 'none', 'edgea', 0)
%     plotGridVolumesBak(G, (states{i}.s(:,3)), 'N', 100, 'padnan', false, 'cmap', @autumn)
    plotGridVolumes(G, (states{i}.s(:,1)), 'N', 100, 'extrudefaces', true, 'cmap', @autumn)
    plotWell(G, W);
%     plotGridVolumes(G, (states{i}.s(:,1)), 'N', 100, 'padnan', false, 'cmap', @winter)
%     plotGridVolumes(G, (states{i}.pressure), 'N', 100, 'padnan', false, 'cmap', @winter)
%     plotGridVolumes(G, log10(rock.perm(:,3)), 'N', 100, 'padnan', false, 'cmap', @winter, 'basealpha', .05)
%     plotGridVolumes(G, log10(rock.poro), 'N', 100, 'padnan', false, 'cmap', @winter, 'basealpha', .05)
    view(35, 40)
    axis tight off
    drawnow
%     pause(.1)
end
%%
for i = 1:numel(states)
    clf;
    plotCellData(G, states{i}.pressure);
    if i == 1
        cbar = caxis();
    else
        caxis(cbar);
    end
    colorbar
    pause(.1)
end

%%

%% Set up plotting
inj = find([wellSols{1}.sign] == 1);
prod = find([wellSols{1}.sign] == -1);

% Put the well solution data into a format more suitable for plotting
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);

d = 2:29;

% Get timesteps for both the reference and the MRST run
T = convertTo(cumsum(schedule.step.val), year);
Tcomp =  smry.get(':+:+:+:+', 'YEARS', d);


%% Plot Producer Gas/Oil ratio
clf
ecl = convertFrom(smry.get('INJE1', 'WGOR', d), 1000*ft^3/stb)';
mrst = qGs(:,prod)./qOs(:,prod);

hold on
plot(T, mrst)
plot(Tcomp, ecl, 'r');
legend({'MRST', 'ECL'})

%% Plot Injector Bottom Hole Pressure
clf
ecl = convertFrom(smry.get('PROD10', 'WOPR', d), 1000*stb)';
mrst = qOs(:,10);
hold on
plot(T, mrst)
plot(Tcomp, ecl, 'r');
legend({'MRST', 'ECL'})

%% Plot Injector Bottom Hole Pressure
clf
ecl = convertFrom(smry.get('INJECTOR', 'WBHP', d), psia)';
mrst = bhp(:,inj);
hold on
plot(T, mrst)
plot(Tcomp, ecl, 'r');
legend({'MRST', 'ECL'})
%%

d = 2:29;

kw = {'WGPR', 'WBHP','WOPR'};
mr = {qGs, bhp, qOs};
un = [100*ft^3, psia, stb];
for i = 2
    figure(i)
    i
    for j = 1:numel(un)
        subplot(1,numel(un),j)
        ecl = convertFrom(smry.get(W(i).name, kw(j), d), un(j))';
        mrst = mr{j}(:,i);
        hold on;
        plot(T, mrst)
        plot(Tcomp, ecl, 'r');
        legend({'MRST', 'ECL'})
        title(kw(j))
    end
end
end
