function setup = ltates_geothermal(varargin)
%Test case describing low-enthalpy aquifer thermal energy storage

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    % Step 1: Test case description and options
    %---------------------------------------------------------------------%
    description = 'Low-enthalpy aquifer thermal energy storage';
    options = struct();
    % Process optinal input arguments
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
                                        options, description, varargin{:});
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%

    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    gravity reset on

    out = importSVG('~/Desktop/outline.svg', 'scale', 1, 'n', 50);

    n = 20;
    arg = out.unpack();
    
    rng(20210821);
    xw = {[0.52, 0.36], [0.15, 0.26], [0.88, 0.1]};
    G2D = pebiGrid2D(1/n, [1,1], arg{:}, 'cellConstraints', xw, 'CCRefinement', true, 'CCEps', 0.1, 'CCFactor', 0.1);
    
    
    [y,x,z]  = peaks(30);
    xmax = max(max(x));
    xmin = min(min(x));
    x = (x - xmin)./(xmax - xmin);
    ymax = max(max(y));
    ymin = min(min(y));
    y = (y - ymin)./(ymax - ymin);
%     z = z - 500;
    z = z./100;
    zmax = max(max(z));
    zmin = min(min(z));
    horizons = {struct('x',x,'y',y,'z',z.*0 + zmin - 0.1), ...
                struct('x',x,'y',y,'z',z), ...
                struct('x',x,'y',2*y+1,'z',z+0.05-0.075*x), ...
                struct('x',2*x,'y',y,'z',z+0.1), ...
                struct('x',x,'y',y,'z',z.*0 + zmax + 0.1 + 0.1)};
    mrstModule add upr
%     G2D = pebiGrid2D(6/25, [6,6], 'polyBdr', [-3,-3; 3,-3; 3,3; -3,3]);
    layers = [3,16,16,3];
    GL = computeGeometry(makeLayeredGrid(G2D, sum(layers)));
    
    
    rock = getSPE10rock(1:sum(layers));
    rock.poro = max(rock.poro, 0.01);
%     [~, model] = setupSPE10_AD('layers', 1:20);
    
    perm = reshape(rock.perm(:,1), [60, 220, sum(layers)]);
    perm = permute(perm, [2,1,3]);
    perm = sampleFromBox(GL, perm);
    
    poro = reshape(rock.poro, [60, 220, sum(layers)]);
    poro = permute(poro, [2,1,3]);
    poro = sampleFromBox(GL, poro);

    
    [G, cellmap]   = makeLayeredHorizonGrid(G2D, horizons,'layers', layers, 'repairFunction2', @min);    
    
    G.nodes.coords = G.nodes.coords*850;
    G.nodes.coords(:,3) = G.nodes.coords(:,3) + 800;
    
    G = mcomputeGeometry(G);
    
    perm = perm(cellmap);
    poro = poro(cellmap);
    
    perm(G.cells.layer == 1) = 0.1*milli*darcy;
    poro(G.cells.layer == 1) = 0.01;
    
    nl = max(G.cells.layer);
    perm(G.cells.layer == nl) = 0.1*milli*darcy;
    poro(G.cells.layer == nl) = 0.01;
    
    perm(G.cells.layer == 2) = perm(G.cells.layer == 2)/10;
    perm(G.cells.layer == 3) = perm(G.cells.layer == 3)*10;
    
    
    rock = makeRock(G, perm, poro);
    rock = addThermalRockProps(rock);
    
    fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'rho', 1, 'mu', 1);
    fluid = addThermalFluidProps(fluid, 'useEOS', true);
    
    model = GeothermalModel(G, rock, fluid);
    
    
    wno = zeros(G2D.cells.num); wno(G2D.cells.tag) = 1:nnz(G2D.cells.tag);
    wno = repmat(wno, sum(layers), 1);
    wno = wno(cellmap);
    
    aquifer = G.cells.layer == 2 | G.cells.layer == 3;
    
    time = 20*year;
    rate = 100*litre/second;
    
    W = [];
    cells = wno == 1 & aquifer;
    zmin = min(G.cells.centroids(cells,3));
    p = 1000*norm(gravity)*zmin;
    W = addWell(W, G, rock, find(cells), 'Name', 'C1', 'type', 'bhp', 'val', p);

    cells = wno == 3 & aquifer;
    zmin = min(G.cells.centroids(cells,3));
    p = 1000*norm(gravity)*zmin;
    W = addWell(W, G, rock, find(cells), 'Name', 'C2', 'type', 'bhp', 'val', p);
    
    cells = wno == 2 & aquifer;
    zmin = min(G.cells.centroids(cells,3));
    p = 1000*norm(gravity)*zmin;
    W = addWell(W, G, rock, find(cells), 'Name', 'H', 'type', 'rate', 'val', rate);
    W(end).lims = struct('bhp', 2*p);

    K0 = 273.15*Kelvin;
    W = addThermalWellProps(W, 'T', [10, 10, 80]*Kelvin + K0);
    
    
    schedule = readSchedule(fullfile(mrstPath('geothermal'), 'example-suite', 'data', 'temperature.xlsx'));
    schedule.control(1).W = W;
    schedule.control(1).type = 'dummy';
    schedule = addLoadPhase(schedule, rate);
    schedule = addRestPhase(schedule);
    schedule = addUnloadPhase(schedule, rate);
    
    schedule.control = arrayfun(@(ctrl) rmfield(ctrl, 'type'), schedule.control);
    
    x = model.G.cells.centroids;
    
    p = @(x) norm(gravity)*1000*x(:,3);
    T = @(x) (K0 + 10)*Kelvin + 30*Kelvin/(kilo*meter)*x(:,3);
    state0 = initResSol(G, p(x), 1);
    state0.T = T(x);
    
    x = model.G.cells.centroids;
    wc = find(cells);
    d = pdist2(x(:,1:2), x(wc(1),1:2));
    monitorCells = d > 150 & d < 160 & aquifer;
    
    
    ctrlFn = @(state, schedule, report, i) ...
            controlLogicFn(state, schedule, report, i, monitorCells);
    
    options.ctrlFn = ctrlFn;
    options.monitorCells = monitorCells;
        
    plotOptions = {'View', [180, 25], ...
                   'Size', [800, 600], ...
                   'PlotBoxAspectRatio', [1.9065    1.2256    1.0000], ...
                   'Box', true};
    
    %---------------------------------------------------------------------%
    
    % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    setup = packTestCaseSetup(mfilename,          ...
                      'description', description, ...
                      'options'    , options    , ...
                      'state0'     , state0     , ...
                      'model'      , model      , ...
                      'schedule'   , schedule   , ...
                      'plotOptions', plotOptions);
    %---------------------------------------------------------------------%
    
end

%-------------------------------------------------------------------------%
function schedule = readSchedule(fn, opt)
    % Read input sheet
    [v, n, a] = xlsread(fn , 'Sheet');
    dt = day * diff(datenum(a(2:end-1, 3), 'dd.mm.yyyy')); 
    %dt = seconds(diff(datetime(a(2:end-1,3)))); % replaced with line over
                                                 % for octave compatibiltiy
    thres = 10;
    temp    = v(:,1);
    % Three controls: 1 = load, 2 = rest, 3 = unload
    % Default control is rest
    ctrl = repmat(2, numel(dt), 1);
    % Unload for all days when temperature is less than threshold
    ctrl(temp < thres) = 3;
    month = year/12;
    % Always start with three months of loading
    ctrl(cumsum(dt) <= 3*month) = 1;
    tmp = reshape(repmat(1:ceil(numel(dt)/5), 5, 1), [], 1);
    tmp = tmp(1:numel(dt));
    [tmp, n] = rlencode([ctrl, tmp], 1);
    
    
    ctrl = tmp(:,1);
    timestep = rldecode((1:numel(ctrl))', n, 1);
    
    dt = accumarray(timestep, dt);
    
    % Make step
    step = struct('val', dt, 'control', ctrl);
    schedule = struct('step', step);
end

%-------------------------------------------------------------------------%
function schedule = addLoadPhase(schedule, rate)
    % Define well controls for loading phase
    ctrl    = 1;
    hasload = strcmpi({schedule.control(ctrl).type}, 'load');
    if ~hasload
        W = schedule.control(ctrl).W;
        % Set deep well rates
        % Cold support wells operate at fixed BHP
        pRef = W(1).refDepth*norm(gravity)*1000;
        [W(1:2).type] = deal('bhp');
        [W(1:2).val ] = deal(pRef);
        [W(1:2).sign] = deal(-1);
        % Hot wells operate at fixed production rate
        [W(3).type] = deal('rate');
        [W(3).val ] = deal(rate);
        [W(3).lims] = deal(struct('bhp', 2*pRef));
        [W(3).sign] = deal(1);
        % Insert control
        schedule.control(ctrl).W = W;
        schedule.control(ctrl).type = 'load';
    end
end

%-------------------------------------------------------------------------%
function schedule = addRestPhase(schedule)
    ctrl = 2;
    hasrest = numel(schedule.control) > ctrl-1 && ...
             strcmpi({schedule.control(ctrl).type}, 'rest');
    if ~hasrest
        W = schedule.control(1).W;
        % Deep wells are inactive
        [W.status] = deal(false);
        % Insert control
        schedule.control(ctrl) = schedule.control(1);
        schedule.control(ctrl).W = W;
        schedule.control(ctrl).type = 'rest';
    end
end

%-------------------------------------------------------------------------%
function schedule = addUnloadPhase(schedule, rate)
    % Define well controls for unloading phase
    ctrl = 3;
    hasunload = numel(schedule.control) > ctrl-1 && ...
                strcmpi({schedule.control(ctrl).type}, 'unload');
    if ~hasunload
        W = schedule.control(1).W;
        pRef = W(1).refDepth*norm(gravity)*1000;
        % Deep wells are inactive
        % Set shallow well controls
        % Cold support wells operate at fixed injection rate
        [W(1:2).type] = deal('rate');
        [W(1:2).val ] = deal(rate/2);
        [W(1:2).lim ] = deal(struct('bhp', 2*pRef));
        [W(1:2).sign] = deal(1);
        % Hot wells operate at fixed production rate
        [W(3).type ] = deal('rate');
        [W(3).val  ] = deal(-rate);
        [W(3).lims ] = deal(struct('bhp', 1*atm));
        [W(3).sign ] = deal(-1);
        % Insert control
        schedule.control(ctrl) = schedule.control(1);
        schedule.control(ctrl).W = W;
        schedule.control(ctrl).type = 'unload';
    end
end


%-------------------------------------------------------------------------%
function [schedule, report, isAltered] = controlLogicFn(state, schedule, report, i, monitorCells)
    if ~isfield(schedule.step, 'control0')
        schedule.step.control0 = schedule.step.control;
    end
    
    if ~isfield(schedule, 'T0')
        schedule.T0 = state.T(monitorCells);
    end
    
    isAltered = false;
    
    ctrl     = schedule.step.control(i);
    loadTime = 14*day;
   
    
    switch ctrl
        case 1
            % Load phase. Adjust injection temperature in shallow reservoir
            % based on production flow rate and temperature from the deep
            % reservoir. If the rate is below the threshold, we switch to
            % the rest phase.
            T = state.T(monitorCells);
            if any(T > schedule.T0 + 5*Kelvin) % negative since we are producing
                stop = find(schedule.step.control(i:end) ~= 1, 1, 'first') + i - 1;
                if isempty(stop), stop = numel(schedule.step.val); end
                schedule.step.control(i+1:stop) = schedule.step.control0(i+1:stop);
                isAltered = true;
            end
        case 3
            % Unload phase. We produce water from the shallow reservoir
            % until the produced temperature drops below the threshold.
            T = getWellOutput({state.wellSol}, 'T');
            if T(3) < (45 + 273.15)*Kelvin
                time = cumsum(schedule.step.val(i+1:end));
                stop = find(time >= loadTime, 1, 'first') + i - 1;
                if isempty(stop), stop = numel(schedule.step.val); end
                schedule.step.control(i+1:stop) = 1;
                isAltered = true;
            end
    end
end