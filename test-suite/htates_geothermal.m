function setup = htates_geothermal(varargin)
%High-temperature aquifer thermal energy storage in synthetic reservoir

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

    % Step 1: Test case description and options
    %---------------------------------------------------------------------%
    description = ['High-temperature aquifer thermal energy storage ', ...
        'in synthetic reservoir. See "Simulation of geothermal '     , ...
        'systems using MRST", DOI: 10.1017/9781009019781.018.'       ];
    options = struct('numCycles', 4); % Cycle = charge-rest-discharge-rest
    % Process optinal input arguments
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
                                        options, description, varargin{:});
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%
    
    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    % Module dependencies
    require ad-core ad-props geothermal compositional upr linearsolvers
    gravity reset on
    data  = load(fullfile(getDatasetPath('geothermal'), 'htates-reservoir.mat'));
    G     = data.G;
    rock  = data.rock;
    % Viscosity and density are p/T-dependent, and will be set later
    fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'mu', 1, 'rho', 1);
    % Assign thermal properties of the fluid, with equation of state from
    % Spivey et. al (2004)
    fluid = addThermalFluidProps(fluid, 'useEOS', true);
    rock  = addThermalRockProps(rock);            % Thermal rock properties
    W     = data.W;
    layer = data.layer;
    % The model has two groups of four wells each, which we refer to as hot and
    % cold. The hot group will be used to inject hot water for storage during
    % the summer months, and extract the hot water for heating in the winter
    % months. The cold well is used to provide pressure support, with a low BHP
    % during storage, and a hight BHP during extraction. In between storage and
    % extraction, there will be a period of rest.
    p0    = 70*barsa;         % Reference reservoir pressure
    K0    = 273.15*Kelvin;    % Zero degrees Celcius
    rate  = 1000*meter^3/day; % Injection rate
    bhpSt = p0;               % Storage BHP
    bhpEx = 85*barsa;         % Extraction BHP
    Thot  = K0 + 100*Kelvin;  % Temperature of stored water
    Tcold = K0 + 30*Kelvin;   % Temperature of injected water
    
    [W.compi]  = deal(1);                                % Single-phase
    W          = addThermalWellProps(W, G, rock, fluid); % Add thermal props
    WSt        = W;
    [hNo, cNo] = deal(0);
    xmid       = mean(G.cells.centroids);                % Reservoir center
    for wNo = 1:numel(W)
        w     = W(wNo);
        cells = w.cells;
        xw    = mean(G.cells.centroids(cells,:));
        w.refDepth = min(xw(:,3));
        if xw(1) < xmid(1)
            % Hot group: set injection rate
            hNo = hNo + 1;
            w.type = 'rate';
            w.val  = rate;
            w.sign = 1;
            w.name = ['H', num2str(hNo)];
            w.T    = Thot;
        else
            % Cold group: set bhp
            cNo    = cNo + 1;
            w.type = 'bhp';
            w.val  = bhpSt;
            w.sign = -1;
            w.name = ['C', num2str(cNo)];
            w.T    = Tcold;
        end
        WSt(wNo) = w;
    end
    [WSt.components] = deal(1);
    % Make well structures for the extraction and rest periods
    [WEx, WRe] = deal(WSt);
    isHot = cellfun(@(n) strcmpi(n(1), 'H'), {WSt.name}); % Hot group
    [WEx(isHot).val ] = deal(-rate); % Reverse rate
    [WEx(~isHot).val] = deal(bhpEx); % Set extraction bhp
    for wNo = 1:numel(WEx)
        WEx(wNo).sign = -WEx(wNo).sign; % Reverse sign
    end
    [WRe.status] = deal(false); % Rest wells are not active (status = false)
    % Storage and extraction periods are four months long, whereas the rest
    % periods are two months
    month  = year/12;
    % Storage period
    [dtSt, dtEx] = deal(rampupTimesteps(4*month, 20*day, 5));
    [nSt , nEx]  = deal(numel(dtSt));
    % Rest period
    dtRe = rampupTimesteps(2*month, 20*day, 2);
    nRe  = numel(dtRe);
    % Controls: storage = 1, rest = 2, extraction = 3
    control   = struct('W', {WSt, WRe, WEx}, 'src', [], 'bc', []);
    controlNo = [1*ones(nSt,1);  % Storage
                 2*ones(nRe,1);  % Rest
                 3*ones(nEx,1);  % Extraction
                 2*ones(nRe,1)]; % Rest
    dt        = [dtSt; dtRe; dtEx; dtRe]; % Timesteps
    % We simulate eight cycles (= eight years) of HTATES
    numCycles   = options.numCycles;
    controlNo = repmat(controlNo, numCycles, 1);
    dt        = repmat(dt, numCycles, 1);
    schedule  = struct('step'   , struct('val', dt, 'control', controlNo), ...
                       'control', control                                );
    model = GeothermalModel(G, rock, fluid); % Make model
    % The EOS is valid for pressure/temperature within a given range. We
    % provide these to the model so that pressure/temperature are within these
    % during the nonlinear solution step
    model.maximumPressure    = 200e6;    % Maximum pressure
    model.minimumTemperature = K0;       % Minimum temperature 
    model.maximumTemperature = K0 + 275; % Maximum temperature 
    % To speed up the simulation, we use mex-accelereated AD backend ...
    model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
    % The reserois should be in thermal and hydrostatic equilibrium
    state0   = initResSol(G, p0, 1);
    dTdz     = 30*Kelvin/(kilo*meter); % Geothermal gradient
    T0       = @(x) K0 + 20*Kelvin + dTdz*x(:,3);
    state0.T = T0(G.cells.centroids);
    % state0.T = K0 + 20*Kelvin;
    % To find the equilibrium state, we simply assign a fixed pressure BC to
    % the 100 topmost faces, and simulate for 50 years with looong timesteps
    f       = boundaryFaces(G);
    [~, ix] = sort((G.faces.centroids(f,3)));
    nf      = 1;
    bc      = addBC([], f(ix(1:nf)), 'pressure', p0, 'sat', 1);
    bc      = addThermalBCProps(bc, 'T', T0(G.faces.centroids(f(ix(1:nf)),:)));
    % Make schedule
    dtPre       = rampupTimesteps(50*year, 10*year, 3);
    schedulePre = simpleSchedule(dtPre, 'bc', bc, 'W', schedule.control(1).W);
    [schedulePre.control(1).W.status] = deal(false);
    schedule.control      = [schedulePre.control, schedule.control];
    schedule.step.control = [ones(numel(dtPre),1); schedule.step.control + 1];
    schedule.step.val     = [dtPre; schedule.step.val];
     % Plotting
    plotOptions = {'View'              , [-10, 45]         , ...
                   'Size'              , [1300, 650]       , ...
                   'PlotBoxAspectRatio', [4.96, 2.10, 1.00]};
    %---------------------------------------------------------------------%
    
    % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    setup = packTestCaseSetup(mfilename,              ...
                          'description', description, ...
                          'options'    , options    , ...
                          'state0'     , state0     , ...
                          'model'      , model      , ...
                          'schedule'   , schedule   , ...
                          'plotOptions', plotOptions);
    %---------------------------------------------------------------------%
    
end