function setup = qfs_foam(varargin)
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
    % One-line description
    description = '';
    
    options = struct( ...
        'cartDims'  , [  20,   20,  20], ...
        'physDims'  , [1400, 1400, 100], ...
        'partition' , true             , ...
        'usePermDep', true             , ...
        'useAds'    , false            , ...
        'noShear'   , true               ...
    );
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
    options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil co2-foam ad-eor
    % ad-eor needed because of computeVelocTPFA
    dataPath = getDatasetPath('co2_foam');
    fn = fullfile(dataPath, 'field', 'RUN', 'QBOX2_GAS.DATA');
    [state0, model, schedule] = initEclipseProblemAD(fn, 'useMex', true, 'rowMajorAD', true);
    
    % Make grid
    G = cartGrid(options.cartDims(1:2), options.physDims(1:2));
    % Five fine cells at the top/bottom
    [nTop, nBot] = deal(min(floor(options.cartDims(3)/2), 5));
    nMid = options.cartDims(3) - (nTop+nBot)';
    dzTop = repmat(2 , nTop, 1);
    dzBot = repmat(10, nBot, 1);
    w = linspace(0,1,nMid)'; % Transition from fine to coarse in-between
    dzMid = dzTop(1).*(1-w) + dzBot(1).*w;
    dz = [dzTop; dzMid; dzBot];
    dz = dz./sum(dz).*options.physDims(3);
    % Make a layered grid
    G = makeLayeredGrid(G, dz);
    G.nodes.coords(:,3) = G.nodes.coords(:,3) + 800;
    G = computeGeometry(G);
    G.cartDims = [options.cartDims(1:2), numel(dz)];
    surfactantSystem = 15; % Brij L23
    % Mobility reduction function parameters
    foam = getFoam(surfactantSystem,'noShear',true);  % Foam properties

    % Surfactant properties
    surf = getSurfactant(surfactantSystem); % Dissolution properties
    surf.Cpart = 0.1; % @@ Does not run with Cpart = 1

    fluid = addSimpleFoamProperties(model.fluid, ...
        'foam', foam, ...
        'cCrit', 2e-3, ...                    % Critical concentration by weight
        'cDecl', 2.5e3, ...                   % Foam effect decline interval
        'adsMax',surf.adsMax, ...             % Maximum adsorption 
        'adsSat',surf.adsSat, ...             % Threshold concentration for adsorption
        'surfingas', false, ...               % Gas-soluble surfactant
        'mobMultIsLinear',true, ...           % Is the mobility multiplier a linear function of gas velocity?
        'usePermDep', options.usePermDep, ... % Is foam strength dependent on absolute permeability?
        'noShear', options.noShear      );    % Velocity dependence
    
    % Turn off adsorption
    if ~options.useAds, fluid.adsMax = 0; end

    rock = makeRock(G, [1,1,0.1]*260*milli*darcy, model.rock.poro(1));
    model.rock.perm = rock.perm;
    model.rock.poro = rock.poro;
    
    fluid.foam.PermDep = @(k) (k/(260*milli*darcy)).^2; 

    model = GasWaterSurfactantFoamModel(G, model.rock, fluid); 
    
    if options.partition
        model.fluid.Cpart = surf.Cpart;
    else
        model.fluid.surfingas = true;
    end
    
    state0 = initResSol(G, mean(state0.pressure), state0.s(1,:));
    % Adsorbed surfactant
    state0.cA = zeros(G.cells.num,1);
    % Initial volume concentration
    [state0.cs, state0.csmax] = deal(zeros(G.cells.num,1));

    % Stored surfactant concentrations in water, gas and adsorbed
    %state0.cW = ones(G.cells.num,1)*initialcW;
    state0.cW = zeros(G.cells.num,1);
    state0.cG = zeros(G.cells.num,1);
    
    
    [ii, jj, kk] = gridLogicalIndices(G);
    W = [];
    lims = struct('bhp', 150*barsa);
    W0 = schedule.control(1).W;
    W = verticalWell(W, G, model.rock, 1, 1, max(kk)-3:max(kk), ...
        'type', W0(1).type, 'val', W0(1).val, 'compi', W0(1).compi, 'lims', lims);
    W = verticalWell(W, G, model.rock, max(ii), max(jj), max(kk)-3:max(kk), ...
        'type', W0(2).type, 'val', W0(2).val, 'compi', W0(2).compi);
    schedule.control(1).W = W;
    [schedule.control.W.cs] = deal(0); % Default to no surfactant injection
    schedule.control(1).W(1).cs = 0.01; % 1.0 wt% surfactant.
    schedule.control(2) = schedule.control(1);
    schedule.control(2).W(1).cs = 0; % Stop surfactant injection.

    % Time step vector
    % Initial period with surfactant injection (control 1) followed by period
    % without surfactant (control 2).
    tvec = rampupTimesteps(4*year, 90*day);
    n1 = length(tvec);
    tvec = [tvec; rampupTimesteps(50*year, 180*day, 3)];
    schedule.step.control = ones(size(tvec))*2;
    schedule.step.control(1:n1) = 1;
    schedule.step.val = tvec;

    plotOptions =  {'View', [70,20], 'PlotBoxAspectRatio', [1,1,0.25]};
        
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
    
end