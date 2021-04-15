function [description, options, state0, model, schedule, plotOptions] = norne_field_bo(varargin)
%Example from the example suite, see description below. Warning: This
%example will likely not run without tweaking the solvers (see
%`fieldModelNorneExample`).
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

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
    % One-line description
    description ...
        = ['The full Norne field model. Made available through the ', ...
           'Open Porous Media Project, https://opm-project.org/'    ];
    % Optional input arguments
    options = struct();
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-blackoil ad-props deckformat
    gravity reset on
    mrstVerbose on
    useMex = true;
    [deck, output] = getDeckOPMData('norne', 'NORNE_ATW2013');    
    opm = mrstPath('opm-tests');
    assert(~isempty(opm), 'You must register https://github.com/opm/opm-tests as a module!');
    % Build model from EGRID/INIT files
    egrid = readEclipseOutputFileUnFmt([output.opm.location, '.EGRID']);
    init  = readEclipseOutputFileUnFmt([output.opm.location, '.INIT']);
    [Ge, rock_ecl, Ne, Te] = initGridFromEclipseOutput(init, egrid, 'outputSimGrid', true);
    Gviz = computeGeometry(Ge{1}); % Visualization grid
    Gsim = Ge{2};                  % Simulation grid
    % Rock
    rock = initEclipseRock(deck);
    rock = compressRock(rock, Gsim.cells.indexMap);
    % Fluid
    fluid = initDeckADIFluid(deck, 'useMex', useMex);
    % Setup model, but skip setting up the operators since we do not have a
    % proper grid
    model = GenericBlackOilModel(Gsim, [], fluid, 'disgas', true, 'vapoil', true, 'inputdata', deck);
    % Finally set up the connectivity graph from output
    model.rock = rock;
    model.operators = setupOperatorsTPFA(Gsim, rock, 'neighbors', Ne, 'trans', Te);
    % Set up everything
    [state0, model, schedule] = initEclipseProblemAD(deck                      , ...
                                                     'model'           , model , ...
                                                     'TimestepStrategy', 'ds'  , ...
                                                     'useCPR'          , true  , ...
                                                     'useMex'          , useMex);
    % Set tolerances
    model.toleranceCNV = 1e-2;
    model.toleranceMB = 1e-7;
    % Set well tolerances
    model.FacilityModel = GenericFacilityModel(model);
    model.FacilityModel.toleranceWellBHP = 1e-3;
    model.FacilityModel.toleranceWellRate = 5e-3;
    % Reset just in case
    model.FlowDiscretization = [];
    model.FlowPropertyFunctions = [];
    model.PVTPropertyFunctions = [];
    % Set autodiff backend stuff
    model.AutoDiffBackend.useMex = true;
    model.AutoDiffBackend.rowMajor = true;

    model = model.validateModel();
    % Use the alternative more rigorous crossflow definition for component
    % fluxes
    xflow = WellComponentTotalVolumeBalanceCrossflow(model);
    xflow.onlyLocalDerivatives = true;
    model.FacilityModel.FacilityFlowDiscretization.ComponentTotalFlux = xflow;

    useFlag = false;
    model.PVTPropertyFunctions.Viscosity.useSaturatedFlag = useFlag;
    model.PVTPropertyFunctions.ShrinkageFactors.useSaturatedFlag = useFlag;

    RvMax = model.PVTPropertyFunctions.getStateFunction('RvMax');
    RvMax.rvReduction = 0.5; % Not tested, but in the deck
    model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('RvMax', RvMax);
    % PLotting
    plotOptions = {'View'              , [85, 30]  , ...
                   'PlotBoxAspectRatio', [1,1,0.3] , ...
                   'Size'              , [700, 600], ...
                   'visualizationGrid' , Gviz      };
end
