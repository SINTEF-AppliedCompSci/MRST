function setup = qfs_geothermal(varargin)
%Test case describing geothermal heat storage in quarter five-spot pattern

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
    description = 'Geothermal heat storage in quarter five-spot pattern';
    options = struct('NaCl', false);
        [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
                                        options, description, varargin{:});
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%
    
     % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    % Module dependencies
    require ad-core ad-props ad-blackoil geothermal compositional upr
    % Set up qfs base case
    setup     = qfs_wo(options.extra{:});
    G        = setup.model.G;
    rock     = setup.model.rock;
    schedule = setup.schedule;
    state0   = setup.state0;
    % Make fluid
    fluid = initSimpleADIFluid('phases', 'W'                  , ...
                               'n'     , 1                    , ...
                               'rho'   , 1000*kilogram/meter^3, ...
                               'mu'    , 0.1*centi*poise      , ...
                               'c'     , 1e-10/Pascal         , ...
                               'pRef'  , 1*atm                );
    % Add thermal properties
    fluid = addThermalFluidProps(fluid, 'useEOS', false               , ...
                                        'cT'    , 1e-6/Kelvin         , ...
                                        'cX'    , log(2)              , ...
                                        'TRef'  , (273.15 + 20)*Kelvin);
    % Make rock
    rock  = addThermalRockProps(rock);
    % Make compositional fluid (optional)
    compFluid = {};
    if options.NaCl
        compFluid = CompositionalBrineFluid(          ...
            {'H2O'             , 'NaCl'            }, ... % Names
            [18.015281*gram/mol, 58.442800*gram/mol], ... % Molar masses
            [0                 , 1e-6              ]);    % Molecular diffusivities
        compFluid = {compFluid};
    end
    % Set up model
    model = GeothermalModel(G, rock, fluid, compFluid{:});
    % Make initial state
    state0.s = ones(G.cells.num, 1);
    state0.T = repmat((273.15 + 20)*Kelvin, G.cells.num, 1);
    state0.components = ones(G.cells.num,1);
    if options.NaCl
        x0 = repmat(0.1, G.cells.num, 1);
        state0.components = [state0.components - x0, x0];
    end
    % Set schedule properties
    [schedule.control(1).W] ...
        = deal(addThermalWellProps(schedule.control(1).W, G, rock, fluid, ...
                                                'T', 273.15*Kelvin + 200));
    [schedule.control(1).W.components] = deal(1);
    [schedule.control.W.compi] = deal(1);
    %---------------------------------------------------------------------%
    
    % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    setup = packTestCaseSetup(mfilename,              ...
                          'description', description, ...
                          'options'    , options    , ...
                          'state0'     , state0     , ...
                          'model'      , model      , ...
                          'schedule'   , schedule   );
    %---------------------------------------------------------------------%
    
end