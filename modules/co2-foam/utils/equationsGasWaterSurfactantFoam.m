function [problem, state] = equationsGasWaterSurfactantFoam(state0, state, model, dt, drivingForces, varargin)
%
% SYNOPSIS:
%   function [problem, state] = equationsGasWaterSurfactantFoam(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: 
%   Assemble the linearized equations for an gas-water-surfactant
%   system with foam generation, computing both the residuals and the Jacobians. 
%   Returns the result as an instance of the class LinearizedProblem which 
%   can be solved using instances of LinearSolverAD.
%
%   This implementation uses the statefunctions formalism.
%
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD
%

%{ 
Copyright 2009-2023 SINTEF

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

    opt = struct( ...
        'Verbose'        , mrstVerbose, ...
        'reverseMode'    , false      , ...
        'velocCompMethod', 'square'   , ...
        'resOnly'        , false      , ...
        'iteration'      , -1           ...
    );
    opt = merge_options(opt, varargin{:});

    op    = model.operators;
    fluid = model.fluid;

    % Properties at current timestep
    [p, sW, cs, wellSol] = model.getProps(state, ...
        'pressure', 'water', 'foam', 'wellsol');

    % Properties at previous timestep
    [p0, sW0, cs0, wellSol0] = model.getProps(state0, ...
        'pressure', 'water', 'foam', 'wellsol');

    [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
    % Typically the primary well variables are:
    %  - the phase well rates (qWell)
    %  - the bottom hole pressures
    %  - the surfactant concentrations, at injection and production wells,
    %    contained in the wellVars, wellVarNames, wellMap structures

    % Initialize independent variables.
    if ~opt.resOnly
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode
            [p, sW, cs, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, sW, cs, wellVars{:});
        else
            wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
            [p0, sW0, cs0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, sW0, cs0, wellVars0{:}); %#ok
        end
    end

    % We will solve for pressure, water saturation (oil saturation follows via
    % the definition of saturations), surfactant concentration and well rates +
    % bhp.
    primaryVars = ['pressure', 'sW', 'foam', wellVarNames];

    sG   = 1 - sW;
    sG0  = 1 - sW0;
    sat  = {sW, sG};
    sat0 = {sW0, sG0};
    
    % Update state with AD-variables
    state  = model.setProps(state , {'s', 'pressure', 'foam'}, {sat , p , cs });
    state0 = model.setProps(state0, {'s', 'pressure', 'foam'}, {sat0, p0, cs0});
    % Set up properties
    state = model.initStateFunctionContainers(state);
    
    % Get properties
    %---------------------------------------------------------------------%
    [b, pv]            = model.getProps(state , 'ShrinkageFactors','PoreVolume');
    [b0, pv0]          = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
    [Cwga]             = model.getProps(state , 'ConcentrationsPartitioning');
    [cW0, cG0, cads0]  = model.getProps(state0, 'sfwat', 'sfgas', 'sfads'); 
    [pressures,  rho]  = model.getProps(state , 'PhasePressures', 'Density');
    mob                = model.getProps(state , 'Mobility'); 
    [phaseFlux, flags] = model.getProps(state , 'PhaseFlux', 'PhaseUpwindFlag');

    [cW, cG, cads, ~, ~, mSrf] = deal(Cwga{:});
    [bW, bG]                   = deal(b{:});
    [bW0, bG0]                 = deal(b0{:});
    mobG                       = mob{2};

    [vW, vG]      = deal(phaseFlux{:});
    [upW, upG]  = deal(flags{:});    
    
    % Upstream weight b factors and multiply by interface fluxes to obtain
    % the fluxes at standard conditions.
    bGvG = op.faceUpstr(upG, bG).*vG;
    bWvW = op.faceUpstr(upW, bW).*vW;
    %---------------------------------------------------------------------%

    % Equations
    %---------------------------------------------------------------------%
    % Conservation of mass for water
    water = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0 );
    divWater = op.Div(bWvW);
    water = water + divWater;
    
    % Conservation of mass for gas
    gas = (1/dt).*( pv.*bG.*sG - pv0.*bG0.*sG0 );
    divGas = op.Div(bGvG);
    gas = gas + divGas;
    
    % Computation of adsoprtion term and poro
    poro = model.rock.poro;
    ads_term = fluid.rhoRSft.*((1-poro)./poro).*(cads - cads0);
    
    % Determine surfactant partitioning
    [SW, SG] = model.surfactantPartitioning();

    % Initialize to all-zero ADI variable of same dimensions as water and gas.
    acc = p.*0;
    tra = acc;
    % Add accumulation and transport terms. 
    % Unit: mass of surfactant in each cell
    if SG
        % Add accumulation and transport term for surfactant in gas phase
        acc   = acc + fluid.rhoGS*(pv.*bG.*sG.*cG - pv0.*bG0.*sG0.*cG0);
        vSft  = op.faceUpstr(upG, cG).*vG;
        bvSft =  op.faceUpstr(upG, bG).*vSft;        
        tra   = tra + fluid.rhoGS*op.Div(bvSft);
    end
    if SW
        % Add accumulation and transport term for surfactant in water phase
        acc   = acc + fluid.rhoWS*(pv.*bW.*sW.*cW - pv0.*bW0.*sW0.*cW0);
        vSft  = op.faceUpstr(upW, cW).*vW;
        bvSft = op.faceUpstr(upW, bW).*vSft;
        tra   = tra + fluid.rhoWS*op.Div(bvSft);
    end
    
    % The combined mass conservation equation 
    foam = (op.pv/dt).*ads_term+1.0/dt.*acc + tra;

    % Safeguarding surfactant concentration for elements with small
    % Jacobian
    if ~isreal(p)
        small   = 1e-6;
        typical = sqrt(small)*mean(abs(diag(foam.jac{3})));
        bad     = abs(diag(foam.jac{3}))<=typical;
        if any(bad)
            foam(bad) = cs(bad);
        end
    end
    
    % Postprocess equations (not supported with diagonal autodiff backend)
    if ~opt.resOnly ...
            && ~isa(model.AutoDiffBackend, 'DiagonalAutodiffBackend')
        epsilon = 1.e-8;
        % The first way is based on the diagonal values of the resulting
        % Jacobian matrix
        eps = sqrt(epsilon)*mean(abs(diag(foam.jac{3})));
        % Sometimes there is no water in the whole domain
        if (eps == 0.)
            eps = epsilon;
        end
        % bad marks the cells prolematic in evaluating Jacobian
        bad = abs(diag(foam.jac{3})) < eps;
        % the other way is to choose based on the water saturation
        foam(bad) = cs(bad);
    end
    
    eqs        = {water   , gas   , foam  };
    names      = {'water' , 'gas' , 'foam'};
    types      = {'cell'  , 'cell', 'cell'};
    components = {{cW, cG}};
    %---------------------------------------------------------------------%
    
    % Add well contributions and boundary conditions
    %---------------------------------------------------------------------%
    % Add non-well boundary conditions
    [eqs, state] = addBoundaryConditionsAndSources( ...
        model, eqs, names, types, state, pressures, sat, mob, rho, {}, ...
        components, drivingForces);
    % Finally, add in and setup well equations
    model.FacilityModel.ReservoirModel.sfwat = cW;
    model.FacilityModel.ReservoirModel.sfgas = cG;
    [eqs, names, types, state.wellSol] ...
        = model.insertWellEquations(eqs, names, types, ...
        wellSol0, wellSol, wellVars, wellMap, ...
        p, mob, rho, {}, components, dt, opt);
    %---------------------------------------------------------------------%
    
    % Set extra state output
    %---------------------------------------------------------------------%
    % Save surfactant concentrations in state variable
    state = model.setProp(state, 'sfwat', value(cW));
    state = model.setProp(state, 'sfgas', value(cG));
    state = model.setProp(state, 'sfads', value(cads));                                               
    
    state.mobG = value(mobG);
    state.mSrf = value(mSrf);
    
        % @@ Do we need BOM?
%     BOM  = model.getProps(state,'BlackoilMobility');
%     state.mob0 = value(BOM{2});
    %---------------------------------------------------------------------%

    % Pack linearized problem
    %---------------------------------------------------------------------%
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    %---------------------------------------------------------------------%

end


