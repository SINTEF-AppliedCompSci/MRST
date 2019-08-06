function [problem, state] = equationsThreePhaseSurfactantPolymer(state0, state, model, ...
    dt, drivingForces, varargin)
    %
    % SYNOPSIS:
    %   function [problem, state] = equationsThreePhaseSurfactantPolymer(state0, state, model, dt, drivingForces, varargin)
    %
    % DESCRIPTION:
    %   Assemble the linearized equations for a blackoil surfactant-polymer system,
    %   computing both the residuals and the Jacobians. Returns the result as
    %   an instance of the class LinearizedProblem which can be solved using
    %   instances of LinearSolverAD.
    %
    %   A description of the modeling equations can be found in the directory
    %   ad-eor/docs.
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
    % SEE ALSO: LinearizedProblem, LinearSolverAD, OilWaterSurfactantModel
    %
    %{
    Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('Verbose'        , mrstVerbose, ...
                 'reverseMode'    , false      , ...
                 'velocCompMethod', 'square'   , ... 
                 'resOnly'        , false      , ...
                 'iteration'      , -1 );
    opt = merge_options(opt, varargin{:});

    % Shorter names for some commonly used parts of the model and forces.
    W     = drivingForces.W;
    fluid = model.fluid;
    s    = model.operators;
    G = model.G;

    % Properties at current timestep
    [p, sW, sG, rs, rv, cp, cpmax, cs, csmax, wellSol] = model.getProps(state, ...
        'pressure', 'water', 'gas', 'rs', 'rv', 'polymer', 'polymermax', ...
        'surfactant', 'surfactantmax', 'wellSol');

    % Properties at previous timestep
    [p0, sW0, sG0, rs0, rv0, cp0, cpmax0, cs0, csmax0, wellSol0] = model.getProps(state0, ...
        'pressure', 'water', 'gas', 'rs', 'rv', 'polymer', 'polymermax', ...
        'surfactant', 'surfactantmax', 'wellSol');

    [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
    % Typically the primary well variables are :
    %  - the phase well rates (qWell)
    %  - the bottom hole pressures
    %  - the surfactant concentrations, at injection and production wells,
    %    contained in the wellVars, wellVarNames, wellMap structures

    % Initialization of primary variables ----------------------------------
    st  = model.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
    st0 = model.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);

    if model.disgas || model.vapoil
        % X is either Rs, Rv or Sg, depending on each cell's saturation status
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
        gvar = 'x';
    else
        x = sG;
        gvar = 'sG';
    end

    % Initialize independent variables.
    if ~opt.resOnly
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode
            [p, sW, x, cp, cs, wellVars{:}] = initVariablesADI(p, sW, x, cp, cs, wellVars{:});
        else
            x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
            % Set initial gradient to zero
            wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
            [p0, sW0, x0, cp0, cs0, wellVars0{:}] = initVariablesADI(p0, sW0, x0, cp0, cs0, wellVars0{:});
            %         clear zw;
            [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1-sW0, x0, rs0, rv0, p0);
        end
    end

    if ~opt.reverseMode
        % Compute values from status flags. If we are in reverse mode, these
        % values have already converged in the forward simulation.
        [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, 1-sW, x, rs, rv, p);
    end

    % We will solve for pressure, water and gas saturation (oil saturation follows via
    % the definition of saturations), surfactant concentration and well rates +
    % bhp.
    primaryVars = {'pressure', 'sW', gvar, 'polymer', 'surfactant', wellVarNames{:}};

    % Evaluate relative permeability
    
    sO  = 1 - sW  - sG;
    sO0 = 1 - sW0 - sG0;
    
    if model.water
        sat = {sW, sO, sG};
        sat0 = {sW0, sO0, sG0};
    else
        sat = {sO, sG};
        sat0 = {sO0, sG0};
    end

    % Update state with AD-variables
    state = model.setProps(state  , {'s', 'pressure', 'rs', 'rv', 'polymer', 'surfactant'}, {sat , p , rs , rv, cp, cs});
    state0 = model.setProps(state0, {'s', 'pressure', 'rs', 'rv', 'polymer', 'surfactant'}, {sat0, p0, rs0, rv0, cp0, cs0});
    % Set up properties
    state = model.initPropertyContainers(state);

    % EQUATIONS ---------------------------------------------------------------
    pBH = wellVars{wellMap.isBHP};
    Nc = computeCapillaryNumber(p, cs, pBH, W, fluid, G, s, 'velocCompMethod', ...
                                opt.velocCompMethod);
    state.CapillaryNumber = Nc;

    [b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
    [b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
    [phaseFlux, flags] = model.getProps(state, 'PolymerPhaseFlux',  'PhaseUpwindFlag');
    [pressures, mob, rho] = model.getProps(state, 'PhasePressures', 'Mobility', 'Density');

    [bW, bO, bG]       = deal(b{:});
    [bW0, bO0, bG0]    = deal(b0{:});
    [vW, vO, vG, vP]       = deal(phaseFlux{:});
    [upcw, upco, upcg] = deal(flags{:});
    [mobW, mobO, mobG] = deal(mob{:});
    
    muWeffMult = model.getProp(state, 'EffectiveMixturePolymerViscMultiplier');


    % Store fluxes / properties for debugging / plotting, if requested.
    if model.outputFluxes
        state = model.storeFluxes(state, vW, vO, vG);
    end
    if model.extraStateOutput
        state = model.storebfactors(state, bW, bO, bG);
        state = model.storeMobilities(state, mob{:});
        state = model.storeUpstreamIndices(state, upcw, upco, upcg);
    end

    % EQUATIONS -----------------------------------------------------------

    % Upstream weight b factors and multiply by interface fluxes to obtain the
    % fluxes at standard conditions.
    
    bWvW = s.faceUpstr(upcw, bW).*vW;
    bOvO = s.faceUpstr(upco, bO).*vO;
    bGvG = s.faceUpstr(upcg, bG).*vG;
    bWvP = s.faceUpstr(upcw, bW).*vP;

    % Conservation of mass for water
    water = (1/dt).*(pv.*bW.*sW - pv0.*bW0.*sW0);
    divWater = s.Div(bWvW);

    % Conservation of mass for oil
    if model.vapoil
        % The model allows oil to vaporize into the gas phase. The conservation
        % equation for oil must then include the fraction present in the gas
        % phase.
        rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
        % Final equation
        oil = (1/dt).*( pv .*(bO.* sO  + rv.* bG.* sG) - ...
                        pv0.*(bO0.*sO0 + rv0.*bG0.*sG0));
        divOil = s.Div(bOvO + rvbGvG);
    else
        oil = (1/dt).*(pv.*bO.*sO - pv0.*bO0.*sO0 );
        divOil = s.Div(bOvO);
    end

    % Conservation of mass for gas
    if model.disgas
        % The gas transported in the oil phase.
        rsbOvO = s.faceUpstr(upco, rs).*bOvO;

        gas = (1/dt).*(pv.* (bG.* sG  + rs.* bO.* sO) - ...
                       pv0.*(bG0.*sG0 + rs0.*bO0.*sO0 ));
        divGas = s.Div(bGvG + rsbOvO);
    else
        gas = (1/dt).*(pv.*bG.*sG - pv0.*bG0.*sG0 );
        divGas = s.Div(bGvG);
    end

    % Computation of adsoprtion term
    vSft   = s.faceUpstr(upcw, cs).*vW;
    bWvSft = s.faceUpstr(upcw, bW).*vSft;
    poro = model.rock.poro;
    adsp  = model.getProp(state, 'PolymerAdsorption');
    adsp0 = model.getProp(state0, 'PolymerAdsorption');
    adss  = model.getProp(state, 'SurfactantAdsorption');
    adss0 = model.getProp(state0, 'SurfactantAdsorption');
    adss_term = fluid.rhoRSft.*((1-poro)./poro).*(adss - adss0);

    % Conservation of polymer in polymer phase:
    polymer = ((1-fluid.dps)/dt).*(pv.*bW.*sW.*cp - ...
                                     pv0.*fluid.bW(p0).*sW0.*cp0) + (s.pv/dt).* ...
              ( fluid.rhoR.*((1-poro)./poro).*(adsp - adsp0));
    divPolymer = s.Div(bWvP);
    
    % Conservation of surfactant in water:
    surfactant = (1/dt)*((pv.*bW.*sW.*cs - pv0.*bW0.*sW0.*cs0) + adss_term);
    divSurfactant = s.Div(bWvSft);
    
    % Applying correction to the surfactant and polymer equation when the Jacobian is
    % prolematic for some cells.
    % Typically it is due to totally and almost non-existence of water.    
    if ~opt.resOnly
        epsilon = 1.e-8;
        % the first way is based on the diagonal values of the resulting
        % Jacobian matrix
        epsP = sqrt(epsilon)*mean(abs(diag(polymer.jac{4})));
        epsS = sqrt(epsilon)*mean(abs(diag(surfactant.jac{5})));
        % sometimes there is no water in the whole domain
        if (epsP == 0.)
            epsP = epsilon;
        end
        if (epsS == 0.)
            epsS = epsilon;
        end
        % bad marks the cells prolematic in evaluating Jacobian
        badP = abs(diag(surfactant.jac{4})) < epsP;
        badS = abs(diag(surfactant.jac{5})) < epsS;
        % the other way is to choose based on the water saturation
        polymer(badP) = cp(badP);
        surfactant(badS) = cs(badS);
    end

    eqs   = {water, oil, gas, polymer, surfactant};
    divTerms = {divWater, divOil, divGas, divPolymer, divSurfactant};
    names = {'water', 'oil', 'gas', 'polymer', 'surfactant'};
    types = {'cell', 'cell', 'cell', 'cell', 'cell'};
    components = {cp, cs};

    % Add in any fluxes / source terms prescribed as boundary conditions.
    dissolved = model.getDissolutionMatrix(rs, rv);
    [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                    pressures, sat, mob, rho, ...
                                                    dissolved, components, ...
                                                    drivingForces);
    
    % Finally, add in and setup well equations
    wc    = vertcat(W.cells);
    perf2well = getPerforationToWellMapping(W);
    sgn = vertcat(W.sign);

    wc_inj = wc(sgn(perf2well) > 0);
    cpw     = cp(wc_inj);

    % remove the old viscosity and applying the fully mixed viscosity
    %     muWMultW = muWMult(wc_inj);
    muWMultW = muWeffMult(wc_inj);
    
    muWFullyMixed = model.fluid.muWMult(cpw);

    mob{1}(wc_inj) = mob{1}(wc_inj) ./ muWFullyMixed .* muWMultW;

    % Finally, add in and setup well equations
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                    names, types, wellSol0, ...
                                                    wellSol, ...
                                                    wellVars, wellMap, ...
                                                    p, mob, rho, dissolved, ...
                                                    components, dt, opt);
    % Finally, adding divergence terms to equations
    for i = 1:numel(divTerms)
        eqs{i} = eqs{i} + divTerms{i};
    end
    
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);




end
