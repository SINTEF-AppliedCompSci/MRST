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
    Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
            [p, sW, x, cp, cs, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, sW, x, cp, cs, wellVars{:});
        else
            x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
            % Set initial gradient to zero
            wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
            [p0, sW0, x0, cp0, cs0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, sW0, x0, cp0, cs0, wellVars0{:});
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
    state = model.initStateFunctionContainers(state);

    % EQUATIONS ---------------------------------------------------------------

    [b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
    [b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
    [phaseFlux, flags] = model.getProps(state, 'PhaseFlux',  'PhaseUpwindFlag');
    [pressures, mob, rho] = model.getProps(state, 'PhasePressures', 'Mobility', 'Density');

    [bW, bO, bG]       = deal(b{:});
    [bW0, bO0, bG0]    = deal(b0{:});
    [vW, vO, vG]       = deal(phaseFlux{:});
    [upcw, upco, upcg] = deal(flags{:});
    [mobW, mobO, mobG] = deal(mob{:});
    vP = model.getProps(state, 'PolymerPhaseFlux');

    muWeffMult = model.getProp(state, 'EffectiveMixturePolymerViscMultiplier');

    if model.usingShear || model.usingShearLog || model.usingShearLogshrate
        % calculate well perforation rates :
        if ~isempty(W)
            if ~opt.reverseMode
                wc    = vertcat(W.cells);
                perf2well = getPerforationToWellMapping(W);
                sgn = vertcat(W.sign);

                wc_inj = wc(sgn(perf2well) > 0);
                cpw        = cp(wc_inj);

                muWMultW = muWeffMult(wc_inj);
                muWFullyMixed = model.fluid.muWMult(cpw);

                mob{1}(wc_inj) = mob{1}(wc_inj) ./ muWFullyMixed .* muWMultW;

                dissolved = model.getDissolutionMatrix(rs, rv);


                [src, wellsys, state.wellSol] = ...
                model.FacilityModel.getWellContributions(wellSol0, wellSol, wellVars, ...
                                        wellMap, p, mob, rho, dissolved, {cp,cs}, ...
                                        dt, opt.iteration);

            else
                error('not supported yet!');
            end
        else
            error('The model with polymer does not support scenarios without wells now!');
        end

        % s = model.operators;  % The previous s was overwritten with saturations.
        poro =  s.pv./G.cells.volumes;
        poroFace = s.faceAvg(poro);
        faceA = G.faces.areas(s.internalConn);

        % Bw * Fw should be flux
        Vw = vW./(poroFace .* faceA);

        % Using the upstreamed viscosity multiplier due to PLYVISC
        muWMultf = s.faceUpstr(upcw, muWeffMult);
        wc = vertcat(W.cells);
        muWMultW = muWeffMult(wc);

        % We assume the viscosity multiplier should be consistent with current
        % way in handling the injection mobility, while the assumption is not
        % verfied with any tests yet due to lack of the reference result.
        [~, wciPoly, iInxW] = getWellPolymer(W);
        cpw = cp(wc);
        muWMultW(iInxW) = model.fluid.muWMult(cpw(iInxW));

        % Maybe should also apply this for PRODUCTION wells.
        muWMultW((iInxW(wciPoly==0))) = 1;

        % The water flux for the wells.
        cqs = vertcat(state.wellSol.cqs);
        fluxWaterWell = value(cqs(:, 1));
        poroW = poro(wc);

        % the thickness of the well perforations in the cell
        welldir = { W.dir };
        i = cellfun('prodofsize', welldir) == 1;
        welldir(i) = arrayfun(@(w) repmat(w.dir, [ numel(w.cells), 1 ]), ...
                              W(i), 'UniformOutput', false);
        welldir = vertcat(welldir{:});
        [dx, dy, dz] = cellDims(G, wc);
        thicknessWell = dz;
        thicknessWell(welldir == 'Y') = dy(welldir == 'Y');
        thicknessWell(welldir == 'X') = dx(welldir == 'X');

        % For the wells
        % The water velocity is computed at the reprensentative radius rR.
        if ~isfield(W, 'rR')
            error('The representative radius of the well is not initialized');
        end
        rR = vertcat(W.rR);
        VwW = bW(wc).*fluxWaterWell./(poroW .* rR .* thicknessWell * 2 * pi);
        muWMultW = value(muWMultW);
        VwW = value(VwW);
        muWMultf = value(muWMultf);
        Vw = value(Vw);

        if model.usingShearLogshrate
            % calculating the shear rate based on the velocity
            if ~opt.resOnly
                krwF = s.faceUpstr(upcw, krW.val);
                swF = s.faceUpstr(upcw, sW.val);
            end

            if opt.resOnly
                krwF = s.faceUpstr(upcw, krW);
                swF = s.faceUpstr(upcw, sW);
            end

            permF = s.T ./faceA;
            temp = permF.*swF.*krwF;
            index = find(abs(Vw) > 0.);
            Vw(index) = 4.8 * Vw(index).*sqrt(poroFace(index) ./temp(index));

            % calculating the shear rate for the wells
            rW = vertcat(W.r);
            VwW = 4.8 * VwW ./(2*rW);
        end

        if model.usingShear
            shearMultf = computeShearMult(model.fluid, abs(Vw), muWMultf);
            shearMultW = computeShearMult(model.fluid, abs(VwW), muWMultW);
        end

        if model.usingShearLog || model.usingShearLogshrate
            shearMultf = computeShearMultLog(model.fluid, abs(Vw), muWMultf);
            shearMultW = computeShearMultLog(model.fluid, abs(VwW), muWMultW);
        end

        vW = vW ./ shearMultf;
        vP = vP ./ shearMultf;
    end

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
    wflux = bWvW;

    % Conservation of mass for oil
    if model.vapoil
        % The model allows oil to vaporize into the gas phase. The conservation
        % equation for oil must then include the fraction present in the gas
        % phase.
        rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
        % Final equation
        oil = (1/dt).*( pv .*(bO.* sO  + rv.* bG.* sG) - ...
                        pv0.*(bO0.*sO0 + rv0.*bG0.*sG0));
        oflux = bOvO + rvbGvG;
    else
        oil = (1/dt).*(pv.*bO.*sO - pv0.*bO0.*sO0 );
        oflux = bOvO;
    end

    % Conservation of mass for gas
    if model.disgas
        % The gas transported in the oil phase.
        rsbOvO = s.faceUpstr(upco, rs).*bOvO;

        gas = (1/dt).*(pv.* (bG.* sG  + rs.* bO.* sO) - ...
                       pv0.*(bG0.*sG0 + rs0.*bO0.*sO0 ));
        gflux = bGvG + rsbOvO;
    else
        gas = (1/dt).*(pv.*bG.*sG - pv0.*bG0.*sG0 );
        gflux = bGvG;
    end

    % Computation of adsoprtion term
    vSft   = s.faceUpstr(upcw, cs).*vW;
    bWvSft = s.faceUpstr(upcw, bW).*vSft;
    poro = model.rock.poro;
    adsp  = model.getProp(state, 'PolymerAdsorption');
    adsp0 = model.getProp(state0, 'PolymerAdsorption');
    adss  = model.getProp(state, 'SurfactantAdsorption');
    adss0 = model.getProp(state0, 'SurfactantAdsorption');
    adsp_term = fluid.rhoR.*((1-poro)./poro).*(adsp - adsp0);
    adss_term = fluid.rhoRSft.*((1-poro)./poro).*(adss - adss0);

    % Conservation of polymer in polymer phase:
    polymer = ((1-fluid.dps)/dt).*(pv.*bW.*sW.*cp - ...
                                     pv0.*fluid.bW(p0).*sW0.*cp0) + (s.pv/dt).* adsp_term;
    pflux = bWvP;
    
    % Conservation of surfactant in water:
    surfactant = (1/dt).*(pv.*bW.*sW.*cs - pv0.*bW0.*sW0.*cs0) + (s.pv/dt).*adss_term;
    sflux = bWvSft;
    
    % Applying correction to the surfactant and polymer equation when the Jacobian is
    % prolematic for some cells.
    % Typically it is due to totally and almost non-existence of water.    
    if ~opt.resOnly
        epsilon = 1.e-8;
        % the first way is based on the diagonal values of the resulting
        % Jacobian matrix
        [dp, ds] = getDiags(polymer, surfactant);
        epsP = sqrt(epsilon)*mean(abs(dp));
        epsS = sqrt(epsilon)*mean(abs(ds));
        % sometimes there is no water in the whole domain
        if (epsP == 0.)
            epsP = epsilon;
        end
        if (epsS == 0.)
            epsS = epsilon;
        end
        % bad marks the cells prolematic in evaluating Jacobian
        badP = abs(dp) < epsP;
        badS = abs(ds) < epsS;
        % the other way is to choose based on the water saturation
        polymer(badP) = cp(badP);
        surfactant(badS) = cs(badS);
    end

    eqs   = {water, oil, gas, polymer, surfactant};
    fluxes = {wflux, oflux, gflux, pflux, sflux};
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
    muWMultW = muWeffMult(wc_inj);
    
    muWFullyMixed = model.fluid.muWMult(cpw);

    mob{1}(wc_inj) = mob{1}(wc_inj) ./ muWFullyMixed .* muWMultW;

    if model.usingShear || model.usingShearLog || model.usingShearLogshrate
        % applying the shear effects
        mob{1}(wc) = mob{1}(wc)./shearMultW;
    end
    
    % Finally, add in and setup well equations
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
                                                    names, types, wellSol0, ...
                                                    wellSol, ...
                                                    wellVars, wellMap, ...
                                                    p, mob, rho, dissolved, ...
                                                    components, dt, opt);
    % Finally, adding divergence terms to equations
    for i = 1:numel(fluxes)
        eqs{i} = s.AccDiv(eqs{i}, fluxes{i});
    end
    
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

%--------------------------------------------------------------------------

function [wPoly, wciPoly, iInxW] = getWellPolymer(W)
    if isempty(W)
        wPoly = [];
        wciPoly = [];
        iInxW = [];
        return
    end
    inj   = vertcat(W.sign)==1;
    polInj = cellfun(@(x)~isempty(x), {W(inj).cp});
    wPoly = zeros(nnz(inj), 1);
    W_inj = W(inj);
    wPoly(polInj) = vertcat(W_inj(polInj).cp);
    wciPoly = rldecode(wPoly, cellfun(@numel, {W_inj.cells}));

    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw    = numel(W);
    perf2well = rldecode((1:nw)', nPerf);
    compi = vertcat(W.compi);
    iInx  = rldecode(inj, nPerf);
    iInx  = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end

%--------------------------------------------------------------------------

function [dx, dy, dz] = cellDims(G, ix)
% cellDims -- Compute physical dimensions of all cells in single well
%
% SYNOPSIS:
%   [dx, dy, dz] = cellDims(G, ix)
%
% PARAMETERS:
%   G  - Grid data structure.
%   ix - Cells for which to compute the physical dimensions
%
% RETURNS:
%   dx, dy, dz -- [dx(k) dy(k)] is bounding box in xy-plane, while dz(k) =
%                 V(k)/dx(k)*dy(k)

    n = numel(ix);
    [dx, dy, dz] = deal(zeros([n, 1]));

    ixc = G.cells.facePos;
    ixf = G.faces.nodePos;

    for k = 1 : n
       c = ix(k);                                     % Current cell
       f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1); % Faces on cell
       e = mcolon(ixf(f), ixf(f + 1) - 1);            % Edges on cell

       nodes  = unique(G.faces.nodes(e, 1));          % Unique nodes...
       coords = G.nodes.coords(nodes,:);            % ... and coordinates

       % Compute bounding box
       m = min(coords);
       M = max(coords);

       % Size of bounding box
       dx(k) = M(1) - m(1);
       if size(G.nodes.coords, 2) > 1
          dy(k) = M(2) - m(2);
       else
          dy(k) = 1;
       end

       if size(G.nodes.coords, 2) > 2
          dz(k) = G.cells.volumes(ix(k))/(dx(k)*dy(k));
       else
          dz(k) = 0;
       end
    end
end

%--------------------------------------------------------------------------

% Computer the Shear Mutliplier by solving EQ 52.12 in TD
% Vw should be the absolute value?
function v = computeShearMult(fluid, Vw, muWMultf)
    f = fluid;
    % The solution of the shear multipler will be performed through an
    % iterative non-linear solution of the EQ. 52.12 in TD.
    % V ( 1+(P-1)M(V) ) / P = Vw;
    % P is the muWmultf, which is from PLYVISC

    % give the initial guess of the Vsh
    Vsh = Vw;

    Vsh = initVariablesADI(Vsh);

    plyshearMult = f.plyshearMult;

    shFunc = @(x) x.*(1+(muWMultf-1.).*plyshearMult(x))-muWMultf.*Vw;
    eqs = shFunc(Vsh);

    resnorm = norm(value(eqs), 'inf');
    iter = 0;
    maxit = 30;
    abstol = 1.e-15;

    while (resnorm > abstol) && (iter <= maxit)

      J = eqs.jac{1};
      dVsh = -(J \ eqs.val);
      Vsh.val = Vsh.val + dVsh;

      eqs = shFunc(Vsh);
      resnorm = norm(value(eqs), 'inf');

      iter = iter + 1;

    end

    if (iter >= maxit) && (resnorm > abstol)
        error('Convergence failure within %d iterations\nFinal residual = %.8e', maxit, resnorm);
    end

    if(resnorm <= abstol)
        M = plyshearMult(Vsh.val);
        v = (1 + (muWMultf - 1.).* M) ./ muWMultf;
    end

end

%--------------------------------------------------------------------------

function zSh = computeShearMultLog(fluid, vW, muWMultf)
    % The current version handles all the values one by one
    % It needs to be improved for better performance.

    refConcentration = fluid.plyshlog.refcondition(1);
    refViscMult = fluid.muWMult(refConcentration);

    plyshlogTable = fluid.plyshlog.data{1};

    % the water velocity/shear rate in the PLYSHLOG table
    waterVel = plyshlogTable(:,1);
    % the viscosity factors in the PLYSHLOG table
    refM = plyshlogTable(:,2);

    % converting the table using the reference condition
    refM = (refViscMult * refM -1)/(refViscMult-1);

    % the minimum velocity/shear rate value specified in the table
    % it should be the first entry value
    minWaterVel = min(waterVel);

    % only calcuate shear factors when polymer exists and velocity/shear rate
    % is big enough
    iShear = find((vW > minWaterVel) & (muWMultf > 1.));

    P = muWMultf(iShear);

    % convert the velocity/shear rate to their logarithms
    V = log(waterVel);
    vW0 = log(vW(iShear));

    % function to decide the relative location of a point to a line
    f = @(x,y, x0) x+y -x0;

    % the shear factors remain to be calculated
    z = ones(numel(iShear),1);

    for i = 1:numel(iShear)

       % for each value from P, we need to generate a table to calculate
       % the shear factor
       Z = (1 + (P(i) - 1) * refM)/ P(i);
       % covert the shear factors to their logarithms
       Z = log(Z);

       % to flag the relative location of the line and line segments
       % to help to determine the intersection point faster
       sign = f(V, Z, vW0(i));

       % to check if there is sign change
       temp = sign(1:end-1).*sign(2:end);

       % to find the index that sign changed
       j = find(temp <= 0);

       % make sure one and only one intersection point is found.
       assert(numel(j) <= 1);

       if (numel(j) == 1)
           [~,z(i)] = findIntersection([V(j), Z(j); V(j+1), Z(j+1)], [0, vW0(i); vW0(i), 0]);
       end

       if (numel(j) == 0)
           % out of the table range
           % since we handled the small value already, this must be a big
           % value.
           assert(vW0(i) - Z(end) > V(end));
           z(i) = Z(end);
       end

    end

    z = exp(z);

    % the shear factors
    zSh = ones(numel(vW), 1);
    zSh(iShear) = z;

end

function [dp, ds] = getDiags(polymer, surfactant)
    if isa(polymer, 'GenericAD') && isa(polymer.jac{1}, 'DiagonalJacobian');
        dp = polymer.jac{1}.diagonal(:, 4);
        ds = surfactant.jac{1}.diagonal(:, 5);
    else
        dp = diag(polymer.jac{4});
        ds = diag(surfactant.jac{5});
    end
end

%--------------------------------------------------------------------------

% finding the intersection of one line segment l1 and one straight line l2
% l1 is defined with the beginning and ending points.
% l2 is defined with two points along the line.
function [x, y] = findIntersection(l1, l2)

    assert(all(size(l1) == [2 2]));
    assert(all(size(l2) == [2 2]));

    x1 = l1(1,1);
    y1 = l1(1,2);
    x2 = l1(2,1);
    y2 = l1(2,2);

    x3 = l2(1,1);
    y3 = l2(1,2);
    x4 = l2(2,1);
    y4 = l2(2,2);

    d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);

    assert(d ~= 0);

    x = ((x3-x4)*(x1*y2-y1*x2)-(x1-x2)*(x3*y4-y3*x4))/d;
    y = ((y3-y4)*(x1*y2-y1*x2)-(y1-y2)*(x3*y4-y3*x4))/d;

    assert(x >= min(x1,x2) && x <= max(x1,x2));
end
