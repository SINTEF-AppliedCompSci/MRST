function output = mandelrun(params)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}

    
    mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui mpsaw mpfa

    output.params = params;
    
    % Setup grid

    % discretization parameters
    nx = params.nx;
    ny = params.ny;
    dt = params.dt;
    totime = params.totime;
    rampup = params.rampup;
    fixedtsteps = params.fixedtsteps;
    
    % Flow parameters
    perm = params.perm;
    muW  = params.muW;
    poro = params.poro;
    % fluid compressibility
    cW = params.cW;
    % Mechanical parameters
    % Bulk's modulus
    K = params.K;
    % Poisson's ratio
    nu = params.nu;    
    
    resolution = [nx, ny];
    physdim = [1, 1]*meter;
    G = cartGrid(resolution, physdim);
    G = computeGeometry(G);
    output.G = G;
    
    iM = cW; % iM = 1/M (iM = 0 if fluid is incompressible)

    % Second Lamé parameter (also called Shear modulus and denoted G)
    mu = 3/2*((1 - 2*nu)/(1 + nu))*K;
    Gm = mu;

    % First Lamé parameter
    lambda = K - 2/3*Gm;
    Gm = mu;

    % Consolidation coefficient
    cv = perm/muW*(K + 4/3*Gm); % See reference verruijt2013theory (3.13)
    output.cv = cv;
    
    % Biot's coefficient
    alpha = 1;

    % Coussy value for diffusivity coefficient
    % we have checked that we obtain the same values
    docoussy = false;
    if docoussy
        iMKu = K*iM + alpha^2; % page 86 (4.66) with iMKu = Ku/M
        cv = perm/muW*(K + 4/3*mu)/(iMKu + 4/3*mu*iM);
    end

    % force at top
    top_force = 1;

    % Value of normalized pressure (Coussy), which corresponds to pressure at t = 0
    if docoussy
        B = alpha/(iMKu); % page 86 (4.68)
        a = alpha*B*(1 - 2*nu)/3;
        nuu = (a + nu)/(1 - a); % page 87 (4.71)
        pnorm = 1/3*B*(1 + nuu)*top_force; % page 143 (5.168)
    else
        pnorm = top_force/2;
    end
    output.pnorm = pnorm;
    
    % setup rock

    rock.poro  = poro * ones(G.cells.num, 1);
    rock.perm  = (perm/muW) * ones(G.cells.num, 1);
    rock.alpha = alpha * ones(G.cells.num, 1);

    % reference pressure on the side
    pref = 0*barsa;

    % setup mechanics mech structure (with field prop and loadstruct)

    lambda = lambda*ones(G.cells.num, 1);
    mu = mu*ones(G.cells.num, 1);
    mechprop = struct('lambda', lambda, 'mu', mu);

    [tbls, mappings] = setupStandardTables(G);

    % We recover the top, bottom and lateral faces using function pside
    dummy = 0;
    bc = pside([], G, 'Ymax', dummy); 
    topfaces = bc.face;
    bc = pside([], G, 'Xmin', dummy); 
    leftfaces = bc.face;
    bc = pside([], G, 'Xmax', dummy); 
    rightfaces = bc.face;
    bc = pside([], G, 'Ymin', dummy); 
    bottomfaces = bc.face;

    lateralfaces = leftfaces;

    % At the bottom, we have rolling condition in x-direction
    bottomlinform  = repmat([0, 1], numel(bottomfaces), 1);
    % On the lateral walls, we have rolling condition in y-direction
    laterallinform = repmat([1, 0], numel(lateralfaces), 1);
    % At the top we have an unknown y-displacement but which is constant in x-direction. We incorporate it in the Dirichlet
    % condition and add in the equations an extra variable, see mandelEquations
    toplinform = repmat([0, 1], numel(topfaces), 1);

    linform  = [bottomlinform; laterallinform; toplinform];
    extfaces = [bottomfaces; lateralfaces; topfaces];
    bcvals   = zeros(numel(extfaces), 1);

    bc = struct('linform'    , linform , ...
                'extfaces'   , extfaces, ...
                'linformvals', bcvals);

    bc = setupFaceBC(bc, G, tbls);

    loadstruct.bc = bc;

    nodefacecoltbl = tbls.nodefacecoltbl;
    cellcoltbl = tbls.cellcoltbl;
    loadstruct.extforce = zeros(nodefacecoltbl.num, 1); 
    loadstruct.force = zeros(cellcoltbl.num, 1);

    % Setup mech structure 
    mech.prop = mechprop;
    mech.loadstruct = loadstruct;

    % Setup flow parameters (with field c and bcstruct)

    fluid.c = cW;
    fluid.src = [];

    % Setup boundary conditions for flow

    bcfaces = rightfaces;
    bcvals = pref*ones(numel(bcfaces));

    nodefacetbl = tbls.nodefacetbl;

    bcfacetbl.faces = bcfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    bcnodefacetbl = crossIndexArray(nodefacetbl, bcfacetbl, {'faces'});

    map = TensorMap();
    map.fromTbl = bcfacetbl;
    map.toTbl = bcnodefacetbl;
    map.mergefds = {'faces'};
    map = map.setup();

    bcvals = map.eval(bcvals);

    bcdirichlet = struct('bcnodefacetbl', bcnodefacetbl, ...
                         'bcvals', bcvals);
    bcneumann = [];

    bcstruct = struct('bcdirichlet', bcdirichlet, ...
                      'bcneumann'  , bcneumann);

    fluid.bcstruct = bcstruct;

    % Setup Biot model

    model = MandelModel(G, rock, fluid, mech, topfaces);
    model = model.validateModel();

    % Setup schedule

    dt = rampupTimesteps(totime, dt, rampup);
    dtinit = 1; % length of initial phase (when top force is equal to zero). This phase is necessary for proper initialization
    output.dtinit = dtinit;
    tt = cumsum(dt);
    tt = [tt; fixedtsteps];
    tt = uniquetol(tt, 1e-12);
    tt = [dtinit; dtinit + tt];
    t = diff([0; tt]);
    schedule.step.val = t;

    %

    schedule.step.control = 2*ones(numel(schedule.step.val), 1);
    schedule.step.control(1) = 1;
    schedule.control(1) = struct('W', [], 'avgtopforce', 0);
    schedule.control(2) = struct('W', [], 'avgtopforce', top_force);

    % Setup initial state
    clear initState;
    % fluid
    initState.pressure = zeros(G.cells.num, 1);
    nlf = size(bcstruct.bcdirichlet.bcvals, 1);
    initState.lambdafluid = zeros(nlf, 1);
    % mech
    cellcoltbl = tbls.cellcoltbl;
    initState.u = zeros(cellcoltbl.num, 1);
    nlm = size(loadstruct.bc.linformvals, 1);
    initState.lambdamech = zeros(nlm, 1);
    initState.avgtopforce = 0;
    initState.vd = 0;

    solver = NonLinearSolver('maxIterations', 100);
    [~, states] = simulateScheduleAD(initState, model, schedule, 'nonlinearsolver', solver);

    % Capture the solutions at the given fixed time steps.
    tt = cumsum(schedule.step.val);
    [lia, tinds] = ismembertol(fixedtsteps, tt - dtinit);
    assert(all(lia), 'did not find one of the given time steps');
    
    nfts = numel(fixedtsteps);
    pressures = cell(nfts, 1);
    for i = 1 : nfts
        tind = tinds(i);
        pressures{i} = states{tind}.pressure;
    end

    output.pressures = pressures;
end
