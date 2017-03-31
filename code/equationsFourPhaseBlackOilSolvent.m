function equationsFourPhaseBlackOilSolvent(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;
G = model.G;
W = drivingForces.W;

% Currently we do not support senario without wells.
assert(isempty(drivingForces.bc) && isempty(drivingForces.src));

% Properties at current timestep
[p, sW, sG, rs, rv, c, cmax, wellSol] = model.getProps(state, ...
   'pressure', 'water', 'gas', 'solvent', 'rs', 'rv', 'wellsol');

% Properties at previous timestep
[p0, sW0, sG0, rs0, rv0, c0, cmax0] = model.getProps(state0, ...
   'pressure', 'water', 'gas', 'solvent', 'rs', 'rv');

bhp    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qGs    = vertcat(wellSol.qGs);
qSs    = vertcat(wellSol.qSs);
qWPoly = vertcat(wellSol.qWPoly);

%Initialization of primary variables ----------------------------------
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

if ~opt.resOnly
    if ~opt.reverseMode
        % define primary varible x and initialize
        [p, sW, x, c, qWs, qOs, qGs, qWPoly, bhp] = ...
            initVariablesADI(p, sW, x, c, qWs, qOs, qGs, qWPoly, bhp);
    else
        x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
        % Set initial gradient to zero
        zw = zeros(size(bhp));
        [p0, sW0, x0, c0, zw, zw, zw, zw, zw] = ...
            initVariablesADI(p0, sW0, x0, c0, zw, zw, zw, zw, zw); %#ok
        clear zw;
        [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1-sW, x0, rs0, rv0, p0);
    end
end
