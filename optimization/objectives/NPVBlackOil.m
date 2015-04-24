function obj = NPVBlackOil(G, wellSols, schedule, varargin)
% Compute net present value of a schedule with well solutions

opt     = struct('OilPrice',             1.0 , ...
                 'GasPrice',             0.1 , ...
                 'GasInjectionCost',     0.1 , ...
                 'WaterProductionCost',  0.1 , ...
                 'WaterInjectionCost',   0.1 , ...
                 'DiscountFactor',       0.0 , ...
                 'ComputePartials',      false, ...
                 'tStep' ,               []);
opt     = merge_options(opt, varargin{:});

ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
riw  = opt.WaterInjectionCost / stb;
rg  = opt.GasPrice   / stb;
rig  = opt.GasInjectionCost   / stb;


d   = opt.DiscountFactor;


% pressure and saturaton vectors just used for place-holding
p  = zeros(G.cells.num, 1);
sW = zeros(G.cells.num, 1);
sG = zeros(G.cells.num, 1);
rs = zeros(G.cells.num, 1);

dts   = schedule.step.val;

tSteps = opt.tStep;
if isempty(tSteps) %do all
    time = 0;
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    time = sum(dts(1:(opt.tStep-1)));
    numSteps = 1;
    dts = dts(opt.tStep);
end

obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
    sol = wellSols{tSteps(step)};
    nW  = numel(sol);
    pBHP = zeros(nW, 1); %place-holder
    qWs  = vertcat(sol.qWs);
    qOs  = vertcat(sol.qOs);
    qGs  = vertcat(sol.qGs);

    if opt.ComputePartials
        [qWs, qWs, qWs, qWs, qWs, qOs, qGs, ignore] = ...
           initVariablesADI(p, sW, sG, rs, qWs, qOs, qGs, pBHP);       %#ok

        clear ignore
    end

    dt = dts(step);
    time = time + dt;

    injInx  = (vertcat(sol.sign) > 0);
    prodInx = ~injInx;
    obj{step} = ( dt*(1+d)^(-time/year) )*...
                spones(ones(1, nW))*( (-ro*prodInx).*qOs +....
                              (-rg*prodInx - rig*injInx).*qGs ...
                             +(rw*prodInx - riw*injInx).*qWs );
end
