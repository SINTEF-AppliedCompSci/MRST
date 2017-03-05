function [state, model, schedule]  = setupSPE10_AD(varargin)
    opt = struct('layers', 1:85, ...
                 'dt',      30*day, ...
                 'T',       2000*day, ...
                 'minporo', 0.01);
    opt = merge_options(opt, varargin{:});
    
    mrstModule add spe10 ad-props ad-blackoil ad-core
    
    srw = 0.2;
    sro = 0.2;
    pRef = 6000*psia;
    
    % Fluid relative permeabilities
    fluid.krW = coreyPhaseRelpermAD(2, srw, 1, srw + sro);
    fluid.krO = coreyPhaseRelpermAD(2, sro, 1, srw + sro);
    
    % Water props
    bW0 = 1./(1.01);
    fluid.bW = @(p) bW0*exp((p - pRef)*3e-6/psia);
    fluid.muW = @(p) 0*p + 0.3*centi*poise;
    
    % Oil props
    p = [300; 800; 8000; 8001]*psia;
    b = 1./[1.05; 1.02; 1.01; 1.01];
    mu = [2.85; 2.99; 3; 3]*centi*poise;
    [fluid.bO, fluid.muO] = tableByPressureLinearAD(p, b, mu);
    
    
    fluid.rhoWS = 64*pound/(ft^3);
    fluid.rhoOS = 53*pound/(ft^3);

    % Rock compressibility
    cR = 1e-6/psia;
    fluid.cR = cR;
    fluid.pvMultR = @(p)(1 + cR.*(p-pRef));

    
    rock = getSPE10rock(opt.layers);

    % Compute pore volume fraction of the full model
    volFrac = sum(rock.poro)/1.9141e+05;
    rock.poro(rock.poro < opt.minporo) = opt.minporo;
    
    
    % Grid
    cartDims = [  60,  220,   numel(opt.layers)];
    physDims = [1200, 2200, 2*cartDims(end)] .* ft();

    G = cartGrid(cartDims, physDims);
    try
        mrstModule add libgeometry
        G = mcomputeGeometry(G);
    catch
        G = computeGeometry(G);
    end
    model = TwoPhaseOilWaterModel(G, rock, fluid, 'gravity', [0, 0, 9.80665]);
    model.minimumPressure = 1000*psia;
    
    state = initResSol(G, pRef, [srw, 1-srw]);
    
    % Wells
    makeProd = @(W, name, I, J) verticalWell(W, G, rock, I, J, [],...
        'Name', name, 'radius', 5*inch, 'sign', -1, 'Type', 'bhp',...
        'Val', 4000*psia, 'comp_i', [.5, .5]);
    W = [];
    W = makeProd(W, 'P1', 1, 1);
    W = makeProd(W, 'P2', 60, 1);
    W = makeProd(W, 'P3', 60, 220);
    W = makeProd(W, 'P4', 1, 220);
    W = verticalWell(W, G, rock, 30, 110, [], 'Name', 'I1', 'radius', 5*inch, ...
        'Type', 'rate', 'Val', volFrac*5000*stb/day, 'comp_i', [1, 0], 'Sign', 1);
    
    dt = rampupTimesteps(opt.T, opt.dt);
    
    schedule = simpleSchedule(dt, 'W', W);
end
