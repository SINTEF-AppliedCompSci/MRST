function [state, model, schedule]  = setupSPE10_AD(varargin)
%Undocumented Utility Function

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

    opt = struct('layers', 1:85, ...
                 'dt',      30*day, ...
                 'T',       2000*day, ...
                 'rampup',  8, ...
                 'I',       1:60, ...
                 'J',       1:220, ...
                 'make2D',  false, ...
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
    fluid.bW = @(p, varargin) bW0*exp((p - pRef)*3e-6/psia);
    fluid.muW = @(p, varargin) 0*p + 0.3*centi*poise;
    
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

    
    rock = getSPE10rock(opt.I, opt.J, opt.layers);

    % Compute pore volume fraction of the full model
    volFrac = sum(rock.poro)/1.9141e+05;
    rock.poro(rock.poro < opt.minporo) = opt.minporo;
    
    
    % Grid
    nl = numel(opt.layers);
    i = numel(opt.I);
    j = numel(opt.J);
    assert(i > 0 && i <= 60);
    assert(j > 0 && j <= 220);

    cartDims = [ i, j, nl ];
    physDims = cartDims .* [ 20, 10, 2 ]*ft;

    if opt.make2D
        assert(nl == 1);
        cartDims = cartDims(1:2);
        physDims = physDims(1:2);
        rock.perm = rock.perm(:, 1);
    end

    G = cartGrid(cartDims, physDims);
    try
        mrstModule add libgeometry
        G = mcomputeGeometry(G);
    catch
        G = computeGeometry(G);
    end
    model = GenericBlackOilModel(G, rock, fluid, 'gravity', [0, 0, 9.80665],...
                                    'water', true, 'oil', true, 'gas', false);
    model.minimumPressure = 1000*psia;
    
    state = initResSol(G, pRef, [srw, 1-srw]);
    
    % Wells
    makeProd = @(W, name, I, J) verticalWell(W, G, rock, I, J, [],...
        'Name', name, 'radius', 5*inch, 'sign', -1, 'Type', 'bhp',...
        'Val', 4000*psia, 'comp_i', [.5, .5]);
    I = G.cartDims(1);
    J = G.cartDims(2);
    W = [];
    W = makeProd(W, 'P1', 1, 1);
    W = makeProd(W, 'P2', I, 1);
    W = makeProd(W, 'P3', I, J);
    W = makeProd(W, 'P4', 1, J);
    W = verticalWell(W, G, rock, ceil(I/2), ceil(J/2), [], 'Name', 'I1', 'radius', 5*inch, ...
        'Type', 'rate', 'Val', volFrac*5000*stb/day, 'comp_i', [1, 0], 'Sign', 1);
    
    dt = rampupTimesteps(opt.T, opt.dt, opt.rampup);
    
    schedule = simpleSchedule(dt, 'W', W);
end
