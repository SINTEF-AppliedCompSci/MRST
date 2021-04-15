function obj = getObjectiveDiagnostics(G, rock, type, varargin)
%Get an objective function for the diagnostics optimization routines
%
% SYNOPSIS:
%  obj = getObjectiveDiagnostics(G, rock, 'minlorenz')
%
% DESCRIPTION:
% This function produces a function handle that is a realization of some
% objective function based on quantities used in flow diagnostics, that is:
%    - Pressure
%    - Reservoir fluxes
%    - Well fluxes and bottom hole pressure
%    - Forward / backward time of flight
%    - Well tracers (production / injection are treated seperately)
%
%   The returned function handle will accept input in the form
%
%      value = obj(state, D)
%
%   where state is the reservoir state and D the flow diagnostics object
%   containing time of flight and tracer values. The 'D' struct can be
%   returned from solveStationaryPressure or computeTOFandTracer and the
%   state may come from any pressure solver contained in MRST. Most of the
%   optimization routines get both these quantities from
%   solveStationaryPressure, which is the recommended way of evaluating the
%   objective functions if gradients are required.
%
%   The value return will be an instance of the ADI class and contains both
%   the objective function value as well as gradients with respect to time
%   of flight etc. The treatment is consistent with the ordering used in
%   solveStationaryPressure.
%
%
% REQUIRED PARAMETERS:
%   G    - The grid which will be used to evaluate the objective function
%
%   rock - Rock with valid <poro> field containing the porosity. Several of
%          the included objective functions use the pore volumes.
%
%   type - Either:
%
%          A string indicating the type of objective functions. Currently the
%          available functions include,
%          'minlorenz' - Minimize the Lorenz' coefficient,
%
%          or
%
%          A function handle which accepts a single struct argument
%          containing the fields
%             'pressure':    Reservoir pressure
%             'state':       Full reservoir state
%             'pv' :         Pore volume
%             'fluxes':      Well fluxes
%             'bhp':         Well bottom hole pressures
%             'forward_tof': Time of flight originating from injectors
%             'backward_tof':Time of flight originating from producers
%
%         This option is primarily intended for rapid experimentation, i.e.
%         allowing someone who knows what they are doing to define
%         objective functions on the fly in the form
%
%         tmp = @(packed) sqrt(sum(packed.forward_tof.^2))
%
%         obj = getObjectiveDiagnostics(G, rock, tmp)
%
%         which would make obj an objective function that is now the norm
%         of the forward time of flight, and can be used to optimize well
%         rates and placement.
%
% OPTIONAL PARAMETERS:
%
%   Any optional parameters are passed onto the objective function.
%   Different objective functions handle different parameters.
%
% RETURNS:
%   obj - Function handle accepting state and D struct.
%


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


    if isa(type, 'function_handle')
        obj = @(state, D) genericFunction(state, D, type, poreVolume(G, rock));
    else
        switch type
            case 'mintof'
                obj = @(state, D) minimizeTOF(state, D, varargin{:});
            case 'minpvdiff'
                obj = @(state, D) minimizeRegionVolume(state, D, poreVolume(G, rock), varargin{:});
            case 'minlorenz'
                obj = @(state, D) minimizeLorenz(state, D, poreVolume(G, rock), varargin{:});
            case 'minlorenzlayers'
                obj = @(state, D) minimizeLayeredLorenz(G, state, D, poreVolume(G, rock), varargin{:});
            case 'mintargettof'
                obj = @(state, D) minTargetTOF(state, D, varargin{:});
            case 'oilvoldisc'
                obj = @(state, D) phaseDiscount(state, D, varargin{:});
            case 'maximizeNPV'
                obj = @(state, D) maximizeNPV(state, D, varargin{:});
            otherwise
                error('I do not know that objective function...')
        end
    end
end

function obj = genericFunction(state, D, func, pv)
    [p, fluxes, bhp, forward_tof, forward_tracer, backward_tof, backward_tracer] = getValues(state, D);
    packed = struct(  'pressure', p,...
                        'state', state, ...
                        'pv', pv, ...
                        'fluxes', fluxes, ...
                        'bhp', bhp, ...
                        'forward_tof', forward_tof, ...
                        'backward_tof', backward_tof);

    packed.forward_tracer = forward_tracer;
    packed.backward_tracer = backward_tracer;
    obj = func(packed);
end

function obj = minimizeTOF(state, D, varargin)
    [p, fluxes, bhp, forward_tof, forward_tracer, backward_tof, backward_tracer] = getValues(state, D);
    obj = sum(forward_tof + backward_tof);
end

function obj = minimizeRegionVolume(state, D, pv, varargin)
    [p, fluxes, bhp, forward_tof, forward_tracer, backward_tof, backward_tracer] = getValues(state, D);


    avgVol = sum(pv)./numel(forward_tracer);
    % Idea: Try to let each injector sweep the same total volume for a
    % balanced reservoir sweep optimization loop
    for i = 1:numel(forward_tracer)
        if i == 1
            obj = abs(sum(forward_tracer{i}.*pv) - avgVol);
        else
            obj = obj + abs(sum(forward_tracer{i}.*pv) - avgVol);
        end
    end

    avgVol = sum(pv)./numel(backward_tracer);
    for i = 1:numel(backward_tracer)
        obj = obj + abs(sum(backward_tracer{i}.*pv) - avgVol);
    end
end

function obj = minTargetTOF(state, D, targets)
    [p, fluxes, bhp, forward_tof, forward_tracer, backward_tof, backward_tracer] = getValues(state, D);
    tmp = zeros(numel(forward_tof.val), 1);
    tmp(targets) = 1;



    obj = sum(forward_tof.*tmp);
end

function obj = phaseDiscount(state, D, pv, r, tscale, varargin)
    [p, fluxes, bhp, forward_tof, forward_tracer, backward_tof, backward_tracer] = getValues(state, D);
    if nargin == 4
        tscale = year;
    end
    if nargin > 5
        maxTof = varargin{1};
        forward_tof(forward_tof.val > maxTof) = maxTof;
        backward_tof(backward_tof.val > maxTof) = maxTof;
    end

%     obj = state.s(:,1).*(1-r).^(forward_tof + backward_tof);
    oil = state.s(:,2);
    wat = state.s(:,1);
    if numel(r) == 1;
        r = [r; r];
    end

    t = forward_tof + backward_tof;
%     t = backward_tof;

    if 1
        a = 1; b = 1/sum(pv);
    else
        a = sum(pv); b = 1;
    end
    if 1
%         tscale = year;
%         tscale = 0.001*year;
%         tscale = max(t.val);
%         sum(pv);
%         tscale = sum(pv)*day;
%         tscale = (log10(0.01)/log10(1+r(1)))/tscale;
%         tscale
        if 1
            tscale = 1/year;
            t = t*tscale;
        else
            t = log(t);
        end
%         f = gcf;
%         figure(4); hist(log10(t.val));
%         figure(f)
        discount_oil = exp(-t.*log(1+r(1)));
        discount_wat = exp(-t.*log(1+r(2)));
%         obj = a - sum(pv.*oil.*discount)*b;
        obj = - (sum(pv.*oil.*discount_oil) + sum(pv.*wat.*discount_wat))*b;
    end

    f = gcf;
    figure(3); clf; plot(sort(discount_oil.val, 'descend'))
    figure(f)
end

function obj = minimizeLorenz(state, D, pv, varargin)
    nc = size(D.tof, 1);
    opt = struct('removeOutliers', false, ...
                'oilModifier', false, ...
                'perPartition', false, ...
                'computeSweep', false, ...
                'mask', true(nc, 1));

    opt = merge_options(opt, varargin{:});

    cntr = 0;
    if opt.perPartition
        assert(~(isempty(D.ipart) || isempty(D.ppart)), ...
            'For the per partition option to work, tracers must be computed!')
        D.inj = [1];
        Ni = numel(D.inj); Np = numel(D.prod);
        for ii = 1:Ni
            for ip = 1:Np
                mask = opt.mask & D.ipart == ii & D.ppart == ip;
                if ~any(mask); continue; end

                tmp = minimizeLorenz(state, D, pv, ...
                                    'oilModifier', opt.oilModifier, ...
                                    'removeOutliers', opt.removeOutliers, ...
                                    'mask'       , mask);
                if cntr == 0
                    obj = tmp;
                else
                    obj = obj + tmp;
                end
                cntr = cntr + 1;
            end
        end
        obj = obj/cntr;
        return
    end

    mask = opt.mask;


    [p, fluxes, bhp, forward_tof, forward_tracer, backward_tof, backward_tracer] = getValues(state, D);

    t = forward_tof + backward_tof;

    if opt.removeOutliers
        vals = t.val;
        vals = log10(vals);

        mv = mean(vals);
        stdv = std(vals);
        inclusion = vals < mv + 1*stdv;

        mask = mask & inclusion;
        fprintf('Excluding %d of %d\n', sum(~inclusion), numel(vals));
    end

    t(~mask) = inf;
    pv = pv.*(mask);
    if opt.oilModifier
        % Scale pore volume by oil to increase sweep of *oil* specifically
        pv = pv.*state.s(:,2);
    end

    [~,order] = sort(t.val);    % sort cells base on travel time

    [~,invorder] = sort(order);

    %Write Lc as
    % Lc ~ 2*[(w'*pvt)/(1'*pvt)]-1
    % where pvt = pv./t and
    % w = [sumcum(dPhi)](invorder)
    % sumcum = flipud*cumsum*flipud
    pvs  = pv(order);
    dPhi = pvs/sum(pvs);
    %trapezoidal rule:
    dPhi(1:end-1) = dPhi(1:end-1) + dPhi(2:end);
    dPhi = .5*dPhi;

    w = flipud(cumsum(flipud(dPhi)));
    w = w(invorder);

    pvt = pv./t;

    obj = 2*((w'*pvt)./(ones(size(w))'*pvt)) - 1;

end

function obj = minimizeLayeredLorenz(G, state, D, pv, removeOutliers)
    ijk = gridLogicalIndices(G);
    for i = 1:G.cartDims(3)
        lc = minimizeLorenz(state, D, pv, removeOutliers, ijk{3} == i);
        lc = lc./G.cartDims(3);
        if i == 1
            obj = lc;
        else
            obj = obj + lc;
        end
    end

end

function obj = maximizeNPV(state, D, varargin)
    % diagnistics NPV:
    %  * values   : value (in dollars) for each grid-cell (or constant value)
    %  * discount : discount factor
    %  * rateCost : cost (in dollars) for injecting one unit volume for
    %               each well
    opt = struct('values',            1, ...
                 'discount',        0.1, ...
                 'rateCost',         [], ...
                 'endTime',      1*year, ...
                 'cutoffSmoothing', 0.1);
    opt = merge_options(opt, varargin{:});
    
    ws    = state.wellSol;
    nw    = numel(ws);
    nperf = arrayfun(@(x)numel(x.flux), ws);
    isInj = arrayfun(@(x)(sum(x.flux) > 0), ws);
    if ~isempty(opt.rateCost)
        rateCost = opt.rateCost;
    else
        rateCost = zeros(nw,1);
    end
    
    fac = rldecode(isInj(:).*rateCost(:), nperf(:));
    
    [p, fluxes, bhp, forward_tof, forward_tracer, backward_tof, backward_tracer] = getValues(state, D);
    tau    = backward_tof;
    tauMax = opt.endTime;
    
    values = opt.values;
    if numel(values) == 1
        values = ones(numel(state.pressure), 1)*values;
    end
    values = values(:);
    % use a smooth cut-off of at endTime
    if opt.cutoffSmoothing > 0
        dTau  = tauMax-tau;
        epsi  = opt.cutoffSmoothing*tauMax;
        cutoff = .5*dTau.*(dTau.^2+epsi^2).^(-.5)+.5;
    else
        cutoff = ones(size(values));
    end
    
    d = opt.discount; 
    costFac = (1-(1+d)^(-opt.endTime))/log(1+d);
    obj = values' * (cutoff.*(1+opt.discount).^(-tau)) - (fac'*fluxes)*costFac;
end

function c = expandMatrix(data)
    N = size(data, 2);
    c = cell(1, N);
    for i = 1:N
        c{i} = data(:, i);
    end
end

function [p, fluxes, bhp, forward_tof, forward_tracer, backward_tof, backward_tracer] = getValues(state, D)

    p = state.pressure;

    fluxes = vertcat(state.wellSol.flux);
    bhp = vertcat(state.wellSol.pressure);

    forward_tof = D.tof(:,1);
    backward_tof = D.tof(:,2);

    forward_tracer = expandMatrix(D.itracer);
    backward_tracer = expandMatrix(D.ptracer);

    if isempty(D.itracer) || isempty(D.ptracer)
        % No tracer solution provided, skip
        [p, fluxes, bhp, forward_tof, backward_tof] = initVariablesADI(...
        p, fluxes, bhp, forward_tof, backward_tof);
    else
        [p, fluxes, bhp, forward_tof, forward_tracer{:}, backward_tof, backward_tracer{:}] = initVariablesADI(...
        p, fluxes, bhp, forward_tof, forward_tracer{:}, backward_tof, backward_tracer{:});
    end

end

function v = adiDiagCumsum(v, order, invorder)
% Cumsum for large adi objects is impossible, but the structure is here is
% that of a (permuted) diagonal, so we exploit this
    v.val = cumsum(v.val);
    for i = 1:numel(v.jac)
        J = v.jac{i};
        if ~(size(J, 1) == size(J, 2))
            assert(nnz(J) == 0);
            continue
        end
        tmp = J(:, invorder).*sparse(1:numel(order), 1:numel(order), order);
        v.jac{i} = tmp(:, order);
    end
end
