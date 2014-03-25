function obj = getObjectiveDiagnosticsStein(G, rock, type, varargin)
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
%     if isempty(varargin) || varargin{1} == 1
%         target = forward_tracer;
%     else
%         target = backward_tracer;
%     end

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

    [~,order]    = sort(t.val);    % sort cells base on travel time
    %tmp       = (1:numel(t.val)).';
    %invorder  = tmp(order);
    [~,invorder] = sort(order);

    %Write Lc as
    % Lc ~ 2*[(w'*pvt)/(1'*pvt)]-.5
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
    1;
%    ts = t(order);
%
%    tmp = (1:numel(t.val)).';
%    invorder = tmp(order);


%     v      = pv(order);  % put volumes in order
%     Phi    = cumsum(v);      % cumulative sum
%     vt     = full(Phi(end)); % total volume of region
%     v      = v./vt;          % normalize to units of pore volumes
%     Phi    = [0; Phi/vt];    % normalize to units of pore volumes
%     flux   = v./ts/vt;       % back out flux based on incompressible flow
%
%     ff  = adiDiagCumsum(flux, order, invorder);
%     for i = 1:numel(ff.jac)
%         % Weighting values by ordering is sort of like a sign change in the
%         % derivative
%         ff.jac{i} = -ff.jac{i};
%     end
%     F      = [0*ff(1); ff./ff.val(end)];     % normalize and store flux
%
%     v      = diff(Phi,1);
%
%     obj = ADI();
%     if opt.computeSweep
%         obj.val =   [0; (Phi(2:end)-Phi(1:end-1))./(F.val(2:end) - F.val(1:end-1))];
%
%
%         obj.jac = {};
%         for i = 1:numel(F.jac)
%             J = F.jac{i};
%             [n m] = size(J);
%
%             J = J(2:end, :) - J(1:end-1, :);
%             [ij,jj,sj] = find(J);
%             obj.jac{i} = sparse(ij, jj, Phi(ij)./sj, n, m);
%         end
%         tD = obj.val;
%         obj = sum(1 - Phi - (1-F).*obj);
%         % Reference for non-ad version
%         %[ev, td] = computeSweep(F.val, Phi);
%     else
%         obj.val = F.val(1:end-1) + F.val(2:end);
%         obj.jac = {};
%         for i = 1:numel(F.jac)
%             J = F.jac{i};
%             obj.jac{i} = J(1:end-1, :) + J(2:end, :);
%         end
%
%         obj = 2*sum(obj.*v./2) - 1;
%     end
%
%     if 0
%         figure(1); clf
%         plot(Phi, F.val);
%         title(num2str(computeLorenz(F.val,Phi)))
%         drawnow
%         figure(3);
%         hist(log10(t.val(order)));
%
%         drawnow
%     end
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
