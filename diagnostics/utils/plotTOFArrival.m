function plotTOFArrival(state, W, pv, fluid, prod, D, useMob)
    tof = D.tof(:,2)./year;
    localtof = tof(D.ppart == prod);

    nbins = 250;
    nPh = size(state.s, 2);


    % Compute fluid mobilities
    if useMob
        % Use mobility as data based on fluid.
        if nPh == 3
            [mu rho] = fluid.properties(state);
        else
            [mu rho] = fluid.properties();
        end
        data = fluid.relperm(state.s)./repmat(mu, size(state.s, 1), 1);
    else
        data = state.s;
    end
    data = data.*repmat(pv, 1, size(state.s, 2));

    % Bin selection
    ltof = log10(localtof);
    tof_sub = localtof(ltof < 3*std(ltof));
    [N, X] = hist(tof_sub, nbins);

    h = X(2) - X(1);
    s = zeros(nbins, nPh);
    for i = 1:nbins
        cind = abs(tof - X(i)) < h/2 & D.ppart == prod;
        s(i, :) = sum(data(cind, :), 1);
    end
    clear i
    ah = plotNormalizedArrival(X, s, tof_sub, true);

    if isfield(fluid, 'names') && numel(fluid.names) == nPh
        names = fluid.names;
    else
        names = arrayfun(@(x) ['Phase ', num2str(x)], 1:nPh, 'UniformOutput', false);
    end
    for i = 1:nPh
        localVolume = convertTo(sum(data(:,i).*D.ptracer(:, prod)), stb);
        names{i} = [names{i}, sprintf(' (%2.1G stb)', localVolume)];
    end
    legend(names);

    ch = get(ah, 'Children');
    for i = 1:numel(ch);
        set(ch{i}, 'ButtonDownFcn', @(varargin) plotTOFPhase(i));
        set(ch{i}, 'HitTest', 'on');
    end

    function plotTOFPhase(phase)
        figure;
        nInj = max(D.ipart);
        sphase = zeros(nbins, nInj);
        for xi = 1:nbins
            cind = abs(tof - X(xi)) < h/2 & D.ppart == prod;
            for j = 1:nInj
%                 sphase(xi, j) = sum(data(cind & D.ipart == j, phase), 1);
                sphase(xi, j) = sum(data(cind, phase).*D.itracer(cind, j), 1);
            end
        end
        plotNormalizedArrival(X, sphase, tof_sub, false);
        legend(arrayfun(@(x) x.name, W(D.inj), 'uniformoutput', false),'Location', 'NorthWest')
        title([W(D.prod(prod)).name ': Injector distribution for ' names{phase}]);
    end
end

function h = plotNormalizedArrival(X, data, tof_sub, normalize)
    cs = cumsum(data, 1);
    if normalize
       cs = bsxfun(@rdivide, cs, sum(cs, 2));
    end
    cs(isnan(cs)) = 0;

    hold on
    h = area(X, cs);
    %axis([min(tof_sub), max(tof_sub), 0, 1])
    axis tight
    xlabel('TOF distance in years')
end
