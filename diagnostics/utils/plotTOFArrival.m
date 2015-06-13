function plotTOFArrival(state, W, pv, fluid, prod, D, varargin)
    opt = struct('maxTOF',   [], ...
             'inj_ix', 1:numel(D.inj));
    opt = merge_options(opt, varargin{:});
    if isempty(opt.maxTOF)
        opt.maxTOF = inf;
    end
    opt.maxTOF = min(opt.maxTOF, 1.25*min( D.tof( W(D.prod(prod)).cells, 1 )) );
    tof = D.tof(:,2)./year;
    tof_ix  = find(and(D.ppart == prod,  tof < .99*opt.maxTOF/year ));
    [tof, ix] = sort(tof(tof_ix));
    tof_ix = tof_ix(ix);
    
    nPh = size(state.s, 2);
    % Compute fluid mobilities
    if isfield(state, 'mob')
        % Use fraction of total mobility to estimate how much will flow
        data = bsxfun(@rdivide, state.mob, sum(state.mob, 2));
    else
        data = state.s;
    end
    % Weight by pore volumes to get actual volumes that will flow
    data = bsxfun(@times, data, pv);
    data = data(tof_ix,:);

    
    
%     % Bin selection
%     tof_sub = localtof(localtof < opt.maxTOF/year);
%     [N, X] = hist(tof_sub, nbins);
% 
%     h = X(2) - X(1);
%     s = zeros(nbins, nPh);
%     for i = 1:nbins
%         cind = abs(tof - X(i)) < h/2 & D.ppart == prod;
%         s(i, :) = sum(data(cind, :), 1);
%     end
%     clear i
    ah  = plotNormalizedArrival(tof, data, true);
    
    if isfield(fluid, 'names') && numel(fluid.names) == nPh
       names = fluid.names;
    else
       names = arrayfun(@(x) ['Phase ', num2str(x)], 1:nPh, 'UniformOutput', false);
    end
    
    %if (~isnan(D.ptracer))
    %    for i = 1:size(data,2)
    %        localVolume = convertTo(sum(data(:,i).*D.ptracer(:, prod)), stb);
    %        names{i} = [names{i}, sprintf(' (%2.1G stb)', localVolume)];
    %    end
    %    legend(names);
    %end

   % ch = get(ah, 'Children');
    for i = 1:numel(ah);
        set(ah(i), 'ButtonDownFcn', @(varargin) plotTOFPhase(i));
        set(ah(i), 'HitTest', 'on');
    end

    function plotTOFPhase(phase)
        inj_ix = opt.inj_ix;
        figure;
        itr = D.itracer(tof_ix, :);
        itr = [itr 1-itr];
        sphase = bsxfun(@times, data(:, phase), itr(:, inj_ix));
        area(tof, cumsum(sphase));
        nms = {W(D.inj).name, 'reservoir'};     
        legend(nms{inj_ix},'Location', 'NorthWest')
        title([W(D.prod(prod)).name ': Injector distribution for ' names{phase}]);
    end
end

function h = plotNormalizedArrival(X, data, normalize)
    phcol = {[.4 .4 1],[.3 0 0], [.5 1 .5]};
    cs = cumsum(data, 1);
    if normalize
       cs = bsxfun(@rdivide, cs, sum(cs, 2));
    end
    cs(isnan(cs)) = 0;

    hold on
    h = area(X, cs);
    for ph = 1:numel(h)
        set(h(ph), 'FaceColor', phcol{ph});
    end
        
    %axis([min(tof_sub), max(tof_sub), 0, 1])
    axis tight
    xlabel('TOF distance in years')
end
