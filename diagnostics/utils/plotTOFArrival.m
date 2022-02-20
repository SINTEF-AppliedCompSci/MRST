function plotTOFArrival(state, W, pv, fluid, prod, D, varargin)
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

    opt = struct('maxTOF',   [], ...
             'inj_ix', 1:numel(D.inj));
    opt = merge_options(opt, varargin{:});
    if isempty(opt.maxTOF)
        opt.maxTOF = inf;
    end
    % Set default maxTOF to 1.25*minimal arrival time: 
    opt.maxTOF = min(opt.maxTOF, 1.25*min( D.tof( W(D.prod(prod)).cells, 1 )) );
    tof = D.tof(:,2)./year;
    % Don't include maxTOF (avoid jump at maxTOF in plot)
    tof_ix = find(and(D.ptracer(:,prod)>.01,  tof < .99*opt.maxTOF/year ));
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
    % Weight by pore volumes*travervalue to get actual volumes 
    data = bsxfun(@times, data, pv.*D.ptracer(:,prod));
    data = data(tof_ix,:);

    %data = cumsum(data);
    
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

    for i = 1:numel(ah);
        set(ah(i), 'ButtonDownFcn', @(varargin) plotTOFPhase(i));
        set(ah(i), 'HitTest', 'on');
    end

    function plotTOFPhase(phase)
        inj_ix = opt.inj_ix;
        figure;
        itr = D.itracer(tof_ix, :);
        itr = [itr 1-sum(itr,2)];
        sphase = bsxfun(@times, data(:, phase), itr(:, inj_ix));
        sphase = cumsum(sphase);
        % use approx 50 points in plot
        di   = ceil(numel(tof)/50);
        area(tof(1:di:end), sphase(1:di:end,:));
        nms = {W(D.inj).name, 'reservoir'};     
        legend(nms{inj_ix},'Location', 'NorthWest')
        title([W(D.prod(prod)).name ': Injector distribution for ' names{phase}]);
    end
end

function h = plotNormalizedArrival(X, data, normalize)
    phcol = {[.4 .4 1],[.3 0 0], [.5 1 .5]};
    data = cumsum(data);
    % use approx 50 points in plot
    di   = ceil(numel(X)/50);
    X    = X(1:di:end);
    data = data(1:di:end, :);
    if normalize
       data = bsxfun(@rdivide, data, sum(data, 2));
    end
    data(isnan(data)) = 0;

    hold on
    h = area(X, data);
    for ph = 1:numel(h)
        set(h(ph), 'FaceColor', phcol{ph});
    end
        
    %axis([min(tof_sub), max(tof_sub), 0, 1])
    axis tight
    xlabel('TOF distance in years')
end
