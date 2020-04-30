function [data,tof,tof_ix] = getAllPhaseTOFDistributions(state, W, pv, prod, D, varargin)

    opt = struct('maxTOF',   [],...
        'normalize',    'true');
    opt = merge_options(opt, varargin{:});
    if isempty(opt.maxTOF)
        opt.maxTOF = inf;
    end
    
    
    % Set default maxTOF to 1.25*minimal arrival time: 
%     opt.maxTOF = min(opt.maxTOF, 1.25*min( D.tof( W(D.prod(prod)).cells, 1 )) );
%     opt.maxTOF = min(opt.maxTOF, 5*min( D.tof( W(D.prod(prod)).cells, 1 )) );
    tof = D.tof(:,2)./year;
    
    
    opt.maxTOF = 500*year;
    % Don't include maxTOF (avoid jump at maxTOF in plot)
    tof_ix = find(and(D.ptracer(:,prod)>.01,  tof < .99*opt.maxTOF/year ));
    
    % Sort tof
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
    

    % Cumsum data
    data = cumsum(data);
    
    
    % use approx 50 points in plot
    di   = ceil(numel(tof)/50);
    tof    = tof(1:di:end);
    tof_ix    = tof_ix(1:di:end);
    data = data(1:di:end, :);
    if opt.normalize
       data = bsxfun(@rdivide, data, sum(data, 2));
    end
    data(isnan(data)) = 0;
end

