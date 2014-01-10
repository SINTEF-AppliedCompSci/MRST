function plotWellRates(W, data, varargin)
% Simple utility for plotting well rates. Used by optimizeTOF.
    if size(data, 2) ~= numel(W)
        data = data .';
    end
    if size(data, 1) < 2
        return
    end
    assert(size(data, 2) == numel(W), 'Input data and well number must match!');
    subs = repmat(1:size(data, 1), 2, 1);
    subs = subs(:);

    t = 1:size(data, 1);
    % Normalize
    data = bsxfun(@rdivide, data, sum(data(end, :), 2));

    area(t(subs(2:end-1)), data(subs(1:end-2), :));
    axis tight
    l = arrayfun(@(x) x.name, W, 'UniformOutput', false);
    caxis([.5 size(data, 2) - .5]);
    if nargin > 2 && any(varargin{1})
        % We were given indices of steps where improvement to the objective
        % function happened, plot them as well!
        hold on
        improvement = varargin{1};
        if islogical(varargin{1})
            improvement = find(improvement);
        end
        ni = numel(improvement);
        plot(improvement, zeros(ni, 1), 'd', 'markerfacecolor', 'g', 'markersize', 12)

        l = [l; 'Improvement'];
    end
    legend(l);
end