function [description, options, state0, model, schedule, plotOptions] = egg_wo(varargin)
    % One-line description
    description = '';
    % Optional input arguments
    options = struct('realization', 0);
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil
    % Model
    deck = getDeckEGG('realization', options.realization);
    [state0, model, schedule] = initEclipseProblemAD(deck);
    plotOptions = {'View', [-30, 20], ...
                   'PlotBoxAspectRatio', [1,1,0.2], ...
                   'Box', true, ...
                   'Size', [800, 500]};
end