function [description, options, state0, model, schedule, plotOptions] = egg_wo(varargin)
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.
    % Description
    description = ['Egg model, see Jansen, J. D., et al. '       , ...
                   '"The egg model–a geological ensemble for '   , ...
                   'reservoir simulation." '                     , ...
                   'Geoscience Data Journal 1.2 (2014): 192-195.'];
    % Optional input arguments
    options = struct('realization', 0); % Realization [0, 100]
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil
    % Get deck
    deck = getDeckEGG('realization', options.realization);
    % Initialize MRST problem from deck
    [state0, model, schedule] = initEclipseProblemAD(deck);
    % Plotting
    plotOptions = {'View'              , [-30, 20] , ...
                   'PlotBoxAspectRatio', [1,1,0.2] , ...
                   'Box'               , true      , ...
                   'Size'              , [800, 500]};
end