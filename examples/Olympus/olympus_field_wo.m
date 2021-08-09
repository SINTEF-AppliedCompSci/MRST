function [description, options, state0, model, schedule, plotOptions] = olympus_field_wo(varargin)


    description = 'Olympus field model from the ISAPP challenge';
    options = struct('refine'     , [], ...
                     'realization', 1);
    options = merge_options(options, varargin{:});

    if nargout <= 2, return; end
    
    gravity reset on
    
    fldr     = fullfile(mrstDataDirectory, 'olympus', 'ECLIPSE');
    caseName = ['OLYMPUS_', num2str(options.realization)];
    deck     = fullfile(fldr, caseName, [caseName, '.DATA']);
    
    [state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'useMex', true);

    % Plotting
    plotOptions = {'View'              , [10, 50]  , ...
                   'Size'              , [900, 450], ...
                   'PlotBoxAspectRatio', [33.7444   16.8243    5.0000], ...
                   'Box'               , true};
    
end