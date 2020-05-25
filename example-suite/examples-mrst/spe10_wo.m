function [description, state0, model, schedule, options, plotOptions] = spe10_wo(varargin)
    % SPE 10 Model 2
    description = 'SPE10 Model 2, with support for subset of layers';
    if nargout == 1, return; end
    
    options = struct('layers', []);
    require spe10
    gravity reset on
    % Get base model
    [options, spe10Args] = merge_options(options, varargin{:});
    if isempty(options.layers)
        options.layers = 1:85;
    elseif ischar(options.layers)
        switch options.layers
            case 'tarbert'
                options.layers = 1:35;
            case 'upper_ness'
                options.layers = 36:85;
        end
    else
        assert(min(options.layers) >= 1 && max(options.layers) <= 85);
    end
        
    [state0, model, schedule] = setupSPE10_AD(spe10Args{:}, 'layers', options.layers);
    % Plotting
    if model.G.cartDims(3) > 2
        view = [120,26];
        proj = 'perspective';
        size = [800, 400];
    else
        view = [0,90];
        proj = 'orthographic';
        size = [500, 800];
    end
    plotOptions = {'PlotBoxAspectRatio', [1,1.83,0.3], ...
                   'Projection'        , proj        , ...
                   'View'              , view        , ...
                   'Size'              , size        };
end