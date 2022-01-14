function [options, fullSetup, setup] = processTestCaseInput(name, options, description, varargin)
    fullSetup = true;
    if nargin > 3 && islogical(varargin{1})
        fullSetup = varargin{1}; varargin = varargin(2:end);
    end
    [options, extra] = merge_options(options, varargin{:});
    options.extra = extra;
    if fullSetup
        setup = [];
    else
        % Pack setup
        setup = packTestCaseSetup(name,                       ...
                                  'description', description, ...
                                  'options'    , options    );
    end
end