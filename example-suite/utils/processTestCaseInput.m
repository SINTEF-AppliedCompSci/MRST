function [options, optOnly, setup] = processTestCaseInput(name, options, description, varargin)
    optOnly = false;
    if nargin > 3 && islogical(varargin{1})
        optOnly = varargin{1}; varargin = varargin(2:end);
    end
    [options, extra] = merge_options(options, varargin{:});
    options.extra = extra;
    if optOnly
        % Pack setup
        setup = packTestCaseSetup(name,                       ...
                                  'description', description, ...
                                  'options'    , options    );
    else
        setup = [];
    end
end