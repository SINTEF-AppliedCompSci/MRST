function [options, fullSetup, setup] = processMRSTRegressionSetupInput(name, options, description, varargin)
%Utility function for processing regression test input.
%
% SYNOPSIS:
%   setup = processMRSTRegressionSetupInput(name, options, description)
%   setup = processMRSTRegressionSetupInput(fullSetup, name, options, description)
%   setup = processMRSTRegressionSetupInput(... , 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   name       - Name of test case as string
%   options    - Test case default options struct
%   description - One-line test case description
%
% OPTIONAL PARAMETERS:
%   fullSetup - Flag indicating if we intend to set up the entire test
%               case or not. If false, the third return variable is a test
%               case setup with only the name, description and options.
%               Default value is true.
%   Test case options - The function processes optional parameters for the
%               test case. Parameters that are not part of the options
%               struct are stored as a cell array on the form {'pn1', pv1,
%               ...} to options.extra to facilitate passing these to
%               functions called from inside the test case setup function.
%
% RETURNS:
%   options   - Test case options struct with user-defined options
%   fullSetup  - Flag indicating if we should build the full test case
%   setup     - If fullSetup is false, this is a simplified setup struct
%               with only name, description and user-defined options. If
%               fullSetup is true, this is empty.

    fullSetup = true;
    if nargin > 3 && islogical(varargin{1})
        fullSetup = varargin{1}; varargin = varargin(2:end);
    end

    [options, extra] = merge_options(options, varargin{:});
    options.extra = extra;

    if fullSetup
        setup = [];
    else
        setup = packRegressionSetup(name, ...
                                    'description', description, ...
                                    'options'    , options    );
    end

    if ~isempty(extra)
        pl = ''; if numel(extra(1:2:end)) > 1, pl = 's'; end
        warning(['Option%s %s \n not part of default regression test ', ...
                 'options. Storing to options.extra'], ...
                 pl, sprintf('\n * %s', extra{1 : 2 : end}))
    end
end
