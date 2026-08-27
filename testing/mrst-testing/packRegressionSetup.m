function setup = packRegressionSetup(name, varargin)
%Utility function for packing regression test setup.
%
% SYNOPSIS:
%   setup = packRegressionSetup(name)
%   setup = packRegressionSetup(name, 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   name - name of regression test case as string
%
% OPTIONAL PARAMETERS:
%   description             - One-line test case description
%   options                 - Test case options struct
%   state0, model, schedule - Test case setup
%   plotOptions             - Options for plotting
%   extra                   - Field for anything that does not fit into the
%                             categories above (e.g., reference results)
%
% RETURNS:
%   setup - Struct with all parameters described above

    setup = struct('name'       , name, ...
                   'description', []  , ...
                   'options'    , []  , ...
                   'state0'     , []  , ...
                   'model'      , []  , ...
                   'schedule'   , []  , ...
                   'plotOptions', {{}}, ...
                   'extra'      , []  );
    setup = merge_options(setup, varargin{:});
end
