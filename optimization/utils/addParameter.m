function params = addParameter(params, SimulatorSetup, varargin)
% Simple utility function for adding a new parameter to a list of
% parameters.
% See  ModelParameter.m for more information about the parameters class.
%
p = ModelParameter(SimulatorSetup, varargin{:});
if isempty(params)
    params = {p};
elseif iscell(params)
    params = [params, {p}];
elseif isa(params, class(p))
    params = [{params}, {p}];
else
    error('Unknown format of input ''params''');
end
end
        