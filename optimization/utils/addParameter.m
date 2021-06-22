function params = addParameter(params, problem, varargin)
% Simple utility function for adding a new parameter to a list of
% parameters
% 
p = ModelParameter(problem, varargin{:});
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
        