function params = addModelParameter(params, problem, varargin)
p = ModelParameter_V2(problem, varargin{:});
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
        