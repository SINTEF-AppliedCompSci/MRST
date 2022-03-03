function out = getVarName(varargin)
    out = {};
    for i = 1:nargin
        out(end+1) = cellstr(inputname(i));
    end
end