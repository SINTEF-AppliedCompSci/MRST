function varargout = tableByPressureLinearAD(pRef, varargin)
    varargout = cell(1, nargout);
    for i = 1:numel(varargin)
        v = varargin{i};
        v = reshape(v, [], 1);
        assert(numel(v) == numel(pRef), 'Inconsistent table dimensions');
        varargout{i} = @(p) interpTable(pRef, v, p);
    end
end
