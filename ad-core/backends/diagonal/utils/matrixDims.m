function varargout = matrixDims(varargin)
% Overloadable version of size
    varargout = cell(1, nargout);
    [varargout{:}] = size(varargin{:});
end