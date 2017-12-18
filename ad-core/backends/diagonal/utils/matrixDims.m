function varargout = matrixDims(varargin)
    varargout = cell(1, nargout);
    [varargout{:}] = size(varargin{:});
end