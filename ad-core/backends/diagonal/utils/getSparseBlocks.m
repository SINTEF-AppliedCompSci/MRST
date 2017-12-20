function varargout = getSparseBlocks(varargin)
    varargout = cell(nargout, 1);
    [varargout{:}] = getSparseArguments(varargin{:});
end