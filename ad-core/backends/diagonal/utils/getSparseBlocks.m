function varargout = getSparseBlocks(varargin)
% Get sparse blocks
    varargout = cell(nargout, 1);
    [varargout{:}] = getSparseArguments(varargin{:});
end