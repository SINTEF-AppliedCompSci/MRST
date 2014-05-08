function varargout = dealvec(varargin)
% Deal the elements in the vector in 'varargin' to the output arguments.
% (Works like 'deal', but treats vectors in the same way as cell arrays).

if nargin==1 && numel(varargin{1}) == nargout
    for i = 1:nargout
        varargout{i} = varargin{1}(i);
    end
else
    error('Error: Number of outputs should match number of input elements.');
end
