function varargout = initVariablesFastAD(varargin)
   assert (nargin == nargout, ...
          ['Number of output variables must equal the number ', ...
           'of input variables.']);

   n_el = cellfun(@numel, varargin);
   m = n_el(1);
   assert(all(n_el == m));
   
   n         = nargin;
   varargout = cell([1, n]);

   for i = 1 : n,
      jac = zeros(m, n);
      jac(:, i) = 1;
      varargout{i} = FastAD(varargin{i}, jac);
   end
end