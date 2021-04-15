function properties = constantProperties(mu, rho)
%Construct fluid property evaluator from constant data
%
% SYNOPSIS:
%   properties = constantProperties(mu, rho)
%
% PARAMETERS:
%   mu  - Phase viscosities in units of Pa*s.
%
%   rho - Phase densities in units of kilogram/meter^3.
%
% RETURNS:
%   properties -
%         Function handle that supports the following calling convention
%
%            mu       = properties()
%           [mu, rho] = properties()
%
%        The return values are exactly the input parameters 'mu' and 'rho'
%        to function 'constantProperties', respectively.
%
%        For compatiblity with the rest of MRST, the 'properties' function
%        accepts any number of input arguments.  All arguments are ignored.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


   properties = @(varargin) eval_const_prop(mu, rho);
end

%--------------------------------------------------------------------------

function varargout = eval_const_prop(mu, rho)
   varargout{1}                 = mu;
   if nargout > 1, varargout{2} = rho; end
end
