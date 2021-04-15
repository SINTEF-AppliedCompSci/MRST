function varargout = dispif(bool, varargin)
%Produce textual output contingent upon predicate.
%
% SYNOPSIS:
%        dispif(bool, format, arg, ...)
%   nc = dispif(bool, format, arg, ...)
%
% PARAMETERS:
%   bool   - Boolean variable.
%
%   format - SPRINTF format specification.
%
%   arg    - OPTIONAL arguments to complete 'format'.
%
% RETURNS:
%   nc     - Number of characters printed to output device.  If 'bool' is
%            `false`, then `nc=0`.  Only returned if specifcially
%            requested. 
%
% NOTE:
%   Function used for making code cleaner where `verbose` option is used.
%
% SEE ALSO:
%   `sprintf`, `tocif`.

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


%error(nargoutchk(0, 1, nargout, 'struct'));

nc = 0;
if bool, nc = fprintf(varargin{:}); end

if nargout > 0, varargout{1} = nc; end
end
