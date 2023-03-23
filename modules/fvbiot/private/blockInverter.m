function bi = blockInverter(opt)
% Set up operator for inverting block matrices. 
%
% This method is identical to a subroutine in MRST, 
%   see mrst/modules/mpfa/computeMultiPointTrans.m
% It has been copied to a separate file to enable access from other
% funcions.
%
% Copyright statement from the original file is pasted below.
%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

if ~ischar(opt.invertBlocks),
    dispif(opt.verbose, ...
        ['Unsupported option value of type ''%s'' in ', ...
        'option ''invertBlocks''. Reset to default ' , ...
        '(''matlab'')\n'], class(opt.invertBlocks));
    opt.invertBlocks = 'matlab';
end

switch lower(opt.invertBlocks),
    case {'matlab', 'm', 'builtin'},
        bi = @invertDiagonalBlocks;
    case {'mex', 'c', 'accelerated'},
        bi = @invertDiagonalBlocksMex;
    otherwise
        dispif(opt.verbose, ...
            ['Unsupported value ''%s'' in option ', ...
            '''invertBlocks''.\nMust be one of ''matlab'' or ', ...
            '''mex''.\nReset to default (''matlab'').'], ...
            opt.invertBlocks);
        
        bi = @invertDiagonalBlocks;
end