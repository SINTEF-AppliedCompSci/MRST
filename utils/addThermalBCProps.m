function bc = addThermalBCProps(bc, varargin)
%Add thermal boundary conditions to existing bc struct.
%
% SYNOPSIS:
%  fluid = addThermalFluidProps(bc, 'pn1', pv1, ...);
%
% PARAMETERS:
%   bc   - Boundary condition structure created with e.g., addBC
%
% OPTIONAL PARAMETERS:
%   'T'     - Temperature. Either scalar, or one value per face in bc.face.
%             For faces with heat flux thermal BC, set 'T' to nan.
%
%   'Hflux' - Heat flux. Either scalar, or one value per face in bc.face.
%             For faces with temperature thermal BC, set 'Hflux' to nan.
%
% RETURNS:
%   bc - Updated boundary condition struct with thermal properties
%
% NOTE:
%   Each face in bc must have either a given temperature or heat flux

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('T'    , nan, ...
                 'Hflux', nan);
    opt = merge_options(opt, varargin{:});
    
    opt = checkInput(bc, opt, 'T'    );
    opt = checkInput(bc, opt, 'Hflux');
    
    assert(~any(opt.T > 0 & ~isnan(opt.Hflux)), ...
                           'Multiple thermal BCs given for the same face');
    assert(~any(isnan(opt.T) & isnan(opt.Hflux)), ...
           ['Please provide either a temperature or heat flux for ', ...
            'each boundary face']                                  );
    
    % Set thermal BC properties
    bc.T     = opt.T;
    bc.Hflux = opt.Hflux;
    
end

%-------------------------------------------------------------------------%
function opt = checkInput(bc, opt, name)

    nf = numel(bc.face);
    if numel(opt.(name)) == 1
        opt.(name) = repmat(opt.(name), nf, 1);
    end
    assert(numel(opt.(name)) == nf, ...
                  [name, ' must be a scalar, or one per face in bc.face']);

end