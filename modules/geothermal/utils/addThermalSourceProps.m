function src = addThermalSourceProps(src, varargin)
%Add thermal source conditions to existing source struct.
%
% SYNOPSIS:
%  fluid = addThermalSourceProps(bc, 'pn1', pv1, ...);
%
% PARAMETERS:
%   src   - Source structure created with e.g., addSource
%
% OPTIONAL PARAMETERS:
%   'T'     - Temperature. Either scalar, or one value per face in
%             src.cell. For cells with heat flux thermal source, set 'T' to
%             nan.
%
%   'Hflux' - Heat flux. Either scalar, or one value per face in src.cell.
%             For cells with temperature thermal source, set 'Hflux' to
%             nan.
%
% RETURNS:
%   src - Updated source struct with thermal properties
%
% NOTE:
%   Each face in src must have either a given temperature or heat flux

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
    
    opt = checkInput(src, opt, 'T'    );
    opt = checkInput(src, opt, 'Hflux');
    
    assert(~any(opt.T > 0 & ~isnan(opt.Hflux)), ...
                           'Multiple thermal BCs given for the same face');
    assert(~any(isnan(opt.T) & isnan(opt.Hflux)), ...
           ['Please provide either a temperature or heat flux for ', ...
            'each boundary face']                                  );
    
    % Set thermal BC properties
    src.T     = opt.T;
    src.Hflux = opt.Hflux;
    
end

%-------------------------------------------------------------------------%
function opt = checkInput(src, opt, name)

    nc = numel(src.cell);
    if numel(opt.(name)) == 1
        opt.(name) = repmat(opt.(name), nc, 1);
    end
    assert(numel(opt.(name)) == nc, ...
                  [name, ' must be a scalar, or one per face in bc.face']);

end