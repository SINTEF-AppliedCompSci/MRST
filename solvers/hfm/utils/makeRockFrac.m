function G = makeRockFrac(G, K_frac, varargin)
% makeRockFrac can be used to set homogeneous rock properties to the
% fracture grid stored in G.FracGrid given a scalar value for fracture
% permeability and porosity.
%
% SYNOPSIS:
%   G = makeRockFrac(G, K_frac)
%   G = makeRockFrac(G, K_frac, 'pn1', 'pv1')
%   G = makeRockFrac(G, K_frac, 'pn1', 'pv1', 'pn2', pn2)
%
% REQUIRED PARAMETERS:
%
%   G      - Grid data structure containing fracture grids (for each
%            fracture line or plane) in G.FracGrid
%
%   K_frac - Scalar Darcy permeability for homogeneous fractures.
%
% OPTIONAL PARAMETERS:
%
%   permtype - 'homogeneous' or 'heterogeneous'. If 'heterogeneous' is
%               passed as the permtype, this function assigns a random
%               permeability distribution to each fracture grid cell. To
%               manually assign a specific permeability distribution, the
%               user must access G.FracGrid.Frac#.rock.perm as shown in the
%               code written below. In that case, one does not need to call
%               makeRockFrac.
%
%   porosity - Scalar value for rock porosity (0<porosity<1)
%   
% RETURNS:
%   G - Grid data structure with structure rock added to each fracture grid
%       (Frac#) in G.FracGrid.
%
% NOTE:
%   This function uses randi() to generate a random permeability field when
%   'heterogeneous' is specified as 'permtype'.
%
% SEE ALSO:
%   FracTensorGrid2D

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

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


opt = struct('permtype', 'homogeneous', 'porosity', []);
opt = merge_options(opt, varargin{:});
if ~isempty(varargin) && strcmp(opt.permtype,'homogeneous')==0 && strcmp(opt.permtype,'heterogeneous')==0
    fprintf(['Wrong argument string. Valid Arguments: ''Homogeneous'' or ''heterogeneous''.\n',...
        'Setting homogeneous permeability in fractures...\n']);
end
for i = 1:numel(fieldnames(G.FracGrid))
    Gf = G.FracGrid.(['Frac',num2str(i)]);
    if isempty(varargin) || strcmp(opt.permtype,'homogeneous')
        G.FracGrid.(['Frac',num2str(i)]).rock.perm = ones(Gf.cells.num, 1)*darcy*K_frac;
    elseif strcmp(opt.permtype,'heterogeneous')
        G.FracGrid.(['Frac',num2str(i)]).rock.perm = ...
            (randi(100,Gf.cells.num,1)*darcy./randi(100,Gf.cells.num,1))*K_frac;
    else
        G.FracGrid.(['Frac',num2str(i)]).rock.perm = ones(Gf.cells.num, 1)*darcy*K_frac;
    end
    if ~isempty(opt.porosity)
        assert(opt.porosity<1 && opt.porosity>0,...
            'Rock porosity must be a single real number between 0 and 1');
        G.FracGrid.(['Frac',num2str(i)]).rock.poro = ...
            opt.porosity*ones(Gf.cells.num,1);
    end
end