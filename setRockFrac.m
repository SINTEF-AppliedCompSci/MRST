function G = setRockFrac(G, K_star, varargin)
% Sets rock properties inside fracture grid stored in G.FracGrid given a
% permeability ratio between matrix and fractures.
%
% SYNOPSIS:
%   G = setRockFrac(G, K_frac)
%   G = setRockFrac(G, K_frac, 'pn1', 'pv1')
%   G = setRockFrac(G, K_frac, 'pn1', 'pv1', 'pn2', pn2)
%
% REQUIRED PARAMETERS:
%
%   G      - Grid data structure containing fracture grids (for each
%            fracture line) in G.FracGrid as defined by FracTensorGrid2D.
%
%   K_star - Used to scale fracture permeability by K_star orders of
%            magnitude higher than matrix permeability
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   permtype - 'homogeneous' or 'heterogeneous'
%
%   rockporo - Single numeric value for rock porosity (0<rockporo<1)
%   
% RETURNS:
%   G - Grid data structure with structure rock added to each fracture grid
%       (Line#) in G.FracGrid.
%
% NOTE:
%   This function uses randi() to generate a random permeability field.
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


opt = struct('permtype', [], 'rockporo', []);
opt = merge_options(opt, varargin{:});
if ~isempty(varargin) && strcmp(opt.permtype,'homogeneous')==0 && strcmp(opt.permtype,'heterogeneous')==0
    fprintf(['Wrong argument string. Valid Arguments: ''Homogeneous'' or ''heterogeneous''.\n',...
        'Setting homogeneous permeability in fractures...\n']);
end
for i = 1:numel(fieldnames(G.FracGrid))
    Gf = G.FracGrid.(['Line',num2str(i)]);
    if isempty(varargin) || strcmp(opt.permtype,'homogeneous')
        G.FracGrid.(['Line',num2str(i)]).rock.perm = ones(Gf.cells.num, 1)*darcy()*K_star;
    elseif strcmp(opt.permtype,'heterogeneous')
        G.FracGrid.(['Line',num2str(i)]).rock.perm = ...
            (randi(100,Gf.cells.num,1)*darcy()./randi(100,Gf.cells.num,1))*K_star;
    else
        G.FracGrid.(['Line',num2str(i)]).rock.perm = ones(Gf.cells.num, 1)*darcy()*K_star;
    end
    if isnumeric(opt.rockporo)
        assert(opt.rockporo<1 && opt.rockporo>0,...
            'Rock porosity must be a single real number between 0 and 1');
        G.FracGrid.(['Line',num2str(i)]).rock.poro = ...
            opt.rockporo*ones(Gf.cells.num,1);
    end
end