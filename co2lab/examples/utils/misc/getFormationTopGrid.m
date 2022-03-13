function [Gt, rock2D, petrodata] = getFormationTopGrid(formation,coarsening_level)
% Load the formation grid at a given coarsening level
%
% SYNOPSIS:
%   function [Gt rock2D] = getFormationTopGrid('Utsirafm',N)
%
% DESCRIPTION:
%
% PARAMETERS:
%   formation        - Name of formation in CO2Atlas directory
%   coarsening_level - Coarsening factor.  If set to 1, a grid with
%                      approximately one cell per datapoint is produced.
%
% RETURNS:
%   Gt      - top surface grid
%   rock2D  - rock structure. Formed using averages contained in petrodata,
%             unless heterogeneous properties were available.
%   petrodata (optional) - contains avg rock properties as reported by
%                          NDP's Atlas, or computed from any available
%                          heterogeneous properties.

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

    moduleCheck('libgeometry');
    [grdecl dataset petroinfo] = ...
        getAtlasGrid(formation, 'coarsening', coarsening_level);%#ok
    
    assert(~isempty(grdecl), ...
        'Empty grdecl returned. Check the spelling of the formation name. Also ensure you have downloaded the dataset.')
    
    % Computing the top-surface grid
    G = processGRDECL(grdecl{1});
    ncvec=nan(numel(G),1);
    for i=1:numel(ncvec)
       ncvec(i)=G(i).cells.num; 
    end
    [nc,j]=max(ncvec); %#ok
    G=G(j);
    try
       G = mcomputeGeometry(G);
    catch
       G = computeGeometry(G);
    end
    
    Gt = topSurfaceGrid(G);
    
    % Setting up the rock structure
    petrodata   = petroinfo{1};
    rock2D.perm = repmat(petrodata.avgperm, Gt.cells.num, 1);
    rock2D.poro = repmat(petrodata.avgporo, Gt.cells.num, 1);
    if isfield(petrodata,'avgntg')
        rock2D.ntg  = repmat(petrodata.avgntg,  Gt.cells.num, 1);
    end
    
    
    % Update rock structure to include heterogeneous values, if they exist.
    
    % NB: A conditional statement is done to ensure PERMX exists before
    % calling grdecl2Rock as the function specifically assigns the
    % permeability tensor values, but will also assign porosity and
    % net-to-gross if they exist. If PERMX doesn't exist, PORO and NTG are
    % added to rock2D in a separate step, using a hack to construct the 3D
    % rock fields (perm and poro) that are required by averageRock()
    
    g = grdecl{1};
    
    if isfield(g,'PERMX')
        
        % NB: may need to check for PERMY, PERMZ fields if isotropic rock
        % properties are available @@
        rock    = grdecl2Rock(grdecl{1}, G.cells.indexMap);
        
        if ~isfield(rock,'poro')
            rock.poro = rock2D.poro;
        end
        
        rock2D  = averageRock(rock, Gt);
        updateNaNs();
        
    elseif isfield(g,'PORO') || isfield(g,'NTG')

        rockprop = {'PORO', 'NTG'};
        for i = 1 : numel(rockprop),
          prop = rockprop{i};
          if isfield(g, prop),
             rock.(lower(prop)) = reshape(g.(prop)(G.cells.indexMap), [], 1);
          end
        end
        rock.perm   = rock2D.perm;
        if ~isfield(rock,'poro')
            rock.poro = rock2D.poro;
        end
        rock2D  = averageRock(rock, Gt);
        updateNaNs();
        
    end
    
    
    % ---------------------------------------------------------------------
    
    function updateNaNs()
    % Replace any nan values with the property's average value. Since perm,
    % poro, or ntg may not exist, an additional conditional statement is
    % required.

        if isfield(rock2D,'perm') && any(isnan(rock2D.perm))
            warning(['Nan permeability values are being replaced with the '...
                'average value.']);
            rock2D.perm(isnan(rock2D.perm)) = petrodata.avgperm; 
        end
        if isfield(rock2D,'poro') && any(isnan(rock2D.poro))
            warning(['Nan porosity values are being replaced with the '...
                'average value.']);
            rock2D.poro(isnan(rock2D.poro)) = petrodata.avgporo;
        end
        if isfield(rock2D,'ntg') && any(isnan(rock2D.ntg))
            warning(['Nan net-to-gross values are being replaced with the '...
                'average value.']);
            rock2D.ntg(isnan(rock2D.ntg))   = petrodata.avgntg;
        end

    end

end
