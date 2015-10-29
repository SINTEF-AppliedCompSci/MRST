function [Gt rock2D, petrodata] = getFormationTopGrid(formation,coarsening_level)
% Load the formation grid at a given coarsening level
%
% SYNOPSIS:
%   function [Gt rock2D] = getUtsiraTopGrid(coarsening_level)
%
% DESCRIPTION:
%
% PARAMETERS:
%   coarsening_level - Coarsening factor.  If set to 1, a grid with
%                      approximately one cell per datapoint is produced.
%
% RETURNS:
%   Gt     - top surface grid
%   rock2D - rock structure
%
    moduleCheck('libgeometry');
    [grdecl dataset petroinfo] = ...
        getAtlasGrid(formation, 'coarsening', coarsening_level);

    % Computing the Utsira top-surface grid
    G = processGRDECL(grdecl{1});
    ncvec=nan(numel(G),1);
    for i=1:numel(ncvec)
       ncvec(i)=G(i).cells.num; 
    end
    [nc,j]=max(ncvec);
    G=G(j);
    G = mcomputeGeometry(G);
    Gt = topSurfaceGrid(G);
    
    % Setting up the rock structure
    petrodata   = petroinfo{1};
    rock2D.perm = repmat(petrodata.avgperm, Gt.cells.num, 1);
    rock2D.poro = repmat(petrodata.avgporo, Gt.cells.num, 1); 
    
    
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
        
    elseif isfield(g,'PORO') || isfield(g,'NTG')

        rockprop = {'PORO', 'NTG'};
        for i = 1 : numel(rockprop),
          prop = regexprep(rockprop{i}, '\W', '_');
          if isfield(g, prop),
             rock.(lower(prop)) = reshape(g.(prop)(G.cells.indexMap), [], 1);
          end
        end
        rock.perm   = [rock2D.perm, rock2D.perm, rock2D.perm];
        if ~isfield(rock,'poro')
            rock.poro = rock2D.poro;
        end
        rock2D      = averageRock(rock, Gt);
        
    end
    
    
    % Replace any nan values with the property's average value. Since ntg
    % may not exist, an additional conditional statement is required.
    if any(isnan(rock2D.perm))
        dispif(mrstVerbose, 'Nan permeability values are being replaced with the average value.\n')
        rock2D.perm(isnan(rock2D.perm)) = petrodata.avgperm; 
    end
    if any(isnan(rock2D.poro))
        dispif(mrstVerbose, 'Nan porosity values are being replaced with the average value.\n')
        rock2D.poro(isnan(rock2D.poro)) = petrodata.avgporo;
    end
    if isfield(rock2D,'ntg') && any(isnan(rock2D.ntg))
        dispif(mrstVerbose, 'Nan net-to-gross values are being replaced with the average value.\n')
        rock2D.ntg(isnan(rock2D.ntg))   = petrodata.avgntg;
    end
    
end
