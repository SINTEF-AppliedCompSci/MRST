function [Gt_base, rock_base, seainfo] = makeMyGeomodel(varargin)

    opt.modify_base_rock = true;
    opt = merge_options(opt, varargin{:});
    
    fmName = 'Utsirafm';
    c_level = 4;
    rhoCref = 760;
    
    % Get formation top surface and rock model from CO2 Atlas
    [Gt_base,rock_base] = getFormationTopGrid(fmName,c_level);
    
    
    if opt.modify_base_rock
        % Modify base rock according to poro-perm model

        % Function handles for phi-depth, perm-phi model (where
        % coefficients were computed previously using Sto geomodel which
        % contained heterogeneous rock data)
        P_poro = [-0.000070460825418 0.323426386160052];
        P_perm = [34.971830212726744 -7.245959920003443];
        porofun = @(depth) -5e-5 * depth + 0.24; % for Utsira
        permfun = @(poro) exp(P_perm(1) * poro + P_perm(2)); % darcy

        % Modification
        rock_base.poro = porofun(Gt_base.cells.z);
        rock_base.perm = convertFrom(permfun(rock_base.poro),darcy);
    end

    % Optional cutting of grid:
    ind = (max(Gt_base.cells.centroids(:,2)) - Gt_base.cells.centroids(:,2)) ...
        > 200 * kilo * meter;
    [g, cellmap] = removeCells(Gt_base.parent, find(~ind)); clear Gt_base
    Gt_base      = topSurfaceGrid(g); clear g
    rock_base.poro = rock_base.poro(cellmap);
    rock_base.perm = rock_base.perm(cellmap);
    if isfield(rock_base,'ntg')
        rock_base.ntg = rock_base.ntg(cellmap);
    end
    clear cellmap
    
    
    seainfo = getSeaInfo(fmName, rhoCref);
    
end