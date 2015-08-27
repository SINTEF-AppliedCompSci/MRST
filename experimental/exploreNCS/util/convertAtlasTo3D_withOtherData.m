function grdecl = convertAtlasTo3D_withOtherData(meta_thick, meta_top, data_thick, data_top, nz, meta_perm, meta_poro, meta_ntg, data_perm, data_poro, data_ntg)
%Create GRDECL struct from CO2 storage atlas thickness/top data
%
% SYNOPSIS:
%   grdecl = convertAtlasTo3D(m_thick, m_top, d_thick, d_top, 3, ...
%                   m_perm, m_poro, m_ntg, d_perm, d_poro, d_ntg)
%
% DESCRIPTION:
%   Given two datasets (plus extra data) with possibily non-matching nodes,
%   interpolate and combine to form a GRDECL struct suitable for 3D
%   simulations after a call to processGRDECL or for further manipulation.
%
% REQUIRED PARAMETERS:
%   meta_thick - Metainformation for the thickness data as produced from
%                readAAIGrid
%
%   meta_top   - Metainformation for the top data as produced from
%                readAAIGrid
%
%   data_thick - Data for the thickness as produced from readAAIGrid
%
%   data_top   - Data for the top as produced from readAAIGrid
%
%   nz         - 
%
%   meta_perm  - Metainformation for the permeability data as produced from
%                rawdata of getAtlasGrid
%
%   meta_poro  - Metainformation for the porosity data as produced from
%                rawdata of getAtlasGrid
%
%   meta_ntg   - Metainformation for the net-to-gross data as produced from
%                rawdata of getAtlasGrid
%
%   data_perm  - Data for permeability as produced from rawdata of
%                getAtlasGrid
%
%   data_poro  - Data for porosity as produced from rawdata of
%                getAtlasGrid
%
%   data_ntg   - Data for net-to-gross as produced from rawdata of
%                getAtlasGrid
%
%
% RETURNS:
%   grdecl - GRDECL struct suitable for processGRDECL
%
%
% NOTES:
%   This function has the additional capability of adding perm, poro, and
%   ntg data to the grdecl structure
%
% SEE ALSO:
%   getAlasGrid, convertAtlasTo3D_withOtherData

%{
#COPYRIGHT#
%}

    ndims = [meta_top.dims, nz + 1];
    dims = ndims - 1;
    h = meta_top.cellsize;

    
    grdecl.cartDims = reshape(dims, 1, []);
    
    xl = meta_top.xllcorner;
    yl = meta_top.yllcorner;
    
    % Create grids (node coordinates)
    [X, Y, Z]  = ndgrid(linspace(xl, dims(1)*h + xl, ndims(1)), ...
                        linspace(yl, dims(2)*h + yl, ndims(2)), ...
                        linspace(0,  1,              ndims(3)));
                    
    % Create grids (cell-center coordinates)
    [XC, YC]  = ndgrid(linspace(xl+h/2, dims(1)*h + xl - h/2, dims(1)), ...
                        linspace(yl+h/2, dims(2)*h + yl - h/2, dims(2)));
    

    % make function handles that will return data which is aligned to the
    % x,y coordinates passed into function, i.e., alignedData = Fh(x,y)
    
    F_top   = interpolateData(meta_top,   data_top);
    F_thick = interpolateData(meta_thick, data_thick);
    
    F_perm  = interpolateData(meta_perm, data_perm);
    F_poro  = interpolateData(meta_poro, data_poro);
    F_ntg   = interpolateData(meta_ntg, data_ntg);
    
    %if(false)
    %   x = squeeze(X(:,:,1));
    %   y = squeeze(Y(:,:,1));  
    %   thick = F_thick(x, y)/dims(3);
    %   top   = F_top(x,y);
    %else
    
    % get all data variant are appropriate coordinates:
    % perm, poro, and ntg need to be cell-center coordinates (i.e., must
    % match with grdecl.cartDims), whereas top and thick are at pillars
    % coordinates, which are node-centered.
    
    % for node coordinates
    x = reshape(X(:,:,1), [], 1);
    y = reshape(Y(:,:,1), [], 1);
    thick = reshape(F_thick(x, y), ndims(1), ndims(2))/dims(3);
    top   = reshape(F_top(x, y)  , ndims(1), ndims(2));
    
    % for cell-center coordinates
    xc = reshape(XC(:,:,1), [], 1);
    yc = reshape(YC(:,:,1), [], 1);
    perm  = reshape(F_perm(xc, yc) , dims(1), dims(2));
    poro  = reshape(F_poro(xc, yc) , dims(1), dims(2));
    ntg   = reshape(F_ntg(xc, yc)  , dims(1), dims(2));
    
    %end
    
   
    
    %
    
    thick(thick<=0) = NaN;
    %
    
    % Uniformly partition the thickness across the layers
    for i = 1:ndims(3)
        Z(:,:,i) = top + thick*(i-1);
    end
    Z = sort(Z, 3);
    
    
    % Make pillars
    n = prod(ndims(1:2));
    lines = zeros([n, 6]);
    lines(:, [1, 4]) = reshape(X(:,:,[1, end]), [n, 2]);
    lines(:, [2, 5]) = reshape(Y(:,:,[1, end]), [n, 2]);
    lines(:, [3, 6]) = reshape(Z(:,:,[1, end]), [n, 2]);
    grdecl.COORD = reshape(lines.', [], 1);

    % Assign z-coordinates
    % ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
    ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
    z   = Z(ind(1), ind(2), ind(3));

    grdecl.ZCORN = z(:);

    % Assign active cells by removing cells containing NaN points
    
    z = squeeze(isnan(Z(:,:,1)));
    
    tmp = z;
    for i = -1:1
        for j = -1:1
            % mask away any cells touching nan-points
            tmp = tmp | circshift(z, [i,j]);
        end
    end
    
    tmp = ~(tmp(1:end-1, 1:end-1));
    for i = 2:dims(3)
        tmp(:,:,i) = tmp(:,:,1);
    end
    grdecl.ACTNUM = reshape(int32(tmp), [], 1);
    
    grdecl.ZCORN(isnan(grdecl.ZCORN)) = inf;
    grdecl.COORD(isnan(grdecl.COORD)) = inf;
    
    % Assign perm, poro, ntg
    grdecl.PERMX = reshape(perm, [dims(1)*dims(2), 1]);
    grdecl.PERMX(isnan(grdecl.PERMX)) = inf;
    
    [grdecl.PERMY, grdecl.PERMZ] = deal(grdecl.PERMX);
    
    grdecl.PORO = reshape(poro, [dims(1)*dims(2), 1]);
    grdecl.PORO(isnan(grdecl.PORO)) = inf;
    
    grdecl.NTG = reshape(ntg, [dims(1)*dims(2), 1]);
    grdecl.NTG(isnan(grdecl.NTG)) = inf;
    
end

function F = interpolateData(meta, data)
    dims = meta.dims;
    % We have a cell centered grid, so subtract by one
    gdims = dims - 1;
    h = meta.cellsize;
    xl = meta.xllcorner;
    yl = meta.yllcorner;
    [X, Y]  = ndgrid(linspace(xl, gdims(1)*h + xl, dims(1)), ...
                     linspace(yl, gdims(2)*h + yl, dims(2))); 
                 
    % alignedData = INTERP2(X,Y,data,x,y)
    % interpolates to find alignedData, the values of the data at the
    % coordinates given by x and y. Matrices X and Y specify the points at
    % which the data is given.
    F =@(x,y) interp2(X',Y',data', x, y);              
    %F = TriScatteredInterp(X(:), Y(:), data(:));
end
