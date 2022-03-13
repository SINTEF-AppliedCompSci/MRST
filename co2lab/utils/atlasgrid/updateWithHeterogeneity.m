function [ deck, petroinfo ] = updateWithHeterogeneity( deck, petroinfo, meta_top, opt )
% Update deck and petroinfo to include heterogeneous rock properties
%
% SYNOPSIS:
%   [deck, petroinfo] = updateWithHeterogeneity( deck, petroinfo, meta_top, opt);
%
% DESCRIPTION:
%   The CO2Atlas directory is searched and any heterogeneous rock data that
%   exists for a formation matching deck.name is loaded and added to the
%   deck structure. If no rock data files for a formation matching
%   deck.name exist, deck and petroinfo are returned from function
%   unchanged.
%
% PARAMETERS:
%   deck       - cell arrays of data members in ECLIPSE grid structures.
%
%   petroinfo  - average perm / porosity as given in the atlas.
%                Will contain NaN where values are not provided or possible
%                to calculate.
%
%   meta_top   - Metainformation for top data, as returned by
%                readAtlasGrids(), a local function in getAtlasGrid()
%
% RETURNS:
%   deck      - updated cell arrays of data members in ECLIPSE grid
%               structures. If heterogeneous rock properties are available
%               for a formation, these are added to deck structure. If not,
%               deck is returned unmodified.
%
%   petroinfo - heterogeneous rock properties and updated average perm /
%               porosity computed from these heterogeneous properties. If
%               no heterogeneous data is available for a formation,
%               petroinfo is returned unmodified.
%
% SEE ALSO:
%   `getAtlasGrid`, `convertAtlasTo3D`

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

    %% Look in co2atlas directories and load any perm, poro, or ntg data
    % that belong to frm_name.
    if any(strcmpi(deck.name, getNorthSeaNames()))
        atlas_dir     = getDatasetPath('co2atlas');
    elseif any(strcmpi(deck.name, getNorwegianSeaNames()))
        atlas_dir     = getDatasetPath('co2atlasnorwegiansea');
    elseif any(strcmpi(deck.name, getBarentsSeaNames()))
        atlas_dir     = getDatasetPath('co2atlasbarentssea');
    else
        error(['Could not find ' deck.name ' in any Atlas directory.']);
    end
    dataset_info  = dir(atlas_dir);
    dataset_names = {dataset_info(~[dataset_info.isdir]).name};
    [rawdata] = deal({});
    for i = 1:numel(dataset_names)

        name        = dataset_names{i};
        datasetType = regexp(name, '_', 'split');

        % Find datasetType that contains _poro or _perm
        ind = strcmpi(datasetType, 'poro') | strcmpi(datasetType, 'perm') ...
            | strcmpi(datasetType, 'ntg');
        tmp.name = [datasetType{~ind}];

        if ~strcmpi(tmp.name, deck.name)
            continue
        end

        if any(strcmpi(datasetType, 'poro'))
            tmp.variant = 'porosity';
        elseif any(strcmpi(datasetType, 'perm'))
            tmp.variant = 'permeability';
        elseif any(strcmpi(datasetType, 'ntg'))
            tmp.variant = 'net_to_gross';
        else
            % continue thru loop, without proceeding to reading of dataset
            continue
        end

        % reading of dataset
        [meta, data] = readAAIGrid(fullfile(atlas_dir, name));

        % perform coarsening if required
        meta.cellsize   = meta.cellsize*opt.coarsening;
        tmp.meta        = meta;
        tmp.data        = data(1:opt.coarsening:end, 1:opt.coarsening:end);
        tmp.meta.ncols  = size(tmp.data, 2);
        tmp.meta.nrows  = size(tmp.data, 1);
        tmp.meta.dims   = [tmp.meta.nrows, tmp.meta.ncols];

        % perform refinement if called for
        if opt.refining > 1
            tmp = refineInfo(tmp, opt.refining);
        end

        % Add all property varients found to rawdata cell structure
        rawdata{end+1} = tmp; %#ok

    end

    %% Assess the variants of the data in rawdata cell structure
    % The rawdata variants are passed into a function which assess their
    % coordinations and performs interpolation to get aligned data in order
    % to make the grdecl (deck) structure.
    for i = 1:numel(rawdata)
        
        rd_variant = rawdata{i}.variant;
        
        if strcmpi(rd_variant,'net_to_gross')
            meta_ntg = rawdata{i}.meta;
            data_ntg = rawdata{i}.data;
            deck.NTG = convertAtlasHeterogenValsTo3D(meta_top, opt.nz, meta_ntg, data_ntg);
            
        elseif strcmpi(rd_variant,'permeability')
            % NB: data_perm is in units of m2, while rawdata is in mD.
            meta_perm   = rawdata{i}.meta;
            data_perm   = convertFrom(rawdata{i}.data, milli*darcy);
            deck.PERMX  = convertAtlasHeterogenValsTo3D(meta_top, opt.nz, meta_perm, data_perm);
            [deck.PERMY, deck.PERMZ] = deal(deck.PERMX);
            
        elseif strcmpi(rd_variant,'porosity')
            meta_poro = rawdata{i}.meta;
            data_poro = rawdata{i}.data;
            deck.PORO = convertAtlasHeterogenValsTo3D(meta_top, opt.nz, meta_poro, data_poro);
            
        end
        
    end
    
    %% Update petroinfo (if required)
    
    if isfield(deck,'PERMX')
        petroinfo.avgperm = mean(deck.PERMX(~isnan(deck.PERMX)));
    end
    
    if isfield(deck,'PORO')
        petroinfo.avgporo = mean(deck.PORO(~isnan(deck.PORO)));
    end
    
    if isfield(deck,'NTG')
        petroinfo.avgntg  = mean(deck.NTG(~isnan(deck.NTG)));
    end
    

end

function PROP = convertAtlasHeterogenValsTo3D(meta_top, nz, meta_prop, data_prop)
% PROP is interpolated property data aligned to a cell-centered grid set-up
% using corner node and cellsize info found in meta_top.

    ndims   = [meta_top.dims, nz + 1];
    dims    = ndims - 1;
    h       = meta_top.cellsize;
    
    xl = meta_top.xllcorner;
    yl = meta_top.yllcorner;

    % Create grids (cell-center coordinates)
    [XC, YC]  = ndgrid(linspace(xl+h/2, dims(1)*h + xl - h/2, dims(1)), ...
                       linspace(yl+h/2, dims(2)*h + yl - h/2, dims(2)));
    
    % Function handles that will return data which is aligned to the
    % x,y coordinates passed into function, i.e., alignedData = Fh(x,y)
    F_prop  = interpolateData(meta_prop, data_prop);

    % Get all data variant in appropriate coordinates:
    % perm, poro, and ntg need to be cell-center coordinates (i.e., must
    % match with grdecl.cartDims), whereas top and thick are at pillars
    % coordinates, which are node-centered.
    
    % for cell-center coordinates
    xc = reshape(XC(:,:,1), [], 1);
    yc = reshape(YC(:,:,1), [], 1);
    prop  = reshape(F_prop(xc, yc) , dims(1), dims(2));
 
    % Assign
    PROP = reshape(prop, [dims(1)*dims(2), 1]);
    % NB: any NaN values are left as is (rather than converting to Inf)

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