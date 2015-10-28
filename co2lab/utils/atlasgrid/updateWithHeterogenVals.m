function [ deck, petroinfo ] = updateWithHeterogenVals( deck, petroinfo, meta_top, opt )
% Update deck and petroinfo to include heterogeneous rock properties
%
% SYNOPSIS:
%   [deck, petroinfo] = updateWithHeterogenVals( deck, petroinfo, meta_top, opt);
%
%
% PARAMETERS:
%   deck       - cell arrays of data members in ECLIPSE grid structures.
%
%   petroinfo  - average perm / porosity as given in the atlas.
%               Will contain NaN where values are not provided or possible
%               to calculate.
%
% RETURNS:
%   deck       - updated cell arrays of data members in ECLIPSE grid
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
%   getAtlasGrid


    %% Look in co2atlas directory and load any perm, poro, or ntg data
    % that belong to frm_name. These are added to rawdata, and then
    % converted into deck (grdecl) format.
    atlas_dir     = getDatasetPath('co2atlas');
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
            meta_perm   = rawdata{i}.meta;
            data_perm   = rawdata{i}.data;
            deck.PERMX  = convertAtlasHeterogenValsTo3D(meta_top, opt.nz, meta_perm, data_perm);
            [deck.PERMY, deck.PERMZ] = deal(deck.PERMX);
            
        elseif strcmpi(rd_variant,'porosity')
            meta_poro = rawdata{i}.meta;
            data_poro = rawdata{i}.data;
            deck.PORO = convertAtlasHeterogenValsTo3D(meta_top, opt.nz, meta_poro, data_poro);
            
        end
        
    end
    
    % deck can now be returned
    
    %% Update petroinfo (if required)
    
    if isfield(deck,'PERMX')
        petroinfo.perm    = deck.PERMX;
        petroinfo.avgperm = mean(deck.PERMX(~isinf(deck.PERMX)));
    end
    
    if isfield(deck,'PORO')
        petroinfo.poro    = deck.PORO;
        petroinfo.avgporo = mean(deck.PORO(~isinf(deck.PORO)));
    end
    
    if isfield(deck,'NTG')
        petroinfo.ntg     = deck.NTG;
        petroinfo.avgntg  = mean(deck.NTG(~isinf(deck.NTG)));
    end
    
    % petroinfo can now be returned

end

function PROP = convertAtlasHeterogenValsTo3D(meta_top, nz, meta_prop, data_prop)
% PROP is interpolated property data aligned to a cell-centered grid set-up
% using corner node and cellsize info found in meta_top.
% See also convertAltasTo3D

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
    PROP(isnan(PROP)) = inf;

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