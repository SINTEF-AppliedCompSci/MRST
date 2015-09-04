function grdecl = addPermPoroNtgData2grdecl( frm_name, varargin )

% frm_name      formation name
% N             coarsening level (pos for coarsening, neg for refining)

opt = struct('nz',1, 'coarsening', 1, 'refining', 1, 'Verbose', mrstVerbose);
opt = merge_options(opt, varargin{:});

N = opt.coarsening;
R = opt.refining;

    % First get grdecl and rawdata. The grdecl returned by getAtlasGrid
    % contains COORD, cartDims, etc., but not PORO, PERMX, or NTG. Rawdata
    % returned here contains thickness and top data variants for formation.
    if N > 1
        fprintf('\n Dataset will be coarsened. \n')
        [grdecl, rawdata, ~] = getAtlasGrid(frm_name, 'coarsening',N);
    elseif R > 1
        fprintf('\n Dataset will be refined. \n')
        [grdecl, rawdata, ~] = getAtlasGrid(frm_name, 'refining',R);
    else
        fprintf('\n Dataset will not be coarsened nor refined. \n')
        [grdecl, rawdata, ~] = getAtlasGrid(frm_name, 'coarsening',N);
    end

    % Next look in co2atlas directory and load any perm, poro, or ntg data
    % that belong to frm_name
    atlas_dir     = getDatasetPath('co2atlas');
    dataset_info  = dir(atlas_dir);
    dataset_names = {dataset_info(~[dataset_info.isdir]).name};
    coarsening    = N;
    for i = 1:numel(dataset_names)

        dataset     = dataset_names{i};
        datasetStr  = regexp(dataset, '_', 'split');

        % Find datasetStr that contains _poro or _perm
        ind = strcmpi(datasetStr, 'poro') | strcmpi(datasetStr, 'perm') ...
            | strcmpi(datasetStr, 'ntg');
        tmp.name = [datasetStr{~ind}];

        if ~strcmpi(tmp.name,frm_name)
            continue
        end

        if any(strcmpi(datasetStr, 'poro'))
            tmp.variant = 'porosity';
        elseif any(strcmpi(datasetStr, 'perm'))
            tmp.variant = 'permeability';
        elseif any(strcmpi(datasetStr, 'ntg'))
            tmp.variant = 'net_to_gross';
        else
            % continue thru loop, without proceeding to reading of dataset
            continue
        end

        % reading of dataset
        [meta, data] = readAAIGrid(fullfile(atlas_dir, dataset));

        
        % perform coarsening if required
        meta.cellsize = meta.cellsize*coarsening;
        tmp.meta = meta;
        tmp.data = data(1:coarsening:end, 1:coarsening:end);
        tmp.meta.ncols = size(tmp.data, 2);
        tmp.meta.nrows = size(tmp.data, 1);
        tmp.meta.dims = [tmp.meta.nrows, tmp.meta.ncols];

        
        % perform refinement if called for
        if R > 1
            tmp = refineInfo(tmp, R);
        end
        
        % Add poro, perm, or ntg to rawdata cell structure
        rawdata{end+1} = tmp; %#ok

    end
    

    
    % Assess the variants of the data in rawdata cell structure
    for i = 1:numel(rawdata)
        
        rd_variant = rawdata{i}.variant;
        
        if strcmpi(rd_variant,'thickness')
            meta_thick = rawdata{i}.meta;
            data_thick = rawdata{i}.data;
        elseif strcmpi(rd_variant,'top')
            meta_top = rawdata{i}.meta;
            data_top = rawdata{i}.data;
        elseif strcmpi(rd_variant,'net_to_gross')
            meta_ntg = rawdata{i}.meta;
            data_ntg = rawdata{i}.data;
        elseif strcmpi(rd_variant,'permeability')
            meta_perm = rawdata{i}.meta;
            data_perm = rawdata{i}.data;
        elseif strcmpi(rd_variant,'porosity')
            meta_poro = rawdata{i}.meta;
            data_poro = rawdata{i}.data;
        end
        
    end
    

    % The rawdata variants are passed into a function which assess their
    % coordinations and performs interpolation to get aligned data in order
    % to make the grdecl structure.
    % TODO: implement such that ntg data does not exist.
    grdecl = convertAtlasTo3D_withOtherData(meta_thick, meta_top, data_thick, data_top, 1, ...
        meta_perm, meta_poro, meta_ntg, data_perm, data_poro, data_ntg);

    % grdecl can now be returned


end

