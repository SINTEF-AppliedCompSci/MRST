function [ G, Gt, rock, rock2D ] = makeSleipnerModelGrid( varargin )
% this should eventually be placed in co2lab/utils/

% SYNOPSIS:
%  [G, Gt, rock, rock2D] = makeSleipnerModelGrid()
%  [G, Gt, rock, rock2D] = makeSleipnerModelGrid('useIEAGHGmodel', true)
%  [G, Gt, rock, rock2D] = makeSleipnerModelGrid('useIEAGHGmodel',true, 'refineLevel',2)
%
% Default varargin:
%  modelName    = 'IEAGHGmodel' (other options: 'ORIGINALmodel')
%  refineLevel  = 1             (other options: 2,3,4,...) 
%  plotsOn      = false         (other options: true)
%
% PARAMETERS:
%   G      - Data structure for 3D grid
%   Gt     - Data structure for topsurface grid
%   rock   - Data structure for 3D rock parameters
%   rock2D - Data structure for rock parameters for topsurface grid


% SEE ALSO:
%   makeSleipnerVEmodel 


opt = struct('modelName','IEAGHGmodel', 'refineLevel',1, 'plotsOn',false);
opt = merge_options(opt, varargin{:});




if strcmpi(opt.modelName,'IEAGHGmodel')

    % there should be a directory named "co2lab/data/sleipner/" which
    % contains files such as M9X1.grdecl, M9X1_perm_X_mD_.inc, etc.
    try

        disp([' -> Reading SleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat']);
        datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
        load(fullfile(datadir,['SleipnerGlobalCoords_numRef', num2str(opt.refineLevel)])); clear datadir
        return;

    catch %#ok<*CTCH>
        disp(['    SleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat has not yet been created.']);
        disp('    Building G, Gt, and rock2D from grdecl files...')

        % First loading of Sleipner Eclipse grid (to get PERMX, PERMZ,
        % PORO)
        sdir    = fullfile('data', 'sleipner');
        disp([' -> Reading data from: ' sdir]);
        grdecl  = readGRDECL(fullfile(mrstPath('co2lab'), sdir, 'SLEIPNER.DATA'));
        % this grdecl contains: cartDims, COORD, ZCORN, ACTNUM, PERMX,
        % PERMZ, PORO
        clear sdir

        % Reshaping
        lines = reshape(grdecl.COORD,6,[]);
        grdecl.COORD = lines(:); clear lines

        
        % Second loading of Sleipner Eclispe grid, to get MAPAXES
        moduleCheck('deckformat', 'mex');
        sl_file = fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'M9X1.grdecl'); % IEAGHG 
        fn      = fopen(sl_file);
        gr  = readGRID(fn, fileparts(sl_file), initializeDeck());
        % this grdecl contains: GRID, and others. grdecl.GRID contains
        % MAPUNITS, MAPAXES, cartDims, COORD, ZCORN, ACTNUM
        fclose(fn);


        % Add data loaded from first loading of Sleipner Eclispe grid
        grdecl.MAPAXES = gr.GRID.MAPAXES;
        grdecl.MAPUNITS = gr.GRID.MAPUNITS;
        clear gr sl_file
        
        
        %%%%%% First assess grids
        [G, Gt, rock, rock2D] = getGrids( grdecl );
        % Then visualize grids.
        if opt.plotsOn; [ hfig, hax ] = plot3DandTopGrids( G, Gt ); end
        
        
        % required?
        % Then, we de-activate cells (i.e., remove) that correspond to the
        % bottom and top layers that contain shale
        %grdecl.ACTNUM(grdecl.PERMX<200) = 0;
        
        %%%%%% Then perform refinement of grdecl data, which includes a
        %%%%%% step to remove the cell layers corresponding to the caprock
        %%%%%% and bottom shale
        grdecl_refined = getRefinedGrdecl( grdecl );
        
        
        %%%%%% And assess new grids
        [G, Gt, rock, rock2D] = getGrids( grdecl_refined );
        % Then visualize grids.
        if opt.plotsOn; [ hfig, hax ] = plot3DandTopGrids( G, Gt ); end
        

       


        % Store data
        disp(' ')
        disp([' -> Writing SleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat'])
        if ~isdir(datadir)
           mkdir(datadir);
        end
        save(fullfile(datadir,['SleipnerGlobalCoords_numRef', num2str(opt.refineLevel)]), 'G', 'Gt', 'rock', 'rock2D');
        clear datadir

    end

    
elseif strcmpi(opt.modelName,'ORIGINALmodel')
    
    % there should be a directory named "co2lab/data/sleipner/original/"
    % which contains files such as sleipner_prep.data, etc.
    % etc.
    try

        disp([' -> Reading OriginalSleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat']);
        datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
        load(fullfile(datadir,['OriginalSleipnerGlobalCoords_numRef', num2str(opt.refineLevel)])); clear datadir
        return;

    catch %#ok<*CTCH>
        disp(['    OriginalSleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat has not yet been created.']);
        disp('    Building REFINED G, Gt, and rock2D from grdecl files...')

        % Open and read original/sleipner_prep.data
        moduleCheck('deckformat', 'mex');
        sl_file = fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'original', 'sleipner_prep.data'); % IEAGHG original
        fn      = fopen(sl_file);
        grdecl  = readGRID(fn, fileparts(sl_file), initializeDeck());
        % this grdecl contains: GRID, and others. grdecl.GRID contains
        % cartDims, MAPUNITS, MAPAXES, COORD, ZCORN, ACTNUM as well as
        % PERMX, PERMY, PERMZ, and PORO
        fclose(fn);
        
        % Rename
        grdecl = grdecl.GRID;

        % Reshaping
        lines = reshape(grdecl.COORD,6,[]);
        grdecl.COORD = lines(:); clear lines
        
        
        
        %%%%%% First assess grids
        [G, Gt, rock, rock2D] = getGrids( grdecl );
        % Then visualize grids.
        if opt.plotsOn; [ hfig, hax ] = plot3DandTopGrids( G, Gt ); end
        
        % Then, we de-activate cells (i.e., remove) that correspond to the
        % bottom and top layers that contain shale
        %grdecl.ACTNUM(grdecl.PERMX<200) = 0;

        
        %%%%%% Then perform refinement of grdecl data, which includes a
        %%%%%% step to remove the cell layers corresponding to the caprock
        %%%%%% and bottom shale
        grdecl_refined = getRefinedGrdecl( grdecl );
        
        
        %%%%%% And assess new grids
        [G, Gt, rock, rock2D] = getGrids( grdecl_refined );
        % Then visualize grids.
        if opt.plotsOn; [ hfig, hax ] = plot3DandTopGrids( G, Gt ); end
        

        

        % Store data
        disp(' ')
        disp([' -> Writing OriginalSleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat'])
        if ~isdir(datadir)
           mkdir(datadir);
        end
        save( fullfile(datadir,['OriginalSleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat']), ...
            'G', 'Gt', 'rock', 'rock2D', ...
            '-v7.3'); % use flag to save as v7.3 to avoid issues with compression of a lot of data
        clear datadir
        
    end
        
    
% elseif useOriginal_model
%     
%     % there should be a directory named "co2lab/data/sleipner/original/"
%     % which contains files such as sleipner_prep.data, etc.
%     % etc.
%     try
% 
%         disp(' -> Reading OriginalSleipnerGlobalCoords.mat');
%         datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
%         load(fullfile(datadir,'OriginalSleipnerGlobalCoords')); clear datadir
%         %return;
% 
%     catch %#ok<*CTCH>
%         disp('    OriginalSleipnerGlobalCoords.mat has not yet been created.');
%         disp('    Building G, Gt, and rock2D from grdecl files...')
% 
%         % Open and read original/sleipner_prep.data
%         moduleCheck('deckformat', 'mex');
%         sl_file = fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'original', 'sleipner_prep.data'); % IEAGHG original
%         fn      = fopen(sl_file);
%         grdecl  = readGRID(fn, fileparts(sl_file), initializeDeck());
%         % this grdecl contains: GRID, and others. grdecl.GRID contains
%         % cartDims, MAPUNITS, MAPAXES, COORD, ZCORN, ACTNUM as well as
%         % PERMX, PERMY, PERMZ, and PORO
%         fclose(fn);
%         
%         % Rename
%         grdecl = grdecl.GRID;
% 
%         % Reshaping
%         lines = reshape(grdecl.COORD,6,[]);
%         grdecl.COORD = lines(:); clear lines
% 
%         % Then, we remove the bottom and top layers that contain shale
%         grdecl.ACTNUM(grdecl.PERMX<200) = 0;
% 
% 
%         % Recompute X and Y coordinates in terms of the provided axes
%         % (depths, Z, do not require any recomputation)
%         coords        = reshape(grdecl.COORD,3,[])';
%         coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
%         coords        = coords';
%         grdecl.COORD  = coords(:); clear coords
% 
% 
%         % Next, we process the grid and compute geometry
%         mrstModule add libgeometry opm_gridprocessing
%         G = processGRDECL(grdecl); % note: processgrid() didn't work
%         G = mcomputeGeometry(G);
% 
%         % Adding tags needed by topSurfaceGrid
%         G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];
% 
%         % Construct petrophysical model
%         rock = grdecl2Rock(grdecl, G.cells.indexMap);
%         rock.perm = convertFrom(rock.perm, milli*darcy);
%         clear grdecl
% 
%         % Construct top-surface grid
%         disp(' -> Constructing top-surface grid');
%         [Gt, G] = topSurfaceGrid(G);
%         rock2D  = averageRock(rock, Gt);
% 
% 
%         % Store data
%         disp(' ')
%         disp(' -> Writing OriginalSleipnerGlobalCoords.mat')
%         if ~isdir(datadir)
%            mkdir(datadir);
%         end
%         save(fullfile(datadir,'OriginalSleipnerGlobalCoords'), 'G', 'Gt', 'rock', 'rock2D');
%         clear datadir
%          
%     end
    
    
else
    error('Undefined model grid name.')

end

% -------------------------------------------------------------------------


    function [G, Gt, rock, rock2D] = getGrids( grdecl )
        
            % Recompute X and Y coordinates in terms of the provided axes
            % (depths, Z, do not require any recomputation)
            coords        = reshape(grdecl.COORD,3,[])';
            coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
            coords        = coords';
            grdecl.COORD  = coords(:); clear coords


            % Next, we process the grid and compute geometry
            mrstModule add libgeometry opm_gridprocessing
            G = processGRDECL(grdecl); % processgrid() didn't work
            G = mcomputeGeometry(G);

            % Adding tags needed by topSurfaceGrid
            G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];

            % Construct petrophysical model
            rock = grdecl2Rock(grdecl, G.cells.indexMap);
            rock.perm = convertFrom(rock.perm, milli*darcy);
            clear grdecl_refined

            % Construct top-surface grid
            disp(' -> Constructing top-surface grid');
            [Gt, G] = topSurfaceGrid(G);
            rock2D  = averageRock(rock, Gt);
            disp(' ')
    end
    
% -------------------------------------------------------------------------

    function grdecl_refined = getRefinedGrdecl( grdecl )
        
        % We perform a few steps for GRID REFINEMENT:

        % 1- first cut the grdecl by removing the caprock layer(s) and the
        % bottom shale layer(s). Removing the caprock and shale layers must
        % be handling uniquely for each grid model.
        
        if strcmpi(opt.modelName,'IEAGHGmodel')
            % The IEAGHGmodel is comprised of 43 cells in the vertical
            % direction: the caprock is comprised of 2 cell layers, and the
            % bottom shale is comprised of 7 cell layers. Thus, keep layers
            % 3 to 36, and remove layers 1,2,37-43. There will be 34 layers
            % remaining.
            lowerBound = [1 1 3];
            upperBound = [grdecl.cartDims(1) grdecl.cartDims(2) 36];
            ind = [lowerBound; upperBound]';
            numSandLayers = 34;
        
        elseif strcmpi(opt.modelName,'ORIGINALmodel')
            % The ORIGINALmodel is comprised of 7 cells in the vertical
            % direction: the caprock and bottom shale are each 1 cell layer
            % thick. Thus, keep layers 2 to 6, and remove layers 1 and 7.
            % There will be 5 layers remaining.
            lowerBound = [1 1 2];
            upperBound = [grdecl.cartDims(1) grdecl.cartDims(2) 6];
            ind = [lowerBound; upperBound]';
            numSandLayers = 5;
        end
        
        grdecl_cut = cutGrdecl(grdecl, ind);
        % grdecl_new contains fields that correspond to the remaining SAND
        % layers.

        
        % 2- second coarsen the grdecl in the z direction since we don't
        % need to keep the resolution of the SAND layers (i.e., we treat
        % the SANE layers as 1 homogeneous layer for VE). In order to
        % ensure PERMX, PERMY, PERMZ, and PORO are coarsened, 'only_grid'
        % must be set to false. Note that coarseGrdecl() assumes grdecl_cut
        % contains PERM in all 3 directions. In the case that PERMY does
        % not exist in the grdecl structure, we assume it is the same as
        % PERMX and thus make a copy of it so as to pass it into the
        % function
        if ~isfield('grdecl_cut','PERMY')
            grdecl_cut.PERMY = grdecl_cut.PERMX;
        end
        dim = [1 1 numSandLayers];
        grdecl_coarsened = coarseGrdecl(grdecl_cut, dim, 'only_grid',false);
        % grdecl_new now contains fields that correspond to the single
        % layer.

        
        % 3- third refine the grdecl in the X and Y direction (Z already has
        % a resolution of 1 cell) by the refinement level specified
        dim = [opt.refineLevel; opt.refineLevel; 1];
        grdecl_refined = refineGrdecl(grdecl_coarsened, dim);
        
        
        % Add some original fields to the refined grdecl structure:
        grdecl_refined.MAPAXES = grdecl.MAPAXES;
        grdecl_refined.MAPUNITS = grdecl.MAPUNITS;
            
    end


end

