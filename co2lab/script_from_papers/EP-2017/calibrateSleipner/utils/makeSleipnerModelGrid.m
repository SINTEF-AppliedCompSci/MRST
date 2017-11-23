function [ G, Gt, rock, rock2D ] = makeSleipnerModelGrid( varargin )
% gets appropriate grdecl data files and makes MRST-type grid (G, Gt, etc)
% as per user-defined resolution

% SYNOPSIS:
%  [G, Gt, rock, rock2D] = makeSleipnerModelGrid()
%  [G, Gt, rock, rock2D] = makeSleipnerModelGrid('modelName','IEAGHGmodel' )
%  [G, Gt, rock, rock2D] = makeSleipnerModelGrid('modelName','IEAGHGmodel', 'refineLevel',2)
%
% Default varargin:
%  modelName    = 'IEAGHGmodel'
%  refineLevel  = 1             (other options: ...-4,-3,-2, 2, 3, 4,...) 
%  plotsOn      = false         (otherwise true to make grids plots)
%  saveBaseGrids = false        (otherwise true to save base grids to file)
%
% PARAMETERS:
%   G      - Data structure for 3D grid
%   Gt     - Data structure for topsurface grid
%   rock   - Data structure for 3D rock parameters
%   rock2D - Data structure for rock parameters for topsurface grid
%
% OTHER:
%   baseGrids - contains G, Gt, rock, rock2D of: 1) grid made up of all 3
%               layers (caprock, sand, and shale layers), and 2) grid of
%               sand layer only (caprock and bottom shale layers removed)
%               without any vertical-averaging of G.
%             - not passed out of function, only saved to .mat file


% SEE ALSO:
%   `makeSleipnerVEmodel` 


opt = struct(   'modelName'     ,   'IEAGHGmodel'   , ...
                'refineLevel'   ,   1               , ...
                'plotsOn'       ,   false           , ...
                'saveBaseGrids' ,   false           );
opt = merge_options(opt, varargin{:});




if strcmpi(opt.modelName,'IEAGHGmodel')

    % there should be a directory named "co2lab/data/sleipner/" which
    % contains files such as M9X1.grdecl, M9X1_perm_X_mD_.inc, etc.
    try

        fprintf([' -> Reading SleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat']);
        datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
        load(fullfile(datadir,['SleipnerGlobalCoords_numRef', num2str(opt.refineLevel)])); clear datadir
        fprintf('\n   File loaded.\n\n')
        return;

    catch %#ok<*CTCH>
        disp(['    SleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat has not yet been created.']);
        if opt.refineLevel>1
            disp('    Building REFINED G, Gt, and rock2D from grdecl files...')
        elseif opt.refineLevel==1
            disp('    Building G, Gt, and rock2D from grdecl files...')
        elseif opt.refineLevel<0
            disp('    Building COARSENED G, Gt, and rock2D from grdecl files...')
        end

        % First loading of Sleipner Eclipse grid (to get PERMX, PERMZ,
        % PORO)
        sdir    = fullfile('data', 'sleipner');
        fprintf(['\n -> Reading data from: ' sdir '\n']);
        grdecl  = readGRDECL(fullfile(mrstPath('co2lab'), sdir, 'SLEIPNER.DATA'));
        % this grdecl contains: cartDims, COORD, ZCORN, ACTNUM, PERMX,
        % PERMZ, PORO
        clear sdir
        
        % Add name of grid
        grdecl.name = ['Sleipner_' opt.modelName];

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
        
        
        % required?
        % Then, we de-activate cells (i.e., remove) that correspond to the
        % bottom and top layers that contain shale
        %grdecl.ACTNUM(grdecl.PERMX<200) = 0;
        
        
        % Perform refinement of grdecl data, which includes a step to
        % remove the cell layers corresponding to the caprock and bottom
        % shale
        [grdecl_refined, grdecl_cut] = getRefinedGrdecl( grdecl );
        
        % Get refined and vertically-averged grids
        [G, Gt, rock, rock2D] = getGrids( grdecl_refined, true, false );
        

        % Store data
        disp(' ')
        disp([' -> Writing SleipnerGlobalCoords_numRef', num2str(opt.refineLevel), '.mat'])
        if ~isdir(datadir)
           mkdir(datadir);
        end
        save(fullfile(datadir,['SleipnerGlobalCoords_numRef', num2str(opt.refineLevel)]), ...
            'G', 'Gt', 'rock', 'rock2D', ...
            '-v7.3'); % use flag to save as v7.3 to avoid issues with compression of a lot of data
        
        
        % Store baseGrids data
        % Note: constucting and saving the base grids is meant for
        % visualization purposes of permeability and depth plots, when
        % caprock and shale layers are required, or when the full vertical
        % discretization is required. Since the base grids are not refined
        % or coarsened from their original resolution, saving base grids
        % needs to be done only once and not for mulitple numRef cases,
        % since the resulting base grids will be the same each time.
        if opt.saveBaseGrids
            
            % Construct layered, fully discrete grids and model properties
            [layers.G, layers.Gt, layers.rock, layers.rock2D] = getGrids( grdecl, true, false );
            
            % Construct cut grids (that have not yet been
            % vertically-averaged or refined)
            [cut.G, cut.Gt, cut.rock, cut.rock2D] = getGrids( grdecl_cut, true, false );
            
            if ~isdir(datadir)
                mkdir(datadir);
            end
            
            fprintf('\n Writing base grids to SleipnerGlobalCoords_baseGrids.mat file. \n')
            save(fullfile(datadir, 'SleipnerGlobalCoords_baseGrids'), ...
                'layers','cut', '-v7.3');
        end
        

        clear datadir

    end

 
else
    error('Undefined model grid name.')

end

% -------------------------------------------------------------------------


    function [G, Gt, rock, rock2D] = getGrids( grdecl, doCoordTransform, useAvgRock )
        
            if doCoordTransform
                % If required, recompute X and Y coordinates in terms of
                % the provided axes (depths, Z, do not require any
                % recomputation)
                coords        = reshape(grdecl.COORD,3,[])';
                coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
                coords        = coords';
                grdecl.COORD  = coords(:); clear coords
            end


            % Next, we process the grid and compute geometry
            mrstModule add libgeometry opm_gridprocessing
            G = processGRDECL(grdecl); % processgrid() didn't work
            G = mcomputeGeometry(G);

            % Adding tags needed by topSurfaceGrid
            G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];

            % Construct petrophysical model
            % TODO: if no PORO or PERMX, PERMY, etc fields exist in grdecl,
            % the average values of the Sleipner benchmark should be
            % retrieved instead.
            fprintf(' -> Constructing rock data\n\n');
            if useAvgRock
                % rock structure for 3D grid.
                avgrock      = getAvgRock(grdecl.name);
                rock.perm(1:G.cells.num,1) = avgrock.avgperm;
                rock.perm(1:G.cells.num,2) = avgrock.avgperm;
                rock.perm(1:G.cells.num,3) = avgrock.avgperm;
                rock.poro(1:G.cells.num,1) = avgrock.avgporo;
                %rock2D.perm = rock.avgperm;
                %rock2D.poro = rock.avgporo;
            else
                rock        = grdecl2Rock(grdecl, G.cells.indexMap);
                rock.perm   = convertFrom(rock.perm, milli*darcy);
            end
            
            % Construct top-surface grid
            fprintf(' -> Constructing top-surface grid\n\n');
            [Gt, G] = topSurfaceGrid(G);
            rock2D  = averageRock(rock, Gt);
            

            

    end

% -------------------------------------------------------------------------
% getAvgRock() is modified from local function found within getAtlasgrid.m
    function avgrock = getAvgRock(name)
        if strcmpi(name(end-4:end), 'model')
            name = name(1:end-5);
        end
        switch lower(name)

            % Sleipner grids
            % TODO: average values could be obtained from the grdecl data
            % files instead of using the avg val provided in benchmark.
            case 'sleipner_ieaghg'
                tmp = [2000, 0.36]; % Singh et al 2010 (SPE 134891)

            otherwise
                tmp = [NaN NaN];
        end
        avgrock.avgperm = tmp(1)*milli*darcy;
        avgrock.avgporo = tmp(2);
    end
    
% -------------------------------------------------------------------------

    function [grdecl_refined, grdecl_cut] = getRefinedGrdecl( grdecl )
        
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
            
        end
        
        grdecl_cut = cutGrdecl(grdecl, ind);
        % grdecl_cut contains fields that correspond to the remaining SAND
        % layers.
        
        % 2- second coarsen the grdecl in the z direction since we don't
        % need to keep the resolution of the SAND layers (i.e., we treat
        % the SAND layers as 1 homogeneous layer for VE). In order to
        % ensure PERMX, PERMY, PERMZ, and PORO are coarsened, 'only_grid'
        % must be set to false. Note that coarseGrdecl() assumes grdecl_cut
        % contains PERM in all 3 directions. In the case that PERMY does
        % not exist in the grdecl structure, we assume it is the same as
        % PERMX and thus make a copy of it so as to pass it into the
        % function
        if ~isfield(grdecl_cut,'PERMY') && isfield(grdecl_cut,'PERMX')
            grdecl_cut.PERMY = grdecl_cut.PERMX;
        end
        dim = [1 1 numSandLayers];
        
        % perms are coarsened along with grid only if they exist
        if isfield(grdecl_cut,{'PERMX','PERMY','PERMZ'})
            grdecl_coarsened = coarseGrdecl(grdecl_cut, dim, 'only_grid',false);
        else
            grdecl_coarsened = coarseGrdecl(grdecl_cut, dim, 'only_grid',true);
        end
        % grdecl_new now contains fields that correspond to the single
        % layer.

        
        % 3- third refine the grdecl in the X and Y direction (Z already has
        % a resolution of 1 cell) by the refinement level specified
        if(opt.refineLevel>0)
            dim = [opt.refineLevel; opt.refineLevel; 1];
            grdecl_refined = refineGrdecl(grdecl_coarsened, dim);
        else
            dim = abs([opt.refineLevel; opt.refineLevel; 1])';
            if isfield(grdecl_cut,{'PERMX','PERMY','PERMZ'})
                grdecl_refined = coarseGrdecl(grdecl_coarsened, dim,'only_grid',false);
            else
                grdecl_refined = coarseGrdecl(grdecl_coarsened, dim,'only_grid',true);
            end
        end
        
        % Add some original fields to the refined (and cut) grdecl structure:
        if isfield(grdecl,{'MAPAXES','MAPUNITS'})
            grdecl_refined.MAPAXES = grdecl.MAPAXES;
            grdecl_refined.MAPUNITS = grdecl.MAPUNITS;
            grdecl_cut.MAPAXES = grdecl.MAPAXES;
            grdecl_cut.MAPUNITS = grdecl.MAPUNITS;
        end
        grdecl_refined.name = grdecl.name;
        grdecl_cut.name = grdecl.name;
            
    end


end

