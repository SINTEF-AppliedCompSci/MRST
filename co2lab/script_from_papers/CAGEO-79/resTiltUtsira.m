function resTiltUtsira(varargin)
    
    opt.load_previous           = false;
    opt.utsira_coarsening_level = 1;
    opt.save_tilt_filename      = 'utsira_subtrap_function.mat';
    opt = merge_options(opt,varargin{:});
    
    mrstModule('add', 'matlab_bgl', 'coarsegrid', 'opm_gridprocessing', 'libgeometry');
    close all;
    
    %% Load the two datasets (Utsira from CO2 atlas, and Sleipner from IEAGHG)
    % If IEAGHG Sleipner dataset is located elsewhere, modify the following line
    %dir_sleipner = fullfile(VEROOTDIR, 'data', 'sleipner', 'M9X1.grdecl'); % IEAGHG
    dir_sleipner = fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'M9X1.grdecl'); % IEAGHG   
    [GtS, GtU] = prepareGrids(dir_sleipner, opt.utsira_coarsening_level);
    
    %% Establishing grid interpolants
    [FS, FU] = createInterpolants(GtS, GtU);
    
    %% Determining the sample points, based on the domain of Sleipner
    sizeS       = 2 * GtS.cartDims; % twice the resolution of the Sleipner grid
    get_range = @(d, pts, n) linspace(min(pts(:, d)), max(pts(:, d)), n);
    xrange    = get_range(1, GtS.nodes.coords, sizeS(1));
    yrange    = get_range(2, GtS.nodes.coords, sizeS(2));
    xdom   = xrange(floor(sizeS(1)*0.1):floor(sizeS(1)*0.9)); % Shrink domain 
    ydom   = yrange(floor(sizeS(2)*0.1):floor(sizeS(2)*0.9)); % slightly
    [x, y]    = meshgrid(xdom, ydom); % This will be the sampling grid
    
    %% Sampling the two grids and recentering them vertically at zero depth
    smpl_U      = FU(x, y); 
    smpl_S      = FS(x, y); 
    smpl_U_mean = mean(smpl_U(isfinite(smpl_U(:))));
    smpl_S_mean = mean(smpl_S(isfinite(smpl_S(:)))); 
    smpl_U      = smpl_U - smpl_U_mean;
    smpl_S      = smpl_S - smpl_S_mean;
    
    %% Define difference surface interpolant
    FD = @(x, y) FS(x,y) - FU(x,y) - (smpl_S_mean - smpl_U_mean);

    %% Display sampled grids, and difference surface
    figure(1); 
    subplot(1,2,1); 
    hold on;
    mesh(x, y, smpl_U); mesh(x, y, smpl_S); title('Recentered surfaces');
    set(gca,'zdir', 'reverse'); hold on; view(-62,22);
    subplot(1,2,2); 
    mesh(x,y, FD(x,y)); title('Difference surface');
    set(gca, 'zdir', 'reverse'); view(-64,60);
    
    %% Compute trapping volume and area of complete Utsira and Sleipner grids
    %  (Used for reference later)
    [trapvolU, areaU, resU] = computeTrapVolAndArea(GtU);%#ok
   
    %% Measuring trap volumes for difference surface, for different tilts
    min_tilt_x  = -0.04;    min_tilt_y = -0.03; % in radians
    max_tilt_x  =  0.05;    max_tilt_y =  0.03; % in radians 
    steps       = 20;                          % same number of steps in x and y
    theta_x_vec = linspace(min_tilt_x, max_tilt_x, steps);
    theta_y_vec = linspace(min_tilt_y, max_tilt_y, steps);

    iplot       = true; % whether to plot each grid during analysis 
    figure(2); % Prepare a figure window to draw plots during loop execution
    
    if opt.load_previous
       load(opt.save_tilt_filename);
    else
       [fine_scale_trap_vols, areaS] = ...
           computeTrapsAllTilts(theta_x_vec, theta_y_vec, xdom, ydom, FD, iplot);
    
       % Saving result (the output is used e.g. by the UtsiraLongTerm script)
       tv_per_area = fine_scale_trap_vols/areaS;%#ok
       save(opt.save_tilt_filename, 'theta_x_vec', 'theta_y_vec', 'tv_per_area', ...
            'fine_scale_trap_vols', 'areaS');
    end
    %load utsira_subtrap_function;
    
    %% Finalize analysis, and report results
    figure(3);
    fprintf('REPORT FOR ANALYSIS COMPARING SLEIPNER AND UTSIRA DATA\n');
    reportAnalysis(GtU, resU, theta_x_vec, theta_y_vec, fine_scale_trap_vols/areaS, iplot);

    fprintf('Push enter to continue.\n');
    pause;
    
    %% Redo fine-scale/coarse-scale analysis. This time, use a smoothed
    %% version of the Sleipner grid as 'coarse scale'
    
    % We make a smoothed version of the Sleipner grid, so that the level of
    % detail should be comparable to the same geographical section of the
    % Utsira grid.  For that, we create a Gaussian kernel of appropriate
    % size, and convolute the top surface with it.
    k        = makeKernel(max(GtU.cells.volumes), max(GtS.cells.volumes));
    FSsmooth = makeSmoothedGridInterpolator(GtS, k);

    % Surface difference interpolator
    FDsmooth = @(x, y) FS(x, y) - FSsmooth(x,y) + 1000;

    % To avoid boundary effect of smoothing operator, we shrink the domain
    % sufficiently. (One might think that we should subtract half the size of
    % k at each side, but since the intervals between sampling points in
    % 'xrange' and 'yrange' are only half of the spacing in the GtS grid used
    % when constructing the kernel, we achieve the desired cutoff by removing
    % a number of intervals equal to the number of intervals in k).
    xdom = xrange(ceil(size(k,1)) : end - ceil(size(k,1)));
    ydom = yrange(ceil(size(k,2)) : end - ceil(size(k,2)));
    
    % Repeating the measuring of trap volumes for different tilts, this time
    % for the diffence surface between the smoothed and the original Slepiner
    % dataset. 
    figure(2);
    [smoothed_trap_vols, areaSsmooth] = ...
        computeTrapsAllTilts(theta_x_vec, theta_y_vec, xdom, ydom, FDsmooth, iplot);
    fprintf('REPORT FOR ANALYSIS USING FINE-SCALE ESTIMATE FROM SMOOTHED SLEIPNER\n');
    figure(3);
    reportAnalysis(GtU, resU, theta_x_vec, theta_y_vec, smoothed_trap_vols/areaSsmooth, iplot);
                   
end

% ----------------------------------------------------------------------------
function F = makeSmoothedGridInterpolator(Gt, kernel)
    nx = Gt.cartDims(1)+1;
    ny = Gt.cartDims(2)+1;

    % X and Y coordinates remain as before
    X  = reshape(Gt.nodes.coords(:,1), nx, ny);
    Y  = reshape(Gt.nodes.coords(:,2), nx, ny);
    Z  = reshape(Gt.nodes.z          , nx, ny);
    
    Z = filter2(kernel, Z);
    
    F = TriScatteredInterp(X(:), Y(:), Z(:));
end

% ----------------------------------------------------------------------------
function k = makeKernel(large, small)
    ksz = 2 * sqrt(large/small);
    sigma = (1/2).^2;
    gridsize = floor(ksz/2);
    [x, y] = meshgrid(-gridsize:gridsize, -gridsize:gridsize);
    x = x/gridsize; y = y/gridsize;
    k = exp(-sigma * (x.^2 + y.^2));
    k = k/sum(k(:));
end

% ----------------------------------------------------------------------------
function reportAnalysis(coarse_grid, coarsetraps, theta_x_vec, theta_y_vec, ...
                        tilt_trap_table, do_plot)

    % Define result grid interpolant, optionally visualize a 2D plot of trapping
    % volume per area as function of tilt
    FTVol = @(x, y) interp2(theta_x_vec, theta_y_vec, tilt_trap_table', x, y, 'spline');
    tilt_samples_x = linspace(theta_x_vec(1), theta_x_vec(end), 100); % 100 samples
    tilt_samples_y = linspace(theta_y_vec(1), theta_y_vec(end), 100); % ditto
    ftvol_sampled  = computeSamples(FTVol, tilt_samples_x, tilt_samples_y);

    % Plot it
    if do_plot
        plotSubscaleTrapVolPerArea(tilt_samples_x, tilt_samples_y, ftvol_sampled);
    end
    
    % Compute estimates of fine-scale trapping as a percentage of the coarse
    % trap volume.
    
    % To avoid extrapolation issues, we will limit the analysis to those
    % cells of the coarse grid that are not tilted _more_ than the tilt used 
    % for analysis
    
    % First, we compute the angle between the normals and the z axis 
    % (really, the tangent of the angle, but considered a sufficiently good
    % approximation here)
    min_tilt_x = theta_x_vec(1); max_tilt_x = theta_x_vec(end);
    min_tilt_y = theta_y_vec(1); max_tilt_y = theta_y_vec(end);
    coarse_area  = sum(coarse_grid.cells.volumes);
    z_norm = bsxfun(@rdivide, coarse_grid.cells.normals, coarse_grid.cells.normals(:,3));
    tan_xy = -z_norm(:, 1:2);
    
    % Then, we use this information to identify cells in the coarse grid that
    % are not tilted beyond what we have sampled above:
    ind = ( (tan_xy(:,1) < max_tilt_x) & (tan_xy(:,1) > min_tilt_x) & ...
            (tan_xy(:,2) < max_tilt_y) & (tan_xy(:,2) > min_tilt_y) );
    % Also, disregard cells that are traps anyway on the coarse grid
    ind = ind & (coarsetraps.traps == 0); 
    
    trapvol_coarse  = sum(volumesOfTraps(coarse_grid, coarsetraps, []));
    
    % Identify the unique tilt that was found to lead to maximum subscale
    % trapping capacity, using our subsampled grid:
    [max_ftvol, max_ix]  = max(ftvol_sampled(:));
    [max_x, max_y]       = ind2sub(size(ftvol_sampled), max_ix);
    
    % As a maximum estimate, we ignore local tilt, and choose the maximum
    % subscale trapping capacity per area as a basis for computing the total
    est_max = max_ftvol * coarse_area;
    
    % We compute another estimate, this time using the local tilt in each
    % case to determine the subscale trapping volume for a given cell
    est_tiltadjusted = ...
        sum(FTVol(tan_xy(ind, 1), tan_xy(ind,2)) .* coarse_grid.cells.volumes(ind));
    
    % In the third estimate, we recenter the tilt table so that maximum value
    % occurs for zero tilt, and recompute
    top_angle = [tilt_samples_x(max_x), tilt_samples_y(max_y)];
    est_tiltadjusted_recentered = ...
        sum(FTVol(tan_xy(ind,1) - top_angle(1), tan_xy(ind,2) - top_angle(2)) ...
            .* coarse_grid.cells.volumes(ind));
    
    % Report results
    fprintf('Fine-scale trapping as a percentage of coarse scale\n');
    fprintf('Upper bound estimate:    %2.5g\n', est_max / trapvol_coarse * 100);
    fprintf('Tilt-adjusted estimate:  %2.5g\n', est_tiltadjusted / trapvol_coarse * 100);
    fprintf('With recentered tilt:    %2.5g\n', ...
            est_tiltadjusted_recentered / trapvol_coarse * 100); 
end


% ----------------------------------------------------------------------------
function [tiltTraps, area] = ...
        computeTrapsAllTilts(theta_x_vec, theta_y_vec, xdom, ydom, FD, do_plot)

    steps = numel(theta_x_vec);  assert(steps == numel(theta_y_vec));
    tiltTraps = nan(steps, steps);
    
    count = 1;  % solely for progress reporting purposes
    for ix = 1:steps
        for iy = 1:steps
            %% Reporting progress
            fprintf('Running inclination case %i out of %i.\n', count, steps^2);
            count = count + 1;
            
            %% Constructing tilted fine-scale grid and computing its trap volume

            % Setting inclinations
            theta_x     = theta_x_vec(ix);
            theta_y     = theta_y_vec(iy);
            
            Gt_tilt           = makeTiltedTopSurface(theta_x, theta_y, xdom, ydom, FD);
            [vol, area, res]  = computeTrapVolAndArea(Gt_tilt);
            % (Remark: although 'area' is here recomputed at each iteration,
            %          its value remains the same, regardless of degree of tilt.)

            % storing results
            tiltTraps(ix, iy) = vol;
            
            %% Plotting current (inclined) grid, with detected traps (if any)
            if do_plot
                clf;
                plotCellData(Gt_tilt, res.traps, res.traps ~=0);
                plotGrid(Gt_tilt, 'FaceColor', 'none');
                view(3); drawnow;
            end
        end
    end
end

% ----------------------------------------------------------------------------
function fun_sampled = computeSamples(fun, samples_x, samples_y)
% Sample 2D function 'fun' at a cartesian grid of sample points, defined by
% 'samples_x' and 'samples_y'
    [x, y] = meshgrid(samples_x, samples_y);
    fun_sampled = fun(x, y);
end
% ----------------------------------------------------------------------------

function plotSubscaleTrapVolPerArea(xvec, yvec, f_sampled)
    % Display result
    mesh(xvec, yvec, f_sampled);
    view([-1 -1 2]); grid off; axis tight;
    xlabel('\theta_x');
    ylabel('\theta_y');
    zlabel('Vol/Area (m)');
    set(gca,'FontSize',16);
end

% ----------------------------------------------------------------------------
function Gt = makeTiltedTopSurface(theta_x, theta_y, xdom, ydom, interp)

    % Starting by constructing a grid in the shape of a regular prism
    G = tensorGrid(xdom, ydom, [0, 30]);

    % computing the interpolant surface
    isurf = interp(G.nodes.coords(:,1), G.nodes.coords(:,2));
    
    % Computing tilt
    tilt = G.nodes.coords(:,1) * sin(theta_x) + G.nodes.coords(:,2) * sin(theta_y);
           
    % adjusting z-coordinates of the grid by adding interpolant and tilt           
    G.nodes.coords(:,3) = G.nodes.coords(:,3) + isurf + tilt;

    % repositioning to a depth of 1000m
    G.nodes.coords(:,3) = G.nodes.coords(:,3) - min(G.nodes.coords(:,3)) + 1000;
    
    % computing geometry and extracting top surface grid 
    try
       Gt = topSurfaceGrid(mcomputeGeometry(G));
    catch
       Gt = topSurfaceGrid(computeGeometry(G));       
    end
end

% ----------------------------------------------------------------------------
function [vol, area, resS] = computeTrapVolAndArea(Gt)
    resS  = trapAnalysis(Gt, true); % use 'false' for node-based algorithm
    vol  = sum(volumesOfTraps(Gt, resS, []));
    area = sum(Gt.cells.volumes);
end

% ----------------------------------------------------------------------------
function [Gt_sleipner, Gt_utsira] = prepareGrids(sl_file, utsira_coarsening_level)
%% Load Sleipner and Utsira datasets
    moduleCheck('deckformat', 'mex');
    
    %% Load Sleipner Eclipse grid and convert it to MRST format
    fn     = fopen(sl_file);
    grdecl = readGRID(fn, fileparts(sl_file), initializeDeck());
    grdecl = grdecl.GRID;

    fclose(fn);
    
    % Recompute coordinates in terms of the provided axes
    coords        = reshape(grdecl.COORD,3,[])';
    coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
    coords        = coords';
    grdecl.COORD  = coords(:);
    
    % Converting to MRST format, computing geometry and extracting the top
    % surface grid
    try
       Gt_sleipner = topSurfaceGrid(mcomputeGeometry(mprocessGRDECL(grdecl)));
    catch
       Gt_sleipner = topSurfaceGrid(computeGeometry(processGRDECL(grdecl)));
    end
       
    %% Load Utsira atlas grid
    grdecl = getAtlasGrid('Utsirafm', 'coarsening', utsira_coarsening_level);
    
    % Converting to MRST format, computing geometry and extracting the top
    % surface grid
    try
       Gt_utsira = topSurfaceGrid(mcomputeGeometry(mprocessGRDECL(grdecl{1})));
    catch
       Gt_utsira = topSurfaceGrid(computeGeometry(processGRDECL(grdecl{1})));
    end
    
end

% ----------------------------------------------------------------------------
function [FS, FU] = createInterpolants(GtS, GtU)
%% Prepare and return grid interpolants
%  (to facilitate comparison, since grids are of different resolution) A
%  scattered-point interpolant is useful here, given the irregular shape of the
%  Utsira domain (many points in the regular grid remain undefined).
    
    %% Interpolant for Sleipner grid
    [nx, ny]  = deal(GtS.cartDims(1)+1, GtS.cartDims(2)+1);

    X         = reshape(GtS.nodes.coords(:,1), nx, ny); 
    Y         = reshape(GtS.nodes.coords(:,2), nx, ny);
    Z         = reshape(GtS.nodes.z          , nx, ny);

    FS        = TriScatteredInterp(X(:), Y(:), Z(:));
    
    %% Interpolant for Utsira grid
    FU = TriScatteredInterp(GtU.nodes.coords(:,1), ...
                            GtU.nodes.coords(:,2), ...
                            GtU.nodes.z);
end
