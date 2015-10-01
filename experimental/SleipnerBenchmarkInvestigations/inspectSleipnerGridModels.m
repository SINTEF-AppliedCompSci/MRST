% Comparison of all three Sleipner grid models, given specified refinement
% levels. The refinement level(s) (i.e., -2, -3, etc for coarsening; 2, 3,
% etc for refinement) are defined, and each refinement case is analyzed
% separately. Results are reported in table displayed in command window and
% with plots.

% Refinement cases:
numRefCases = [-4]; %-4 -2 1 2 4];

% Do you want top surface grids to be re-centered when making surface
% comparisons? (true for recentering, false for no recentering)
recenterGrids = false;

% pre-allocate report cell
report = cell( numel(numRefCases), 3);
for k = 1:numel(numRefCases)
    
    % get current refinement case
    numRef = numRefCases(k);
    
    %% Construct base grids (or load if already written to .mat file)
    
    % IEAGHG grid
    try
        load SleipnerGlobalCoords_baseGrids.mat
        bg_ieaghg.layers = layers; bg_ieaghg.cut = cut;
        return;
    catch
        [ ~, ~, ~, ~ ] = makeSleipnerModelGrid( 'modelName','IEAGHGmodel', 'saveBaseGrids',true);
        load SleipnerGlobalCoords_baseGrids.mat
        bg_ieaghg.layers = layers; bg_ieaghg.cut = cut;
    end
    
    layers  = bg_ieaghg.layers;
    cut     = bg_ieaghg.cut;
    clear bg_ieaghg
    
    % visualize base grids and permeabilities
    [ hfig, hax ] = plot3DandTopGrids( layers.G, layers.Gt );
    [ hfig ]      = plotGridPerms( layers.G, layers.rock ); 
    [ hfig ]      = plotGridPerms( layers.Gt, layers.rock2D );
    [ hfig, hax ] = plot3DandTopGrids( cut.G, cut.Gt );
    [ hfig ]      = plotGridPerms( cut.G, cut.rock ); 
    [ hfig ]      = plotGridPerms( cut.Gt, cut.rock2D );
    
    % some analysis
    zmax_3ls = max(layers.G.nodes.coords(:,3)); %
    zmax_cut = max(cut.G.nodes.coords(:,3));
    zmin_3ls = min(layers.G.nodes.coords(:,3));
    zmin_cut = min(cut.G.nodes.coords(:,3));
    fprintf('\n Removed a %d meter thickness from top of grid. \n', abs(zmin_3ls-zmin_cut) )
    fprintf('\n Removed a %d meter thickness from bottom of grid. \n', abs(zmax_3ls-zmax_cut) )
    % above message assumes uniform cell size in vertical direction
    
    
    % ORIGINAL grid
    try
        load OriginalSleipnerGlobalCoords_baseGrids.mat
        bg_original.layers = layers; bg_original.cut = cut;
        return;
    catch
        [ ~, ~, ~, ~ ] = makeSleipnerModelGrid( 'modelName','ORIGINALmodel', 'saveBaseGrids',true);
        load OriginalSleipnerGlobalCoords_baseGrids.mat
        bg_original.layers = layers; bg_original.cut = cut;
    end
    
    layers  = bg_original.layers;
    cut     = bg_original.cut;
    clear bg_original
    
    % visualize base grids and permeabilities
    [ hfig, hax ] = plot3DandTopGrids( layers.G, layers.Gt );
    [ hfig ]      = plotGridPerms( layers.G, layers.rock ); 
    [ hfig ]      = plotGridPerms( layers.Gt, layers.rock2D );
    [ hfig, hax ] = plot3DandTopGrids( cut.G, cut.Gt );
    [ hfig ]      = plotGridPerms( cut.G, cut.rock ); 
    [ hfig ]      = plotGridPerms( cut.Gt, cut.rock2D );
    
    % some analysis
    zmax_3ls = max(layers.G.nodes.coords(:,3)); %
    zmax_cut = max(cut.G.nodes.coords(:,3));
    zmin_3ls = min(layers.G.nodes.coords(:,3));
    zmin_cut = min(cut.G.nodes.coords(:,3));
    fprintf('\n Removed a %d meter thickness from top of grid. \n', abs(zmin_3ls-zmin_cut) )
    fprintf('\n Removed a %d meter thickness from bottom of grid. \n', abs(zmax_3ls-zmax_cut) )
    % above message assumes uniform cell size in vertical direction
    
    
    %% Create Sleipner model grids for each of the three types
    [ G_ieaghg, Gt_ieaghg, rock_ieaghg, rock2D_ieaghg ]          = makeSleipnerModelGrid( 'modelName','IEAGHGmodel','refineLevel', numRef);
    [ G_original, Gt_original, rock_original, rock2D_original ]  = makeSleipnerModelGrid( 'modelName','ORIGINALmodel','refineLevel', numRef );
    [ G_inhouse, Gt_inhouse, rock_inhouse, rock2D_inhouse ]      = makeSleipnerModelGrid( 'modelName','INHOUSEmodel','refineLevel', numRef );

    
    gts = [Gt_ieaghg; Gt_original; Gt_inhouse];
    gs  = [G_ieaghg; G_original; G_inhouse];

    names = {'IEAGHG'; 'ORIGINAL'; 'INHOUSE'};
    mymap = [ 0 0 1; 0 1 0; 1 0 0 ]; % for coloring grids

    %% Inspect Top Surface Grid and H
    figure; set(gcf,'Position',[1 1 1200 800])

    % pre-allocate res cell (for current numRef)
    res = cell( 1, numel(gts));
    for i = 1:numel(gts)
        
        Gt = gts(i);
        G  = gs(i);
        name = names{1};
        avthk = mean(Gt.cells.H);
        numCells = Gt.cells.num;

        plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', mymap(i,:));

        % Comparisons of grid details
        res{i}.name      = names{i};
        res{i}.refLevel  = numRef;
        res{i}.cells     = Gt.cells.num;
        res{i}.zmin      = min(Gt.cells.z);
        res{i}.zmax      = max(Gt.cells.z);
        res{i}.volume    = sum(G.cells.volumes);
        res{i}.surfarea  = sum(Gt.cells.volumes);
        res{i}.avgthk    = mean(Gt.cells.H);
        
        % trapping analysis
        tan     = trapAnalysis(Gt, false);
        tac     = trapAnalysis(Gt, true);
        
        res{i}.ctrapvols = volumesOfTraps(Gt,tac);
        res{i}.ccapacity = sum(res{i}.ctrapvols);
        res{i}.ntrapvols = volumesOfTraps(Gt,tan);
        res{i}.ncapacity = sum(res{i}.ntrapvols);
        
        % store for later use
        report{k,i} = res{i};
    end

    legend(cellfun(@(x) x.name, res, 'UniformOutput', false), 'Location', 'EastOutside')

    view(3)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Top Surfaces elevation')

    %% Create table to compare grids:
   
   fprintf('\n\n--------------Num Ref. =  %-2d-----------------------\n', numRef);
   fprintf('%-20s&   Cells  & Min Surf Ele. & Max Surf Ele. &  Volume  &  Surf. Area  & Avg. Thk. \\\\\n', 'Name');

   for i = 1:3
   fprintf('%-20s&  %6d  &     %4.0f      &     %4.0f      & %4.2e &   %4.2e   & %5.2f     \\\\\n',...
      res{i}.name, res{i}.cells, res{i}.zmin, res{i}.zmax, res{i}.volume, res{i}.surfarea, ...
      res{i}.avgthk );
   end
   fprintf('------------------------------------------------\n');

    %% Compare three top surface elevation grids (taken from resTiltUtsira.m):
    %% Establishing grid interpolants
    [FS_ieaghg, FS_original, FS_inhouse] = createGridInterpolants(Gt_ieaghg, Gt_original, Gt_inhouse);
    
    %% Determining the sample points, based on the domain of INHOUSE grid
    GtS = Gt_inhouse;
    
    sizeS     = GtS.cartDims;
    get_range = @(d, pts, n) linspace(min(pts(:, d)), max(pts(:, d)), n);
    xrange    = get_range(1, GtS.nodes.coords, sizeS(1));
    yrange    = get_range(2, GtS.nodes.coords, sizeS(2));
    xdom      = xrange(floor(sizeS(1)*0.1):floor(sizeS(1)*0.9)); % Shrink domain 
    ydom      = yrange(floor(sizeS(2)*0.1):floor(sizeS(2)*0.9)); % slightly
    [x, y]    = meshgrid(xdom, ydom); % This will be the sampling grid
    
    %% Sampling the grids, and recentering them vertically at zero depth (optional)
    smpl_S_ieaghg      = FS_ieaghg(x, y); 
    smpl_S_original    = FS_original(x,y);
    smpl_S_inhouse     = FS_inhouse(x, y);
    
    smpl_S_ieaghg_mean   = mean(smpl_S_ieaghg(isfinite(smpl_S_ieaghg(:))));
    smpl_S_original_mean = mean(smpl_S_original(isfinite(smpl_S_original(:)))); 
    smpl_S_inhouse_mean  = mean(smpl_S_inhouse(isfinite(smpl_S_inhouse(:)))); 

    if recenterGrids
        smpl_S_ieaghg      = smpl_S_ieaghg - smpl_S_ieaghg_mean;
        smpl_S_original    = smpl_S_original - smpl_S_original_mean;
        smpl_S_inhouse     = smpl_S_inhouse - smpl_S_inhouse_mean;
        title_str = 'Recentered Surfaces';
    else
        smpl_S_ieaghg_mean   = 0;
        smpl_S_original_mean = 0; 
        smpl_S_inhouse_mean  = 0; 
        title_str = 'Surfaces';
    end
    
    %% Define difference surface interpolant (wrt INHOUSE grid)
    FD_ieaghg_inhouse = @(x, y) FS_ieaghg(x,y) - FS_inhouse(x,y) - (smpl_S_ieaghg_mean - smpl_S_inhouse_mean);
    FD_original_inhouse = @(x, y) FS_original(x,y) - FS_inhouse(x,y) - (smpl_S_original_mean - smpl_S_inhouse_mean);
    
    %% Display sampled grids, and difference surface
    figure; set(gcf,'Position',[1 1 2000 1000]);
    
    % IEAGHG wrt INHOUSE:
    subplot(2,2,1); 
    hold on;
    mesh(x, y, smpl_S_ieaghg, 'EdgeColor', mymap(:,1));
    mesh(x, y, smpl_S_inhouse, 'EdgeColor', mymap(:,3));
    title([title_str ': IEAGHG/INHOUSE']);
    legend('IEAGHG','INHOUSE')
    set(gca,'zdir', 'reverse'); view(-62,22);
    
    subplot(2,2,2); 
    mesh(x,y, FD_ieaghg_inhouse(x,y));
    title('Difference surface, IEAGHG/INHOUSE');
    set(gca, 'zdir', 'reverse'); view(-64,60);

    
    % ORIGINAL wrt INHOUSE:
    subplot(2,2,3); 
    hold on;
    mesh(x, y, smpl_S_original, 'EdgeColor', mymap(:,2));
    mesh(x, y, smpl_S_inhouse, 'EdgeColor', mymap(:,3));
    title([title_str ': ORIGINAL/INHOUSE']);
    legend('ORIGINAL','INHOUSE')
    set(gca,'zdir', 'reverse'); view(-62,22);
    
    subplot(2,2,4); 
    mesh(x,y, FD_original_inhouse(x,y));
    title('Difference surface, ORIGINAL/INHOUSE');
    set(gca, 'zdir', 'reverse'); view(-64,60);

    %% Get ready for next refinement case 
    clearvars -except numRefCases k recenterGrids report

end



 %% Create table to compare grids:
 
fprintf('\n\n--------------SUMMARY-----------------------\n');
fprintf('%-10s&   Cells  & Min Surf Ele. & Max Surf Ele. &  Volume  &  Surf. Area  & Avg. Thk. \\\\\n', 'Name');
for i = 1:3
    for k = 1:numel(numRefCases)
       fprintf('%-10s&  %6d  &     %4.0f      &     %4.0f      & %4.2e &   %4.2e   & %5.2f     \\\\\n',...
       report{k,i}.name, report{k,i}.cells, report{k,i}.zmin, report{k,i}.zmax, ...
       report{k,i}.volume, report{k,i}.surfarea, report{k,i}.avgthk );
    end
end
fprintf('------------------------------------------------\n');


fprintf('\n\n--------------SUMMARY----------------------------------- Node-based ------------ Cell-based ------------------\n');
fprintf('%-10s& Refined & Cells  & Min  & Max  & Volume   & Capacity  & Percent &  Capacity & Percent\\\\\n', 'Name');
for i = 1:3
    for k = 1:numel(numRefCases)
        fprintf('%-10s&   %2d    & %6d & %4.0f & %4.0f & %4.2e & %4.2e  & %5.2f   & %4.2e  & %5.2f \\\\\n',...
          report{k,i}.name, report{k,i}.refLevel, report{k,i}.cells, ...
          report{k,i}.zmin, report{k,i}.zmax, report{k,i}.volume, ...
          report{k,i}.ncapacity, report{k,i}.ncapacity/report{k,i}.volume*100, ...
          report{k,i}.ccapacity, report{k,i}.ccapacity/report{k,i}.volume*100);
    end
end
fprintf('------------------------------------------------\n');
  
fprintf('\n\n--------------SUMMARY--------------------------------- Node-based --------------------------------------- Cell-based ---------------------\n');
fprintf('%-10s& Refined & Num. global traps & Tot. trap vol. (m3) & Avg. global trap vol. (m3) & Num. global traps & Tot. trap vol. (m3) & Avg. global trap vol. (m3)\\\\\n', 'Name');
for i = 1:3
    for k = 1:numel(numRefCases)
        fprintf('%-10s&   %2d    &     %6d        &       %4.2e      &         %4.2e           &     %6d        &       %4.2e      &    %4.2e       \\\\\n',...
          report{k,i}.name, report{k,i}.refLevel, ...
          numel(report{k,i}.ntrapvols), ...
          sum(report{k,i}.ntrapvols), ...
          mean(report{k,i}.ntrapvols), ...
          numel(report{k,i}.ctrapvols), ...
          sum(report{k,i}.ctrapvols), ...
          mean(report{k,i}.ctrapvols) );
    end
end
fprintf('------------------------------------------------\n');


%% Create table to study grid specs:
% TODO: put this into a table...

% bounds of IEAGHG geological model:
fprintf('\n IEAGHG geological model ... \n')
xmin3ls = min(grid3ls_ieaghg.G.nodes.coords(:,1));
xmax3ls = max(grid3ls_ieaghg.G.nodes.coords(:,1));
ymin3ls = min(grid3ls_ieaghg.G.nodes.coords(:,2));
ymax3ls = max(grid3ls_ieaghg.G.nodes.coords(:,2));
zmin3ls = min(grid3ls_ieaghg.G.nodes.coords(:,3));
zmax3ls = max(grid3ls_ieaghg.G.nodes.coords(:,3));
avthk3ls = mean(grid3ls_ieaghg.Gt.cells.H);
fprintf('\n The bounds of the grid are (x1,y1)=( %d , %d ) to (x2,y2)=( %d, %d ). \n', ...
    xmin3ls, ymin3ls, xmax3ls, ymax3ls);
fprintf('\n The average thickness is %d, ranging from %d to %d meters. \n', ...
    avthk3ls, zmin3ls, zmax3ls);
fprintf('\n The areal coverage is %d by %d = %d km^2. \n', ...
    (xmax3ls-xmin3ls)/1000, (ymax3ls-ymin3ls)/1000, ((xmax3ls-xmin3ls)*(ymax3ls-ymin3ls))/(1000*1000));

zmin = min(G_ieaghg.nodes.coords(:,3));
zmax = max(G_ieaghg.nodes.coords(:,3));
avthk = mean(Gt_ieaghg.cells.H);
fprintf('\n After removing caprock and bottom shale:\n')
fprintf('\n The average thickness is %d, ranging from %d to %d meters. \n', ...
    avthk, zmin, zmax);

avpermx = mean(rock_ieaghg.perm(:,1));
avpermy = mean(rock_ieaghg.perm(:,2));
avpermz = mean(rock_ieaghg.perm(:,3));

permxmin = min(rock_ieaghg.perm(:,1));
permxmax = max(rock_ieaghg.perm(:,1));
permymin = min(rock_ieaghg.perm(:,2));
permymax = max(rock_ieaghg.perm(:,2));
permzmin = min(rock_ieaghg.perm(:,3));
permzmax = max(rock_ieaghg.perm(:,3));
fprintf('\n The average horizontal perm is %d meters^2 (%d Darcy), \n ranging from %d (%d) to %d (%d). \n', ...
    avpermx, avpermx/convertFrom(1, darcy()), ...
    permxmin, permxmin/convertFrom(1, darcy()), ...
    permxmax, permxmax/convertFrom(1, darcy()) );
fprintf('\n The average vertical perm is %d meters^2 (%d Darcy), \n ranging from %d (%d) to %d (%d). \n', ...
    avpermz, avpermz/convertFrom(1, darcy()), ...
    permzmin, permzmin/convertFrom(1, darcy()), ...
    permzmax, permzmax/convertFrom(1, darcy()) );


% bounds of ORIGINAL geological model:
fprintf('\n ORIGINAL geological model ... \n')
xmin3ls = min(grid3ls_original.G.nodes.coords(:,1));
xmax3ls = max(grid3ls_original.G.nodes.coords(:,1));
ymin3ls = min(grid3ls_original.G.nodes.coords(:,2));
ymax3ls = max(grid3ls_original.G.nodes.coords(:,2));
zmin3ls = min(grid3ls_original.G.nodes.coords(:,3));
zmax3ls = max(grid3ls_original.G.nodes.coords(:,3));
avthk3ls = mean(grid3ls_original.Gt.cells.H);
fprintf('\n The bounds of the grid are (x1,y1)=( %d , %d ) to (x2,y2)=( %d, %d ). \n', ...
    xmin3ls, ymin3ls, xmax3ls, ymax3ls);
fprintf('\n The average thickness is %d, ranging from %d to %d meters. \n', ...
    avthk3ls, zmin3ls, zmax3ls);
fprintf('\n The areal coverage is %d by %d = %d km^2. \n', ...
    (xmax3ls-xmin3ls)/1000, (ymax3ls-ymin3ls)/1000, ((xmax3ls-xmin3ls)*(ymax3ls-ymin3ls))/(1000*1000));

zmin = min(G_original.nodes.coords(:,3));
zmax = max(G_original.nodes.coords(:,3));
avthk = mean(Gt_original.cells.H);
fprintf('\n After removing caprock and bottom shale:\n')
fprintf('\n The average thickness is %d, ranging from %d to %d meters. \n', ...
    avthk, zmin, zmax);

avpermx = mean(rock_original.perm(:,1));
avpermy = mean(rock_original.perm(:,2));
avpermz = mean(rock_original.perm(:,3));

permxmin = min(rock_original.perm(:,1));
permxmax = max(rock_original.perm(:,1));
permymin = min(rock_original.perm(:,2));
permymax = max(rock_original.perm(:,2));
permzmin = min(rock_original.perm(:,3));
permzmax = max(rock_original.perm(:,3));
fprintf('\n The average horizontal perm is %d meters^2 (%d Darcy), \n ranging from %d (%d) to %d (%d). \n', ...
    avpermx, avpermx/convertFrom(1, darcy()), ...
    permxmin, permxmin/convertFrom(1, darcy()), ...
    permxmax, permxmax/convertFrom(1, darcy()) );
fprintf('\n The average vertical perm is %d meters^2 (%d Darcy), \n ranging from %d (%d) to %d (%d). \n', ...
    avpermz, avpermz/convertFrom(1, darcy()), ...
    permzmin, permzmin/convertFrom(1, darcy()), ...
    permzmax, permzmax/convertFrom(1, darcy()) );
