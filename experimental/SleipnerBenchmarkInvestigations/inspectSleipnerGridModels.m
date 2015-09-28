% Comparison of all three Sleipner grid models, given specified refinement
% levels. The refinement level(s) (i.e., -2, -3, etc for coarsening; 2, 3,
% etc for refinement) are defined, and each refinement case is analyzed
% separately. Results are reported in table displayed in command window and
% with plots.

% Refinement cases:
numRefCases = [-6 ]; %-4 -2 1 2 4];

% Do you want top surface grids to be re-centered when making surface
% comparisons? (true for recentering, false for no recentering)
recenterGrids = false;

% pre-allocate report cell
report = cell( numel(numRefCases), 3);
for k = 1:numel(numRefCases)
    
    % get current refinement case
    numRef = numRefCases(k);
    
    %% Create Sleipner model grids for each of the three types
    [ G_ieaghg, Gt_ieaghg, rock_ieaghg, rock2D_ieaghg ]         = makeSleipnerModelGrid( 'modelName','IEAGHGmodel','refineLevel', numRef);
    [ G_original, Gt_original, rock_original, rock2D_original ] = makeSleipnerModelGrid( 'modelName','ORIGINALmodel','refineLevel', numRef );
    [ G_inhouse, Gt_inhouse, rock_inhouse, rock2D_inhouse ]     = makeSleipnerModelGrid( 'modelName','INHOUSEmodel','refineLevel', numRef );

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
  
