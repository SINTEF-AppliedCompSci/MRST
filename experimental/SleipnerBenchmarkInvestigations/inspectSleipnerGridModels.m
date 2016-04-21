function [ hfig1, hfig2, hfig3, hfigA, hfigB, hfigC, hfigA2, hfigB2, hfigC2  ] = inspectSleipnerGridModels( varargin )
% Comparison of all three Sleipner grid models, given specified refinement
% levels. The refinement level(s) (i.e., -2, -3, etc for coarsening; 2, 3,
% etc for refinement) are defined, and each refinement case is analyzed
% separately. Results are reported in table displayed in command window and
% with plots.

opt.refineLevels    = -6;
opt.add2008Plume    = true;
opt.addlegend       = true;
opt = merge_options(opt, varargin{:});

% Refinement cases:
numRefCases = opt.refineLevels; %[-4, -3, -2, 1, 2, 3, 4];

% Do you want to do a lot of analysis of the top-surfaces comparison?
moreAnalysis = false;

% pre-allocate report cell
report = cell( numel(numRefCases), 3);
for k = 1:numel(numRefCases)
    
    % get current refinement case
    numRef = numRefCases(k); 

    %% Create Sleipner model grids for each of the three types
    [ G_ieaghg, Gt_ieaghg, rock_ieaghg, rock2D_ieaghg ]          = makeSleipnerModelGrid( 'modelName','IEAGHGmodel','refineLevel', numRef);
    [ G_original, Gt_original, rock_original, rock2D_original ]  = makeSleipnerModelGrid( 'modelName','ORIGINALmodel','refineLevel', numRef );
    %[ G_inhouse, Gt_inhouse, rock_inhouse, rock2D_inhouse ]      = makeSleipnerModelGrid( 'modelName','INHOUSEmodel','refineLevel', numRef );
    [ G_seismic, Gt_seismic, rock_seismic, rock2D_seismic ]      = makeSeismicModelGrid( numRef );
    
    % add grid names
    Gt_ieaghg.name = 'IEAGHG';
    Gt_original.name = 'GHGT';
    Gt_seismic.name = 'SEISMIC';
    
    
    %% Visualize Grids
    % use bounds of GHGT grid for consistency of plotting limits
    G = G_original;
    bounds.xmin = min(G.nodes.coords(:,1)); bounds.xmax = max(G.nodes.coords(:,1));
    bounds.ymin = min(G.nodes.coords(:,2)); bounds.ymax = max(G.nodes.coords(:,2));
    bounds.zmin = min(G.nodes.coords(:,3)); bounds.zmax = max(G.nodes.coords(:,3));
    
    [ hfig1 ] = makeGridPlot(G_original, bounds);
    [ hfig2 ] = makeGridPlot(G_ieaghg,   bounds);
    [ hfig3 ] = makeGridPlot(G_seismic,  bounds);
    
    %% Compare top surfaces of grids
    [ hfigA, hfigB, hfigC ] = compareSurfaces(Gt_ieaghg, Gt_original, opt.add2008Plume, opt.addlegend);
    [ hfigA2, hfigB2, hfigC2 ] = compareSurfaces(Gt_seismic, Gt_original, opt.add2008Plume, opt.addlegend);
    %[ hfigA2, hfigB2, hfigC2 ] = compareSurfaces(Gt_ieaghg, Gt_seismic, opt.add2008Plume, opt.addlegend);

    
    %% Other (todo)
    %gts = [Gt_ieaghg; Gt_original; Gt_inhouse];
    %gs  = [G_ieaghg; G_original; G_inhouse];

    %names = {'IEAGHG'; 'GHGT'; 'INHOUSE'};

    
%     %% Inspect Top Surface Grid and H
%     figure; set(gcf,'Position',[1 1 1200 800])
% 
%     % pre-allocate res cell (for current numRef)
%     res = cell( 1, numel(gts));
%     for i = 1:numel(gts)
%         
%         Gt = gts(i);
%         G  = gs(i);
%         name = names{1};
%         avthk = mean(Gt.cells.H);
%         numCells = Gt.cells.num;
% 
%         plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', mymap(i,:));
% 
%         % Comparisons of grid details
%         res{i}.name      = names{i};
%         res{i}.refLevel  = numRef;
%         res{i}.cells     = Gt.cells.num;
%         res{i}.zmin      = min(Gt.cells.z);
%         res{i}.zmax      = max(Gt.cells.z);
%         res{i}.volume    = sum(G.cells.volumes);
%         res{i}.surfarea  = sum(Gt.cells.volumes);
%         res{i}.avgthk    = mean(Gt.cells.H);
%         
%         % trapping analysis
%         tan     = trapAnalysis(Gt, false);
%         tac     = trapAnalysis(Gt, true);
%         
%         res{i}.ctrapvols = volumesOfTraps(Gt,tac);
%         res{i}.ccapacity = sum(res{i}.ctrapvols);
%         res{i}.ntrapvols = volumesOfTraps(Gt,tan);
%         res{i}.ncapacity = sum(res{i}.ntrapvols);
%         
%         % store for later use
%         report{k,i} = res{i};
%     end
% 
%     legend(cellfun(@(x) x.name, res, 'UniformOutput', false), 'Location', 'EastOutside')
% 
%     view(3)
%     set(gca,'DataAspect',[1 1 0.02]); grid
%     xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
%     title('Top Surfaces elevation')

    %% Create table to compare grids:
    
   % TODO: fix error when a particular grid's results are empty.
%    fprintf('\n\n--------------Num Ref. =  %-2d-----------------------\n', numRef);
%    fprintf('%-20s&   Cells  & Min Surf Ele. & Max Surf Ele. &  Volume  &  Surf. Area  & Avg. Thk. \\\\\n', 'Name');
% 
%    for i = 1:3
%    fprintf('%-20s&  %6d  &     %4.0f      &     %4.0f      & %4.2e &   %4.2e   & %5.2f     \\\\\n',...
%       res{i}.name, res{i}.cells, res{i}.zmin, res{i}.zmax, res{i}.volume, res{i}.surfarea, ...
%       res{i}.avgthk );
%    end
%    fprintf('------------------------------------------------\n');


    %% Get ready for next refinement case 
    clearvars -except numRefCases k recenterGrids report moreAnalysis hfig1 hfig2 hfig3 hfigA hfigB hfigC hfigA2 hfigB2 hfigC2
end


%% Do more analysis if required.
if moreAnalysis

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


    %% Some more analysis on base grids.
    % TODO: put this into a table...
    clear layers cut

    % bounds of IEAGHG geological model:
    load SleipnerGlobalCoords_baseGrids.mat
    fprintf('\n IEAGHG geological model ... \n')
    xmin3ls = min(layers.G.nodes.coords(:,1));
    xmax3ls = max(layers.G.nodes.coords(:,1));
    ymin3ls = min(layers.G.nodes.coords(:,2));
    ymax3ls = max(layers.G.nodes.coords(:,2));
    zmin3ls = min(layers.G.nodes.coords(:,3));
    zmax3ls = max(layers.G.nodes.coords(:,3));
    avthk3ls = mean(layers.Gt.cells.H);
    fprintf('\n The bounds of the grid are (x1,y1)=( %d , %d ) to (x2,y2)=( %d, %d ). \n', ...
        xmin3ls, ymin3ls, xmax3ls, ymax3ls);
    fprintf('\n The average thickness is %d, ranging from %d to %d meters. \n', ...
        avthk3ls, zmin3ls, zmax3ls);
    fprintf('\n The areal coverage is %d by %d = %d km^2. \n', ...
        (xmax3ls-xmin3ls)/1000, (ymax3ls-ymin3ls)/1000, ((xmax3ls-xmin3ls)*(ymax3ls-ymin3ls))/(1000*1000));

    zmin = min(cut.G.nodes.coords(:,3));
    zmax = max(cut.G.nodes.coords(:,3));
    avthk = mean(cut.Gt.cells.H);
    fprintf('\n After removing caprock and bottom shale:\n')
    fprintf('\n The average thickness is %d, ranging from %d to %d meters. \n', ...
        avthk, zmin, zmax);

    avpermx = mean(cut.rock.perm(:,1));
    avpermy = mean(cut.rock.perm(:,2));
    avpermz = mean(cut.rock.perm(:,3));

    permxmin = min(cut.rock.perm(:,1));
    permxmax = max(cut.rock.perm(:,1));
    permymin = min(cut.rock.perm(:,2));
    permymax = max(cut.rock.perm(:,2));
    permzmin = min(cut.rock.perm(:,3));
    permzmax = max(cut.rock.perm(:,3));
    fprintf('\n The average horizontal perm is %d meters^2 (%d Darcy), \n ranging from %d (%d) to %d (%d). \n', ...
        avpermx, avpermx/convertFrom(1, darcy()), ...
        permxmin, permxmin/convertFrom(1, darcy()), ...
        permxmax, permxmax/convertFrom(1, darcy()) );
    fprintf('\n The average vertical perm is %d meters^2 (%d Darcy), \n ranging from %d (%d) to %d (%d). \n', ...
        avpermz, avpermz/convertFrom(1, darcy()), ...
        permzmin, permzmin/convertFrom(1, darcy()), ...
        permzmax, permzmax/convertFrom(1, darcy()) );


    % bounds of ORIGINAL geological model:
    load OriginalSleipnerGlobalCoords_baseGrids.mat
    fprintf('\n ORIGINAL geological model ... \n')
    xmin3ls = min(layers.G.nodes.coords(:,1));
    xmax3ls = max(layers.G.nodes.coords(:,1));
    ymin3ls = min(layers.G.nodes.coords(:,2));
    ymax3ls = max(layers.G.nodes.coords(:,2));
    zmin3ls = min(layers.G.nodes.coords(:,3));
    zmax3ls = max(layers.G.nodes.coords(:,3));
    avthk3ls = mean(layers.Gt.cells.H);
    fprintf('\n The bounds of the grid are (x1,y1)=( %d , %d ) to (x2,y2)=( %d, %d ). \n', ...
        xmin3ls, ymin3ls, xmax3ls, ymax3ls);
    fprintf('\n The average thickness is %d, ranging from %d to %d meters. \n', ...
        avthk3ls, zmin3ls, zmax3ls);
    fprintf('\n The areal coverage is %d by %d = %d km^2. \n', ...
        (xmax3ls-xmin3ls)/1000, (ymax3ls-ymin3ls)/1000, ((xmax3ls-xmin3ls)*(ymax3ls-ymin3ls))/(1000*1000));

    zmin = min(cut.G.nodes.coords(:,3));
    zmax = max(cut.G.nodes.coords(:,3));
    avthk = mean(cut.Gt.cells.H);
    fprintf('\n After removing caprock and bottom shale:\n')
    fprintf('\n The average thickness is %d, ranging from %d to %d meters. \n', ...
        avthk, zmin, zmax);

    avpermx = mean(cut.rock.perm(:,1));
    avpermy = mean(cut.rock.perm(:,2));
    avpermz = mean(cut.rock.perm(:,3));

    permxmin = min(cut.rock.perm(:,1));
    permxmax = max(cut.rock.perm(:,1));
    permymin = min(cut.rock.perm(:,2));
    permymax = max(cut.rock.perm(:,2));
    permzmin = min(cut.rock.perm(:,3));
    permzmax = max(cut.rock.perm(:,3));
    fprintf('\n The average horizontal perm is %d meters^2 (%d Darcy), \n ranging from %d (%d) to %d (%d). \n', ...
        avpermx, avpermx/convertFrom(1, darcy()), ...
        permxmin, permxmin/convertFrom(1, darcy()), ...
        permxmax, permxmax/convertFrom(1, darcy()) );
    fprintf('\n The average vertical perm is %d meters^2 (%d Darcy), \n ranging from %d (%d) to %d (%d). \n', ...
        avpermz, avpermz/convertFrom(1, darcy()), ...
        permzmin, permzmin/convertFrom(1, darcy()), ...
        permzmax, permzmax/convertFrom(1, darcy()) );

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%      Local functions - for plotting    %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ hfig ] = makeGridPlot(G, bounds)


    figure; set(gcf,'Position',[1 1 1000 800])
    %plotGrid(G, 'FaceColor', 'y', 'EdgeColor', [.8 .8 .8]);
    plotGrid(G, 'FaceColor', [1 .9 .9], 'EdgeColor', 'none');% 'EdgeAlpha', .05);
    %plotCellData(G, G.cells.centroids(:,3))
    %plotGrid(Gt, 'FaceColor', 'y', 'EdgeColor', 'k');
    %plotCellData(Gt, Gt.cells.z, 'EdgeColor', [.8 .8 .8]);
    light('Position',[-1 -1 1],'Style','infinite');lighting phong

    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x','Fontweight','bold'); ylabel('y','Fontweight','bold'); zlabel('depth (meters)','Fontweight','bold');
    % Create textarrow
    view(44,22)
    %annotation(gcf,'textarrow',[0.4734 0.5391],[0.7825 0.81],'String',{'North'});
    set(gca,'FontSize',22); box
    xlim([bounds.xmin bounds.xmax]);
    ylim([bounds.ymin bounds.ymax]);
    zlim([bounds.zmin bounds.zmax]);

    % get lateral boundary faces:
    % see gridFactoryExamples.m for grid details.
    % see also boundaryFaceIndices.m
    boundary        = any(G.faces.neighbors==0,2);
    facelist        = 1:G.faces.num;
    bdryfaces       = facelist( boundary);
    %internalfaces   = facelist(~boundary);
    ef = G.cells.faces(:,2)==1;
    wf = G.cells.faces(:,2)==2;
    sf = G.cells.faces(:,2)==3;
    nf = G.cells.faces(:,2)==4;
    keepfaces = ef + wf + sf + nf;
    tmp = G.cells.faces(:,1);
    keepfaces = unique( tmp(logical(keepfaces)) );

    latbdryfaces = intersect(bdryfaces,keepfaces); % excludes top and bottom bdry faces
    plotFaces(G, latbdryfaces, 'FaceColor',[.8 .8 .8], 'EdgeColor',[.8 .5 0], 'EdgeAlpha', .05)

    hfig = gcf;
end


function [ hfigA, hfigB, hfigC ] = compareSurfaces(Gt1, Gt2, add2008Plume, addlegend)
    % Compare top surface elevation grids (taken from resTiltUtsira.m):
    % Gt1 is the grid used for sampling points
    mymap = [ 1 0 0; 0 1 0 ]; % for coloring grids
    %mymap = [ 0.5 0 0.9; 0 0.8 0.8; 0 0 1 ];

    
    %% Establishing grid interpolants
    [FS1, FS2] = createGridInterpolants(Gt1, Gt2);
    
    %% Determining the sample points (using smaller grid)
    GtS = Gt1; % todo: ensure this is the smaller grid.
    
    sizeS     = GtS.cartDims;
    get_range = @(d, pts, n) linspace(min(pts(:, d)), max(pts(:, d)), n);
    xrange    = get_range(1, GtS.nodes.coords, sizeS(1));
    yrange    = get_range(2, GtS.nodes.coords, sizeS(2));
    xdom      = xrange(floor(sizeS(1)*0.1):floor(sizeS(1)*0.9)); % Shrink domain 
    ydom      = yrange(floor(sizeS(2)*0.1):floor(sizeS(2)*0.9)); % slightly
    [x, y]    = meshgrid(xdom, ydom); % This will be the sampling grid
    
    
    %% Compute and Display sampled grids, and difference surface
    % using difference surface interpolant FD
    hfigA = figure; %set(gcf,'Position',[1 1 2000 500]);
    
    % Sampling the grids:
    smpl_S1   = FS1(x,y); 
    smpl_S2   = FS2(x,y);

    % First, plot original surface elevations:
    smpl_S1_mean = 0;
    smpl_S2_mean = 0; 
    title_str = {'Surfaces';' '};
    FD = @(x, y) FS1(x,y) - FS2(x,y) - (smpl_S1_mean - smpl_S2_mean);
    
    %subplot(1,3,1);
    hold on;
    mesh(x, y, smpl_S1, 'EdgeColor', mymap(1,:)); strname1 = Gt1.name;
    mesh(x, y, smpl_S2, 'EdgeColor', mymap(2,:)); strname2 = Gt2.name;
    
    % adjust plot
    set(gca,'zdir', 'reverse'); view(24, 34); %view(-40,40);
    axis equal tight; set(gca,'DataAspect',[1 1 0.02]);
    grid; box
    
    % add title
    %title(title_str, 'FontSize',22);
    
    % add legend
    hl = legend(strname1,strname2,'Location','nw'); set(hl,'FontSize',22);

    % adjust axis ticks and fontsize
    hax = gca;
    set(hax,'YTick',[6470000 6472000 6474000])
    set(hax,'XTick',[437000 438000 439000])
    set(hax,'ZTick',[800 815 830])
    set(hax,'FontSize',16);
    
    
    % Second, get recenterGrids:
    smpl_S1_mean = mean(smpl_S1(isfinite(smpl_S1(:))));
    smpl_S2_mean = mean(smpl_S2(isfinite(smpl_S2(:)))); 
    
    smpl_S1    = smpl_S1 - smpl_S1_mean;
    smpl_S2    = smpl_S2 - smpl_S2_mean;

    title_str = {'Recentered Surfaces';' '};
    FD = @(x, y) FS2(x,y) - FS1(x,y) - (smpl_S2_mean - smpl_S1_mean);
    
    hfigB = figure;
    %subplot(1,3,2);
    hold on;
    mesh(x, y, smpl_S1, 'EdgeColor', mymap(1,:)); strname1 = Gt1.name;
    mesh(x, y, smpl_S2, 'EdgeColor', mymap(2,:)); strname2 = Gt2.name;
    
    % adjust plot
    set(gca,'zdir', 'reverse'); view(24, 34); %view(-40,40);
    axis equal tight; set(gca,'DataAspect',[1 1 0.02]);
    grid; box
    
    % add title
    %title(title_str, 'FontSize',22);
    
    % add legend
    hl = legend(strname1,strname2,'Location','nw'); set(hl,'FontSize',22);
    
    % adjust axis ticks and fontsize
    hax = gca;
    set(hax,'YTick',[6470000 6472000 6474000])
    set(hax,'XTick',[437000 438000 439000])
    set(hax,'ZTick',[-10 0 10 20])
    set(hax,'FontSize',16);
    
    
    % Lastly, plot recentered grids elevation difference:
    hfigC = figure;
    %subplot(1,3,3);
    surf(x,y, FD(x,y), 'Edgecolor','none');
    
    % adjust plot
    set(gca, 'zdir', 'reverse'); %view(-64,60);
    view(2); axis equal tight;
    
    % add title
    %title({'Difference';' '},'FontSize',22);
    
    % adjust axis ticks and fontsize
    hax = gca;
    ymin = min(min(y)); ymax = max(max(y));
    xmin = min(min(x)); xmax = max(max(x));
    set(hax,'YTick',[ymin (ymax-ymin)/2+ymin ymax])
    set(hax,'XTick',[xmin xmax])
    set(hax,'FontSize',16);
    
    % add colorbar. adjust label and fontsize
    hcb = setColorbarHandle(gcf, 'fontSize',22, 'LabelName','meters');
    
    if add2008Plume
        % (optional): superimposed 2008 plume outline onto last subplot
        plumes = getLayer9CO2plumeOutlines();
        Year2plot = 2008;
        for j = 1:numel(plumes)
            if plumes{j}.year == Year2plot
                disp('Plotting Observed CO2 plume outline...')
                hp = line(plumes{j}.outline(:,1), plumes{j}.outline(:,2), -150*ones(numel(plumes{j}.outline(:,2)),1), 'LineWidth',2, 'Color','r');
            end
        end
        % add legend (optional)
        if addlegend
            hl = legend([hp],{'2008 plume'}, 'Location','nw');
            set(hl,'FontSize',22)
        end
    end
    
end

function [FS1, FS2] = createGridInterpolants(Gt1, Gt2)
% Prepare and return grid interpolants
%  (to facilitate comparison, since grids are of different resolution) A
%  scattered-point interpolant is useful when the domain of the formation
%  has an irregular shape.

    FS1 = TriScatteredInterp(Gt1.nodes.coords(:,1), Gt1.nodes.coords(:,2), ...
                            Gt1.nodes.z);
    
    FS2 = TriScatteredInterp(Gt2.nodes.coords(:,1), Gt2.nodes.coords(:,2), ...
                            Gt2.nodes.z);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%          End of Local Functions        %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%