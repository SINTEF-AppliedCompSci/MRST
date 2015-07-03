%% CO2 Storage Atlas: the Barents Sea and Norwegian Sea
% Similar to showCO2atlas.m, this script analyzes the formations present in
% the CO2 Atlas Seas (Barents and Norwegian), and estimates its CO2 storage
% capacity.

moduleCheck('co2lab', 'opm_gridprocessing');


%% 1. Load CO2 atlas formations of a specified Sea.
% This step assumes that the formation datasets of the Barents Sea and/or
% Norwegian Sea have already been downloaded (from NPD) and are stored in
% mrst-core/examples/data/CO2Atlas (with file extensions removed). Note
% that the filenames 'bcu' stands for base cretaceous unconformity, and
% they should be appropriately renamed 'bcu_BarentsSea' and
% 'bcu_NorwegianSea'.

% TODO: put Barents Sea and Norwegian Sea datasets on SINTEF server and
% modify downloadDataSets() function in order to download these
% formations in addition to the current North Sea formation datasets.

% Since some formations are only given in terms of top data (not
% thickness), we get the rawdata of the altas formations for analysis and
% visualization. First, select the sea you want to study, i.e., 'Barents',
% 'Norwegian', or 'NorwegianNorth'.

studySea = 'Barents';
makeFigures = false;
saveFigures = true;

if strcmpi(studySea,'Barents')
    names = {'tubЖen', 'stЫ', 'nordmela', 'knurr', 'fruholmen', ...
    'bjarmeland', 'bcuBarentsSea'};

elseif strcmpi(studySea,'Norwegian')
    names = {'Пre', 'ror', 'not', 'ile', 'garn', 'bcuNorwegianSea', ...
    'baseПre', 'tilje'};

elseif strcmpi(studySea,'NorwegianNorth')
    %TODO: or simply run showCO2atlas.m
    fprintf('Use showCO2atlas.m to study North Sea formations.\n')
    return
    
else
    fprintf('Sea not specified.\n')
    return
end

fprintf('Loading %-14s Sea atlas data ...', studySea);
[grdecls, rawdata] = getAtlasGrid(names);
fprintf('done\n');


%% Description of raw data
% Show the raw data. Each dataset contains four fields:
% - Name, which is the name of the formation
% - Variant: Either thickness or height, indicating wether the dataset
% represents height data or thickness data.
% - Data: The actual datasets as a matrix.
% - Meta: Metadata. The most interesting field here is the
%   xllcorner/yllcorner variable which indicates the position in ED50
%   datum space.

fprintf('\nRaw data of %-14s Sea:\n', studySea)
fprintf('----------------------------------------------------------------\n');
for i=1:numel(rawdata);
    rd = rawdata{i};
    fprintf('Dataset %-2i is %-12s (%-9s). Resolution: %4i meters\n', ...
            i, rd.name, rd.variant,  rd.meta.cellsize)
end
fprintf('----------------------------------------------------------------\n\n');

% Store names for convenience (required for correct name ordering)
names = cellfun(@(x) x.name, rawdata, 'UniformOutput', false)';


%%
% Visualize formations and save into folder, only if makeFigures is true
% TODO: remove this part and visualize formation top surfaces in next part.
if makeFigures
    
    folderName = 'CO2AtlasFormations_Plots';
    if ~exist(folderName, 'dir')
        mkdir(folderName)
    end

    for i = 1:numel(names)

        frd = rawdata{1,i};
        clf; set(gcf,'Position', [100 100 1000 1000]);
        surf(frd.data'); shading interp;

        % make adjustments to figure, such as axes labels
        title([frd.name ' ' frd.variant])
        view(0,90); %axis tight
        set(gca,'DataAspect',[1 1 frd.meta.cellsize/10])
        ax = gca;
        ax.XTickLabel = ax.XTick*frd.meta.cellsize/1000;
        ax.YTickLabel = ax.YTick*frd.meta.cellsize/1000;
        ax.ZTickLabel = ax.ZTick/1000;
        hold off
        xlabel('x, km'); ylabel('y, km'); zlabel('z, km');
        
        hcb = colorbar;
        hcb.TickLabels = hcb.Ticks/-1000;
        hcb.Label.String = 'Depth below sealevel, km';
        
        if saveFigures
            % save figure
            saveas(gcf, [folderName '/' frd.name], 'png')
        end
        clear frd

    end

end


%% Hammerfest Basin Aquifer:
% Next, we visualize the "Hammerfest Basin Aquifer", which is made up of
% three formations present in the Barents Sea: Tubaen, Sto, and Nordmela
% (see NPD's CO2 Storage Atlas report, chp 6, found at 
% <http://www.npd.no/en/Publications/Reports/Compiled-CO2-atlas/>).

names_Ham = {'tubЖen', 'stЫ', 'nordmela'};

if ~strcmpi(studySea,'Barents')
    fprintf('Loading formations of Hammerfest Basin Aquifer ...');
    [grdecls, rawdata] = getAtlasGrid(names_Ham);
    fprintf('done\n');
    
    fprintf('\nRaw data of Hammerfest Basin Aquifer:\n')
    fprintf('----------------------------------------------------------------\n');
    for i=1:numel(rawdata);
        rd = rawdata{i};
        fprintf('Dataset %-2i is %-12s (%-9s). Resolution: %4i meters\n', ...
                i, rd.name, rd.variant,  rd.meta.cellsize)
    end
    fprintf('----------------------------------------------------------------\n\n');
    
    % Store names for convenience (required for correct name ordering)
    names = cellfun(@(x) x.name, rawdata, 'UniformOutput', false)';
end

if makeFigures
    
    %mymap = [0 0 0; 0.9 0.9 0.9; 0.5 0.5 0.5];
    mymap = ['k'; 'r'; 'b'];
    clf; set(gcf,'Position', [100 100 2000 1000]); hold on

    [cellsize, xllcorner, yllcorner] = deal( zeros(1,numel(names_Ham)) );
    nstr = cell(1,numel(names_Ham));

    for i = 1:numel(names_Ham)

        formation_rd = rawdata(strcmpi(names, names_Ham{i}));
        frd = formation_rd{1,1};
        if ~strcmpi(frd.variant, 'top')
            continue
        end

        % store important info
        cellsize(i)  = frd.meta.cellsize;
        xllcorner(i) = frd.meta.xllcorner; % TODO: check corner alignment
        yllcorner(i) = frd.meta.yllcorner;
        nstr{i}      = frd.name;

        % make surface plot
        surf(frd.data', 'FaceColor', 'none', 'EdgeColor', mymap(i,:)); %shading interp

        clear frd formation_rd

    end

    % make adjustments to figure, such as axes labels
    view(3); axis tight; legend(nstr, 'Location', 'EastOutside')
    title('Hammerfest Aquifer in Barents Sea')

    if all(cellsize == cellsize(1))
        set(gca,'DataAspect',[1 1 cellsize(1)/10])
        ax = gca;
        [ax.XTick, ax.YTick] = deal([0 100 200 300]);
        ax.XTickLabel = ax.XTick*cellsize(1)/1000;
        ax.YTickLabel = ax.YTick*cellsize(1)/1000;
        ax.ZTickLabel = -ax.ZTick/1000;

        hold off
        xlabel('x, km'); ylabel('y, km'); zlabel('z, km');
        grid
    end

    if saveFigures
        % save plot into folder
        folderName = 'CO2AtlasFormations_Plots';
        if ~exist(folderName, 'dir')
            mkdir(folderName)
        end
        saveas(gcf, [folderName '/' 'HammerfestAquifer'], 'png')
    end
    
end


% Plot of x-section through aquifer to show three layers of formations:
% TODO
close all;




%% 3. Analyze traps with spill-point analysis, and plot structural traps.
% To create 3D grids, the thickness information that is described in the
% CO2 Storage Atlas may be used, or alternatively a hypothetical thickness
% could be set assuming CO2 plume is thin and under the formation top.


% Select formation to analyze.
i = 2;
fthk = 2; % km, an approximation (should get data from CO2 Atlas)

if strcmpi(rawdata{1,i}.variant,'top') % re-write to check for only top surface data and no thickness data

    % Set level of coarsening cases:
    Res = [1 2 3 4 5 6];
    numRes = numel(Res);
    [Grids, res] = deal(cell(numel(numRes),1));

for N = 1:numRes
    
    % coarsen the data
    coarsening = Res(N);
    
    rd = rawdata{1,i};
    
    rd.meta.cellsize = rd.meta.cellsize*coarsening;
    rd.data = rd.data(1:coarsening:end, 1:coarsening:end);  
    rd.meta.ncols = size(rd.data, 2);
    rd.meta.nrows = size(rd.data, 1);   
    rd.meta.dims = [rd.meta.nrows, rd.meta.ncols];
    

    % make formation grid
    nx = rd.meta.dims(1)-1;
    ny = rd.meta.dims(2)-1;
    nz = 1;
    Lx = (rd.meta.dims(1)-1)*rd.meta.cellsize/1000; % km
    Ly = (rd.meta.dims(2)-1)*rd.meta.cellsize/1000;
    Lz = -min(min(rd.data))/1000 + fthk; 

    G = cartGrid([nx,ny,nz], [Lx,Ly,Lz]);
    %figure; plotGrid(G); view(3); title(['Grid of ' rd.name])
    %xlabel('x, km'); ylabel('y, km'); zlabel('z (depth below sealevel), km');

    % replace z coordinate with top surface data
    z = reshape(-rd.data/1000,rd.meta.dims(1)*rd.meta.dims(2),1);
    G.nodes.coords(1:rd.meta.dims(1)*rd.meta.dims(2),3) = z;

    % project top surface data down to make bottom of formation
    G.nodes.coords(rd.meta.dims(1)*rd.meta.dims(2)+1:end,3) = z + fthk;

    % first computeGeometry, then remove NaNs from G.cells.centroids
    G = computeGeometry(G);
    G = removeCells(G, isnan(G.cells.centroids(:,3)) );
    %figure; plotGrid(G); view(3); title(['Grid of ' rd.name]); grid
    %xlabel('x, km'); ylabel('y, km'); zlabel('z, km');

    % compute topSurfaceGrid
    Gt = topSurfaceGrid(G);
    
    % compute traps with spill-point analysis, using both node and cell
    % methods:
    tan     = trapAnalysis(Gt, false);
    tac     = trapAnalysis(Gt, true);
    
    res{N}.name      = rawdata{1,i}.name;
    res{N}.cells     = Gt.cells.num;
    res{N}.zmin      = min(Gt.cells.z);
    res{N}.zmax      = max(Gt.cells.z);
    res{N}.volume    = sum(G.cells.volumes);
    res{N}.ctrapvols = volumesOfTraps(Gt,tac);
    res{N}.ccapacity = sum(res{N}.ctrapvols);
    res{N}.ntrapvols = volumesOfTraps(Gt,tan);
    res{N}.ncapacity = sum(res{N}.ntrapvols);
    
    res{N}.volumes = volumesOfTraps(Gt, tac); % using cell-based

    % visualization
    clf; set(gcf,'Position', [100 100 1000 1000]);
    
    %data = repmat(i, Gt.cells.num, 1); plotCellData(Gt, data, 'facea', .8, 'edgea', .01, 'edgec', 'k');
    plotCellData(Gt,Gt.cells.z,'FaceAlpha',.95,'EdgeColor','none');  
    plotGrid(Gt, tac.traps>0, 'FaceColor', 'red', 'EdgeColor', 'r');
    
    view(2);
    set(gca,'DataAspect',[1 1 rd.meta.cellsize/1000/10])
    title(['Grid of ' rd.name ', Coarsening level ' num2str(Res(N))]);
    xlabel('x, km'); ylabel('y, km'); zlabel('z, km');
    axis equal
    
    hcb = colorbar;
    hcb.TickLabels = hcb.Ticks;
    hcb.Label.String = 'Depth below sealevel, km';
    
    if saveFigures
        % save plot into folder
        folderName = 'CO2AtlasFormations_Plots';
        if ~exist(folderName, 'dir')
            mkdir(folderName)
        end
        saveas(gcf, [folderName '/' rd.name '_traps_cl' num2str(Res(N))], 'png')
    end
    
    
    Grids{N} = Gt;
    clear G Gt tan tac
end



end

%% (from resJohansen.m script in CAGEO-75): Number and size distribution of global traps
% Output number of global traps, as well as their average volumes, for each
% degree of coarsening.  Also plot their size distribution.
for i = 1:numRes
    fprintf('Coarsening level %d:\n', i);
    fprintf('  Num. global traps: %d\n', numel(res{i}.volumes));
    fprintf('  Total trap volume: %e m3\n', sum(res{i}.volumes));
    fprintf('  Avg. global trap size: %e m3\n', mean(res{i}.volumes));
end

figure;
defaultpos = get(0, 'DefaultFigurePosition');
set(gcf, 'Position', [defaultpos(1:2) - [300 100], 1200 800],...
   'PaperPositionMode','auto');
colorize = 'rgcmyb';

% Plot the total volume as subplot
subplot(2,2,1); cla
hold on
vol = cellfun(@(x) sum(x.volumes), res);
for i = 1:numRes
    bar(i, vol(i), colorize(i))
end
title('Total trap volume [m^3]', 'FontSize', 14)
set(gca, 'Color', get(gcf, 'Color'), 'FontSize', 14, ...
   'XTickLabel', regexp(num2str((1:numRes)*500), '  ', 'split'))
xlabel('Lateral resolution [m]', 'FontSize', 12);
axis tight

% Plot cumulative colume covered
subplot(2, 2, 2); cla
hold on;
cum_max = 0;
for i = 1:numRes
    cumul_vol = cumsum(sort(res{i}.volumes, 'descend'));
    plot([0, cumul_vol], colorize(i),'LineWidth', 1, 'Marker', '.', 'MarkerSize', 16);
    cum_max = max(cum_max, cumul_vol(end));
end
axis([0 60 0 cum_max*1.05]);
set(gca, 'Color', get(gcf,'Color'), 'FontSize', 14);
h=legend(regexp(num2str((1:numRes)*500), '  ', 'split'), 4);
title('Cumulative trap volume [m^3]', 'FontSize', 14);
xlabel('Number of traps counted', 'FontSize', 12);

% Bar plot: number of traps at each coarsening level
subplot(2, 2, 3);
hold on
for i = 1:numRes
    bar(i, numel(res{i}.volumes), colorize(i));
end
title('Total number of traps', 'FontSize', 14)
axis tight;
set(gca, 'Color', get(gcf, 'Color'), 'FontSize', 14, ...
   'XTick', [1 2 3 4 5 6], ...
   'XTickLabel', regexp(num2str((1:numRes)*500), '  ', 'split'))
xlabel('Lateral resolution [m]', 'FontSize', 12);

% Plot the average volume as subplot
subplot(2,2,4); cla
hold on
mvol = cellfun(@(x) mean(x.volumes), res);
for i = 1:numRes
    bar(i, mvol(i), colorize(i))
end
title('Average trap volume [m^3]', 'FontSize', 14)
set(gca, 'Color', get(gcf, 'Color'), 'FontSize', 14, ...
   'XTickLabel', regexp(num2str((1:numRes)*500), '  ', 'split'))
axis tight
xlabel('Lateral resolution [m]', 'FontSize', 12);

if saveFigures
    % save plot into folder
    folderName = 'CO2AtlasFormations_Plots';
    if ~exist(folderName, 'dir')
        mkdir(folderName)
    end
    saveas(gcf, [folderName '/' rd.name '_traps_compare'], 'png')
end






%% 4. Use method to place wells:
% a) manual
% b) use optimization (objective function) to determine well placement of N
% wells.
% c) use array of wells




%% 5. Given placed wells, compute their injection rates:
% a) compute rate such that an injection period of X years will fill the
% volume of traps along the spill-paths over a period of Y years.
% (user-specified X and Y)




%% *************************************************************************
% 6. Run a VE simulation, and obtain CO2 inventory breakdown:
% a) compare with and without dissolution
% b) compare bdry types (open, closed, semi)




%% 7. Use objective function to get optimized injection rates at wells.
