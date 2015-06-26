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

studySea = 'Norwegian';

if strcmpi(studySea,'Barents')
    names = {'tubЖen', 'stЫ', 'nordmela', 'knurr', 'fruholmen', ...
    'bjarmeland', 'bcuBarentsSea'};

elseif strcmpi(studySea,'Norwegian')
    names = {'Пre', 'ror', 'not', 'ile', 'garn', 'bcuNorwegianSea', ...
    'baseПre', 'tilje'};

elseif strcmpi(studySea,'NorwegianNorth')
    %TODO: or simply run showCO2atlas.m
    disp('Use showCO2atlas.m to study North Sea formations.')
    return
    
else
    disp('Sea not specified.')
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
% Visualize formations and save into folder
folderName = 'CO2AtlasFormations_Plots';
if ~exist(folderName, 'dir')
    mkdir(folderName)
end

for i = 1:numel(names)
    
    frd = rawdata{1,i};
    clf; set(gcf,'Position', [100 100 2000 1000]);
    surf(frd.data); shading interp
    
    % make adjustments to figure, such as axes labels
    title([frd.name ' ' frd.variant])
    view(3); axis tight
    set(gca,'DataAspect',[1 1 frd.meta.cellsize/10])
    ax = gca;
    ax.XTickLabel = ax.XTick*frd.meta.cellsize/1000;
    ax.YTickLabel = ax.YTick*frd.meta.cellsize/1000;
    ax.ZTickLabel = ax.ZTick/1000;
    hold off
    xlabel('km'); ylabel('km'); zlabel('km');

    % save figure
    saveas(gcf, [folderName '/' frd.name], 'png')
    clear frd
    
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

mymap = [0 0 0; 0.94 0.94 0.94; 0.5 0.5 0.5];
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
    surf(frd.data, 'FaceColor', 'none', 'EdgeColor', mymap(i,:)); %shading interp

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
    ax.ZTickLabel = ax.ZTick/1000;
   
    hold off
    xlabel('km'); ylabel('km'); zlabel('km');
end

% save plot into folder
folderName = 'CO2AtlasFormations_Plots';
if ~exist(folderName, 'dir')
    mkdir(folderName)
end
saveas(gcf, [folderName '/' 'HammerfestAquifer'], 'png')


% Plot of x-section through aquifer to show three layers of formations:






%% 3. Analyze traps with spill-point analysis, and plot structural traps.
% Create trap analysis




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
