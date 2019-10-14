pth     = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'spe10', 'fig');
savepng = @(name) print(fullfile(pth, name), '-dpng', '-r300');

%%

m = [500,  1, 11;
     11,  21,  1;
       5, 11,  1];

type = {'METIS', 'cake-cutter', 'cake-cutter', 'cake-cutter'};

partitions = makeNestedPartitions(modelFI, m, 'type', type);

close all
grids = cell(numel(m),1);
pos = [0, 0, 500, 800];
mdl = modelFI;
for pNo = 1:numel(partitions)
    figure('Position', pos)
    
    mdl = upscaleModelTPFA(modelFI, partitions{pNo});
    plotCellData(mdl.G, log10(mdl.rock.perm(:,1)), 'edgec', 'none');
%     grids{pNo} = generateCoarseGrid(modelFI.G, partitions{pNo});
    plotGrid(mdl.G, 'facec', 'none');
    colormap(pink)
    axis equal tight
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    savepng(['spe10-p', num2str(pNo)]);
end

%%

m = getSequentialModelFromFI(modelFI);
m.transportModel = AdaptiveTransportModel(m.transportModel, partitions);

level1 = 2*ones(modelFI.G.cells.num,1);

level2 = 2*ones(modelFI.G.cells.num,1);
cNo = 7601;
level2(cNo) = 1;

levels = [level1, level2];

close all

lperm = log10(modelFI.rock.perm(:,1));
for lNo = 1:2

    figure
    
    level = levels(:, lNo);
    partition = m.transportModel.makePartition(level);

    mdl = upscaleModelTPFA(modelFI, partition);
    cNo_c = mdl.G.partition(cNo);
    if lNo == 1
        x_c = mdl.G.cells.centroids(cNo_c,:);
    end
    
    ix = pdist2(mdl.G.cells.centroids, x_c) < 40;
    
    plotCellData(mdl.G, log10(mdl.rock.perm(ix,1)), ix, 'edgec', 'none');
    plotGrid(mdl.G, ix, 'facec', 'none');

    colormap(pink)

    caxis([min(lperm), max(lperm)]);
    axis equal off
    
    savepng(['spe10-zoom-', num2str(lNo)]);
    
end