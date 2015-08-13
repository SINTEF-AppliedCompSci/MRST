% for comparing all three Sleipner grid models

numRefCases = [-6 -4 -2 1 2 4];

for k = 1:numel(numRefCases)
    
    numRef = numRefCases(k);
    
    [ G_ieaghg, Gt_ieaghg, rock_ieaghg, rock2D_ieaghg ] = makeSleipnerModelGrid( 'modelName','IEAGHGmodel','refineLevel', numRef);

    [ G_original, Gt_original, rock_original, rock2D_original ] = makeSleipnerModelGrid( 'modelName','ORIGINALmodel','refineLevel', numRef );

    [ G_inhouse, Gt_inhouse, rock_inhouse, rock2D_inhouse ] = makeSleipnerModelGrid( 'modelName','INHOUSEmodel','refineLevel', numRef );

    % take the top and bottom of the SAND layer to be the same 2Dgrids from the
    % first two models. --> check top and bottom surfaces.


    gts = [Gt_ieaghg; Gt_original; Gt_inhouse];
    gs  = [G_ieaghg; G_original; G_inhouse];

    names = {'IEAGHG'; 'ORIGINAL'; 'INHOUSE'};
    mymap = [
             0         0    1.0000
             0    1.0000         0
        1.0000         0         0 ];

    % Inspect Top Surface Grid and H
    figure; set(gcf,'Position',[1 1 1000 1000])

    for i = 1:numel(gts)
        Gt = gts(i);
        G  = gs(i);
        name = names{1};
        avthk = mean(Gt.cells.H);
        numCells = Gt.cells.num;

        plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', mymap(i,:));
        legstr(i) = {[name avthk numCells]};
        %plotCellData(Gt, Gt.cells.z)

        % Compare three model grids
        res{i}.name      = names{i};
        res{i}.cells     = Gt.cells.num;
        res{i}.zmin      = min(Gt.cells.z);
        res{i}.zmax      = max(Gt.cells.z);
        res{i}.volume    = sum(G.cells.volumes);
        res{i}.surfarea  = sum(Gt.cells.volumes);
        res{i}.avgthk    = mean(Gt.cells.H);
    end

    legend(cellfun(@(x) x.name, res, 'UniformOutput', false), 'Location', 'EastOutside')

    view(3)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Top Surfaces elevation')



    % create table:
       fprintf('\n\n--------------Num Ref. =  %-2d-----------------------\n', numRef);
       fprintf('%-20s&   Cells  & Min Surf Ele. & Max Surf Ele. &  Volume  &  Surf. Area  & Avg. Thk. \\\\\n', 'Name');

       for i = 1:3
       fprintf('%-20s&  %6d  &     %4.0f      &     %4.0f      & %4.2e &   %4.2e   & %5.2f     \\\\\n',...
          res{i}.name, res{i}.cells, res{i}.zmin, res{i}.zmax, res{i}.volume, res{i}.surfarea, ...
          res{i}.avgthk );
       end
       fprintf('------------------------------------------------\n');

       clearvars -except numRefCases k

end


% Check top-surface elevation differences:
% relative to ORIGINAL:
%eleDiff_a = Gt_original - Gt_ieaghg.cells;
