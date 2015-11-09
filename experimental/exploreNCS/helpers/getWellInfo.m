function [ wellinfo ] = getWellInfo( Gt, capOutput, varargin )
% Set-up wells for a generic injection scenario.

% Wells: injectors and/or producers
%
% Layout of wells could be:
%   - array of wells
%   - a well based on spill-paths, structural traps, reachable capacity,
%   - ...

% DESCRIPTION:
%   Several options are available for injector and producer well placement.
%   
%   Firstly, an array of injectors can be placed such that they cover the
%   entire top-surface, or only in a limited area. Or, the array of
%   injectors can be placed in a trap region, ideally the trap region with
%   the highest structural trap capacity. Info on structural trap capacity
%   is passed in by the capOutput argument.
%
%   Producers are placed by specifying a percentage of injector wells to
%   steal (i.e., an injector becomes a producer). The deepest elevation of
%   the top-surface will contain the producers.
%
%   A buffer along the formation boundary can be set such that no wells
%   exist within this area.

    opt.inj             = true;
    opt.inj_layout      = 'array_of_wells';
    
    opt.prod            = true;
    opt.prod_layout     = 'steal_some_inj_wells';
    
    
    % spacing for array of inj wells
    opt.DX      = 3*5000; % meters
    opt.DY      = 3*5000; % meters
    
    % spacing for array of prod wells
    opt.DXp     = 17500; % meters
    opt.DYp     = 17500; % meters
    
    % coverage of wells
    %opt.limits = 'portion';
    opt.minX = 4.6e5;
    opt.minY = 7.25e6;
    opt.maxX = 5.4e5;
    opt.maxY = 7.31e6;
    
    opt.limits          = 'bestReachCap'; % wells are placed in first N reachable traps of highest capacity.
    opt.numTopReachTrap = 3;
    
    
    
    % amount of inj wells to steal
    opt.steal        = 30; % percent
    opt.steal_layout = 'deepest'; % 'centralgroup', 'northgroup', 'nearbdry', 'deepest', 'shallowest', 'line', 'random', ...
    
    % distance from bdry where no well will be placed
    opt.buffer  = 1500; % meters
    
    opt = merge_options(opt, varargin{:});
    
    
    %% Initialize all output as empty;
    wellCoords_inj  = []; cinx_inj  = [];
    wellCoords_prod = []; cinx_prod = [];
    ta = trapAnalysis(Gt,'false'); % for plotting
    
    
    %% Get injection well coords according to specified layout:
    if opt.inj
        fprintf('Setting up injector well(s).\n')
        
        if strcmpi(opt.inj_layout,'array_of_wells')
            
            % default: place array of wells covering entire top-surface
            [l.minX, l.minY, l.maxX, l.maxY] = getLimitsOfGrid(Gt);
            
            if strcmpi(opt.limits,'portion')
                [l.minX, l.minY, l.maxX, l.maxY] = deal(opt.minX, opt.minY, opt.maxX, opt.maxY);
            end
            [ Xcoord, Ycoord, cinx ] = getArrayCoords(Gt, opt.DX, opt.DY, l);
            
            % keep wells inside specified area
            if strcmpi(opt.limits,'bestReachCap')
                
%                 % trap_region with highest structural reachable capacity
%                 [val,ind] = max(capOutput.structural_mass_reached);
%                 bestTrapRegion = ta.trap_regions(ind);
%                 
%                 cells = 1:Gt.cells.num;
%                 cinx2keep = cells( logical(ta.trap_regions==bestTrapRegion) ); 
                
                % OR:
                %cellsOfMaxReachCap = find(capOutput.structural_mass_reached == max(capOutput.structural_mass_reached));
                %cinx2keep = cellsOfMaxReachCap;
                
                
                % cell index pertaining to top N trap_regions of highest
                % structural reachable capacity.
                topReachCap  = unique(capOutput.structural_mass_reached);
                topReachCap  = sort(topReachCap,'descend');
                topReachCap  = topReachCap(1:opt.numTopReachTrap);
                cinx2keep = [];
                for i = 1:numel(topReachCap)
                    [ind]     = find(capOutput.structural_mass_reached == topReachCap(i));
                    cinx2keep = [cinx2keep; ind];
                end
                

                % any wells located in a cell inx that does not match
                % inx2keep is discarded, as it's outside the specified area
                [~,inx2keep,~] = intersect(cinx,cinx2keep);
                if isempty(inx2keep)
                   error(['No injectors have been placed inside '...
                       'specified region. Try reducing array spacing.']) 
                end
                cinx = cinx(inx2keep);
                
            end

        else
            error('Unknown type of injector well set-up.')

        end
        
        cinx_inj            = cinx;
        [ Xcoord, Ycoord ]  = getXYcoords(Gt, cinx_inj);
        wellCoords_inj      = [ Xcoord, Ycoord ];
        
    end
    
    
    %% Get producer well coords according to specified layout:
    if opt.prod
        fprintf('Setting up producer well(s).\n')
        
        if strcmpi(opt.prod_layout,'array_of_wells')
            [ Xcoord, Ycoord, cinx ] = getArrayCoords(Gt, opt.DXp, opt.DYp);

        elseif strcmpi(opt.prod_layout,'single')
            % TODO
            
        elseif strcmpi(opt.prod_layout,'steal_some_inj_wells')
            assert(~isempty(wellCoords_inj),'No injector wells have been set-up.')
            
            % number of inj wells to steal
            N = round(size(wellCoords_inj,1) * opt.steal/100);
            assert(N~=0, ['%d percent of injectors equates to ' ...
                '< 1 well, thus no producers will be placed. Try '...
                'increasing percentage to steal.'], opt.steal)
            
            % location of inj wells to steal
            if strcmpi(opt.steal_layout,'deepest')
                [ cinx ]     = getDeepestCellIndex(Gt, cinx_inj, N);
                [ Xcoord, Ycoord ]  = getXYcoords(Gt, cinx);
            else
                error('Steal layout unknown.')
            end

        else
            error('Unknown type of producer well set-up.')

        end
        
        wellCoords_prod = [ Xcoord, Ycoord ];
        cinx_prod       = cinx;
        
    end
    
    
    %% Checks:
    % check if any injectors and producers are in same cell, too close
    % together, or have been repeated
    
    % remove any well cell indexes on the bdry cells
    cinx_inj  = removeAnyBdryCellIndex(Gt, cinx_inj, opt.buffer);
    cinx_prod = removeAnyBdryCellIndex(Gt, cinx_prod, opt.buffer);
    
    
    % remove any zeros from well cell indexes
    cinx_inj  = cinx_inj( logical(cinx_inj) );
    cinx_prod = cinx_prod( logical(cinx_prod) );
    
    % remove any repeated well cell indexes, and then get updated
    % coords
    cinx_prod = unique(cinx_prod);
    cinx_inj  = unique(cinx_inj);

    if ~isempty(cinx_inj) && ~isempty(cinx_prod)
        
        tmp_inj = zeros(Gt.cells.num,1);
        tmp_inj(cinx_inj) = 1;
        tmp_prod = zeros(Gt.cells.num,1);
        tmp_prod(cinx_prod) = 1;
        conflicts = tmp_inj.*tmp_prod; % 1 means conflict at inx
        fprintf('%d injector(s) is/are being assigned as producer(s).\n', sum(conflicts)')
        
        % squeeze out conflicting inx from cinx_inj
        tmp_inj( logical(conflicts) ) = 0;
        cinx_inj = find(tmp_inj);
        
    end
    
    
    [ Xcoord, Ycoord ] = getXYcoords(Gt, cinx_prod);
    wellCoords_prod = [ Xcoord, Ycoord ];

    [ Xcoord, Ycoord ] = getXYcoords(Gt, cinx_inj);
    wellCoords_inj = [ Xcoord, Ycoord ];
    
    % final check is to visualize injector/producer well locations
    plotWells(Gt, ta, cinx_inj, cinx_prod);

    
    
    %% Return well info:
    % re-structure so it's wellinfo.prod... and wellinfo.inj... fields

    wellinfo.wellCoords_inj     = wellCoords_inj;   % physical coords (X,Y)
    wellinfo.cinx_inj           = cinx_inj;         % cell index of inject.
    wellinfo.wellCoords_prod    = wellCoords_prod;
    wellinfo.cinx_prod          = cinx_prod;
% 
%     wellinfo.inj_rate_MtperYr;      % Mt/yr
%     wellinfo.inj_time;
%     wellinfo.inj_steps;
%     
%     wellinfo.prod_rate_MtperYr;
%     wellinfo.prod_time;
%     wellinfo.prod_steps;
%     
%     wellinfo.mig_time;
%     wellinfo.mig_steps;


end

% -------------------------------------------------------------------------

function cinx = removeAnyBdryCellIndex(Gt, cinx, buffer)
% NB: remove on-bdry wells and in-buffer-zone wells separately
    
    % remove cell inx that is a bdry cell inx
    bdrycinx = sum(Gt.faces.neighbors(find( prod(Gt.faces.neighbors,2) == 0 ), :),2);
    [~,ia,~] = intersect(cinx, bdrycinx);
    fprintf('Found %d wells located on bdry cells.\n', numel(ia))
    cinx(ia) = 0;
    cinx = unique(cinx);          % to remove repeated cell indexes
    cinx = cinx( logical(cinx) ); % to remove any 0.
    
    
    % remove cell inx that is within 'buffer' meters from bdry
    cinx2remove = [];
    [ Xcoord, Ycoord ] = getXYcoords(Gt, cinx);
    [ Xcoord_bdry, Ycoord_bdry ] = getXYcoords(Gt, bdrycinx);
    for i = 1:numel(Xcoord)
        dv        = bsxfun(@minus, [Xcoord_bdry, Ycoord_bdry], [Xcoord(i), Ycoord(i)]);
        [v, ind]  = min(sum(dv.^2, 2));
        if v <= buffer
            cinx2remove = [cinx2remove; cinx(i)];
        end
    end
    [~,ia,~] = intersect(cinx, cinx2remove);
    fprintf('Found %d wells within %d meters of bdry.\n', numel(ia), buffer)
    cinx(ia) = 0;
    cinx = unique(cinx);          % to remove repeated cell indexes
    cinx = cinx( logical(cinx) ); % to remove any 0.

    
end

function plotWells(Gt, ta, cinx_inj, cinx_prod)

    figure; set(gcf,'Position',[1 1 1242 825]);

    subplot(1,2,1)
    title('Injectors')
    mapPlot(gcf, Gt, 'traps', ta.traps, 'trapalpha', 0.2, ...
        'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20, 'wellcells', cinx_inj, 'well_numbering', true);

    subplot(1,2,2)
    title('Producers')
    mapPlot(gcf, Gt, 'traps', ta.traps, 'trapalpha', 0.2, ...
        'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20, 'wellcells', cinx_prod, 'well_numbering', true);
    
    % adjust plots
    hfig = gcf;
    set(findobj(hfig.Children,'Type','axes'),'Fontsize',16,'box','on')
    axis(findobj(hfig.Children,'Type','axes'),'equal','tight')

end

% -------------------------------------------------------------------------

function [ Xcoord, Ycoord, cinx ] = getArrayCoords(Gt, DX, DY, limits)
% Array of wells covering top surface
% first well can be offset from cinx=1 by a specified distance @@

        % first get (x,y) coords of wells which span the grid limits
        [minX, minY, maxX, maxY] = deal(limits.minX, limits.minY, limits.maxX, limits.maxY);
        wX = [minX:DX:maxX]';
        wY = [minY:DY:maxY]';
        
        Xcoord = repmat(wX, numel(wY),1);
        Ycoord = [];
        for i = 1:numel(wY)
            Ycoord = [Ycoord; repmat(wY(i), numel(wX),1)];
        end
        
        cinx = getCellIndex(Gt, Xcoord, Ycoord);

%         % get corresponding cellIndex of physical location on top surface
%         for i = 1:numel(wX)
%             for j = 1:numel(wY)
%                
%                 cinx(i,j) = getCellIndex(Gt, wX(i), wY(j));
%                 
%                 % but then remove it if wX and wY were not on top surface
%                 if Gt.cells.centroids(i) ~= wX(i) && Gt.cells.centroids(j) ~= wY(j)
%                     cinx(i,j) = 0;
%                 end
%                 
%             end
%         end
%         cinx = reshape(cinx, numWells,1);
        
        % remove any repeated well cell indexes and zeros
        %cinx = unique(cinx(logical(cinx)));
        
        % update number of wells, and return coords
        [ Xcoord, Ycoord ]  = getXYcoords(Gt, cinx);

end

% -------------------------------------------------------------------------

function cellIndex = getCellIndex(Gt, Xcoord, Ycoord)
% Closest cell index of grid Gt corresponding to physical coordinate (X,Y)

    assert(numel(Xcoord) == numel(Ycoord));
    cellIndex = zeros(numel(Xcoord),1);
    
    for i = 1:numel(Xcoord)
        
        dv        = bsxfun(@minus, Gt.cells.centroids(:,1:2), [Xcoord(i), Ycoord(i)]);
        [v, ind]  = min(sum(dv.^2, 2));
        cellIndex(i) = ind; 

    end
    
end

% -------------------------------------------------------------------------

function [ Xcoord, Ycoord ] = getXYcoords(Gt, cellIndex)

    Xcoord = Gt.cells.centroids(cellIndex,1);
    Ycoord = Gt.cells.centroids(cellIndex,2);

end

% -------------------------------------------------------------------------

function [minX, minY, maxX, maxY] = getLimitsOfGrid(Gt)
% get max and min (X,Y) physical coordinates

    minX = min(Gt.cells.centroids(:,1));
    minY = min(Gt.cells.centroids(:,2));
    maxX = max(Gt.cells.centroids(:,1));
    maxY = max(Gt.cells.centroids(:,2));
    
end

% -------------------------------------------------------------------------

function [ deepCellInx ] = getDeepestCellIndex(Gt, cinx, N)
% z is array of all possible depths (i.e., top surface elevation where a
% well is located)
% N is the number of deepest cell indexes to find and return

    z       = zeros(Gt.cells.num,1);
    z(cinx) = Gt.cells.z(cinx);
    assert(numel(z) == Gt.cells.num);

    [sortedz, sortedzInx] = sort(z,'descend');
    %maxVals = sortedX(1:N);
    maxValsInx = sortedzInx(1:N);
    
    deepCellInx = maxValsInx;

end