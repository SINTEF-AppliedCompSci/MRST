function [G,fracplanes]=EDFMgrid(G,fracplanes,varargin)

if nargout~=2
    error('Output strictly needs to be both G and fracplanes')
end

t1=clock;

opt = struct('Tolerance'   ,    1e-5, ...
             'plotgrid', false, ...
             'fracturelist',1:length(fracplanes), ...
             'Verbose', true,...
             'lowesttol',1e-15);
         
opt = merge_options(opt, varargin{:});

verbose = opt.Verbose;

fracturelist=opt.fracturelist;

%% Rearrange G if hasn't been done before
if ~isfield(G,{'Matrix','FracGrid'})
    % If G has never been preprocessed before, structure will be
    % adjusted.
    Gtemp=G;
    Gtemp.Matrix=G;
    Gtemp.FracGrid=struct;
    Gtemp.nnc=[];
    G=Gtemp;    
end

%% If fracplanes is being preprocessed for the first time, add new field
if ~isfield(fracplanes,'processed')
    for i=1:length(fracplanes)
        fracplanes(i).processed=false;
    end
end

%% Remove from list fractures that have been preprocessed
templist=[];
for i=fracturelist
    if fracplanes(i).processed % if preprocessed before, remove from list
        dispif(verbose,['Fracplane #',num2str(i),' has been preprocessed before.']);
        continue;
    end
    templist=[templist,i];
end
fracturelist=templist;

%% Calculate plane normals

% Plane normals
if ~isfield(fracplanes,'normal')
    for i = 1:numel(fracplanes)
        points = [fracplanes(i).points;fracplanes(i).points(1,:)];
        diffp = diff(points,1);
        normal = cross(diffp(1,:), diffp(2,:));
        fracplanes(i).normal = normal/norm(normal);
    end
end

%% Generate fracture grids
Gm=G.Matrix;
nummatcells=Gm.cells.num;

Fstorecell=cell(length(fracturelist),1);

[cn,cpos]=gridCellNodes(Gm,1:nummatcells);
mcellnodes_par=cell(nummatcells,1);
for i=1:nummatcells
    mcellnodeind=cn(cpos(i):(cpos(i+1)-1));
    mcellnodes=Gm.nodes.coords(mcellnodeind,:);
    mcellnodes_par{i}=mcellnodes;
end

% gaugetol calculation for circular filter
[~,mindex]=max(Gm.cells.volumes);
mcellnodes=mcellnodes_par{mindex};
gaugetol=norm(max(mcellnodes)-min(mcellnodes)); % for out of plane matrix cell filtering
% end of gaugetol calculation

lenfraclist=num2str(length(fracturelist)); % for verbose

for k=1:length(fracturelist)
    tol=opt.Tolerance; % reset to default tolerance for subsequent fractures
    
    i=fracturelist(k);
    points=fracplanes(i).points;
    planenormal=fracplanes(i).normal;
    planecenter=sum(points)/size(points,1);
    aperture=fracplanes(i).aperture;
    
    tightentolerance=true; % set to true to start while loop
    while tightentolerance
        % Ensure tolerance is reasonable
        assert(tol>opt.lowesttol,'Tolerance too tight.');
        
        % Narrow down matrices that intersect fracture
        planepoint=points(1,:);
        temp=bsxfun(@minus,points,planecenter);
        gaugeradius=realsqrt(max(sum(temp.^2,2))); % for in plane matrix cell filtering
        
        
        mcellsintersect=false(nummatcells,1); % preallocate max number of cells, remove excess later
        
        parfor j=1:nummatcells % check for plane intersection. PARFOR HERE.
            % find out if cell 'j' intersects infinite plane
            mcellnodes=mcellnodes_par{j};
            [xtruth,xtype]=checkplaneAABBintersect(mcellnodes,planepoint,planenormal,tol);
            
            if xtruth && (strcmp(xtype,'interior') || strcmp(xtype,'face'))
                circtruth=circularfilter(mcellnodes,planenormal,points,tol,'gaugetol',gaugetol,'gaugeradius',gaugeradius,'center',planecenter);
                if circtruth
                    mcellsintersect(j)=true; % if cell j possibly intersects fracture, assign logical true in row
                end
            end
        end
        
        possiblemcells=find(mcellsintersect);
        dispif(verbose,['Processing ',num2str(k),' out of ',lenfraclist,' fractures.\n']);
        %     disp([num2str(length(possiblemcells)),'/',num2str(nummatcells),' cells were selected by quick filter.']);
        
        % at this point, we have an array possiblecells which contains indices of
        % cells that could potentially intersect the fracture. We proceed to a more
        % rigorous intersection detection in the next step.
        
        % RIGOROUS STEP
        Nm=length(possiblemcells);
        
        mcells=-1*ones(Nm,1);area=-1*ones(Nm,1); % preallocate, remove excess later
        typenum=-1*ones(Nm,1); % preallocate, remove excess later
        fraccellpoints=cell(Nm,1);
        centroidside=false(Nm,1);
        
        % m related data.
        mcellnodes_par_i=cell(Nm,1);
        for j=1:Nm
            mcellnodeind=cn(cpos(possiblemcells(j)):(cpos(possiblemcells(j)+1)-1));
            mcellnodes=Gm.nodes.coords(mcellnodeind,:); % 8x3 matrix, each row containing xyz coords of nodes
            mcellnodes_par_i(j)={mcellnodes};
        end
        
        
        parfor j=1:Nm % PARFOR HERE.
            
            
            % matrix cell data
            mcellnodes=mcellnodes_par_i{j};
            
            % check intersection
            [xtruth,xarea,xtype,xlocation,fraccellpoints{j}]=pebiAABBintersect(points,mcellnodes,tol);
            % [~,~,~,~,fraccellpoints{j}]=pebiAABBintersect(points,mcellnodes,tol);
            
            if xtruth && strcmp(xtype,'polygon')
                % if fracture plane intersects a matrix cell, save details to
                % cells, area and CI.
                mcells(j)=possiblemcells(j);
                
                area(j)=xarea;
                
                switch xlocation
                    case 'none'
                        xlocation=0;
                    case 'boundary'
                        xlocation=1;
                        centroid=sum(mcellnodes)/size(mcellnodes,1);
                        dir=dot(centroid-points(1,:),planenormal);
                        dir=dir/abs(dir);
                        centroidside(j)=dir>0;
                    case 'interior'
                        xlocation=2;
                        centroidside(j)=true;
                end
                
                typenum(j)=xlocation;
            end
        end
        
        % remove excess data in mcells, area, typenum
        intersected=~cellfun('isempty',fraccellpoints);
        mcells=[mcells(intersected & centroidside);mcells(intersected & ~centroidside)]; %
        area=[area(intersected & centroidside);area(intersected & ~centroidside)];
        typenum=[typenum(intersected & centroidside);typenum(intersected & ~centroidside)];
        
        % Construct V
        fraccellpoints=fraccellpoints(intersected & centroidside);
        V=vertcat(fraccellpoints{:});
        
        % Construct C
        C=cellfun(@(c) 1:size(c,1),fraccellpoints,'UniformOutput',false);
        for j=2:size(C,1)
            addTo=C{j-1}(end);
            C{j}=C{j}+addTo;
        end
        
        % How to read the above data
        % Cell f has nodes C{f}, which have coordinates V(C{f},:)
        % Cell f intersects matrix cell mcells(f) and the intersection area is
        % area(f) and the intersection type is typenum(f)
        
        Fgrid=fractureplanegeneralgrid(V,C,points,planenormal,aperture,tol);
        
        % Save matrix-fracture connection information
        Fgrid.matrix_connection.cells=mcells;
        Fgrid.matrix_connection.area=area;
        Fgrid.matrix_connection.type=typenum;
        
        % Check volume. If not equal, tightentolerance
        targetvol = convexpolygonarea(points,opt.Tolerance)*aperture;
        totcellvol = sum(Fgrid.cells.volumes);
        tightentolerance=~(abs(targetvol-totcellvol)<opt.Tolerance);
        tol=tol/10;
        
    end
    Fstorecell{k}=Fgrid;
end


%% Append new fracture grids to G
% Create FracGrid which contains fracture grids for every fracture in
% fracplanes. FracGrid.Frac1,....,N will have a grid structure just like
% the matrix grid. Additionally, they will also contain global starting
% indices for the cells, nodes and faces, e.g. if starting global index for
% Frac1 is 1001, then the second cell in Frac1 will be the 1002-th cell in
% the global grid. The porosity and permeability data from fracplanes will
% be transferred into FracGrid. 
FracGrid_new=G.FracGrid;
cstart = G.cells.num + 1; % new fractures processed will be appended to previous global grid
fstart = G.faces.num + 1;
nstart = G.nodes.num + 1;
lastfrac=length(fieldnames(G.FracGrid));
for j=1:length(fracturelist)
    i=fracturelist(j);
    
    % fieldname 
    fieldname=['Frac',num2str(j+lastfrac)];
    
    % Create new field in FracGrid
    FracGrid_new.(fieldname)=Fstorecell{j};
    
    % Save global starting indices and compute next one
    FracGrid_new.(fieldname).cells.start = cstart;
    FracGrid_new.(fieldname).faces.start = fstart;
    FracGrid_new.(fieldname).nodes.start = nstart;
    cstart = cstart + FracGrid_new.(fieldname).cells.num;
    fstart = fstart + FracGrid_new.(fieldname).faces.num;
    nstart = nstart + FracGrid_new.(fieldname).nodes.num;
    
    % Poroperm
    FracGrid_new.(fieldname).rock.perm = ones(FracGrid_new.(fieldname).cells.num,1)*fracplanes(i).perm;
    FracGrid_new.(fieldname).rock.poro = ones(FracGrid_new.(fieldname).cells.num,1)*fracplanes(i).poro;
    
    % Vertices of fracture plane
    FracGrid_new.(fieldname).points=fracplanes(i).points;
    
    % Mapping back to fracplanes
    FracGrid_new.(fieldname).fracplanenumber=i; % mapping fracgrid to fracplanes
    
    % Set default NNC processing status
    FracGrid_new.(fieldname).matrixnnc=false; % newly processed fracture has not been checked for fracture-matrix nnc
    FracGrid_new.(fieldname).fracgridnnc=[]; % newly processed fracture has not been checked against other fracture for fracture-fracture nnc
    
    % Save information on grid type (PEBI vs Triangle)
    FracGrid_new.(fieldname).gridtype=FracGrid_new.(fieldname).type; 
    FracGrid_new.(fieldname)=rmfield(FracGrid_new.(fieldname),'type');
    
    % Save mapping and processing status in fracplanes
    fracplanes(i).fracgrid=j+lastfrac; % mapping fracplanes to fracgrid
    fracplanes(i).processed=true; % record fracplanes that have been processed
    

    
end
G.FracGrid=FracGrid_new;

%% Plot Grid
if opt.plotgrid
    figure; plotGrid(G.Matrix,'facealpha',0);
    for i = 1:numel(fieldnames(G.FracGrid))
        fieldname=['Frac',num2str(i)];
        if ismember(G.FracGrid.(fieldname).fracplanenumber,opt.fracturelist)
            facecolor='y';
        else
            facecolor='r';
        end
        plotGrid(G.FracGrid.(['Frac',num2str(i)]),'facealpha',0.75,'Facecolor',facecolor);
    end
    view(15,20);
    axis equal tight;
end

%% Regenerate global grid. Fracture and matrix grids are saved under
% G_global.FracGrid and G_global.Matrix (use assembleGlobalGrid)
% G_out=G_in; G_out.FracGrid=FracGrid_new; G_out.nnc=[]; % this line is a trick to fit the input for assembleGlobalGrid
Gtemp=G.Matrix; Gtemp.FracGrid=G.FracGrid; Gtemp.nnc=G.nnc;
G=assembleGlobalGrid(Gtemp);
% at this point, G_global will
% contain G_global.Matrix and G_global.FracGrid which will each contain the
% original matrix and fracture grids.

%% Calculate elapsed time
t2=clock;
e=etime(t2,t1);
dispif(verbose,['EDFM grid assembly took ',num2str(e),' seconds.']);

end