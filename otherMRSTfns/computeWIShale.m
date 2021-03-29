function WI = computeEDFMWIshale(G,radius, cells, varargin)
    nc = size(cells,1);
    NumFracs=length(fieldnames(G.FracGrid));
    opt = struct('Subset', [], ...
                 'Dir', 'z', ...
                 're',[],...
                 'Skin', zeros(nc, 1));
    opt = merge_options(opt, varargin{:});

    %Find Fracture parameter from G.FracGrid
    WI=zeros(nc,1);
    for i=1:NumFracs
        mask=find(cells(:,2)==i);
        if isempty(mask), continue, end
        
        frac_data=G.FracGrid.(['Frac' num2str(i)]);
        startID=frac_data.cells.start;
        localFracCells=cells(mask)-startID+1;

        Kf=frac_data.rock.perm(localFracCells);
        wf=frac_data.aperture;
        
        re_estimate=getWellRe(G,cells(mask));
        if(re_estimate~=-1)%Check if fracture cell is Cartesian
            opt.re(mask)=re_estimate;
        end
        
        WI(mask)=2 * pi * Kf.*wf ./ (log(opt.re(mask) ./ radius) + opt.Skin(mask));
        
        %Ignore Well-NF NNC for now, BinWang@06/20/2020
        if(isfield(G,'numHFplanes'))%NF enabled
            if(i>G.numHFplanes)
                WI(mask)=0.0;
            end
        end
        
    end

    if any(WI < 0)
       if any(re < radius)
          error(id('WellRadius'), ...
               ['Equivalent radius in well model smaller than well ', ...
                'radius causing negative well index'].');
       else
          error(id('SkinFactor'), ...
                'Large negative skin factor causing negative well index.');
       end
    end

    if ~isempty(opt.Subset)
        % Only return calculated WI for requested cells
        WI = WI(opt.Subset);
    end
end

function re=getWellRe(G,cellIDs)
%Compute EDFM well representative radius

%Bin modification
%Code to extract the nodes of a list of cells (lines 68-74 in EDFMgrid.m)
%[cn,cpos]=gridCellNodes(G,cellIDs);
%numFcells = numel(cellIDs); % mcells {matrix_cell,frac_cell}
%mcellnodes_par=cell(numFcells,1);
%for i=1:numFcells
%    mcellnodeind=cn(cpos(i):(cpos(i+1)-1)); %idx of all nodes of each cell
%    mcellnodes=G.nodes.coords(mcellnodeind,:);
%    mcellnodes_par{i}=mcellnodes;
%end

[cn,cpos]=gridCellNodes(G,cellIDs);
mcellnodeind=cn(cpos(1):(cpos(1+1)-1)); %idx of all nodes of each cell
poly_box=G.nodes.coords(mcellnodeind,:);
numVert=size(poly_box,1);
poly_plane=(poly_box(1:numVert/2,:)+poly_box(numVert/2+1:end,:))/2.0;
poly_plane=poly3Dccw(poly_plane);

%Case 1: Fracture intersected by a well is a vertical rectangle fracture
if(size(poly_plane,1)==4)
    tol=1e-3;
    Length_1=norm(poly_plane(1,:)-poly_plane(2,:));
    Length_2=norm(poly_plane(2,:)-poly_plane(3,:));
    Length_3=norm(poly_plane(3,:)-poly_plane(4,:));
    Length_4=norm(poly_plane(4,:)-poly_plane(1,:));
    %Check if polygon is a rectangle
    if abs(Length_1 - Length_3)<=tol && abs(Length_2 - Length_4)<=tol
        
        %OMO modification
        %dy=30.4800; dz=9.1440;
        dy=Length_1; dz=Length_2;

        % Table look-up (interpolation) for mimetic or 0.14 for tpf
        wc  = 0.14; %wellConstant(d1, d2, opt.InnerProduct);
        re  = wc .* sqrt((dy.^2) + (dz.^2));
    else %Non rectangle case
        re=-1;
    end
else %Other case, abritary shape and oriented fracture 
    re=-1;
end


end

function [xyz_order] = poly3Dccw(xyz)
%https://www.mathworks.com/matlabcentral/answers/429265-how-to-
%order-vertices-of-a-flat-convex-polygon-in-3d-space-along-the-edge
    xyzc = mean(xyz,1);
    P = xyz - xyzc;
    [~,~,V] = svd(P,0);
    [~,is] = sort(atan2(P*V(:,1),P*V(:,2)));
    xyz_order = xyz(is([1:end]),:);
end
