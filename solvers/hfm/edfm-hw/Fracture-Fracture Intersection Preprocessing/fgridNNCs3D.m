function [G,fracplanes]=fgridNNCs3D(G,fracplanes,cloudlist,fullline,i,j,fracplanei,fracplanej,tol)
% This function takes in the global grid which contains a subfield 
% FracGrid (which itself contains subfields Frac1,2,3...N which follow the 
% standard grid_structure).
%
% The function also takes in indices i and j which signify that we wish to
% generate NNCs between Frac(i) and Frac(j).
%
% The NNCs will be generated and appended to G's subfield G.nnc which 
% itself contains the following subfields:
% -G.nnc.cells - each row has two global indices of connected cells.
% -G.nnc.T - each row has a conductivity index between corresponding
%             connected cells in G.nnc.cells. Refer to Moinfar et al
%             (2013).
%
% The intersection relationship between planes i and j will be saved in
% fracplanes(i) and frac fracplanes(j).
%
% Warning! Ensure that fracture grids are done with triangles only! This
% code only works with triangles.
%
% Future improvements: Include G.nnc.type.

% Extract Fracture Grid 'i' from G
Gfi=G.FracGrid.(['Frac',num2str(i)]);
Nfi=Gfi.cells.num;
cstart_i=Gfi.cells.start;
[cni,cposi] = gridCellNodes(Gfi,1:Gfi.cells.num);
aperture_i=fracplanes(fracplanei).aperture;
mcells_i=Gfi.matrix_connection.cells;

% Extract Fracture Grid 'i' from G
Gfj=G.FracGrid.(['Frac',num2str(j)]);
Nfj=Gfj.cells.num;
cstart_j=Gfj.cells.start;
[cnj,cposj] = gridCellNodes(Gfj,1:Gfj.cells.num);
aperture_j=fracplanes(fracplanej).aperture;
mcells_j=Gfj.matrix_connection.cells;

%% Filter using mcells
% For every mcells_i, we find a cloud of cells (26 max) around the cell and
% check if Frac_j intersects any of the cells (using mcells_j). If there
% are, save the pairs into a list. Repeat this for all mcells_i. The list
% will contain candidate pairs that might intersect.
pairlist=cell(size(mcells_i));
for i=1:length(mcells_i)
    cloud=cloudlist{i};
    jcells=find(ismember(mcells_j,cloud)); % Find any Frac_j cells that intersect this cloud 
    jcells=mod(jcells,Nfj); jcells(jcells==0)=Nfj;
    icell=mod(i,Nfi); icell(icell==0)=Nfi;
    pairlist{i}=[repmat(icell,length(jcells),1),jcells];
end
candidatelist=vertcat(pairlist{:});



%% Rigorous step
% We proceed from matrix cell to matrix cell, i.e. along mcells.

maxallocate=size(candidatelist,1);
nnc_cells=-1*ones(maxallocate,2); nnc_T=-1*ones(maxallocate,1); % preallocate, remove excess later

% i related data.
icellnodes_par=cell(maxallocate,1);
icellperm_par=zeros(maxallocate,1);
for i=1:maxallocate
    k=candidatelist(i,1);
    icellnodeind=cni(cposi(k):(cposi(k+1)-1));
    Nicellnodes=size(icellnodeind,1); % this is always even
    icellnodes=Gfi.nodes2D.coords(icellnodeind(1:Nicellnodes/2),:); % Nx3 matrix, each row containing xyz coords of nodes
%     icellnodes=0.5*(icellnodes(1:(Nicellnodes/2),:)+icellnodes(((Nicellnodes/2)+1):Nicellnodes,:));
    icellnodes_par(i)={icellnodes};
    icellperm_par(i)=Gfi.rock.perm(k);
end

% j related data. 
jcellnodes_par=cell(maxallocate,1); 
jcellperm_par=zeros(maxallocate,1);
for i=1:maxallocate
    l=candidatelist(i,2);
    jcellnodeind=cnj(cposj(l):(cposj(l+1)-1));
    Njcellnodes=size(jcellnodeind,1); % this is always even
    jcellnodes=Gfj.nodes2D.coords(jcellnodeind(1:Njcellnodes/2),:); % Nx3 matrix, each row containing xyz coords of nodes 
%     jcellnodes=0.5*(jcellnodes(1:(Njcellnodes/2),:)+jcellnodes(((Njcellnodes/2)+1):Njcellnodes,:));
    jcellnodes_par(i)={jcellnodes};
    jcellperm_par(i)=Gfj.rock.perm(l);
end

parfor m=1:maxallocate
    % i cell data
    icellnodes=icellnodes_par{m};
%     ixsect=linesegmentpebiintersect(fullline,icellnodes,tol);
    
    % j cell data
    jcellnodes=jcellnodes_par{m};
%     jxsect=linesegmentpebiintersect(fullline,jcellnodes,tol);
    
%     if ~(ixsect && jxsect)
%         continue;
%     end

    % Intersection calculation
    [xornot,~,xlength,df]=PEBIPEBIintersection3D(icellnodes,jcellnodes,tol);
    
    % if the cells Tik and Tjl intersect (xornot==true). This is redundant
    % and for double checking.
    if xornot
        % save global indices into G.nnc.cells
        indices=candidatelist(m,:);
        Findex_i=indices(1);
        Findex_j=indices(2);
        globindex_i=Findex_i+cstart_i-1;
        globindex_j=Findex_j+cstart_j-1;        
        nnc_cells(m,:)=sort([globindex_i,globindex_j]);
        
        % calculate transmissibility (Moinfar et al, 2013)
        Trans_i=icellperm_par(m)*aperture_i*xlength/df(1);
        Trans_j=jcellperm_par(m)*aperture_j*xlength/df(2);
        Trans_ij=(Trans_i*Trans_j)/(Trans_i+Trans_j);
        nnc_T(m)=Trans_ij;
    end
end

nnc_cells=removeexcess(nnc_cells,-1);
nnc_T=removeexcess(nnc_T,-1);
[~,ind]=unique(nnc_cells,'rows');
nnc_cells = nnc_cells(ind,:);
nnc_T = nnc_T(ind,:);
G.nnc.cells=[G.nnc.cells;nnc_cells];
G.nnc.T=[G.nnc.T;nnc_T];

% Save intersection data
if ~isempty(nnc_T)
    if isfield(fracplanes(fracplanei),'intersects')
        fracplanes(fracplanei).intersects=[fracplanes(fracplanei).intersects;fracplanej];
    else
        fracplanes(fracplanei).intersects=fracplanej;
    end
    
    if isfield(fracplanes(fracplanej),'intersects')
        fracplanes(fracplanej).intersects=[fracplanes(fracplanej).intersects;fracplanei];
    else
        fracplanes(fracplanej).intersects=fracplanei;
    end
end

end
