function [G,fracplanes]=fgridNNCs3D_serial(G,fracplanes,i,j,tol)
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
typei=Gfi.gridtype;
cstart_i=Gfi.cells.start;
Ni=Gfi.cells.num;
[cni,cposi] = gridCellNodes(Gfi,1:Ni);
aperture_i=fracplanes(i).aperture;

% Extract Fracture Grid 'j' from G
Gfj=G.FracGrid.(['Frac',num2str(j)]);
typej=Gfj.gridtype;
cstart_j=Gfj.cells.start;
Nj=Gfj.cells.num;
[cnj,cposj] = gridCellNodes(Gfj,1:Nj);
aperture_j=fracplanes(j).aperture;

% check that both grid types are the same, if not, default to pebi/pebi
% comparison
if ~strcmp(typei,typej)
    type='Pebi';
else
    type=typei;
end

%% Quick filters
% Since before calling this function, we have pre-determined that the
% fracture planes are not parallel and do intersect, we use the
% intersection line as a quick filter to skip cells that definitely do not
% intersect.

pointsi=fracplanes(i).points;
pointsj=fracplanes(j).points;
[~,fullline,~,~]=PEBIPEBIintersection3D(pointsi,pointsj,tol);

%% Rigorous step
% use combined indexing. For frac_i_cell ik and frac_j cell jl,
% the combined index is (ik-1)*Nj + jl (proceed along jl first, then ik)
% Given combined index C, the fraccell number ik = floor(C/Nj)+1
% Given combined index C, the possiblemcell number jl = C - (ik-1)*Nj

maxallocate=Ni*Nj;
nnc_cells=-1*ones(maxallocate,2); nnc_T=-1*ones(maxallocate,1); % preallocate, remove excess later

% i related data. These need to be repeated Nj times for each original
% element since we progress along jl first, then ik.
icellnodes_par=cell(maxallocate,1); % combined index
icellperm_par=zeros(maxallocate,1);
for k=1:Ni
    icellnodeind=cni(cposi(k):(cposi(k+1)-1));
    icellnodes=Gfi.nodes.coords(icellnodeind,:); % Nx3 matrix, each row containing xyz coords of nodes
    Nicellnodes=size(icellnodes,1); % this is always even
    icellnodes=0.5*(icellnodes(1:(Nicellnodes/2),:)+icellnodes(((Nicellnodes/2)+1):Nicellnodes,:));
    icellnodes_par((1:Nj)+(k-1)*Nj)={icellnodes};
    icellperm_par((1:Nj)+(k-1)*Nj)=Gfi.rock.perm(k);
end

% j related data. These need to be stacked Ni times on themselves since we
% progress along jl first, then ik.
jcellnodes_par=cell(maxallocate,1); % combined index
jcellperm_par=zeros(maxallocate,1);
for l=1:Nj
    jcellnodeind=cnj(cposj(l):(cposj(l+1)-1));
    jcellnodes=Gfj.nodes.coords(jcellnodeind,:); % Nx3 matrix, each row containing xyz coords of nodes
    Njcellnodes=size(jcellnodes,1); % this is always even
    jcellnodes=0.5*(jcellnodes(1:(Njcellnodes/2),:)+jcellnodes(((Njcellnodes/2)+1):Njcellnodes,:));
    jcellnodes_par((1:Nj:(1+(Ni-1)*Nj))+(l-1))={jcellnodes};
    jcellperm_par((1:Nj:(1+(Ni-1)*Nj))+(l-1))=Gfj.rock.perm(l);
end

for m=1:maxallocate
    % i cell data
    currenticell=floor(m/Nj)+1; % current fraccell number
    if mod(m,Nj)==0, currenticell=currenticell-1; end % correction in case i cell index is at Nj.
    % disp(['Processing pair ',num2str(currentfraccell),'/',num2str(Nf),'.']);
    Tik=icellnodes_par{m};
    if ~linesegmentpebiintersect(fullline,Tik,tol)
        continue; % continue to next iteration if cell does not intersect line
    end
    
    % j cell data
    currentjcell=m-(currenticell-1)*Nj; % current fraccell number
    Tjl=jcellnodes_par{m};
    if ~linesegmentpebiintersect(fullline,Tjl,tol)
        continue; % continue to next iteration if cell does not intersect line
    end
    
    xornot=false; xlength=0; df=[];
    if strcmp(type,'Triangle')
        [xornot,~,xlength,df]=triangletriangleintersection3D(Tik,Tjl,tol);
    elseif strcmp(type,'Pebi')
        [xornot,~,xlength,df]=PEBIPEBIintersection3D(Tik,Tjl,tol);
    end
    
    % if the cells Tik and Tjl intersect (xornot==true)
    if xornot
        % save global indices into G.nnc.cells
        globindex_ik=currenticell+cstart_i-1;
        globindex_jl=currentjcell+cstart_j-1;
        nnc_cells(m,:)=sort([globindex_ik,globindex_jl]);
        
        % calculate transmissibility (Moinfar et al, 2013)
        Trans_ik=icellperm_par(m)*aperture_i*xlength/df(1);
        Trans_jl=jcellperm_par(m)*aperture_j*xlength/df(2);
        Trans_ikjl=(Trans_ik*Trans_jl)/(Trans_ik+Trans_jl);
        nnc_T(m)=Trans_ikjl;
    end
end

nnc_cells=removeexcess(nnc_cells,-1);
nnc_T=removeexcess(nnc_T,-1);
G.nnc.cells=[G.nnc.cells;nnc_cells];
G.nnc.T=[G.nnc.T;nnc_T];

% if nnc list size is not empty, the two fractures intersect. Save this data.
if ~isempty(nnc_T)
    if isfield(fracplanes(i),'intersects')
        fracplanes(i).intersects=[fracplanes(i).intersects;j];
    else
        fracplanes(i).intersects=j;
    end
    
    if isfield(fracplanes(j),'intersects')
        fracplanes(j).intersects=[fracplanes(j).intersects;i];
    else
        fracplanes(j).intersects=i;
    end
end

end
