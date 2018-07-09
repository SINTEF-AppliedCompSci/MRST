function [G,fracplanes]=fgridNNCs3D_serial_obs(G,fracplanes,i,j,tol)
% OBSOLETE AS OF 20161118: Upgraded to include parallel processing.

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
[cni,cposi] = gridCellNodes(Gfi,1:Gfi.cells.num);

% Extract Fracture Grid 'i' from G
Gfj=G.FracGrid.(['Frac',num2str(j)]);
typej=Gfj.gridtype;
cstart_j=Gfj.cells.start;
[cnj,cposj] = gridCellNodes(Gfj,1:Gfj.cells.num);

% check that both grid types are the same
if ~strcmp(typei,typej)
    error('Fracture grid type need to be the same. No capability to compare PEBI/Triangle yet!');
else
    type=typei;
end

count=0; % number of intersecting fracture gridcells per fracture
% outer loop: for each cell in Frac(i)
for k=1:Gfi.cells.num
    % arrange nodes and put in Tik
    nodes_ik=cni(cposi(k):(cposi(k+1)-1)); 
    Gfinodes=Gfi.nodes.coords(nodes_ik,:);
    numnodes_ik=size(Gfinodes,1);
    Tik=0.5*(Gfinodes(1:(numnodes_ik/2),:)+Gfinodes(((numnodes_ik/2)+1):numnodes_ik,:));
    
    % inner loop: for each cell in Frac (j)
    for l=1:Gfj.cells.num
        % arrange nodes and put in Tjl
        nodes_jl=cnj(cposj(l):(cposj(l+1)-1));
        Gfjnodes=Gfj.nodes.coords(nodes_jl,:);
        numnodes_jl=size(Gfjnodes,1);
        Tjl=0.5*(Gfjnodes(1:(numnodes_jl/2),:)+Gfjnodes(((numnodes_jl/2)+1):numnodes_jl,:));
        
        if strcmp(type,'Triangle')
            [xornot,~,xlength,df]=triangletriangleintersection3D(Tik,Tjl,tol);
        elseif strcmp(type,'Pebi')
            [xornot,~,xlength,df]=PEBIPEBIintersection3D(Tik,Tjl,tol);
        else
            error('Only Triangle or Pebi grids can be handled!');
        end
       
        % if the cells Tik and Tjl intersect (xornot==true)
        if xornot
            % save global indices into G.nnc.cells
            globindex_ik=cstart_i-1+k;
            globindex_jl=cstart_j-1+l;
            G.nnc.cells(end+1,:)=sort([globindex_ik,globindex_jl]);
            
            % calculate transmissibility (Moinfar et al, 2013)
            Trans_ik=Gfi.rock.perm(k)*fracplanes(i).aperture*xlength/df(1);
            Trans_jl=Gfj.rock.perm(l)*fracplanes(j).aperture*xlength/df(2);
            Trans_ikjl=(Trans_ik*Trans_jl)/(Trans_ik+Trans_jl);
            G.nnc.T(end+1)=Trans_ikjl;
            
            % increase number of intersecting grids per fracture by one
            count=count+1;
        end
    end
end

% if count is non-zero, the two fractures intersect. Save this data.
if count~=0
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
