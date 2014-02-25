function trap_struct=findTrappingStructure(Gt,varargin)
% Find the trapping structure
%
% SYNOPSIS:
%     trap_struct=findTrappingStructure(Gt)
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%              use_multipoint -- locical to extend neighbours node neigbours 
%
% PARAMETERS:
%   G     - Valid top surface grid
%
% RETURNS:
% trap_stuct  - structure describing the traps of the top surface by
%
%   trap_struct=struct( 'z_spill_loc',z_spill_loc,...
%                       'g_trap',g_trap,...
%                       'trap_level',trap_level,...
%                       'z_spill_level',z_spill_level,...
%                       'z_spill_loc_level',z_spill_loc_level,...
%                       'Gtop',Gtop);
%
%    z_spill_loc - a vector for each cell with value of the trap
%                  maximal trap surface or 0 if no trapping in the cell
%    g_trap      - a vector for each cell with the value of the trapping
%    trap_level  - struct of sparse matrix define traps assosiated with
%                  each top point for each level
%    z_spill_level - struct of vectors of size number of traps with z_values for
%                    each trap surface.
%    z_spill_loc_level - struct with vector as z_spill_loc but for each
%                      - level
%
%
%  NOTE:
%    All of the surface has to have z values > 0.
%
%  Example:
%     trap_struct=findTrappingStructure(Gt)
%     figure(),clf, 
%     plotGrid(trap_str.Gtop,'FaceColor','none','edgeAlpha',0.1)
%     plotGrid(trap_str.Gtop, trap_str.z_spill_loc>0),axis off;axis tight;view(3)
%
% SEE ALSO:
%    topSurfaceGrid, findCellLines
% 
opt = struct('use_multipoint',false,'include_bc',true);
opt = merge_options(opt, varargin{:});
Gtop=Gt;
grid_flat=false;

if(~opt.include_bc)
    C=maxTPFAGravityMatrix(Gtop,'use_multipoint',opt.use_multipoint);
    [ci,ss] = components(C+C');%#ok
else
    [C,N,CC]=maxTPFAGravityMatrix(Gtop,'use_multipoint',opt.use_multipoint);
    [ci,ss] = components(CC+CC');%#ok
    ci=ci(1:Gtop.cells.num);
end
kk=0;
while ~grid_flat;
    kk=kk+1;
    fprintf(1,'Trap level %d: %d regions identified\n', kk, max(ci));
    volume=0;
    Gtop_new=Gtop;
    for i=1:max(ci)
        m=find(ci==i);
        if(numel(m)>-1)
            if(opt.use_multipoint)
                bc = boundaryCellsSubGrid(Gtop, m,'neigbours',N);
                [z_level,j]=min(Gtop.cells.z(bc));%#ok
                 % include spill cell which may be in other domain into region
                %%{
                faces=any(N==bc(j),2);
                sp_cells=N(faces,:);
                sp_cells=reshape(sp_cells,[],1);
                sp_cells=sp_cells(sp_cells>0);
                sp_cells=unique(sp_cells(:));
                sp_cells=sp_cells(Gtop.cells.z(sp_cells)>z_level);
                [z_level_new,j]=min(Gtop.cells.z(sp_cells));
                if(~isempty(j))
                    z_level=z_level_new;
                end
                %}
                
            else
                bc = boundaryCellsSubGrid(Gtop, m);
                [z_level,j]=min(Gtop.cells.z(bc));%#ok
                % include spill cell which may be in other domain into region
                %%{
                faces=Gtop.cells.faces(Gtop.cells.facePos(bc(j)):Gtop.cells.facePos(bc(j)+1)-1,1);
                faces=unique(faces);
                sp_cells=Gtop.faces.neighbors(faces,:);
                sp_cells=reshape(sp_cells,[],1);
                sp_cells=sp_cells(sp_cells>0);
                sp_cells=unique(sp_cells(:));
                sp_cells=sp_cells(Gtop.cells.z(sp_cells)>z_level);
                [z_level_new,j]=min(Gtop.cells.z(sp_cells));
                if(~isempty(j))
                    z_level=z_level_new;
                end
                %}
            end
            
            % include spill cell which may be in other domain into region            
            %[z_level, sind] =min(Gtop.cells.z(ci==i));  %#ok
            cells=find(ci==i & Gtop.cells.z<=z_level);
            %if(numel(cells)>0)
            Gtop_new.cells.z(cells)=z_level;
            volume=volume+sum(z_level-Gtop.cells.z(cells));%#ok
        end
        %end
    end
    fprintf(1,'Total volume %d\n:',volume);
    %clf,plotCellData(Gtop_new,ci);
    %drawnow;
    %pause()
    if(volume==0)
        grid_flat=true;        
    else
        Gtop=Gtop_new;
        if(~opt.include_bc)
            C=maxTPFAGravityMatrix(Gtop,'use_multipoint',opt.use_multipoint);%ok
            [ci,ss] = components(C+C');%#ok            
        else
            [C,DD,CC]=maxTPFAGravityMatrix(Gtop,'use_multipoint',opt.use_multipoint);%#ok
            [ci,ss] = components(CC+CC');%#ok            
            ci=ci(1:Gtop.cells.num);
        end
    end
end
%z_spill_loc=Gtop.cells.z-Gt.cells.z;
z_spill_loc=Gtop.cells.z;
z_spill_loc(Gtop.cells.z-Gt.cells.z<=0)=0;

cells=find(z_spill_loc>0);
% find all nodes for this cells
eIX = Gtop.cells.facePos;
nn  = double(diff([Gtop.cells.facePos(cells), ...
   Gtop.cells.facePos(cells + 1)], [], 2));

if(~isfield(Gtop.cells,'cellNodes'))
 Gtop.cells.cellNodes = getCellNodes(Gtop);
end
cn  = double(Gtop.cells.cellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));
% set all the node values for this cells to the cell center value
zz=rldecode(z_spill_loc(z_spill_loc>0),nn);
Gtop.nodes.z(cn)=zz;



trap_struct=struct('z_spill_loc',z_spill_loc,...
                   'Gtop',Gtop);