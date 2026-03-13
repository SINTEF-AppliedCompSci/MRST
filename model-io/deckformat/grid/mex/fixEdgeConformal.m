function gm = fixEdgeConformal(gm)% on one cell face
%gmold=gm;
[cells,info]=checkConformalEdges(gm);
faces_nodes =gm.faces.nodes;
faces_nodePos=gm.faces.nodePos;
for iii=1:numel(cells)
    cell=cells(iii);
    hfaces = gm.cells.facePos(cell):gm.cells.facePos(cell+1)-1;
    %faces = gm.cells.faces(hfaces,1);
    for edir = 1:4
        %make orientated outer boundary of face this should include all
        %nodes at top and bottom
        myhfaces = find(gm.cells.faces(hfaces,2) == edir);
        if(not(numel(myhfaces)==0))
        myfaces = gm.cells.faces(hfaces(myhfaces));
        edges = [];
        %% find all edges 
        for i = 1:numel(myfaces)
            face = myfaces(i);
            %fnodes= gm.faces.nodes([gm.faces.nodePos(face):gm.faces.nodePos(face+1)-1]);
            %disp(['face ', num2str(face),' numnodes ', num2str(gm.faces.nodePos(face+1)-gm.faces.nodePos(face))])
            new_edges = [gm.faces.nodePos(face):gm.faces.nodePos(face+1)-2;gm.faces.nodePos(face)+1:gm.faces.nodePos(face+1)-1]';
            new_edges = [new_edges;gm.faces.nodePos(face+1)-1, gm.faces.nodePos(face)];
            if(not(cell == gm.faces.neighbors(face,2)))
                assert(gm.faces.neighbors(face,1)==cell)
                a=new_edges';a=reshape(a(end:-1:1),2,[])';
                new_edges = a;
            end
            %new_egdes
            %
            new_edges = gm.faces.nodes(new_edges);
            edges = [edges;new_edges];
        end
        %%
        % remove dublicated i.e. internal edges
        sedges = sort(edges,2);
        iedges=[sedges,[1:size(sedges,1)]'];
        a=sortrows(iedges);
        nind =find(all(diff(a(:,1:2),1) == 0,2));
        nind = [nind; nind+1];
        ind = true(size(a,1),1);
        ind(nind) = false;
        oedges = edges(a(ind,end),:);
        %%
        % make oriented list of edges
        tmpoedges = oedges;
        soedges = [tmpoedges(1,:)];
        while size(tmpoedges,1)>0
            ind = tmpoedges(:,1)== soedges(end,2);
            assert(sum(ind==1));
            soedges = [soedges;tmpoedges(ind,:)];
            tmpoedges=tmpoedges(not(ind),:);
        end
        soedges=soedges(1:end-1,:);
        %%
        % find if top bottom edges need extra nodes
        for tb = 5:6
            %
            bface = find(gm.cells.faces(hfaces,2) == tb);
            bface = gm.cells.faces(hfaces(bface),1);
            %
            bfnodes = gm.faces.nodes([gm.faces.nodePos(bface):gm.faces.nodePos(bface+1)-1]);
            assert(numel(bfnodes)==4)
                odir = [1,3,2,4];
                if(odir(edir) == 1)
                    bedge = [bfnodes(end),bfnodes(1)];
                else
                    bedge = [bfnodes(odir(edir)-1),bfnodes(odir(edir))];
                end
                fsign=1;
                fsigntb=1; 
                if(not(gm.faces.neighbors(bface,2) == cell) )
                    fsigntb=-1;                    
                end

                if(fsigntb == 1)
                    fsign=-1;
                    bedge=bedge(end:-1:1);
                end

                ind2 = find(soedges(:,1)==bedge(2));
                ind1 = find(soedges(:,1)==bedge(1));
                if(ind1==size(soedges,1) && ind2==1)
                    addnode = 0;
                elseif(ind2<ind1)
                    newedge = [soedges(ind1:end,1);soedges(1:ind2,1)];
                    addnode= numel(newedge)-2;
                else
                    newedge = soedges(ind1:ind2,1);
                    addnode= numel(newedge)-2;
                end

                if(addnode>0)
                    %%insert extra nodes
                    disp('Adding nodes')
                    disp(['face ', num2str(face),' numnodes ', num2str(gm.faces.nodePos(face+1)-gm.faces.nodePos(face))])
                    disp(['Face ', num2str(bface),' new nodes ', num2str(addnode)]);
                    if(fsign<0)
                        newedge=newedge(end:-1:1);
                    end
                    if(odir(edir) == 1)
                        %newfaceorg = [bfnodes(2:end-1);newedge];
                        newfaceorg = [bfnodes(1:end);newedge(2:end-1)];
                    else
                        newfaceorg = [bfnodes(1:odir(edir)-2);newedge;bfnodes(odir(edir)+1:end)] ;
                    end
                    %%
                    %% insert new nodes
                    %disp('Somthing new')
                    bfnodesnew = faces_nodes([faces_nodePos(bface):faces_nodePos(bface+1)-1]);
                    ind2 = find(bfnodesnew == newedge(end));
                    ind1 = find(bfnodesnew ==newedge(1));
                    if(ind1==numel(bfnodesnew) && ind2==1)
                       newface=[bfnodesnew;newedge(2:end-1)];
                    elseif((ind1+1)==ind2)
                        newface=[bfnodesnew(1:ind1);newedge(2:end-1);bfnodesnew(ind2:end)];
                    else
                        %assert(false)
                        disp('Nodes allready added')
                        newface = bfnodesnew;
                        addnode=0;
                    end
                    %newface = newfaceorg;
                    if(numel(bfnodesnew) == numel(bfnodes))
                        assert(all(bfnodesnew == bfnodes));
                        assert(all(newface == newfaceorg));
                    end
                    %% check if order of all old points is the same??
%                    if(numel(bfnodesnew) == numel(bfnodes))
                        newnodes = [faces_nodes(1:(faces_nodePos(bface)-1));...
                            newface;...
                            faces_nodes(faces_nodePos(bface+1):faces_nodePos(end)-1)];
                        %gm.faces.nodePos(bface+1:end)= gm.faces.nodePos(bface+1:end)+addnode;
                        faces_nodePos(bface+1:end) = faces_nodePos(bface+1:end)+addnode;
                        faces_nodes = newnodes;                        
%                     else
%                         disp('other')
%                         if(true)
%                         newnodes = [faces_nodes(1:(faces_nodePos(bface)-1));...
%                             newface;...
%                             faces_nodes(faces_nodePos(bface+1):faces_nodePos(end)-1)];
%                         %gm.faces.nodePos(bface+1:end)= gm.faces.nodePos(bface+1:end)+addnode;
%                         assert(numel(newnodes)==faces_nodePos(end)+addnode-1);
%                         faces_nodePos(bface+1:end) = faces_nodePos(bface+1:end)+addnode;                       
%                         faces_nodes = newnodes;
%                         end
%                     end
                    assert(numel(faces_nodes)==faces_nodePos(end)-1);
                    
                    %size(gm.faces.nodes);
                    %size(newnodes);
                    %gm.faces.nodes=newnodes;
                end
        end
        end
    end
end
gm.faces.nodePos = faces_nodePos;
gm.faces.nodes = faces_nodes;
