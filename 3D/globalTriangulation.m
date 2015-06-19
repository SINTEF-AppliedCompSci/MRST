function GlobTri = globalTriangulation(G, GlobTri)
%
[cn,cnmap] = gridCellNodes(G,1:G.cells.num); i = 1;

tx = delaunayTriangulation(G.nodes.coords(cn(cnmap(i):cnmap(i+1)-1),:));
P = tx.Points;
T = tx.ConnectivityList;
map = [];
map(transpose(1:size(T,1))) = 1;
lookup = zeros(G.nodes.num,1);
lookup(cn(cnmap(i):cnmap(i+1)-1)) = transpose(1:size(P,1));
for i = 2:G.cells.num
    nx = cn(cnmap(i):cnmap(i+1)-1);
    tx = delaunayTriangulation(G.nodes.coords(nx,:));
    map(end+1:end+size(tx,1)) = i;
    if any(lookup(nx))
        conn = tx.ConnectivityList;
        loc = lookup(nx);
        local_ind  = find(loc);
        conn2 = zeros(size(conn));
        for j = 1:numel(local_ind)
            conn2(conn == local_ind(j)) = loc(local_ind(j));
        end
        local_add = find(loc==0);
        global_add = size(P,1) + 1:size(P,1) + numel(local_add);
        for j = 1:numel(local_add)
            conn2(conn == local_add(j)) = global_add(j);
        end
        lookup(nx(local_add)) = global_add;
        P = [P;G.nodes.coords(nx(loc==0),:)]; %#ok
        T = [T;conn2]; %#ok
    else
        P = [P;G.nodes.coords(nx,:)]; %#ok
        conn = tx.ConnectivityList + size(P,1);
        T = [T;conn]; %#ok
    end
end
Gtri = triangulation(T,P);        
GlobTri.Tri = Gtri;
GlobTri.map = map;
return

% Gtri = delaunayTriangulation(G.nodes.coords);
% 
% t = struct;
% for i = 1:G.cells.num
%     t(i).triangulation = delaunayTriangulation(G.nodes.coords(gridCellNodes(G,i),:));
% end
% map = zeros(size(Gtri,1),1);
% Ic = incenter(Gtri);
% for i = 1:size(Gtri,1)
% %     qc = mean(Gtri.Points(Gtri.ConnectivityList(i,:),:));
%     for j = 1:G.cells.num
%         if ~isnan(pointLocation(t(j).triangulation,Ic(i,:)))
%            map(i) = j;
%            break;
%         end
%     end
% end