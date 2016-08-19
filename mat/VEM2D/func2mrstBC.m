function bc = func2mrstBC(bc,G)

v = zeros(numel(bc.face),1);
for i = 1:numel(bc.face)
    g = bc.func{i};
    if strcmp(bc.type{i}, 'flux')
        n = G.faces.nodes(G.faces.nodePos(bc.face(i)):G.faces.nodePos(bc.face(i)+1)-1);
        v(i) = G.faces.areas(bc.face(i))*[1/6, 1/6, 2/3]...
            *g([G.nodes.coords(n,:);G.faces.centroids(bc.face(i),:)]);
    else
        v(i) = g(G.faces.centroids(bc.face(i),:));
    end
end

isNeu = strcmp(bc.type, 'flux');
if abs(sum(v(isNeu))) > eps
    ii = find(isNeu);
    v(ii(1)) = v(ii(1)) - sum(v(isNeu));
end

bc.value = v;    