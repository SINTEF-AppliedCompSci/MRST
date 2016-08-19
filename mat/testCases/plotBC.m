function plotBC(bc, G)

isDir = strcmp(bc.type, 'pressure');

plotGrid(G);
hold on

n = G.faces.nodes(mcolon(G.faces.nodePos(bc.face),G.faces.nodePos(bc.face+1)-1));
X = G.nodes.coords(n,:);

for i = 1:numel(bc.face)
    if isDir(i)
        clr = 'b';
    else
        clr = 'r';
    end
    ii = 2*i-1:2*i;
    plot(X(ii,1), X(ii,2) , clr, 'linewidth', 5)
    hold on
end