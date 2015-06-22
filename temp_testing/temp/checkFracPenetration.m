function isit = checkFracPenetration(G,fracplane,cell,sign,...
    foi,faceNodes,fnodes,fnmap)


% Check if the infinite plane of the fracture intersects any
% face of interest by checking the node signs of that face
doesInt = zeros(numel(foi),1);
for i = 1:numel(foi)
    noi = fnodes(fnmap(foi(i)):fnmap(foi(i)+1)-1);
    if all(sign.NodeSign(noi,1)==0)
        isit = true;
        return
    else
        checkif = abs(sum(sign.NodeSign(noi,1))) < ...
            numel(find(faceNodes(foi(i),:)));
        if checkif
            noi = [noi;noi(1)]; %#ok
            edges = [noi(1:end-1),noi(2:end)];
            doesInt(i) = checkFaceIntersection(G,fracplane,...
                sign.NodeSign, edges);
        end
    end
end
if any(doesInt)
    isit = true;
else
    n = gridCellNodes(G,cell);
    if sum(sign.NodeSign(n,1))==0
        isit = true;
    else
        isit = false;
    end
end
return

function isInt = checkFaceIntersection(G,fracplane,nodeSign,edges)
isInt = 0;
in = zeros(size(edges,1),1);
on = zeros(size(edges,1),1);
for i = 1:size(edges,1);
    sums = abs(sum(nodeSign(edges(i,:))));
    switch sums
        case {0,1}
            if any(nodeSign(edges(i,:)))
                [Int,type] = linePlaneIntersection(fracplane.normal,...
                    mean(fracplane.points,1),G.nodes.coords(edges(i,:),:));
            else
                Int = G.nodes.coords(edges(i,:),:);
                type = 2;
            end
            if type == 1 || type == 2
                [tim,tom] = inPolygon3D(Int,fracplane.points,'tolerance',eps*100);
                in(i) = sum(tim); on(i) = sum(tom);
            end
    end
end
if any(in), isInt = true; end
return