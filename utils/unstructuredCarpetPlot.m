function unstructuredCarpetPlot(G,val,varargin)
    
    opt = struct('nodeval', false, 'extrapolation', 'linear');
    [opt, extra] = merge_options(opt, varargin{:});

    f = G.cells.faces(:,1);
    n = sortNodes(G,f);
    
    if opt.nodeval
        z = val;
    else
        x = G.cells.centroids;
        z = scatteredInterpolant(x(:,1), x(:,2), val, 'linear', opt.extrapolation);
        x = G.nodes.coords;
        z = z(x(:,1), x(:,2));
    end
    vertices = [G.nodes.coords, z];
    
    numn = diff(2*G.cells.facePos);
    mnf = max(numn);
    faces = zeros(G.cells.num, mnf);
    for i = 1:G.cells.num
        nn = n(G.cells.facePos(i):G.cells.facePos(i+1)-1)';
        faces(i,:) = [nn, nan(1,mnf-numel(nn))];
    end
    
    gcf;
    patch('faces', faces, 'vertices', vertices, 'facevertexCdata', val, 'facecolor', 'interp', extra{:});
%     patch('faces', faces, 'vertices', vertices, 'facecolor', 'flat', extra{:});
    colormap(parula);
    
end

function n = sortNodes(G, f)

n = G.faces.nodes(mcolon(G.faces.nodePos(f),G.faces.nodePos(f+1)-1));
s = G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', diff(G.cells.facePos),1);

n = reshape(n, 2, []);
n(:,s) = n([2,1], s);
n = n(:);
n = n(1:2:end);

end