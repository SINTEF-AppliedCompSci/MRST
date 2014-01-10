function altZoom(h, G)
    set(h, 'ButtonDownFcn', @onClick)
    bf = boundaryFaces(G);

function onClick(src, event)
    pt = get(gcf, 'CurrentPoint');
    pt = pt(1,:);
    pts = get(gca, 'CurrentPoint');
    [c f] = nearestCellLine(G, bf, pts);
    c = c(1);
    cpos = get(gca, 'CameraPosition');

    set(gca, 'CameraTarget', G.cells.centroids(c, :))
    set(gcf, 'WindowButtonMotionFcn', @(src, event) zoomfn(src, event, pt, gca, gcf, cpos, G.cells.centroids(c, :)))
    set(src, 'ButtonDownFcn', @(src, event) set(gcf, 'WindowButtonMotionFcn', []));
end

function zoomfn(src, event, pt, ax, fig, campos_old, camtarget)

    pt2 = get(fig, 'CurrentPoint');
    pt2 = pt2(1,:);

    x0 = campos_old;
    x1 = camtarget;
    k = 3*(pt2(2)-pt(2))./pt(2);
    pos = x0 - (x0-x1)*k;
%     pos = campos_old + 0.1*campos_old*(pt2(2)-pt(2))./pt(2)
%     campos_old
    set(ax, 'CameraPosition', pos);
end

end
