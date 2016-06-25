function h = colorizeCatchmentRegions( Gt, ta )
% Add catchment areas to current figure window

    % Plot accumulation regions
    if max(ta.traps)<128
       cmap = lines(max(ta.traps) + 1);
    else
       cmap = jet(max(ta.traps) + 1);
    end
    colormap(cmap);
    %map = greyMask(cmap);
    map = (cmap + repmat(get(gcf, 'Color'), size(cmap, 1), 1))./2;
    map(1,:) = get(gcf, 'Color');
    tmp = ta.trap_regions;
    tmp(tmp>max(ta.traps)) = max(ta.traps);

    h = plotCellData(Gt, ones(Gt.cells.num, 1), 'EdgeColor', 'none');
    set(h, 'FaceVertexCData', map(tmp + 1, :))

end

