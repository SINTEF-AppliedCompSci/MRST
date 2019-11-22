function demoPlotPoly(p, marker, mcol, msz)
    plot([p(:,1);p(1,1)], [p(:,2);p(1,2)], marker, 'markerFaceColor', ...
        mcol, 'markerSize', msz)
end