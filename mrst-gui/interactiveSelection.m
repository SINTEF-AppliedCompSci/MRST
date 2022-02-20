function interactiveSelection(seltype, callback, getPatch)
    ph = getPatch();

    ax = get(ph, 'Parent');
    fig = get(ax, 'Parent');

    start = get(ax, 'CurrentPoint');
    stop = NaN;

    set(fig, 'WindowButtonMotionFcn', @movefn)
    shape = nan;

    set(ph, 'ButtonDownFcn', @endClick)



function movefn(src, event)

    stop = get(gcf, 'CurrentPoint');

    figPos = get(gcf, 'Position');
    axPos = get(gca, 'Position');
    axCoords = axis();

    stop = stop./figPos(3:4);
    stop = stop - axPos(1:2);
    stop = min(stop, axPos(3:4));
    stop = stop./axPos(3:4);
    stop = (stop + axCoords([1,3])).*(axCoords([2 4]) - axCoords([1 3]));

    start = start(1,1:2);
%     hold on
%     plot(gca, stop(1), stop(2), 'r*')
    if ishandle(shape); delete(shape); end;
    switch lower(seltype)
        case 'square'
            shape = patch([start(1), stop(1), stop(1), start(1)],...
                          [start(2), start(2), stop(2), stop(2)],...
                          axCoords(3)*ones(4,1), 'FaceColor', 'None', 'LineWidth', 2, 'ButtonDownFcn', @endClick, 'EdgeColor', 'red');
        case 'circle'
            r = norm(stop - start)/2;
            [X,Y,Z] = cylinder(r);
%             h = axCoords([2 4]) - axCoords([1 3]);
            xyz = min(start, stop);
            shape = patch((X - xyz(1)), Y - xyz(2), Z.*abs(diff(zlim)), 'LineWidth', 2, 'FaceColor', 'blue', 'EdgeColor', 'red');
    end

    axis(axCoords);
    view(0, 90)
%     if get(ax, 'CurrentPoint') ~= start
%         set(fig, 'WindowButtonMotionFcn', []);
%         delete(shape)
%     end
end


function endClick(src, event)
    set(fig, 'WindowButtonMotionFcn', [])
    delete(shape)
    callback(start, stop)
    set(getPatch(), 'ButtonDownFcn', []);
end

end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}