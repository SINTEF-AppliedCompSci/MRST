function plotAdjustiblePlane(G, data, varargin)
 %{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
    opt = struct('Callback', []);
    opt = merge_options(opt, varargin{:});

    doSurf = false;
    doPrune = false;

    if isempty(opt.Callback)
        fig0 = gcf;
        bf = boundaryFaces(G);
    end

    if doSurf
        fig1 = figure;
        surfh = NaN;
    end
    planeh = NaN;
    outlineh = NaN;
    ph = figure('Name','Slice plane adjustment','NumberTitle','off','ToolBar','none', 'MenuBar','none');
    grp = uibuttongroup('Parent', ph, 'Units', 'Normalized', 'Position', [0 0 1 1]);
    createSlider = @(data, loc, val) uicontrol('Parent', grp, 'Style', 'slider', 'min', min(data), 'max', max(data), 'value', val, 'Units', 'normalized', 'Position', loc, 'Callback', @plotPlane);
    createText = @(text, loc)    uicontrol(grp, 'Style', 'Text', 'Units', 'Normalized', 'Position', loc, 'String', text);

    Mgc = max(G.cells.centroids);
    mgc = min(G.cells.centroids);
    crange = Mgc - mgc;


    pts = bsxfun(@rdivide, bsxfun(@minus, G.cells.centroids, mgc), crange);


    x = createSlider([-.5 1.5], [.05  0 .90 .1], .4);
    createText('Translate', [.3  .1 .4 .05])

    phi = createSlider([0 180],  [.05  .15 .90 .1], 0);
    createText('Rotate (x-axis)', [.3  .25 .4 .05])

    theta = createSlider([0 180], [.05  .30 .90 .1], 0);
    createText('Rotate (y-axis)', [.3  .4 .4 .05])

    thickness = createSlider([0 1], [.05  .45 .90 .1], 0.05);
    createText('Selection thickness', [.3  .55 .4 .05])

    selecttype = uicontrol(grp, 'Style', 'Listbox', ...
                                'String', {'Near', 'Above', 'Below'},...
                                'Units', 'normalized', ...
                                'Position', [.01 .6, .2 .3]);


    addlistener(x, 'Value', 'PostSet', @plotPlane);
    addlistener(phi, 'Value', 'PostSet', @plotPlane);
    addlistener(theta, 'Value', 'PostSet', @plotPlane);
    addlistener(thickness, 'Value', 'PostSet', @plotPlane);
    addlistener(selecttype, 'Value', 'PostSet', @plotPlane);

    plotPlane(0, [])
    function plotPlane(src, event)
        if isempty(opt.Callback)
            figure(fig0)
        end
        deleteHandle([planeh, outlineh])
        xv = get(x, 'value');
        pt0 = repmat(xv, 1, 3);
        th = get(theta, 'value')*2*pi/360;
        phv = get(phi, 'value')*2*pi/360;
        pt1 = [pt0(1) , pt0(2:3) + [sin(th),  cos(th)]];
        pt2 = [pt0(1:2) + [sin(phv),  cos(phv)], pt0(3)];
        N = cross(pt1 - pt0, pt2 - pt0);
        N = N./norm(N);

        tmp = bsxfun(@minus, pts, pt0);
        switch get(selecttype, 'Value')
            case 1
                filterfun = @(x,y) abs(x) < y;
            case 2
                filterfun = @(x,y) x <= 0;
            case 3
                filterfun = @(x,y) x >= 0;
        end
        subset = filterfun(sum(bsxfun(@(x,y) (x.*y), tmp, N), 2), get(thickness, 'Value'));
        % Outline plane
        outlinepts = [0 0; 0 1; 1 1; 1 0];
        outlinez = (sum(repmat(N(1:2), 4, 1).*(repmat(pt0(1:2), 4, 1) - outlinepts(:, 1:2)), 2) + N(3)*pt0(3))/N(3);
        outlinepts = [outlinepts, outlinez];
        outlinepts = bsxfun(@plus, bsxfun(@times, outlinepts, crange), mgc);

        if isempty(opt.Callback)
            planeh = plotCellData(G, data, subset);
            outlineh = patch(outlinepts(:, 1), outlinepts(:, 2), outlinepts(:, 3), 'k', 'FaceAlpha', 0);
        else
            opt.Callback(subset, outlinepts);
        end
        if doSurf
            plotSurf(src, event, subset, N, outlinepts, pt0.*crange + mgc, pt1.*crange + mgc, pt2.*crange + mgc);
        end
    end
    function plotSurf(src, event, subset, N, outlinepts, pt0, pt1, pt2)
        set(0, 'CurrentFigure', fig1);
        deleteHandle(surfh);

        ptsd = G.cells.centroids(subset, :);
        f     = intersect(boundaryFaces(G, find(subset)), bf);
        if doPrune
            ptsd = [ptsd; G.faces.centroids(f, :)];
        end
        proj = @(x, y) sum(bsxfun(@(x,y) (x.*y), x, y), 2)./norm(y);

        v1 = pt0-pt1;
        v2 = pt0-pt2;
        v1 = v1/norm(v1);
        v2 = v2/norm(v2);
        v2 = v2 - proj(v2, v1)*v1;

        ptsd1 = proj(ptsd, v1);
        ptsd2 = proj(ptsd, v2);

        ptsd = [ptsd1, ptsd2];
        d = data(subset);
        if doPrune
            d = [d; nan(numel(f), 1)];
        end
        T = TriScatteredInterp(ptsd(:,1), ptsd(:,2),  d);

        mptsd = min(ptsd); Mptsd = max(ptsd);
        getRange = @(i) mptsd(i):(Mptsd(i)-mptsd(i))/G.cartDims(i):Mptsd(i);
        [mx, my] = meshgrid(getRange(1), getRange(2));
        surfh = surf(mx, my, T(mx, my), 'edgec', 'none');
        set(gca, 'zdir', 'reverse')
        view(0,90)
        axis tight off

    end
end
