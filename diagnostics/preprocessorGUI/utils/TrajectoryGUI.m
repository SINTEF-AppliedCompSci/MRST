classdef TrajectoryGUI < handle
    properties
        Figure
        Axes1
        Axes2
        WellPlotOrig
        WellPlotNew = [];
        Data
        Menu
        model
        W
        ccol
        Patch3D
        Patch2D
        Line3D
        Line2D
        WNew
        layoutParams   = struct('menuWidth',        300, ...
                                'itemSpace',         40, ...
                                'includeToolTips', true);
    end
    properties (Dependent)
        wellNo
    end
    
    methods
        function d = TrajectoryGUI(model, W, varargin)
            require wellpaths
            opt = struct('style', 'default', 'wellNo', []);
            opt = merge_options(opt, varargin{:});
            
            fprintf('\n\n')
            fprintf('*  press new to create new well\n')
            fprintf('*  right-click inside axis for drawing options\n')
            fprintf('*  left-click edge of points to drag\n')
            fprintf('*  right-click points for options\n')
            fprintf('*  click save to save well\n')
            fprintf('*  click launch diagnsotics to view diagnostics for saved well\n')
            fprintf('\n\n')
            
            
            % ------ Create figure window ---------------------------------
            screensize = get( groot, 'Screensize' );
            wsize = .75*screensize(3:4);
            wsize = [screensize(3:4)-wsize*1.1, wsize];
            d.Figure = limitedToolbarFigure('Position', wsize);
            
            d.model = model;
            d.W = W;
            d.W = addTrajectories(d.W, d.model.G, 10);
            % add fields for faster comp
             d.model.G = addBoundingBoxFields(d.model.G);
            %--------------------------------------------------------------
            itemOpts = {'Parent', d.Figure, 'Visible','off', 'style', opt.style};
            wsel = WellEditSelector(W(end), 'wellNames', {W(end).name}, itemOpts{:});
            
            
            d.Menu = UIMenu('Title', 'Menu', 'Parent', d.Figure, ...
                             itemOpts{:}, 'items', {wsel});
            
            % setup axes
            d.Axes1 = axes('Parent', d.Figure, 'Units', 'pixels', 'ZDir', 'reverse');
            %d.Axes1 = addAxesContextMenu(d.Axes1);
            d.Axes2 = axes('Parent', d.Figure, 'Units', 'pixels', 'YDir', 'reverse');
            %d.Axes2 = addAxesContextMenu(d.Axes2);
            d.ccol = model.rock.perm(:,1);
            d.Patch3D = CellDataPatch(model.G, d.ccol, ...
                                     'Parent', d.Axes1, 'EdgeColor', [.3 .3 .3], ...
                                     'FaceAlpha', .2, 'EdgeAlpha', .3, ...
                                     'BackFaceLighting', 'lit', ...
                                     'Hittest', 'off');
            
            d.WellPlotOrig = WellPlotHandle(model.G, W, 'Parent', d.Axes1);                    
            d.WNew = [];
            
            wsel.newButton.Callback = @d.putXYLine;
            wsel.viewButton.Callback = @d.viewTraj;
            wsel.saveButton.Callback = @d.saveWell;
            wsel.launchButton.Callback = @d.launchDiagnostics;
            d.Figure.SizeChangedFcn        = @d.layout;
            d.layout();
            wsel.popupCallback([], []);
            wsel.typePopup.Enable = 'off';
        end
        
        function val = get.wellNo(d)
            val = 1;
        end
    
        %%-----------------------------------------------------------------
        function xyLine(d, ~, ~)
            w = d.WNew(d.wellNo);
            pnts = w.trajectory(:, 1:2);
            pnts = uniquetol(pnts, 'ByRows', true);
            if size(pnts, 1) == 1
                dx = d.Axes1.XLim;
                pnts = [pnts(1)+[dx/5; -dx/5], pnts(2)*ones(2,1)];
            end
            if ~isempty(d.Line3D)&&isvalid(d.Line3D)
                delete(d.Line3D);
            end
            d.Line3D = InteractiveLine('XData', pnts(:,1), 'YData', pnts(:,2), 'Parent', d.Axes1);
            d.Line3D.Callback = @d.crossSection;
        end
        
        %%-----------------------------------------------------------------
        function putXYLine(d, ~, ~)
            xx = d.Axes1.XLim;
            yy = d.Axes1.YLim;
            pnts = [mean(xx) yy(1); mean(xx) yy(2)];
            if ~isempty(d.Line3D)&&isvalid(d.Line3D)
                delete(d.Line3D);
            end
            if ~isempty(d.Line2D)&&isvalid(d.Line2D)
                delete(d.Line2D);
            end
            d.Line3D = InteractiveLine('XData', pnts(:,1), 'YData', pnts(:,2), 'Parent', d.Axes1);
            d.crossSection();
            d.Line3D.Callback = @d.crossSection;
            d.empty2DLine();
            d.Line2D.smoothing = 'pchip';
            d.Axes2.XLimMode = 'auto';
            d.Axes2.YLimMode = 'auto';
        end
        
        %%-----------------------------------------------------------------
        function empty2DLine(d, ~, ~)
            if ~isempty(d.Line2D)&&isvalid(d.Line2D)
                delete(d.Line2D);
            end
            d.Line2D = InteractiveLine('Parent', d.Axes2, 'Callback', '');
            d.Line2D.drawMode = 'along';
        end
            
        function crossSection(d, ~, ~)
            if ~isempty(d.Line3D)&&isvalid(d.Line3D)
                x = d.Line3D.trajectory.XData;
                y = d.Line3D.trajectory.YData;
                if numel(x) >1
                    fac = computeVerticalGridIntersection(d.model.G, [x(:),y(:)]);
                    if ~isempty(d.Patch2D)&&isvalid(d.Patch2D)
                        set(d.Patch2D, 'Vertices', fac.coords2D, ...
                            'Faces', fac.nodes, ...
                            'FaceVertexCData', d.ccol(fac.cellIx));
                    else
                        %if ~isempty(d.Line3D)&&isvalid(d.Line3D)
                        %    delete(d.Patch2D);
                        %end
                        d.Patch2D = patch('Vertices', fac.coords2D, 'Faces', fac.nodes, 'FaceColor', 'flat', 'EdgeColor', [.7 .7 .7], ...
                            'EdgeAlpha', .5, 'FaceVertexCData', d.ccol(fac.cellIx), 'Parent', d.Axes2, 'HitTest', 'off');
                        uistack(d.Patch2D, 'bottom')
                    end
                    d.Axes2.XLimMode = 'auto';
                    d.Axes2.YLimMode = 'auto';
                end
            end
        end
        
        function viewTraj(d, ~, ~)
            traj = getTrajectory(d);
            % exit/entry points only given in unprocessed struct (second output)
            [~, t] = computeTraversedCells(d.model.G, traj); 
            
            fe = figure;hold on
            %caxis([166.9325  336.7537])
            plot3(traj(:,1), traj(:,2), traj(:,3), 'Color', [.4 .4 .4],'LineWidth', 5, 'MarkerSize', 14);
            % cell entry points
            plot3(t.p1(:,1), t.p1(:,2), t.p1(:,3), '.k','LineWidth', 2, 'MarkerSize', 3);
            % cell exit points
            plot3(t.p2(:,1), t.p2(:,2), t.p2(:,3), '.k','LineWidth', 2, 'MarkerSize', 3);
            plotCellData(d.model.G, d.ccol, t.cell, 'FaceAlpha', .5, 'EdgeAlpha', .5, 'EdgeColor', 'r')
            view(3),daspect([1, 1, .2])
            
            % Add extra cell layer
            N = getNeighbourship(d.model.G);
            A = getConnectivityMatrix(N,true,d.model.G.cells.num);
            show = zeros(d.model.G.cells.num,1); show(t.cell)=1;
            K=false(d.model.G.cartDims(3),1); K(1:2:end)=true;
            [ijk{1:3}] = ind2sub(d.model.G.cartDims, d.model.G.cells.indexMap);
            h = plotCellData(d.model.G,d.ccol,(A*A*A*show>0)&(show==0)&(K(ijk{3})),'FaceAlpha',.8);
            if numel(d.WNew)>0
                plotWellData(d.model.G, d.WNew, 'lineplot', true);
            end
        end
        
        function saveWell(d, ~, ~)
             traj = getTrajectory(d);
             %t = computeTraversedCells(d.model.G, traj);
             nnew = numel(d.WNew);
             nm   = ['case_', num2str(nnew+1)];
             info = d.Menu.items{1}.wellInfo;
             %seg = t.vec;
             ncompi = numel(d.W(1).compi);
             if info.sign > 0
                 compi = [1, ones(1, ncompi-1)];
             else
                 compi = ones(1, ncompi)/ncompi;
             end
             
             w = addWellFromTrajectory([], d.model.G, d.model.rock, traj, ...
                    'name', nm, 'type', info.type, 'sign', info.sign, ...
                    'val', info.val, 'compi', compi);

             w.trajectory = traj;
             w.cell_origin = ones(size(w.cells));
             d.WNew = [d.WNew; w];
             d.WellPlotNew{nnew+1} = WellPlotHandle(d.model.G, w, 'Parent', d.Axes1);      
             d.Menu.items{1}.updateWells(d.WNew);
             d.Menu.items{1}.wellPopup.Value = numel(d.WNew);
        end
        
        function layout(d, ~, ~)
            [mw, sp] = deal(d.layoutParams.menuWidth, d.layoutParams.itemSpace);
            fip    = d.Figure.Position;
            mPos   = d.Menu.Position;
            mPos   = [5, fip(4)-mPos(4)-1, mw, mPos(4)];
            aw = floor((fip(3)-mw-4*sp)/2);
            aPos1   = [mw+2*sp, 2*sp, max(10, aw), max(10, fip(4)-3*sp)];
            aPos2   = [mw+3*sp+aw, 2*sp, max(10, aw), max(10, fip(4)-3*sp)];
            d.Menu.Position = mPos;
            d.Axes1.Position = aPos1;
            d.Axes2.Position = aPos2;
        end
        
        function launchDiagnostics(d, ~, ~)
            n = numel(d.WNew)+1;
            [models, wells, names] = deal(cell(1,n));
            if n>1
                % add closed for first well-list
                newIx = numel(d.W) + (1:(n-1));
                W = [d.W; d.WNew];
                [W(newIx).status] = deal(false);
                [W(newIx).cstatus] = deal(false(size(W(end).cstatus)));
                [W(newIx).type] = deal('rate');
                [W(newIx).val] = deal(0);
                [models{1}, wells{1}]  = deal(d.model, W);
                names{1} = 'Base';
                for k = 2:n
                    tmp = W;
                    tmp(newIx(k-1)) = d.WNew(k-1);
                    %W = [d.W; d.WNew(k-1)];
                    [models{k}, wells{k}]  = deal(d.model, tmp);
                    names{k} = ['case ', num2str(k-1)];
                end
                DiagnosticsViewer(models,wells, 'modelNames', names);
            end
        end
        
    end
end

function traj = getTrajectory(d)
pxy = [d.Line3D.trajectory.XData(:), d.Line3D.trajectory.YData(:)];
p   = [d.Line2D.trajectory.XData(:), d.Line2D.trajectory.YData(:)];
lenxy = [0; sqrt(dot(diff(pxy), diff(pxy), 2))];
%traj = [interp1(fac.cumlength, fac.segments, p(:,1)), p(:,2)];
traj = [interp1(cumsum(lenxy), pxy, p(:,1), 'linear', 'extrap'), p(:,2)];
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
