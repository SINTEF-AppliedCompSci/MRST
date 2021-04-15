function showAllocation(d, src, ax, s2, s3)
%Undocumented Utility Function

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

ts = s3.tsel.ix;
[iIx, pIx] = deal(s3.wsel.injectorIx, s3.wsel.producerIx);

if s2.asel.panelNo == 1
    src   = s2.asel.leftPopup;
    doAvg = s2.asel.leftSwitch && numel(ts)>1;
else
    src   = s2.asel.rightPopup;
    doAvg = s2.asel.rightSwitch && numel(ts)>1;
end

if doAvg
    WP_avg = timeAverageWellPairDiagnostics(d, ts);
end

switch src.Value
    case 1  % None
        axis(ax,'off');

    case 2
        showWellCommunication(d, ax, s3.wsel.getCommunicationStrength());

    case 3  % Injector volumes
        if numel(iIx)~=1
            text(0,.1,'Injector volumes:','Parent', ax);
            text(0,0,'Please select{\bf one} injector only','Parent',ax);
            axis(ax,'off'); return
        end
        if numel(ts)==1 || doAvg
            if doAvg
                wp = WP_avg;
            else
                wp = d.Data.diagnostics(ts).WP;
            end
%             d.plotPie(ax, wp.vols(wp.pairIx(:,1)==iIx), ...
%                 arrayfun(@(x) x.label.String, ...
%                 d.WellPlot.producers, 'UniformOutput',false), ...
%                 d.Data.prodColors);
              d.plotPie(ax, wp.vols(wp.pairIx(:,1)==iIx), ...
                d.WellPlot.producerNames, d.Data.prodColors);
        else
            %names= arrayfun(@(x) x.label.String, d.WellPlot.producers,'UniformOutput',false);   
            names= d.WellPlot.producerNames;    
            vols = zeros(numel(names),numel(ts));
            for i=1:numel(ts)
                wp = d.Data.diagnostics(ts(i)).WP;
                vols(:,i) = wp.vols(wp.pairIx(:,1)==iIx);
            end
            h = area(ax,vols');
            for i=1:numel(h)
                set(h(i),'FaceColor', d.Data.prodColors(i,:));
            end
            if numel(ts)>10
                ds = round(numel(ts)/10);
                ts = ts(1:ds:numel(ts));
            end
            set(ax,'XTick', ts, 'XTickLabel', ...
                {datestr(d.Data.time.cur(ts) , 'mm.dd.yy')});
            set(ax,'XTickLabelRotation',30,'FontSize',8);
            axis(ax,'tight')
            set(h,'HitTest','off');        
        end
        %title(ax,d.WellPlot.injectors(iIx).label.String, 'Interpreter','none');
        title(ax,d.WellPlot.injectorNames(iIx), 'Interpreter','none');
    case 4  % Injector allocation
        if numel(iIx)~=1
            text(0,.1,'Injector allocation:','Parent', ax);
            text(0,0,'Please select{\bf one} injector only','Parent',ax);
            axis(ax,'off'); return
        end
        names = d.WellPlot.producerNames;    
        %arrayfun(@(x) x.label.String, d.WellPlot.producers,'UniformOutput',false);
        names{end+1}='reservoir';
        if numel(ts)==1 || doAvg
            if doAvg
                inj = WP_avg.inj(iIx);
            else
                inj = d.Data.diagnostics(ts).WP.inj(iIx);
            end
            d.plotPie(ax, abs(sum([inj.alloc, inj.ralloc],1)), ...
                names, d.Data.prodColors);
        else
            alloc = zeros(numel(names),numel(ts));
            for i=1:numel(ts)
                inj = d.Data.diagnostics(ts(i)).WP.inj(iIx);
                alloc(:,i) = abs(sum([inj.alloc, inj.ralloc],1))';
            end
            h = area(ax,alloc');
            for i=1:numel(h)
                set(h(i),'FaceColor', d.Data.prodColors(i,:));
            end
            if numel(ts)>10
                ds = round(numel(ts)/10);
                ts = ts(1:ds:numel(ts));
            end
            set(ax,'XTick', ts, 'XTickLabel', ...
                {datestr(d.Data.time.cur(ts) , 'mm.dd.yy')});
            set(ax,'XTickLabelRotation',30,'FontSize',8);
            axis(ax,'tight')
            set(h,'HitTest','off');
        end
        %title(ax,d.WellPlot.injectors(iIx).label.String);
        title(ax,d.WellPlot.injectorNames(iIx));
    case 5  % PLT plot (production logging tool)
        if numel(iIx)~=1
            text(0,.1,'Injector profile:','Parent', ax);
            text(0,0,'Please select{\bf one} injector','Parent',ax);
            axis(ax,'off'); return
        elseif numel(ts)==1 || doAvg
            if doAvg
                injAlloc = WP_avg.inj(iIx);
            else
                injAlloc = d.Data.diagnostics(ts).WP.inj(iIx);
            end
            d.plotPLT(ax, injAlloc, d.Data.prodColors);
        else
            d.plotPLT3D(ax, arrayfun(@(x) x.WP.inj(iIx),...
                d.Data.diagnostics(ts),'UniformOutput',false), ...
                false, d.Data.prodColors);
            if numel(ts)>10
                ds = round(numel(ts)/10);
                ts = ts(1:ds:numel(ts));
            end
            set(ax,'XTick', ts, 'XTickLabel', ...
                {datestr(d.Data.time.cur(ts) , 'mm.dd.yy')});
            set(ax,'XTickLabelRotation',-30,'FontSize',8);
        end
        %title(ax,d.WellPlot.injectors(iIx).label.String);
        title(ax,d.WellPlot.injectorNames(iIx));
        if ~isempty(pIx)
%             legend(ax, ax.Children(numel(d.WellPlot.producers)+1-pIx), ...
%                 arrayfun(@(x) x.label.String, d.WellPlot.producers(pIx),'UniformOutput',false));
              legend(ax, ax.Children(d.WellPlot.nProd+1-pIx), d.WellPlot.producerNames(pIx));
        end

    case 6  % Producer volumes
        if numel(pIx)~=1
            text(0,.1,'Producer volumes:','Parent', ax);
            text(0,0,'Please select{\bf one} producer only','Parent', ax);
            axis(ax,'off'); return
        end
        if numel(ts)==1 || doAvg
            if doAvg
                wp = WP_avg;
            else
                wp = d.Data.diagnostics(ts).WP;
            end
%             d.plotPie(ax, wp.vols(wp.pairIx(:,2)==pIx), ...
%                 arrayfun(@(x) x.label.String, ...
%                 d.WellPlot.injectors,'UniformOutput',false), ...
%                 d.Data.injColors);
            d.plotPie(ax, wp.vols(wp.pairIx(:,2)==pIx), ...
                d.WellPlot.injectorNames, d.Data.injColors);
        else
            %names= arrayfun(@(x) x.label.String, d.WellPlot.injectors,'UniformOutput',false);    
            names = d.WellPlot.injectorNames;
            vols = zeros(numel(names),numel(ts));
            for i=1:numel(ts)
                wp = d.Data.diagnostics(ts(i)).WP;
                vols(:,i) = wp.vols(wp.pairIx(:,2)==pIx);
            end
            h = area(ax,vols');
            for i=1:numel(h)
                set(h(i),'FaceColor', d.Data.injColors(i,:));
            end
            if numel(ts)>10
                ds = round(numel(ts)/10);
                ts = ts(1:ds:numel(ts));
            end
            set(ax,'XTick', ts, 'XTickLabel', ...
                {datestr(d.Data.time.cur(ts) , 'mm.dd.yy')});
            set(ax,'XTickLabelRotation',30,'FontSize',8);
            axis(ax,'tight')
            set(h,'HitTest','off');        
        end
            
        %title(ax,d.WellPlot.producers(pIx).label.String);
        title(ax,d.WellPlot.producerNames(pIx));
    case 7  % Producer allocation
        if numel(pIx)~=1
            text(0,.1,'Producer allocation:','Parent', ax);
            text(0,0,'Please select{\bf one} producer only','Parent', ax);
            axis(ax,'off'); return
        end
        %names = arrayfun(@(x) x.label.String, d.WellPlot.injectors,'UniformOutput',false);
        names = d.WellPlot.injectorNames;
        names{end+1}='reservoir';
        if numel(ts)==1 || doAvg
            if doAvg
                prod = WP_avg.prod(pIx);
            else
                prod = d.Data.diagnostics(ts).WP.prod(pIx);
            end
            d.plotPie(ax, abs(sum([prod.alloc, prod.ralloc],1)), ...
                names, d.Data.injColors);
        else
            alloc = zeros(numel(names),numel(ts));
            for i=1:numel(ts)
                prod = d.Data.diagnostics(ts(i)).WP.prod(pIx);
                alloc(:,i) = abs(sum([prod.alloc, prod.ralloc],1))';
            end
            h = area(ax,alloc');
            for i=1:numel(h)
                set(h(i),'FaceColor', d.Data.injColors(i,:));
            end
            if numel(ts)>10
                ds = round(numel(ts)/10);
                ts = ts(1:ds:numel(ts));
            end
            set(ax,'XTick', ts, 'XTickLabel', ...
                {datestr(d.Data.time.cur(ts) , 'mm.dd.yy')});
            set(ax,'XTickLabelRotation',30,'FontSize',8);
            axis(ax,'tight')
            set(h,'HitTest','off');
       end
       %title(ax,d.WellPlot.producers(pIx).label.String);
       title(ax,d.WellPlot.producerNames(pIx));
    case 8  % PLT plot (production logging tool)
        if numel(pIx)~=1
            text(0,.1,'Producer profile:','Parent', ax);
            text(0,0,'Please select{\bf one} producer','Parent', ax);
            axis(ax,'off'); return
        elseif numel(ts)==1 || doAvg
            if doAvg
                prodAlloc = WP_avg.prod(pIx);
            else
                prodAlloc = d.Data.diagnostics(ts).WP.prod(pIx);
            end
            d.plotPLT(ax, prodAlloc, d.Data.injColors);
        else
            d.plotPLT3D(ax,arrayfun(@(x) x.WP.prod(pIx),...
                d.Data.diagnostics(ts),'UniformOutput',false), ...
                true, d.Data.injColors);
            if numel(ts)>10
                ds = round(numel(ts)/10);
                ts = ts(1:ds:numel(ts));
            end
            set(ax,'XTick', ts, 'XTickLabel', ...
                {datestr(d.Data.time.cur(ts) , 'mm.dd.yy')});
            set(ax,'XTickLabelRotation',-30,'FontSize',8);
        end
        %title(ax,d.WellPlot.producers(pIx).label.String);
        title(ax,d.WellPlot.producerNames(pIx));
        if ~isempty(iIx)
%             legend(ax, ax.Children(numel(d.WellPlot.injectors)+1-iIx), ...
%                 arrayfun(@(x) x.label.String, d.WellPlot.injectors(iIx),'UniformOutput',false));
            legend(ax, ax.Children(d.WellPlot.nInj+1-iIx), d.WellPlot.injectorNames(iIx));
        end

    otherwise
        disp('functionality not implemented yet');
end
end

