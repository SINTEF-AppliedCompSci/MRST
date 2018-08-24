function showAllocation(d, src, ax, s2, s3)
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
        
    case 2  % Injector volumes
        if (~doAvg && numel(ts)~=1) || numel(iIx)~=1
            text(0,.1,'Injector volumes:','Parent', ax);
            text(0,0,'Please select{\bf one} injector and{\bf one} time step (or select average)','Parent',ax);
            axis(ax,'off'); return
        end
        if doAvg
            wp = WP_avg;
        else
            wp = d.Data.diagnostics(ts).WP;
        end
        d.plotPie(ax, wp.vols(wp.pairIx(:,1)==iIx), ...
            arrayfun(@(x) x.label.String, ...
            d.WellPlot.producers, 'UniformOutput',false), ...
            d.Data.prodColors);
        title(ax,d.WellPlot.injectors(iIx).label.String, 'Interpreter','none');
        
    case 3  % Injector allocation
        if numel(iIx)~=1
            text(0,.1,'Injector allocation:','Parent', ax);
            text(0,0,'Please select{\bf one} injector only','Parent',ax);
            axis(ax,'off'); return
        end
        names= arrayfun(@(x) x.label.String, d.WellPlot.producers,'UniformOutput',false);
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
            axis tight
        end
        title(ax,d.WellPlot.injectors(iIx).label.String);
        
    case 4  % PLT plot (production logging tool)
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
        end
        title(ax,d.WellPlot.injectors(iIx).label.String);
        legend(ax, ax.Children(numel(d.WellPlot.producers)+1-pIx), ...
            arrayfun(@(x) x.label.String, d.WellPlot.producers(pIx),'UniformOutput',false));
        
    case 5  % Producer volumes
        if (~doAvg && numel(ts)~=1) || numel(pIx)~=1
            text(0,.1,'Producer volumes:','Parent', ax);
            text(0,0,'Please select{\bf one} producer and{\bf one} time step (or select average)','Parent', ax);
            axis(ax,'off'); return
        end
        if doAvg
            wp = WP_avg;
        else
            wp = d.Data.diagnostics(ts).WP;
        end
        d.plotPie(ax, wp.vols(wp.pairIx(:,2)==pIx), ...
            arrayfun(@(x) x.label.String, ...
            d.WellPlot.injectors,'UniformOutput',false), ...
            d.Data.injColors);
        title(ax,d.WellPlot.producers(pIx).label.String);
        
    case 6  % Producer allocation
        if numel(pIx)~=1
            text(0,.1,'Producer allocation:','Parent', ax);
            text(0,0,'Please select{\bf one} producer only','Parent', ax);
            axis(ax,'off'); return
        end
        names = arrayfun(@(x) x.label.String, d.WellPlot.injectors,'UniformOutput',false);
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
            axis tight
        end
        title(ax,d.WellPlot.producers(pIx).label.String);
        
    case 7  % PLT plot (production logging tool)
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
        end
        title(ax,d.WellPlot.producers(pIx).label.String);
        legend(ax, ax.Children(numel(d.WellPlot.injectors)+1-iIx), ...
            arrayfun(@(x) x.label.String, d.WellPlot.injectors(iIx),'UniformOutput',false));
        
    otherwise
        disp('functionality not implemented yet');
end
end