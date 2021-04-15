function showAllocation(d, ax, s2, s3)
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

mIx = s3.modsel.ix;
nM = numel(mIx);
[iIx, pIx] = deal(s3.wsel.injectorIx, s3.wsel.producerIx);

if s2.asel.panelNo == 1
    src   = s2.asel.leftPopup;
    doAvg = s2.asel.leftSwitch && numel(mIx)>1;
else
    src   = s2.asel.rightPopup;
    doAvg = s2.asel.rightSwitch && numel(mIx)>1;
end

if doAvg
    WP_avg = modelAverageWellPairDiagnostics(d, mIx);
end

switch src.Value
    case 1  % None
        axis(ax,'off');

    case 2
        showWellCommunication(d, ax, s3.wsel.getCommunicationStrength());

    case 3  % Injector volumes
        if numel(iIx)~=1
%             text(0,.1,'Injector volumes:','Parent', ax);
%             text(0,0,'Please select{\bf one} injector only','Parent',ax);
%             axis(ax,'off'); return
            d.messageBar.String = 'Please select one injector only';
            return
        else
            d.messageBar.String = '';
        end
        names = d.WellPlot.producerNames;
        if nM==1 || doAvg
            if doAvg
                wp = WP_avg;
            else
                wp = d.Data{mIx}.diagnostics.WP;
            end
            d.plotPie(ax, wp.vols(wp.pairIx(:,1)==iIx), ...
                d.WellPlot.producerNames, ...
                d.Data{mIx(1)}.prodColors);
            title(ax,d.WellPlot.injectorNames{iIx});              
        else
            plotVolsMultipleModels(d,ax,s3,mIx,iIx,pIx,names,'injectors');
        end            

    case 4  % Injector allocation
        if numel(iIx)~=1
%             text(0,.1,'Injector allocation:','Parent', ax);
%             text(0,0,'Please select{\bf one} injector only','Parent',ax);
%             axis(ax,'off'); return
            d.messageBar.String = 'Please select one injector only';  
            return
        else
            d.messageBar.String = '';
        end
        names = d.WellPlot.producerNames;
        names{end+1}='reservoir';
        if nM == 1 || doAvg
            if doAvg
                inj = WP_avg.inj(iIx);
            else
                inj = d.Data{mIx}.diagnostics.WP.inj(iIx);
            end
            d.plotPie(ax, abs(sum([inj.alloc, inj.ralloc],1)), ...
                names, d.Data{mIx(1)}.prodColors);
            title(ax,d.WellPlot.injectorNames{iIx});

        else
            plotAllocsMultipleModels(d,ax,s3,mIx,iIx,pIx,names,'injectors');
        end

    case 5  % PLT plot (production logging tool)
        if numel(iIx)~=1
%             text(0,.1,'Injector allocation:','Parent', ax);
%             text(0,0,'Please select{\bf one} injector only','Parent',ax);
%             axis(ax,'off'); return
            d.messageBar.String = 'Please select one injector only';  
            return
        elseif nM==1 || doAvg
            d.messageBar.String = '';
            if doAvg
                injAlloc = WP_avg.inj(iIx);
            else
                injAlloc = d.Data{mIx}.diagnostics.WP.inj(iIx);
            end
            d.plotPLT(ax, injAlloc, d.Data{mIx(1)}.prodColors);
            wplIx = getWellsForLegend(injAlloc.alloc,0.001);      
            
            legend(ax, ax.Children(d.WellPlot.nProd+1-wplIx), ...
                d.WellPlot.producerNames(wplIx));            
        else
            d.messageBar.String = '';
            [injAlloc, wplIxAll] = deal(cell(nM,1));
            for i = 1:nM
                injAlloc{i} = d.Data{mIx(i)}.diagnostics.WP.inj(iIx);
                wplIxAll{i} = getWellsForLegend(injAlloc{i}.alloc,0.01);                
            end

            d.plotPLT3D(ax, injAlloc, ...
                false, d.Data{mIx(1)}.prodColors);
            
            ll = fliplr([wplIxAll{:}]);
            wplIx = unique([wplIxAll{:}]);
            goIx =arrayfun(@(x) find(ll == x,1),wplIx);

            legend(ax, ax.Children(goIx), ...
                d.WellPlot.producerNames(wplIx), 'Location', 'Best');            

        end
        title(ax,d.WellPlot.injectorNames{iIx});


       

    case 6  % Producer volumes
        if numel(pIx)~=1
%             text(0,.1,'Injector allocation:','Parent', ax);
%             text(0,0,'Please select{\bf one} injector only','Parent',ax);
%             axis(ax,'off'); return
            d.messageBar.String = 'Please select one producer only';  
            return
        else
            d.messageBar.String = '';
        end
        names= d.WellPlot.injectorNames;
        if nM==1 || doAvg
            if doAvg
                wp = WP_avg;
            else
                wp = d.Data{mIx}.diagnostics.WP;
            end
            d.plotPie(ax, wp.vols(wp.pairIx(:,2)==pIx), ...
                names, d.Data{mIx(1)}.injColors);
            title(ax,d.WellPlot.producerNames{pIx});            
        else
            plotVolsMultipleModels(d,ax,s3,mIx,iIx,pIx,names,'producers');
        end
            

    case 7  % Producer allocation
        if numel(pIx)~=1
%             text(0,.1,'Injector allocation:','Parent', ax);
%             text(0,0,'Please select{\bf one} injector only','Parent',ax);
%             axis(ax,'off'); return
            d.messageBar.String = 'Please select one producer only';  
            return
        else   
            d.messageBar.String = '';
        end
        names = d.WellPlot.injectorNames; 
        names{end+1}='reservoir';
        if nM == 1 || doAvg
            if doAvg
                prod = WP_avg.prod(pIx);
            else
                prod = d.Data{mIx}.diagnostics.WP.prod(pIx);
            end
            d.plotPie(ax, abs(sum([prod.alloc, prod.ralloc],1)), ...
                names, d.Data{mIx(1)}.injColors);
            title(ax,d.WellPlot.producerNames{pIx});
        else
            plotAllocsMultipleModels(d,ax,s3,mIx,iIx,pIx,names,'producers');
        end


    case 8  % PLT plot (production logging tool)
        if numel(pIx)~=1
%             text(0,.1,'Injector allocation:','Parent', ax);
%             text(0,0,'Please select{\bf one} injector only','Parent',ax);
%             axis(ax,'off'); return
            d.messageBar.String = 'Please select one producer only';  
            return
        elseif nM==1 || doAvg
            d.messageBar.String = '';
            if doAvg
                prodAlloc = WP_avg.prod(pIx);
            else
                prodAlloc = d.Data{mIx}.diagnostics.WP.prod(pIx);
            end
            d.plotPLT(ax, prodAlloc, d.Data{mIx(1)}.injColors);
            wplIx = getWellsForLegend(prodAlloc.alloc,0.001);
            
            legend(ax, ax.Children(d.WellPlot.nInj+1-wplIx), ...
                d.WellPlot.injectorNames(wplIx));
        else
            d.messageBar.String = '';
            [prodAlloc, wplIxAll] = deal(cell(nM,1));
            for i = 1:nM
                prodAlloc{i} = d.Data{mIx(i)}.diagnostics.WP.prod(pIx);
                wplIxAll{i} = getWellsForLegend(prodAlloc{i}.alloc,0.01);
            end
            
            d.plotPLT3D(ax,prodAlloc, ...
                true, d.Data{mIx(1)}.injColors);
            if nM >10
                ds = round(nM/10);
                 set(ax,'XTick', 1:ds:nM, 'XTickLabel', ...
                     {s3.modsel.modelNames{1:ds:nM}});
            else
                set(ax,'XTick', 1:nM, 'XTickLabel', ...
                     {s3.modsel.modelNames{mIx}});
            end
            set(ax,'XTickLabelRotation',-10,'FontSize',8);
            
            ll = fliplr([wplIxAll{:}]);
            wplIx = unique([wplIxAll{:}]);
            goIx =arrayfun(@(x) find(ll == x,1),wplIx);
            
            
            legend(ax,  ax.Children(goIx), ...
                d.WellPlot.injectorNames(wplIx), 'Location', 'Best');
        end
        title(ax,d.WellPlot.producerNames{pIx});
        
    case 9
        
        percentageConnection =  zeros(size(d.Data{1}.diagnostics.wellCommunication));
        for i = 1:nM
              

           com = d.Data{mIx(i)}.diagnostics.wellCommunication;
           tot = sum(com(:));
           n   = max(size(com));
           com = com > (1/100)*tot/n;

        
            percentageConnection = percentageConnection + (com > 0);
        end
        
        percentageConnection = percentageConnection./nM;   

        showWellCommunicationPercentage(d, ax, percentageConnection);

    otherwise
        disp('functionality not implemented yet');
end
end

function wplIx = getWellsForLegend(data,threshold)
    s = abs(sum(data));
    mm = max(s);
    wplIx = find(s>threshold*mm);
end



function plotVolsMultipleModels(d,ax,s3,mIx,iIx,pIx,names,welltype)

    if strcmp(welltype,'injectors')
        idx = 1;
        t2 = 'producers';
        nWell = d.WellPlot.nProd;
        i1 = iIx;
        nm = d.WellPlot.injectorNames{iIx};
        % i2 = pIx;
        cmap = d.Data{mIx(1)}.prodColors;            
    elseif strcmp(welltype,'producers')
        idx = 2;
        t2 = 'injectors';
        nWell = d.WellPlot.nInj;
        i1 = pIx;
        nm = d.WellPlot.producerNames{pIx};
        % i2 = iIx;
        cmap = d.Data{mIx(1)}.injColors;           
    else
        disp('Unknown welltype')
        return
    end
    
    nM = numel(mIx);
    vols = zeros(nM,numel(names));
    wplIxAll = cell(nM,1);
    for i=1:nM
        wp = d.Data{mIx(i)}.diagnostics.WP; 
        vols(i,:) = wp.vols(wp.pairIx(:,idx)==i1);
        wplIxAll{i} = getWellsForLegend(vols,0.01);           
    end
    
    h=bar(ax, vols,'stacked');
    for i=1:numel(h)
        set(h(i),'FaceColor', cmap(i,:));
    end
    
    set(ax,'XTickLabel', ...
         {s3.modsel.modelNames{mIx}});
    set(ax,'XTickLabelRotation',30,'FontSize',8);
    
    wplIx = unique([wplIxAll{:}]);    
    
%     names= arrayfun(@(x) x.label.String, d.WellPlot.(t2),'UniformOutput',false);
%     names{end+1}='reservoir';
    
    legend(ax,  ax.Children(nWell+1-wplIx), ...
        names{wplIx});            
    axis(ax,'tight')
    set(h,'HitTest','off');   
    title(ax,nm, 'Interpreter','none');            
end

function plotAllocsMultipleModels(d,ax,s3,mIx,iIx,pIx,names,welltype)

    if strcmp(welltype,'injectors')
        t2 = 'producers';
        i1 = iIx;
        % i2 = pIx;
        iflux = 'inj';
        cmap = d.Data{mIx(1)}.prodColors;      
        nwell = d.WellPlot.nProd;
        nm    = d.WellPlot.injectorNames{i1};
    elseif strcmp(welltype,'producers')
        t2 = 'injectors';
        i1 = pIx;
        % i2 = iIx;
        iflux = 'prod';
        cmap = d.Data{mIx(1)}.injColors;    
        nwell = d.WellPlot.nInj;
        nm    = d.WellPlot.producerNames{i1};
    else
        disp('Unknown welltype')
        return
    end

    nM = numel(mIx);
    alloc = zeros(nM,numel(names));
    wplIxAll = cell(nM,1);
    for i=1:nM
        inj = d.Data{mIx(i)}.diagnostics.WP.(iflux)(i1);
        alloc(i,:) = abs(sum([inj.alloc, inj.ralloc],1))';
        wplIxAll{i} = getWellsForLegend(alloc,0.001);        
    end
    
    h=bar(ax, alloc,'stacked');
    for i=1:numel(h)
            set(h(i),'FaceColor', cmap(i,:));
    end
    
    set(ax,'XTickLabel', ...
         {s3.modsel.modelNames{mIx}});
    set(ax,'XTickLabelRotation',30,'FontSize',8);
    
    wplIx = unique([wplIxAll{:}]);
    

    legend(ax,  ax.Children(nwell+2-wplIx), ...
        names{wplIx});     
       
    axis(ax,'tight')    
    title(ax, nm, 'Interpreter','none');            
end
