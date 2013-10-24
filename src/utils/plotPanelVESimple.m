function varargout = plotPanelVESimple(G, Gt, W, sol, t, vol, fluid, fluidADI, varargin)
%Make a panel plot used in the examples of the VE module
%
% SYNOPSIS
%       plotPanelVE(G, Gt, sol, t, vol)
%       plotPanelVE(G, Gt, sol, t, vol, 'pn1', pv1, ...)
%   h = plotPanelVE(...)
%
% PARAMETERS:
%   G       - Grid data structure for 3D grid
%   Gt      - Grid data structure for topsurface grid
%   W       - Well
%   sol     - Solution data structure
%   t       - Time
%   vol     - Vector of (free, trapped, total) volumes
%   opt     - Various options
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.
%
% DESCRIPTION:
%   The function will make a composit plot that consists of several parts
%    - a 3D plot of the plume
%    - a pie chart of trapped versus free volume
%    - a plane view of the plume from above
%    - two cross-sections in the x/y directions through the well
%
%   Supported options for the 3D plot of saturations
%      'maxH'   -- Double, maximal height value used for scaling the
%                  colorbar.  Default value: 10
%      'view'   -- Vector, view angle. Default value: [30,60]
%      'wireS'  -- Bool, if true, plot only positive saturation values
%                  on top of a wireframe grid. Default value: false
%      'Wadd'   -- Double, height added to the well. Default value: 100
%      'Saxis'  -- Vector, scaling of caxis. Default value: [0 1]
%
%   Supported options for the plot of CO2 height:
%      'wireH'  -- Bool, if true, plot only positive saturation values
%                  on top of a wireframe grid. Otherwise, make a color plot
%                  with height shown as contours. Default value: false
%      'ptol'   -- Double, smallest number that is considered as positive
%                  in the wireframe plot. Default value: 1e-3
%
%   Supported options for the slice plots:
%      'slice'  -- Vector, index along which to pick the two 2D slices
%                  Default value: slice=[10 10]
%      'plotPlume' -- Boolean indicating if the plume should be plotted as
%                     a sharp interface.
%
%   In addition, the function may plot the CO2 inventory as a function of
%   time if the user clicks the pie chart or if the option 'plotHist' is
%   set to true.
%
% SEE ALSO:
%   runIGEMS

persistent G_x G_y;      % Grids for 2D slices
persistent h h1 h2 h3 h_area;   % Figure and graphics handles
persistent x y;          % Node coordinates
persistent cells_sub_x;  % Cell numbers in x-slice
persistent cells_sub_y;  % Cell numbers in y-slice
persistent cells_sub_gt_x;  % Cell numbers in x-slice
persistent cells_sub_gt_y;  % Cell numbers in y-slice
persistent volHistory;   % Pie chart historic values
persistent thistory;     % Previous timesteps

prm = struct('maxH',      10, ...
    'view',      [30, 60],...
    'wireS',     false, ...
    'wireH',     false, ...
    'Wadd',      100, ...
    'Saxis',     [0 1], ...
    'ptol',      1e-3, ...
    'slice',     [10 10],...
    'plotPlume', false,  ...
    'plotHist',  true);

prm = merge_options(prm, varargin{:});

%% Prepare for subsequent plotting
if (t==0)
    clear ijk
    [ij{1:2}] = ind2sub(Gt.cartDims, Gt.cells.indexMap);
    [ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap);
    ijk = [ijk{:}];
    ij = [ij{:}];
    
    cells_sub_x = find(ijk(:,2) == prm.slice(2));
    cells_sub_gt_x = find(ij(:,2) == prm.slice(2));
    G_x = extractSubgrid(G, cells_sub_x);
    G_x = mcomputeGeometry(G_x);
    G_x.cartDims = G.cartDims;
    
    cells_sub_y = find(ijk(:,1) == prm.slice(1));
    cells_sub_gt_y = find(ij(:,1) == prm.slice(1));
    G_y = extractSubgrid(G, cells_sub_y);
    G_y = mcomputeGeometry(G_y);
    G_y.cartDims = G.cartDims;
    
    scrsz  = get(0,'ScreenSize');
    fsize = min([scrsz(3:4); 1024 700]);
    h = figure('Position', [scrsz(3)-fsize(1), scrsz(4)-fsize(2)-75, fsize]);
    if nargout > 0; varargout{1} = h; end;
    
    set(h, 'color', 'white');
    
    if prod(Gt.cartDims)==Gt.cells.num
        x = reshape(Gt.nodes.coords(:,1),Gt.cartDims+1);
        y = reshape(Gt.nodes.coords(:,2),Gt.cartDims+1);
        %x = reshape(Gt.cells.centroids(:,1),Gt.cartDims);
        %y = reshape(Gt.cells.centroids(:,2),Gt.cartDims);
        %zt = reshape(Gt.cells.z,Gt.cartDims);
        %zb = reshape(Gt.cells.z+Gt.cells.H,Gt.cartDims);
    else
        h2 = [];
    end
    volHistory = vol;
    thistory   = 0;
    if prm.plotHist
        h3 = figure();
    else
        h3 = [];
    end
    h1 = figure();
    h_area=figure();
end


%% Height of the CO2 column
set(0,'CurrentFigure', h);
%{
subplot(2,3,2);
if t==0
   if prm.wireS
      plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
   else
      h1 = plotCellData(G, sol.s, 'EdgeColor','w','EdgeAlpha', .15);
   end
   plotGrid(G_x, 'FaceColor', 'none', 'EdgeColor', 'c');
   plotGrid(G_y, 'FaceColor', 'none', 'EdgeColor', 'm');
   plotWell(G, W, 'height', prm.Wadd, 'color', 'r');
   axis tight off, view(prm.view);
else
   if ishandle(h1), delete(h1); end
   if prm.wireS
      h1 = plotCellData(G, sol.s, sol.s>prm.ptol, 'EdgeAlpha', 0.1);
   else
      h1 = plotCellData(G, sol.s, 'EdgeColor','w','EdgeAlpha', .15);
   end
end
title(['CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
caxis(prm.Saxis);
%}
%% Pie chart of trapped and free CO2
subplot(2,2,2), cla
if numel(vol)==3
    str1 = ['trapped ' num2str(round(vol(1)/vol(3)*100)) ' %'];
    str2 = ['free ' num2str(round(vol(2)/vol(3)*100)) ' %'];
    ph = pie([max(vol(1),eps) max(vol(2),eps)],{str1, str2});
    title(['Total mass : ', ...
        num2str(round(convertTo(vol(3),mega*meter^3))),' M m^3']);
elseif numel(vol)==7
    vplot = max([vol(1:6) vol(7)-sum(vol(1:6))], eps);
    ph = pie(vplot);
    % trick auto-shrink by first placing legend below and then moving it
    lnames = {'Disolved','Struct. residual', 'Residual', 'Free residual', 'Struct. movable', ...
        'Free movable','Leaked'};
    hl = legend(lnames, 'Location','SouthOutside', 'orientation','horizontal');
    set(hl,'Location','SouthEastOutside','orientation','vertical');
    pl = get(hl, 'Position'); set(hl, 'Position', [.85 pl(2)-.035 pl(3:4)]);
    title(['Total injected mass: ', ...
        num2str(round(convertTo(vol(7),1e3*mega))),' M tonn']);
    
    if t ~= 0
        volHistory = [volHistory; vplot];
        thistory = [thistory; t];
        set(ph, 'ButtonDownFcn', @(src, event) makePlotHistory);
        plotHistory()
        legend(lnames{:})
        set(0,'CurrentFigure', h);
    end
else
    disp('Incompatible volume vector. No pie chart created');
end



%% Height of the CO2 plume plotted on the 2D top-surface grid
%subplot(3,2,[1,3])
subplot(1,2,1)
if prm.wireH
    if t==0
        plotGrid(Gt, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
        axis tight off, set(gca,'zdir','reverse'); view([-90, 90])
        colorbar, caxis([0 prm.maxH]);
    else
        if ishandle(h2), delete(h2); end
        h2 = plotCellData(Gt, sol.h, sol.h(:,1)>prm.ptol,'EdgeAlpha', 0.1);
    end
else
    cla
    cvec=sol.h;
    if(~(prod(Gt.cartDims)==Gt.cells.num))        
        cvec(cvec<0.1)=NaN;
    end
    makepColorPlot(cvec);    
    hold off, 
end
title(['Height of CO2-column at time',num2str(t/year),' year']);

%% Saturation along the two slices
subplot(4,4,11:12); cla % x-slice
if prm.plotPlume
    plotPlume(G, Gt, sol.h, cells_sub_x,'EdgeColor','w','EdgeAlpha',.1);
    plotGrid(G, cells_sub_x,'EdgeColor', 'k','EdgeAlpha',.1, 'FaceColor', 'none');
    caxis([0 prm.maxH]);
    axis tight off, view([0 0] + prm.plotPlume*4), title('CO2-saturation, x-slice (cyan)');
else
    plotSlizeHeight(cells_sub_gt_x,1);
end



subplot(4,4,15:16); cla % yslice
if prm.plotPlume
    plotPlume(G, Gt, sol.h, cells_sub_y,'EdgeColor','w','EdgeAlpha',.1);
    plotGrid(G, cells_sub_y,'EdgeColor','k','EdgeAlpha',.1, 'FaceColor', 'none');
    caxis([0 prm.maxH]);
    axis tight off, view([90 0] + prm.plotPlume*4), title('CO2-saturation, y-slice (magenta)');
else
    plotSlizeHeight(cells_sub_gt_y,2);
end
set(0,'CurrentFigure',h_area)
subplot(1,3,1),cla,
cvec=sol.h;
if(~(prod(Gt.cartDims)==Gt.cells.num))   
    cvec(cvec<0.1)=NaN;
end
makepColorPlot(cvec)
%{
hplot = reshape(sol.h(:,1),Gt.cartDims);
hplot = hplot([1:end end],:);
hplot = hplot(:,[1:end end]);
pcolor(x,-y,hplot); shading flat; caxis([0 prm.maxH])
hold on;
colorbar
contour(x,-y,reshape(Gt.nodes.z,Gt.cartDims+1),5,'w');
%}
subplot(1,3,2),cla
cvec=sol.h_max;
if(~(prod(Gt.cartDims)==Gt.cells.num))   
    cvec(cvec<0.1)=NaN;
end
makepColorPlot(cvec);
if(isfield(fluidADI,'dis_max'))
    rsH=Gt.cells.H.*(1-sol.s).*sol.rs/fluidADI.dis_max;
else
    rsH=zeros(Gt.cells.num,1);
end
subplot(1,3,3),cla
cvec=rsH;
if(~(prod(Gt.cartDims)==Gt.cells.num))   
    cvec(cvec<1)=NaN;
end
makepColorPlot(cvec);



set(0,'CurrentFigure', h1);
head=sol.pressure-fluid.rho(2)*norm(gravity)*(Gt.cells.z-mean(Gt.cells.z));
subplot(1,2,1)
%press = reshape(head,Gt.cartDims)/barsa;
%press = press([1:end end],:);
%press = press(:,[1:end end]);
makepColorPlot(head);
%pcolor(x,-y,press); shading flat; %caxis([min])
%colorbar
%shading flat;
title('Water head at mean depth')


subplot(4,4,11:12); cla % x-slice
xx=Gt.cells.centroids(cells_sub_gt_x,1);
plot(xx,head(cells_sub_gt_x)/barsa);
hold off;
axis tight;
box on;
subplot(4,4,15:16)
xx=Gt.cells.centroids(cells_sub_gt_y,2);
plot(xx,head(cells_sub_gt_y)/barsa);
hold off;
axis tight;
box on;
%%


%if(t>0)
%drawnow
%end
    function makePlotHistory()
        if isempty(h3)
            h3 = NaN;
        end
        
        if ~ishandle(h3)
            h3 = figure;
        end
        plotHistory();
    end

    function plotHistory()
        if ishandle(h3)
            set(0, 'CurrentFigure', h3);
            %area(convertTo(thistory, year), convertTo(volHistory, mega*meter^3));
            plot(convertTo(thistory, year), convertTo(volHistory, mega*meter^3));
            xlabel('Years since simulation start');
            ylabel('Volume (M m^3)')
            axis tight
            %legend(lnames, 'location', 'eastoutside', 'orientation', 'vertical');
        end
    end

    function makepColorPlot(cvec)
        if prod(Gt.cartDims)==Gt.cells.num
         cvec = reshape(cvec,Gt.cartDims);
         cvec = cvec([1:end end],:);
         cvec = cvec(:,[1:end end]);
         pcolor(x,-y,cvec); shading flat; caxis([0 prm.maxH])
         %axis auto;
         hold on;
         colorbar
            contour(x,-y,reshape(Gt.nodes.z,Gt.cartDims+1),5,'w');
         hold off;
         axis on;
         axis tight off; axis equal
         box on;
         
        else
            if(false)
          plotFaces(Gt.parent,Gt.cells.map3DFace,cvec,'EdgeColor','none');
          plotFaces(Gt.parent,Gt.cells.map3DFace,'FaceColor','none','EdgeAlpha',0.1)          
            else
                plotFaces(Gt.parent,Gt.cells.map3DFace,cvec);
          plotFaces(Gt.parent,Gt.cells.map3DFace,'FaceColor','none')  
            end
          %plotCellData(Gt,cvec);          
          axis tight off;
          colorbar;
          view(3)
          
          %view(2)  
        end
        
    end
%}

    function plotSlizeHeight(cells_sub,dir)
        xx=Gt.cells.centroids(cells_sub,dir);
        zt=Gt.cells.z(cells_sub);
        zb=Gt.cells.z(cells_sub)+Gt.cells.H(cells_sub);
        zco=zt+sol.h(cells_sub);
        %if(isfield(sol,'sGmax'))
        %   zco_max=Gt.cells.z(cells_sub)+Gt.cells.H(cells_sub).*sol.sGmax(cells_sub);
        %else
        zco_max=zt+sol.h_max(cells_sub);
        %end
        
        if(isfield(sol,'rsH'))
            %rsH=Gt.cells.H.*(1-sol.s).*sol.rs/fluidADI.dis_max;
            rsH=sol.rsH;
            %assert(all(rsH<=Gt.cells.H.*(1-sol.s)))
            rsH=rsH(cells_sub,1);
        else
            rsH=zco_max-zt;
        end
        zrs = zt+rsH;
        plot(xx,zt,'k','LineWidth',2)
        hold on;
        plot(xx,zb,'k','LineWidth',2)
        if(true)
            plot(xx,zrs,'g','LineWidth',1)
            plot(xx,zco_max,'b','LineWidth',1)
            plot(xx,zco,'r','LineWidth',1)
            %assert(all(zco_max>=zco));
            set(gca,'Ydir','reverse')
            set(gca,'LineWidth',2)
            hold off;
            axis tight;
            box on;
        else
            patch(xx([1:end end:-1:1]), ...
                [zco; zt(end:-1:1)], 'g')
            patch(xx([1:end end:-1:1]), ...
                [zco; zco_max(end:-1:1)],[1 1 0])
            patch(xx([1:end end:-1:1]), ...
                [zco_max; zrs(end:-1:1)], [0 0.1 1])
            %patch(xx([1:end end:-1:1]), ...
            %  [zrs; zb(end:-1:1)], [0 0 0.6])
            set(gca,'Ydir','reverse')
            set(gca,'LineWidth',2)
            axis tight;
            hold off
            box on
        end
    end
end
