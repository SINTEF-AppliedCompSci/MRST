%% Visualisation of high-resolution simulations
% This script produces the plots of high-resolution simulations of the two
% slug-injection scenarios from the EOR book chapter. It assumes that you
% have conduced simulations with the setup from the SP_1D_hires.DATA input
% file and transformed the results as follows:
% x  = model.G.cells.centroids(:,1);
% t  = reports.ReservoirTime/day;
% s  = cellfun(@(x) x.s(:,1),  states, 'UniformOutput',false); S =horzcat(S{:});
% cp = cellfun(@(x) x.cp(:,1), states, 'UniformOutput',false); Cp=horzcat(Cp{:});
% cs = cellfun(@(x) x.cs(:,1), states, 'UniformOutput',false); Cs=horzcat(Cs{:});
% save slug1D.mat x t s cp cs ;

slugs = 2; % For Strategy 2. Set slugs=1 for Strategy 1.
tind  = [750 1250 2000]; % This corresponds to times 7.5, 12.5, and 20

%% Plot solution in (x,t) space
figure('Position',[400 400 640 280])
axes('Position',[.065 0.11 .90 .815]);
cols = lines(6);
i = 1:2:2500;
hs=pcolor(x(i),t(i),s(i,i)'); shading flat, colormap(flipud(winter))

hold on
[is, ip] = findStatesSlug(cs(:,end),cp(:,end), slugs);
hcs=plot(x(is.ix'), is.y','Color',cols(4,:),'LineWidth',2);
hcp=plot(x(ip.ix'), ip.y','Color',cols(3,:),'LineWidth',2);
plot(repmat([0; 50],1,3), t(repmat([tind],2,1)),'--w','LineWidth',1);
hold off
axis([0 50 0 30]);
legend([hs,hcs(1),hcp(1)],'S(x,t)','c_s(x,t)','c_p(x,t)','Location','SouthEast');
args = {'FontSize',14,'FontWeight','bold','VerticalAlignment','top'};
inx = unique(is.ix(:,2));
text(x(inx), repmat(t(end),numel(inx),1), ...
    num2str((numel(inx):-1:1)'), 'Color', cols(4,:),...
    'HorizontalAlignment','left', args{:});
inx = unique(ip.ix(:,2));
text(x(inx), repmat(t(end),numel(inx),1), ...
    num2str((numel(inx):-1:1)'), 'Color', cols(3,:),...
    'HorizontalAlignment','right', args{:});

hb = colorbar('Location','NorthOutside','Direction','reverse');
set(gca,'FontSize',14,'YTick',[])

%% Plot solution at different times in physcal and state space
for i=tind
    figure('Position',[400 400 640 280])
    axes('Position',[.065 .11 .515 .815]);
    hold on
    plot(x,s(:,i),    'Color',cols(1,:),'LineWidth',2);
    plot(x,cs(:,i)/50,'Color',cols(4,:),'LineWidth',2);
    plot(x,cp(:,i)/3, 'Color',cols(3,:),'LineWidth',2);
    [is, ip] = findStatesSlug(cs(:,i), cp(:,i), slugs);
    args = {'FontSize',14,'VerticalAlignment','bottom','FontWeight','bold'};
    if size(is.ix,2)>1
        inx = unique(is.ix(:,2));
        plot(x(inx), cs(inx,i)/50, '.', 'Color', cols(4,:),'MarkerSize', 18)
        text(x(inx), cs(inx,i)/50, ...
            num2str((numel(inx):-1:1)'), 'Color', cols(4,:),...
            'HorizontalAlignment','left', args{:});
    end
    if size(ip.ix,2)>1
        inx = unique(ip.ix(:,2));
        plot(x(inx), cp(inx,i)/3, '.', 'Color', cols(3,:),'MarkerSize', 18)
        text(x(inx), cp(inx,i)/3, ...
            num2str((numel(inx):-1:1)'), 'Color', cols(3,:),...
            'HorizontalAlignment','right', args{:});
    end
    hold off
    axis([0 25 -.05 1.05])
    legend('S(x,t)','c_s(x,t)','c_p(x,t)');
    % title(num2str(t(i)));
    set(gca,'FontSize',14,'YTick',[])
    
    axes('Position',[.65 .11 .31 .815]);
    hold on
    %cplot3(s(:,i),cs(:,i),cp(:,i),x);
    plot3(s(:,i),cs(:,i),cp(:,i),'LineWidth',2);
    args = {'FontSize',14,'VerticalAlignment','bottom'};
    if size(is.ix,2)>1
        inx = unique(is.ix(:,2));
        plot3(s(inx,i), cs(inx,i), cp(inx,i), '.', 'Color', ...
            cols(4,:),'MarkerSize', 24)
        text(s(inx,i), cs(inx,i), cp(inx,i), ...
            num2str((numel(inx):-1:1)'), 'Color', cols(4,:),...
            'HorizontalAlignment','left', args{:});
    end
    if size(ip.ix,2)>1
        inx = unique(ip.ix(:,2));
        plot3(s(inx,i), cs(inx,i), cp(inx,i), ...
            '.', 'Color', cols(3,:),'MarkerSize', 18)
        text(s(inx,i), cs(inx,i), cp(inx,i), ...
            num2str((numel(inx):-1:1)'), 'Color', cols(3,:),...
            'HorizontalAlignment','right', args{:});
    end
    hold off, view(3)
    % title(num2str(t(i)));
    axis([-.05 1.05 -1 51 -.1 3.1]), box on
    set(gca,'FontSize',14,'XTick',[],'YTick',[],'ZTick',[])
    xlabel('S'), ylabel('c_s'), zlabel('c_p','Rotation',0)
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
