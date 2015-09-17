function [ hfig ] = makeSideProfilePlots_CO2heights( G, Gt, wellCellIndex, states, fluid)

% TODO: make this plotting function more generic, so any cell index can be
% passed in which the user wants to make two slices through (in the
% direction parallel to the x and y axis).


% Cross-sectional slices through injection point:
% first get index of injection well:
[ii,jj] = ind2sub(Gt.cartDims, wellCellIndex);
disp(['Well cell index of ',num2str(wellCellIndex),' corresponds to I=',num2str(ii),', J=',num2str(jj)])

ijk = gridLogicalIndices(G);

ii = ijk{1}(wellCellIndex);
jj = ijk{2}(wellCellIndex);
kk = ijk{3}(wellCellIndex);
disp(['Well cell index of ',num2str(wellCellIndex),' corresponds to I=',num2str(ii),', J=',num2str(jj)])

% % plot of Saturations (total CO2 --- included both free and residual):
% figure; set(gcf,'Position',[1 1 2200 500])
% subplot(2,1,1)
% plotCellData(G, states{end}.s(:,2), ijk{1} == ii); view([90 0])
% title({'Saturation';['Slice through i=',num2str(ii)];'South -- North (facing West)'})
% ylabel('y-axis, meters'); zlabel('depth, meters')
% zlim([800, 850])
% set(gca,'DataAspect',[1 1 1/10])
% 
% subplot(2,1,2)
% plotCellData(G, states{end}.s(:,2), ijk{2} == jj); view([0 0])
% title({'Saturation';['Slice through j=',num2str(jj)];'West -- East (facing North)'})
% xlabel('x-axis, meters'); zlabel('depth, meters')
% zlim([800, 850])
% set(gca,'DataAspect',[1 1 1/10])


% get the free plume height, and the max plume height reached in past
states = addCO2HeightData(states, Gt, fluid);



% % plot the free plume height:
% figure; set(gcf,'Position',[1 1 2200 500])
% subplot(2,1,1)
% %plotCellData(G, states{end}.h_free, ijk{1} == ii);
% plotGrid(G, ijk{1} == ii, 'FaceColor', 'none')
% plotPlume(G, Gt, states{end}.h_free,  ijk{1} == ii, 'EdgeColor','w','EdgeAlpha',.1)
% view([90 0])
% title({'Free Plume Height';['Slice through i=',num2str(ii)];'South -- North (facing West)'})
% ylabel('y-axis, meters'); zlabel('depth, meters')
% zlim([800, 850])
% set(gca,'DataAspect',[1 1 1/10])
% 
% subplot(2,1,2)
% %plotCellData(G, states{end}.h_free, ijk{2} == jj);
% plotGrid(G, ijk{2} == jj, 'FaceColor', 'none')
% plotPlume(G, Gt, states{end}.h_free,  ijk{2} == jj, 'EdgeColor','w','EdgeAlpha',.1)
% view([0 0])
% title({'Free Plume Height';['Slice through j=',num2str(jj)];'West -- East (facing North)'})
% xlabel('x-axis, meters'); zlabel('depth, meters')
% zlim([800, 850])
% set(gca,'DataAspect',[1 1 1/10])
% 


% % plot the max plume height:
% figure; set(gcf,'Position',[1 1 2200 500])
% subplot(2,1,1)
% %plotCellData(G, states{end}.h_max, ijk{1} == ii);
% plotGrid(G, ijk{1} == ii, 'FaceColor', 'none')
% plotPlume(G, Gt, states{end}.h_max,  ijk{1} == ii, 'EdgeColor','w','EdgeAlpha',.1)
% view([90 0])
% title({'Max Plume Height';['Slice through i=',num2str(ii)];'South -- North (facing West)'})
% ylabel('y-axis, meters'); zlabel('depth, meters')
% zlim([800, 850])
% set(gca,'DataAspect',[1 1 1/10])
% 
% subplot(2,1,2)
% %plotCellData(G, states{end}.h_max, ijk{2} == jj);
% plotGrid(G, ijk{2} == jj, 'FaceColor', 'none')
% plotPlume(G, Gt, states{end}.h_max,  ijk{2} == jj, 'EdgeColor','w','EdgeAlpha',.1)
% view([0 0])
% title({'Max Plume Height';['Slice through j=',num2str(jj)];'West -- East (facing North)'})
% xlabel('x-axis, meters'); zlabel('depth, meters')
% zlim([800, 850])
% set(gca,'DataAspect',[1 1 1/10])


% Colors for the CO2 trapping categories are found in: getInventoryColors()
%            1 - dissolved
%            2 - residual (traps)
%            3 - residual
%            4 - residual (plume)
%            5 - movable (traps)
%            6 - movable (plume)
%            7 - leaked


% plot of both max plume and free plume height:
figure;
% to xlim and ylim:
xmin = min(Gt.nodes.coords(:,1));
ymin = min(Gt.nodes.coords(:,2));
xmax = max(Gt.nodes.coords(:,1));
ymax = max(Gt.nodes.coords(:,2));
zmax = max(Gt.cells.z+Gt.cells.H); % the deepest depth
zmin = min(Gt.cells.z); % the shallowest depth


% to set up dimensions of subplots according to scales of x,y,z-coordinates
width1 = xmax - xmin;
width2 = ymax - ymin;
height = zmax - zmin;
space = 0.15*width1;
sideSpace = 0.2*width1;

normPlotWidth = 1;
PlotWidth = sideSpace + width1 + space + width2 + sideSpace;
normWidth1 = width1/PlotWidth;
normWidth2 = width2/PlotWidth;
normSpace = space/PlotWidth;
normSideSpace = sideSpace/PlotWidth;

left1 = normSideSpace;
left2 = left1 + normWidth1 + normSpace;
bottom1 = normSideSpace;
bottom2 = normSideSpace;
normPlotHeight = 0.8;



set(gcf,'Position',[1 1 PlotWidth/8 PlotWidth/20]); % make position more generic for other formations, to view two slices such that they are roughly scaled

hsp1 = subplot('Position', [left1 bottom1 normWidth1 normPlotHeight]);
hsp2 = subplot('Position',[left2 bottom2 normWidth2 normPlotHeight]);


%figure; set(gcf,'Position',[1 1 2200 350])
subplot(hsp1)
plotGrid(G, ijk{2} == jj, 'FaceColor', 'none'); view([0 0])
% max depth first...
plotPlume(G, Gt, states{end}.h_max,  ijk{2} == jj, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
% then free plume depth...
plotPlume(G, Gt, states{end}.h_free,  ijk{2} == jj, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)

%title('West <--> East')
%legend(['Grid slice through j=',num2str(jj)],'Residual CO2','Free plume CO2')
%hl = legend(['West to East slice through ',num2str( Gt.cells.centroids(wellCellIndex,2) ),' m'],'Residual CO2','Free plume CO2', 'Location','NorthEast');
xlabel('West to East','FontSize',16); zlabel('depth, meters','FontSize',16)
%zlim([800, 850]); % max/min depths for Sleipner. Make more generic for other formation cases.
set(gca,'DataAspect',[1 1 1/50])
xlim([xmin xmax]);
zlim([zmin zmax]);
box
set(gca,'FontSize',14)

subplot(hsp2)
plotGrid(G, ijk{1} == ii, 'FaceColor', 'none'); view([90 0])
% max depth first...
plotPlume(G, Gt, states{end}.h_max,  ijk{1} == ii, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
% then free plume depth...
plotPlume(G, Gt, states{end}.h_free,  ijk{1} == ii, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)

%title('South <--> North')
%legend(['Grid slice through i=',num2str(ii)],'Residual CO2','Free plume CO2')
%hl = legend(['South to North slice through ',num2str( Gt.cells.centroids(wellCellIndex,1) ),' m'],'Residual CO2','Free plume CO2', 'Location','SouthEast');

ylabel('South to North','FontSize',16); %zlabel('depth, meters','FontSize',16)
%zlim([800, 850]); % max/min depths for Sleipner. Make more generic for other formation cases.
set(gca,'DataAspect',[1 1 1/50])
ylim([ymin ymax]);
zlim([zmin zmax]);
box
set(gca,'FontSize',14)

hl = legend('Grid','Residual CO2','Free plume CO2', 'Location','SouthEast');
set(hl,'FontSize',16)




hfig = gcf;


end

