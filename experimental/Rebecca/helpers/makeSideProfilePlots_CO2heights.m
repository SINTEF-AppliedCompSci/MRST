function [ hfig ] = makeSideProfilePlots_CO2heights( G, Gt, wellCellIndex, states, fluid)

% Cross-sectional slices through injection point:
% first get index of injection well:
[ii,jj] = ind2sub(Gt.cartDims, wellCellIndex);
disp(['Well cell index is I=',num2str(ii),', J=',num2str(jj)])

ijk = gridLogicalIndices(G);


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
figure; set(gcf,'Position',[1 1 2200 350])

subplot(2,1,1)
plotGrid(G, ijk{1} == ii, 'FaceColor', 'none'); view([90 0])
% max depth first...
plotPlume(G, Gt, states{end}.h_max,  ijk{1} == ii, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
% then free plume depth...
plotPlume(G, Gt, states{end}.h_free,  ijk{1} == ii, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)

title('South <--> North')
legend(['Grid slice through i=',num2str(ii)],'Residual CO2','Free plume CO2')
ylabel('y-axis, meters'); zlabel('depth, meters')
zlim([800, 850]); set(gca,'DataAspect',[1 1 1/10])


subplot(2,1,2)
plotGrid(G, ijk{2} == jj, 'FaceColor', 'none'); view([0 0])
% max depth first...
plotPlume(G, Gt, states{end}.h_max,  ijk{2} == jj, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
% then free plume depth...
plotPlume(G, Gt, states{end}.h_free,  ijk{2} == jj, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)

title('West <--> East')
legend(['Grid slice through j=',num2str(jj)],'Residual CO2','Free plume CO2')
xlabel('x-axis, meters'); zlabel('depth, meters')
zlim([800, 850]); set(gca,'DataAspect',[1 1 1/10])


hfig = gcf;


end

