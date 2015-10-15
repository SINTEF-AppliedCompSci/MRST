function [ hfig ] = makeSideProfilePlots_CO2heights( G, Gt, sliceCellIndex, states, fluid, varargin)

% sliceCellIndex is the cell index for which a east-west slice and a
% south-north slice will intersect.

% if states and fluid are empty, they are not plotted.

opt.figname = '';

opt.SleipnerBounded = false;
% if true, subplots will be bounded to set region.

opt.legendWithFreeCO2Only = false;
% if true, what to include in legend will be determined by comparing free
% to residual heights.

opt = merge_options(opt, varargin{:});


% Cross-sectional slices through a point:
% first get index of point:
[ii,jj] = ind2sub(Gt.cartDims, sliceCellIndex);
disp(['Well cell index of ',num2str(sliceCellIndex),' corresponds to I=',num2str(ii),', J=',num2str(jj)])

ijk = gridLogicalIndices(G);

ii = ijk{1}(sliceCellIndex);
jj = ijk{2}(sliceCellIndex);
kk = ijk{3}(sliceCellIndex);
disp(['Cell index of ',num2str(sliceCellIndex),' corresponds to I=',num2str(ii),', J=',num2str(jj)])




% get the free plume height, and the max plume height reached in past
if ~isempty(states)
    states = addCO2HeightData(states, Gt, fluid);
end




% Colors for the CO2 trapping categories are found in: getInventoryColors()
%            1 - dissolved
%            2 - residual (traps)
%            3 - residual
%            4 - residual (plume)
%            5 - movable (traps)
%            6 - movable (plume)
%            7 - leaked


% plot of both max plume and free plume height:
figure('name',opt.figname); set(gcf,'Position',[1 1 1500 350]);

% Set up subplot sizes:

if opt.SleipnerBounded
    xmin = 436914;
    xmax = 440114;
    ymin = 6469150;
    ymax = 6475050;
    zmin = 790; %802.0627;
    zmax = 845; %850.6440;
else
    xmin = min(Gt.nodes.coords(:,1));
    ymin = min(Gt.nodes.coords(:,2));
    xmax = max(Gt.nodes.coords(:,1));
    ymax = max(Gt.nodes.coords(:,2));
    zmax = max(Gt.cells.z+Gt.cells.H);       % the deepest depth
    zmin = min(Gt.cells.z);                  % the shallowest depth
end


width1 = xmax - xmin;
width2 = ymax - ymin;
scaleFactor = width1/width2;
normWidth1 = 0.3;
normWidth2 = normWidth1/scaleFactor;

left1 = 0.07;
bottom = 0.17;
gap = 0.06;
normHeight = 0.75;
hsp1 = subplot('Position', [left1                bottom normWidth1 normHeight]);
hsp2 = subplot('Position', [left1+normWidth1+gap bottom normWidth2 normHeight]);


%figure; set(gcf,'Position',[1 1 2200 350])
subplot(hsp1)
plotGrid(G, ijk{2} == jj, 'FaceColor', 'none'); view([0 0])

if ~isempty(states)
% max depth first...
plotPlume(G, Gt, states{end}.h_max,  ijk{2} == jj, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
% then free plume depth...
plotPlume(G, Gt, states{end}.h_free,  ijk{2} == jj, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)
end

%title('West <--> East')
%legend(['Grid slice through j=',num2str(jj)],'Residual CO2','Free plume CO2')
%hl = legend(['West to East slice through ',num2str( Gt.cells.centroids(wellCellIndex,2) ),' m'],'Residual CO2','Free plume CO2', 'Location','NorthEast');
xlabel('West to East','FontSize',16);
zlabel('depth (meters)','FontSize',16)
%set(gca,'DataAspect',[1 1 1/50])
xlim([xmin xmax]);
zlim([zmin zmax]);
box
set(gca,'FontSize',14)



subplot(hsp2)
plotGrid(G, ijk{1} == ii, 'FaceColor', 'none'); view([90 0])

if ~isempty(states)
% max depth first...
hp1 = plotPlume(G, Gt, states{end}.h_max,  ijk{1} == ii, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1);
% then free plume depth...
hp2 = plotPlume(G, Gt, states{end}.h_free,  ijk{1} == ii, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1);
end

%title('South <--> North')
%legend(['Grid slice through i=',num2str(ii)],'Residual CO2','Free plume CO2')
%hl = legend(['South to North slice through ',num2str( Gt.cells.centroids(wellCellIndex,1) ),' m'],'Residual CO2','Free plume CO2', 'Location','SouthEast');

ylabel('South to North','FontSize',16);
%zlabel('depth, meters','FontSize',16)
%set(gca,'DataAspect',[1 1 1/50])
ylim([ymin ymax]);
zlim([zmin zmax]);
box
set(gca,'FontSize',14)

if ~isempty(states)
if opt.legendWithFreeCO2Only
    hl = legend([hp2],{'Free/Structural plume CO2'}, 'Location','SouthEast');
else
    if states{end}.h_max == states{end}.h_free
        hl = legend([hp2],{'Free/Structural plume CO2'}, 'Location','SouthEast');
    else
        hl = legend([hp1 hp2],{'Residual CO2','Free/Structural plume CO2'}, 'Location','SouthEast');
    end
end
set(hl,'FontSize',16)
end




hfig = gcf;


end

