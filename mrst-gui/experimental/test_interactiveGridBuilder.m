mrstModule add upr mrst-gui
% Run the builder
out = interactiveGridBuilder();
%% Load a pre-existing grid sketch
interactiveGridBuilder(out);
%% Convert to MRST grid using UPR module
G = convertBuilderToPEBI(out, 1000, 'wellRefinement', true, 'wellGridFactor', 0.5);
%% Plot and demonstrate
figure; plotGrid(G)
axis tight
hold on
if ~isempty(out.points)
    plot(out.points(:, 1), out.points(:, 2), 'ok', 'markerfacecolor', 'r')
end
plotLinePath(out.wells,'color','blue');
plotLinePath(out.faults,'color','red');

%% Use a background image for tracing 
out = interactiveGridBuilder('image', 'street1.jpg');
