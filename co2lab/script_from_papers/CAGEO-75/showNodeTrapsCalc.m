G = cartGrid([20 20 1]); 
[x, x] = meshgrid(linspace(0, 1, 21));                                     %#ok<ASGLU>
z = max(peaks(21) + 2 * x, 0); 
G.nodes.coords(1:21 * 21, 3) = -z(:); 
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 10; 
G = computeGeometry(G); 
Gt = topSurfaceGrid(G); 

%% Start analysis
% Set breakpoint on line 93 of computeNodeTraps and call this function
%dbstop in computeNodeTraps at 93;
res = computeNodeTraps(Gt, [], []); 

%% Visualization inside computeNodeTraps
% Assuming that you have come to line 63 of computeNodeTraps, run the
% following lines to visualize various structures of the trap analysis
figure('Position', [900 420 860 460]); 
plotGrid(Gt, 'FaceColor', 'none'); view( - 55, 73); 

x = Gt.nodes.coords(:, 1); y = Gt.nodes.coords(:, 2); z = Gt.nodes.z; 
col = colorcube(max(max(regions) + 1, 8)); 
hold on, 
hr = []; 
for n = 0:max(regions)
   i = regions == n; 
   h = plot3(x(i), y(i), z(i) -.05, '.', 'MarkerSize', 16, 'Color', col(n + 1, :)); 
   hr = [hr; h];                                                           %#ok<AGROW>
end
i = spill_edges(:, 6:7)'; 
hse = plot3(x(i), y(i), z(i), 'Color', [.5 .5 .5], 'LineWidth', 2.5); 

faces = vertcat(spoint_edges{:}); 
eIX = Gt.faces.nodePos; 
nn = double(diff([Gt.faces.nodePos(faces),...
                  Gt.faces.nodePos(faces + 1)], [], 2)); 
fn = double(Gt.faces.nodes(mcolon(eIX(faces), eIX(faces + 1) - 1), 1)); 
i = reshape(fn, 2, []); 
hspe = plot3(x(i), y(i), z(i) -.05, 'k', 'LineWidth', 5); 


hold off
axis off; 
h = legend([hr(2:end); hse(1); hspe(1)], 'spill region 1', 'spill region 2',...
'spill region 3', 'spill edges', 'spill point edges', 'Location', 'SouthEastOutside');
dbclear computeNodeTraps; dbcont;