%% Small test for the simple equilibrium routine
G = cartGrid([10, 1, 10], [1, 1, 1]);
G = computeGeometry(G);

% First contact at .37, second at .8 for a total of three phases present
contacts = [.37, .8];

% Compute equilibrium saturations
s = simpleEquilibrium(G, contacts);

% Plot the different saturations, as well as the different contacts
figure(1); clf
for i = 1:size(s, 2)
    subplot(1, size(s, 2), i);
    plotCellData(G, s(:, i))
    caxis([0, 1])
    view(0, 0)
    hold on
    for j = 1:numel(contacts)
        plot3([0, 1], [-0.01, -0.01], contacts([j, j]), 'r', 'linewidth', 3)
    end
end

%% The equilibriation does not have to happen along the z-axis
% By specifying a vector, we can define the saturations to be aligned with
% the force of gravity in any direction.
vec = [1, 0, 1];
s = simpleEquilibrium(G, contacts, vec);

figure(1); clf
for i = 1:size(s, 2)
    subplot(1, size(s, 2), i);
    plotCellData(G, s(:, i))
    caxis([0, 1])
    view(0, 0)
    hold on
end
