%% Example 1: Show aquifer model
% In this example we consider a 1D antiform aquifer with a top surface
% given by the following expression
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
% We make a detailed plot of the model

aquifer = makeAquiferModel();

Gt = aquifer.Gt;
x  = Gt.cells.centroids(:, 1) / 1e3;
z  = Gt.cells.z;
H  = Gt.cells.H;

% Main plot of the antiform aquifer
figure;
hold on
patch(x([1:end, end:-1:1]), [z; z(end:-1:1) + H(end:-1:1)], 'y');
set(gca, 'YDir', 'reverse');
xlabel('Lateral extent [km]', 'FontSize', 12),
ylabel('Depth [m]', 'FontSize', 12);
h1 = gca;

% Inlet: zoom for x in [5,6]
h2 = axes('Position', get(h1, 'Position') * 0.4 + [0.5, 0.13, 0, 0]);
i = (x>4.9) & (x<6.1); [x2, z2, H2] = deal(x(i), z(i), H(i));
patch(x2([1:end, end:-1:1]), [z2; z2(end:-1:1) + H2(end:-1:1)], 'y');
axis([5 6 2120 2205])
set(h2, 'YDir', 'reverse', 'XTick', [5 5.5 6], 'YTick', [2150 2175 2200])
plot(h1, [5 6 6 5 5], [2120 2120 2205 2205 2120],'k')

% Inlet: zoom for x in [25,26]
h3 = axes('Position', get(h1, 'Position') * 0.35 + [0.11, 0.56, 0, 0]);
i = (x>24.9) & (x<26.1); [x3, z3, H3] = deal(x(i), z(i), H(i));
patch(x3([1:end, end:-1:1]), [z3; z3(end:-1:1) + H3(end:-1:1)], 'y');
axis([25 26 1720 1733])
set(h3, 'YDir', 'reverse', 'XTick', [25 26], 'YAxisLocation', 'right')
plot(h1, [25 26 26 25 25]', [1720 1720 1733 1733 1720], 'k')
set([h1 h2 h3], 'FontSize',12);

% print -depsc2 grid_1D_example.eps;