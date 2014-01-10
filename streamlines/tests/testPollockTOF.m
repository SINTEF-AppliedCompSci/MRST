require mimetic

G = cartGrid([60,60]);
G = computeGeometry(G);

% NB non-unit porosity not working in pollock...

rock.perm = repmat(10*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(1,            [G.cells.num, 1]);

%load Udata
%rock.perm = reshape(KU(1,1:60,1:60,37), [], 1);
rock.perm = convertFrom(logNormLayers([G.cartDims, 1], 1), milli*darcy);

fluid = initSimpleFluid('mu',  [1, 1]*centi*poise', ...
                        'rho', [1000, 1000]*kilogram/meter^3, ...
                        'n',   [2, 2]);

IP = computeMimeticIP(G, rock);
src = addSource([], [1, G.cells.num], [1.1, -1.1], 'sat', [1.0 0.0;1.0,0.0]);

x = initResSol(G, 0, 0);
x = solveIncompFlow(x, G, IP, fluid, 'src', src);

numStreamlines = 500;


subplot(1,3,1);
title('Streamlines');
h = plotGrid(G, 'facea', 0.3, 'edgea',0.1);
axis equal tight;
pos = [ones(numStreamlines,1), rand(numStreamlines,2)];
timer=tic;
[xyz,t,c] = pollock(G, x, pos, 'substeps', 1, 'maxsteps',10000);
thetime=toc(timer);
dispif(true, 'Traced %d steps in %d streamlines in %g second\n', ...
   sum(cellfun(@numel, t)), numStreamlines, thetime);
streamline(xyz);

subplot(1,3,2);
title('Streamline time of flight')

% add time of flight along streamlines
ct = cellfun(@cumsum, t, 'uniformoutput', false);

% rearrange to plain double arrays
t  = vertcat(t{:});
i= (t ~= inf);
t  = t(i);
ct = vertcat(ct{:}); ct = ct(i);
c  = vertcat(c{:});  c  = c(i);

% Compute weithed average of time of flight in each cell that has been
% crossed by one or more streamlines,
%
%         ---               /   ---
% tof_i = \   t_ij*ct_ij   /    \     t_ij
%         /__ j           /     /__ j
%
% where t_ij is the time-of-flight diff through cell i of streamline j

tof = accumarray(c, ct.*t, [G.cells.num, 1], [], inf)./...
   accumarray(c, t,     [G.cells.num, 1], [], 1);


plotCellData(G, tof);
axis equal tight;

subplot(1,3,3);
title('Finite-volume time of flight')
fdTOF = computeTimeOfFlight(x, G, rock, 'src', src);
plotCellData(G, fdTOF);
axis equal tight;

% Use fdTOF for caxis scalling as this does not change randomly.
subplot(1,3,2), caxis([0,fdTOF(src.cell(end))])
subplot(1,3,3), caxis([0,fdTOF(src.cell(end))])

