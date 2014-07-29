function normS = twophaseTransport(doPlot, grav)
% twophaseTransport -- Examples of two-phase transport on small Cartesian
%                      grid with homogeneous, isotropic permeability.
%                      Compares explicit and implicit transport solvers.
%
% SYNOPSIS:
%           twophaseTransport()
%           twophaseTransport(doPlot)
%           twophaseTransport(doPlot, gravity)
%   normS = twophaseTransport(...)
%
% DESCRIPTION:
%   Function twophaseTransport runs the following examples:
%     1) Water flooding in homogeneous domain by Dirichlet boundary
%        conditions on left and right boundaries. The reservoir is
%        initially saturated with oil.
%
%     2) Water flooding in a quarter five-spot well configuration on a
%        10-by-10 Cartesian grid with homogeneous, isotropic permeability
%        and porosity.  The reservoir is initially saturated with oil.
%
%     3) A gravity segregation example on the same 10-by-10 Cartesian grid
%        in which all cells initially have a water (and oil) saturation of
%        0.5.
%
%     4) A vertically separated gravity driven flow example on the same
%        10-by-10 Cartesian grid in which all cells to the left of the
%        separation line have an initial water saturation of 0.8, while all
%        cells to the right of the separation line have an initial water
%        saturation of 0.
%
% PARAMETERS:
%   doPlot  - Whether or not to plot results while the simulation is
%             running.  Logical.  Default value: doPlot = true.
%
%   gravity - Whether or not to consider gravity effects.
%             Examples 2) and 3) are only run when gravity is present.
%             Logical.  Default value: gravity = true.
%
% RETURNS:
%   normS - Euclidian norm of saturation differences, NORM(s_e - s_i), with
%           's_e' being the saturations computed with the explicit
%           transport solver and 's_i' being the saturations computed with
%           the implicit transport solver.  One scalar value for each test.


if nargin == 0,
   doPlot = true;
   grav   = true;
end

verbose = true;

if grav, gravity reset on, end

%--------------------------------------------------------------------------
%- Initialize system ------------------------------------------------------
%
dims      = [10, 1, 10];
G         = cartGrid(dims);
G         = computeGeometry(G, 'Verbose', verbose);
rock.poro = repmat(0.3, [G.cells.num, 1]);
rock.perm = repmat(1*darcy, [G.cells.num, 1]);

if false,
   fn    = 'data/immisciblewateroilquadratic.txt';
   fluid = initEclipseFluid(convertDeckUnits(readEclipseDeck(fn)), 'verbose', true);
   fluid = rmfield(fluid,'pc');
else
   fluid = initSimpleCompFluid('mu', [1, 1], 'rho_ref', [700, 1000], ...
                               'press_ref', 100*barsa, 'n', [2, 2], ...
                               'c', 1*[9e-5, 9e-6]/barsa);
   fluid.names = {'Oil', 'Water'};
end

gravity(grav);

T = computeTrans(G, rock);
trans = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ T, [G.faces.num, 1]);

%--------------------------------------------------------------------------
%- Run 3 cases if grav = true, 1 case otherwise ---------------------------
%
if grav, tests = 1:1; else tests = 1; end
normS = zeros([length(tests), 1]);


fun = {(@(varargin) test_grav(G, fluid, 100*barsa))};

for i = tests,
   t = 0;
   [x, drive, tf, dt] = fun{i}(initResSol(G, 0.0, [0,1]), verbose);

   %-----------------------------------------------------------------------
   %- Compute initial solution of pressure equation -----------------------
   %
   x_tpfa  = x;
   x_impes = x;

   if doPlot,
      % Plot inital saturation
      figure(i);

      subplot(2,1,1),
         h = plotCellData(G, x.s(:,1));
         view([0, 0]);
         axis tight; axis equal;
         caxis([0, 1]), colorbar

      subplot(2,1,2),
         h2 = plotCellData(G, x.s(:,1));
         view([0, 0]);
         axis tight; axis equal;
         caxis([0, 1]), colorbar
   end

   %-----------------------------------------------------------------------
   %- Transport loop - solve saturation equation and update pressures -----
   %
   pv = poreVolume(G, rock);
   while t < tf,
      x_impes = impesTPFA        (x_impes, G, trans, fluid, dt, pv, drive{:});
      x_tpfa  = compTPFA         (x_tpfa, G, rock, T, fluid, dt, drive{:});
      x_tpfa  = implicitTransport(x_tpfa, G, dt, rock, fluid, drive{:});

      if doPlot,
         % Plot saturation
         delete(h);
         subplot(2,1,1)
            h = plotCellData(G, x_impes.s(:,1));
            title(['Test ', int2str(i), ...
                   '. ', fluid.names{1}, ' saturation - IMPES solver'])
            view([0, 0]), axis tight equal

         delete(h2);
         subplot(2,1,2),
            h2 = plotCellData(G, x_tpfa.s(:,1));
            title([fluid.names{1}, ...
                   ' saturation - SS (cTPFA/iTransport) solver'])
            view([0, 0]), axis tight equal

         pause(0.2)
      end

      t = t + dt;
   end

   normS(i) = norm(x_impes.s(:,1) - x_tpfa.s(:,1), inf);
end

%--------------------------------------------------------------------------

function [x, drive, tf, dt] = test_grav(G, fluid, p0, varargin)
%- Gravity column ----------------------------------------
% Test of segregation flow caused by gravity.
% Initial water saturation in all cells:
%           ________
%          |        |
%          | s = 0.5|
%          |________|
%

x     = initResSolComp(G, [], fluid, p0, [0.3, 0.7]);
drive = {};

tf    = 30000 * day;
dt    = 100 * day;
