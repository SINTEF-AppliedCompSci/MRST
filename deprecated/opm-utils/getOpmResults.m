function opmres = getOpmResults(output_dir, steps, deck)
%Retrieve OPM simulation results from disk
%
% SYNOPSIS:
%   opmres = getOpmResults(output_dir, steps, deck)
%
% PARAMETERS:
%   output_dir - Master directory containing all simulation results.
%                Individual result vectors are assumed to be contained in
%                named directories directly below 'output_dir'.
%
%                The following result vectors are currently supported
%
%                    concentration, flux, pressure, reorder_it,
%                    residualcounts, saturation, surfvolume
%
%   steps      - List of time steps for which to load simulation results.
%                Must be numeric and valid indices between zero (0) and the
%                total number of simulation steps less one (inclusive).
%
%   deck       - ECLIPSE/FrontSim input deck.  Needed to identify the
%                active fluid phases in the simulation run.
%
% RETURNS:
%   opmres - Simulation results.  Array of structures containing the result
%            vectors detected during directory scanning of 'output_dir',
%            restricted to those steps identified by 'steps'.
%
%            Specifically, opmres(i) contains the result vectors of
%            simulation step 'steps(i)'.
%
% EXAMPLE:
%   % Load all simulation results in 'test-output' from a three-phase
%   % simulation of a FrontSim simulation model identified by 'deck', and
%   % plot the fifth pressure field as well as the saturation of the second
%   % component of the twentieth time step:
%
%   G     = initEclipseGrid(deck.GRID);
%   nstep = numberOfSteps('test-output');
%   res   = getOpmResults('test-output', 0 : nstep-1, deck);
%
%   figure
%   plotCellData(G, res(5).pressure, ...
%                'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.85)
%   view(3)
%
%   figure
%   plotCellData(G, res(20).saturation(:,2), ...
%                'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.85)
%   view(3)
%
% SEE ALSO:
%   numberOfSteps.

%{
#COPYRIGHT#
%}

% $Date: 2012-09-14 20:40:02 +0200 (Fri, 14 Sep 2012) $
% $Revision: 9694 $

   np  = sum(isfield(deck.RUNSPEC, { 'WATER', 'OIL', 'GAS' }));
   res = loadResults(output_dir, np, 'steps', steps + 1);

   fn  = reshape(fieldnames (res), 1, []);
   res = reshape(struct2cell(res), 1, []);
   v   = [ fn ; res ];

   opmres = struct(v{:});
end
