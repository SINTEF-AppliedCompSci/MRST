%% VE-incomp: C-Acceleration on Large Model
% In this example we consider a synthetic sloping aquifier that has a
% significant variation in the top-surface morphology.
%
% We demonstrate the use of C/C++-accelerated MATLAB, using the function
% mtransportVE to replace explicitTransportVE. Using mtransportVE requires
% that you have built the solver in the VEmex directory.
%
% In addition, the routine makeIGEMSmodel that sets up the data model may
% use the following C-accelerated MATLAB routines
%
% * processgrid (replaces processGRDECL)
% * mcomputegeometry (replaces computeGeometry)

mrstModule add co2lab
%% Write header
clc;
disp('================================================================');
disp('   Vertical averaging applied to a large model');
disp('   using C++ accelleration in the transport solver');
disp('================================================================');
disp(' ');

%% Construct stratigraphic, petrophysical, and VE models
% The 3D model consists of a grid (G) and petrophysical parameters (rock).
% The VE model consists of a top-surface grid (Gt), petrophysical data
% (rock2D), and indices to the boundarcy cells where we will supply
% pressure boundary conditions. Called with a true flag, the routine will
% use C-accelerated MATLAB routines to process the data input and compute
% geometry. Once the models are created, they are stored in a data file for
% faster access at a later time.
[G, Gt, rock, rock2D, bcIxVE] = makeSlopingAquiferBig(true);

%% Set time and fluid parameters
% Inject CO2 for 150 years and study subsequent migration until 750 years
% after injection started. The fluid data are chosen so that they are
% resonable at p = 300 bar
gravity on
T          = 750*year();
stopInject = 150*year();
dT         = 2*year();
dTplot     = 2*dT;
fluidVE    = initVEFluidHForm(Gt, 'mu' , [0.056641 0.30860] .* centi*poise, ...
                         'rho', [686.54 975.86] .* kilogram/meter^3, ...
                        'sr', 0.2, 'sw', 0.1, 'kwm', [0.2142 0.85]);

%% Set well and boundary conditions
% We use one well placed down the flank of the model, perforated in the
% bottom layer. Injection rate is 2.8e4 m^3/day of supercritical CO2.
% Hydrostatic boundary conditions are specified on all outer boundaries.
disp(' -> Setting well and boundary conditions');

% Set well in 3D model
wellIx = [G.cartDims(1:2)/5, G.cartDims([3 3])];
rate   = 2.8e4*meter^3/day;
W      = verticalWell([], G, rock, wellIx(1), wellIx(2), ...
                      wellIx(3):wellIx(4), 'Type', 'rate', 'Val', rate, ...
                      'Radius', 0.1, 'comp_i', [1,0], 'name', 'I', ...
                      'InnerProduct', 'ip_simple');

% Well in 2D model
WVE = convertwellsVE(W, G, Gt, rock2D);

% BC in 2D model
bcVE   = addBC([], bcIxVE, 'pressure', ...
            Gt.faces.z(bcIxVE)*fluidVE.rho(2)*norm(gravity));
bcVE   = rmfield(bcVE,'sat');
bcVE.h = zeros(size(bcVE.face));


%% Prepare simulations
% Compute inner products and instantiate solution structure
disp(' -> Initialising solvers');
SVE = computeMimeticIPVE(Gt, rock2D, 'Innerproduct','ip_simple');
preComp = initTransportVE(Gt, rock2D);
sol = initResSolVE(Gt, 0, 0);
sol.wellSol = initWellSol(W, 300*barsa());
sol.s = height2Sat(sol, Gt, fluidVE);

% Select transport solver
% Use C++ acceleration if it exists - NB: requires the VEmex module
% Notice that the two solvers determine the time steps differently and
% may therefore give slightly different answers.
try
   mtransportVE();
   cpp_accel = true;
catch me
   disp('mex-file for C++ acceleration not found');
   disp(['See ', fullfile(mrstPath('co2lab'),'ve','VEmex','README'), ...
      ' for building instructions']);
   disp('Using matlab ve-transport');
   cpp_accel = false;
end
mrstModule add mimetic

% Find trapping structure in grid. Used for calculation of trapped volumes
ts=findTrappingStructure(Gt);

%% Prepare plotting
% We will make a composite plot that consists of several parts: a 3D plot
% of the plume, a pie chart of trapped versus free volume, a plane view of
% the plume from above, and two cross-sections in the x/y directions
% through the well
opts = {'slice', wellIx, 'Saxis', [0 1-fluidVE.res_water], ...
   'maxH', 200, 'Wadd', 1000};
plotPanelVE(G, Gt, W, sol, 0.0, zeros(1,6), opts{:});

%% Main loop
% Run the simulation using a sequential splitting with pressure and
% transport computed in separate steps. The transport solver is formulated
% with the height of the CO2 plume as the primary unknown and the relative
% height (or saturation) must therefore be reconstructed.
t = 0;
totVol = 0.0;
fprintf(1,'\nSimulating %d years of injection', convertTo(stopInject,year));
fprintf(1,' and %d years of migration\n', convertTo(T-stopInject,year));
fprintf(1,'Time: %4d years', convertTo(t,year));
tic
while t<T
   % Advance solution: compute pressure and then transport
   sol = solveIncompFlowVE(sol, Gt, SVE, rock, fluidVE, ...
      'bc', bcVE, 'wells', WVE);

   if cpp_accel
      [sol.h, sol.h_max] = mtransportVE(sol, Gt, dT, rock, ...
                                          fluidVE, 'bc', bcVE, 'wells', WVE, ...
                                          'gravity', norm(gravity));
   else
      sol = explicitTransportVE(sol, Gt, dT, rock, fluidVE, ...
                               'bc', bcVE, 'wells', WVE,    ...
                               'preComp', preComp,          ...
                               'intVert', false);
   end
   
   % Reconstruct 'saturation' defined as s=h/H, where h is the height of
   % the CO2 plume and H is the total height of the formation
   sol.s = height2Sat(sol, Gt, fluidVE);
   assert( max(sol.s(:,1))<1+eps && min(sol.s(:,1))>-eps );
   t = t + dT;

   % Compute total injected, trapped and free volumes of CO2
   if ~isempty(WVE)
      totVol = totVol + WVE.val*dT;
   end
   vol = volumesVE(Gt, sol, rock2D, fluidVE, ts);

   % Check if we are to stop injecting. If so, increase the time step.
   if t>= stopInject
      WVE  = []; dT = 10*year(); dTplot = dT;
   end

   % Plotting
   fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d years', convertTo(t,year));
   if mod(t,dTplot)~= 0 && t<T
      continue
   else
      plotPanelVE(G, Gt, W, sol, t, [vol totVol], opts{:});
      drawnow
   end
end
fprintf(1,'\n\n');

% delete C++ simulator
if cpp_accel, mtransportVE(); end
etime = toc;
disp(['Elapsed simulation time: ', num2str(etime), ' seconds.']);

displayEndOfDemoMessage(mfilename)

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
