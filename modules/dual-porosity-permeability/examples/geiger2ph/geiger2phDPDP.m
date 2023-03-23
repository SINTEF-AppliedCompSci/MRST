%% Waterflooding of the upscaled geometry due to Geiger et al., 2011.
% The geometry due to Geiger et al., 2011 contains a set of orthogonal 
% fractures in a 1x1 m box. A Cartesian DPDP grid of 11x11 blocks 
% is overlaid on the domain and the upscaled matrix and fractures porosities 
% and permeabilities, as well as the shape factors for each DPDP block are 
% calculated as described in Andrianov and Nick, 2021. The DPDP results are
% compared with the reference fine-scale DFM solution from N. Andrianov
% 2021 (submitted).
%
% References
% S. Geiger, M. Dentz, and I. Neuweiler. A novel multi-rate dual-porosity model
%   for improved simulation of fractured and multi-porosity reservoirs. 
%   Society of Petroleum Engineers, 2011. doi:10.2118/148130-MS.
% N. Andrianov and H. M. Nick. Machine learning of dual porosity model closures 
%   from discrete fracture simulations. Advances in Water Resources, 147:103810, 2021.
%   doi:10.1016/j.advwatres.2020.103810.
% N. Andrianov, Upscaling of Two-Phase Discrete Fracture Simulations Using 
%   a Convolutional Neural Network, submitted (2021).


mrstModule add ad-core ad-props ad-blackoil dual-porosity

%% Set up the problem

% Uniform matrix porosity and permeability
phim = 0.3;
Km = 9.8692e-14;

% The parameters for matrix and fracture Corey relative permeabilities
[krw, kro, nw, no, Swr, Sor] = deal(0.393736, 0.201947, 2.05375, 2.01668, 0.1034, 0.355703);

% Densities and viscosities for water and oil
[rhow, rhoo, muw, muo] = deal(1050, 750, 1.09e-3, 9.2e-4);

% Initial pressure (the same for matrix and fractures)
p0 = 1e5;

% Water mass injection rate in kg/(m*s)
qf = -2.100e-04;

% Define the upscaled geometry
[Nx, Ny] = deal(11);
[xmin, xmax] = deal(0, 1);
[ymin, ymax] = deal(0, 1);

% Create the upscaled geometry
dx = (xmax-xmin) / Nx;
dy = (ymax-ymin) / Ny;
G = cartGrid([Nx Ny],[xmax-xmin ymax-ymin]);
G.nodes.coords(:,1) = G.nodes.coords(:,1) + xmin;
G.nodes.coords(:,2) = G.nodes.coords(:,2) + ymin;
G = computeGeometry(G); 

% Reading the upscaled fracture permeabilities, porosities, and transfer rate parameters
dataset = 'geiger_upscaled_props.csv';
T = dlmread(dataset);
    
% The block indices in Cartesian DPDP geometry and the linear index of
% the blocks, containing fractures
I = T(:, 2);
J = T(:, 3);
ind = sub2ind([Nx Ny], I, J);

% Set up the uniform porosities, diagonal permeabilities, and shape factors
rock.phim = ones(Nx * Ny, 1) * phim;   % Initialize with the constant matrix porosity
rock.phif = zeros(Nx * Ny, 1);        % Initialize with zeros (the default value when no fractures are present)  
rock.Km = ones(Nx * Ny, 3) * Km;      % Initialize with the constant matrix permeability
rock.Km(:, 2) = 0;    
rock.Kf = zeros(Nx * Ny, 3);          % Initialize with zeros 
rock.smin = zeros(Nx * Ny, 1);        % Initialize with zeros (the default value when no fractures are present) 

% Read the non-uniform porosities and permeabilities
rock.phim(ind) = T(:,4); 
Km_xx = T(:, 5);
Km_xy = T(:, 6);
Km_yx = T(:, 7);
Km_yy = T(:, 8);
rock.phif(ind) = T(:,9);
Kf_xx = T(:, 10);
Kf_xy = T(:, 11);
Kf_yx = T(:, 12);
Kf_yy = T(:, 13);  

% Read the shape factors
rock.smin(ind) = T(:, 14);

% Dimensionalize the shape factors
vb = G.cells.volumes;
rock.smin = rock.smin ./ vb;

% Set zero fracture porosities to a small value
izp = find(rock.phif == 0);
if ~isempty(izp)
    poro0 = mean(rock.phif) * 1e-8;
    disp(['Avoid MRST complaints by setting zero fracture porosities to ' num2str(poro0)])
    rock.phif(izp) = poro0;
end
 
% Remove the rows with negative Kxx and Kyy
ink = find(Kf_xx < 0);
mink = 1e-3 * min(min(abs([Kf_xx Kf_yy]))); 
if ~isempty(ink)               
    disp(['Replacing ' num2str(length(ink)) ' rows of negative Kf_xx with ' num2str(mink/darcy) ' D'])
    Kf_xx(ink) = mink;
end        
ink = find(Kf_yy < 0);
if ~isempty(ink)
    disp(['Replacing ' num2str(length(ink)) ' rows of negative Kf_yy with ' num2str(mink/darcy) ' D'])
    Kf_yy(ink) = mink;
end       

% Remove the rows with negative shape factor
insf = find(rock.smin < 0);
if ~isempty(insf)
    disp(['Replacing ' num2str(length(insf)) ' rows of negative parameters of the transfer function with zeros..'])
    rock.smin(insf) = 0;
end    

% Approximate the off-diagonal permeabilities using arithmetic average
rock.Km(ind, :) = [Km_xx 0.5*(Km_xy + Km_yx)  Km_yy];
rock.Kf(ind, :) = [Kf_xx 0.5*(Kf_xy + Kf_yx)  Kf_yy];

% Set zero fracture permeabilities to a small value, ensuring that this
% value is less than the minimal fracture permeability (e.g. for
% hanging fracture traces).
if ~isempty(izp)
    perm0 = 1e-16;
    mink = min(min(abs([Kf_xx Kf_xy Kf_yy])));
    perm0 = min(perm0, mink);
    disp(['Avoid convergence problems by setting zero fracture permeabilities to ' ...
        num2str(perm0) '(=' num2str(perm0/milli/darcy) ' mD)'])
    rock.Kf(izp, :) = perm0; 
end
 

%% Set up the DPDP model

% Create the rock structures for fractures and the matrix
rock_fracture = makeRock(G, rock.Kf, rock.phif);
rock_matrix = makeRock(G, rock.Km, rock.phim);


%% Plot the upscaled fracture permeabilities and shape factors together with the fracture geometry 
plot_rock = true;
if plot_rock
    fracsfile = 'geiger_fracs.csv';
    fid = fopen(fracsfile);
    if (fid == -1)
        disp(['Cannot open ' fracsfile]);
        return;
    end

    fclose(fid);

    T = readtable(fracsfile);

    FracIdx = T{:, 1};
    SegIdx = T{:, 2};

    SegXStart = T{:, 5};
    SegYStart = T{:, 6};
    SegZStart = T{:, 7};

    SegXMid = T{:, 8};
    SegYMid = T{:, 9};
    SegZMid = T{:, 10};

    SegXEnd = T{:, 11};
    SegYEnd = T{:, 12};
    SegZEnd = T{:, 13};

    Nseg = length(SegIdx);

    % Plot the upscaled properties
    for k = 1:3
        subplot(2, 2, k)
        plotCellData(G, rock_fracture.perm(:,k)/darcy,  'EdgeColor', 'none')
    end
    subplot(2, 2, 4)
    % Plot the dimensionless shape factor
    plotCellData(G, rock.smin.*vb,  'EdgeColor', 'none')

    % Plot the fracture geometry
    tit = {'Upscaled fracture K_{xx} (D)', ...
           'Upscaled fracture K_{xy} (D)', ...
           'Upscaled fracture K_{yy} (D)', ...
           'Shape factor'};
    for k = 1:4
        subplot(2, 2, k)
        hold on
        for i = 1:Nseg
           hh = line([SegXStart(i) SegXMid(i) SegXEnd(i)], ...
                     [SegYStart(i) SegYMid(i) SegYEnd(i)], ...
                     'LineWidth', 1, 'Color', 'r'); 
        end
        colormap parula
        axis equal
        axis([ xmin, xmax, ymin, ymax ]);
        title(tit{k})
        colorbar
    end

end


%% Read the reference fine-scale DFM solution
refsolution = 'geiger_reference.csv';
disp(['Reading ' refsolution]);
data = readtable(refsolution);

t = data{:, 1};     % Time (sec)
Qw_i = data{:, 4};  % Injected water volume (m3)
Qo_p = data{:, 6};  % Produced oil volume (m3)
Qw_p = data{:, 7};  % Produced water volume (m3)
Nt = length(t);

% Eventually limit the number of time steps in DPDP simulations 
Ntlim = 1e9;
if Nt > Ntlim
    disp(['Limiting the number of time steps in DPDP run to ' num2str(Ntlim)])
    dN = floor(Nt / (Ntlim - 1));
    ind = [1:dN:Nt]';
    t = t(ind);
    Qw_i = Qw_i(ind);
    Qo_p = Qo_p(ind);
    Qw_p = Qw_p(ind);
    Nt = length(t);
end

% Keep only unique time instants (identical time instants may occur due to
% rounding error while outputting a csv)
[~, ind] = unique(t);
if length(ind) ~= Nt
    disp(['Keep only ' num2str(length(ind)) ' unique time instants'])
    t = t(ind);
    Qw_i = Qw_i(ind);
    Qo_p = Qo_p(ind);
    Qw_p = Qw_p(ind);
end

% Reference solutions for cumulative produced oil and water volumes
timescale = day;
figure(2)
subplot(2, 1, 1)
plot(t/timescale, Qo_p, '-b')
hold on
xlabel('Time');
ylabel('Produced oil volume (m3)')
title('Produced oil volume')

subplot(2, 1, 2)
plot(t/timescale, Qw_p, '-b')
hold on
xlabel('Time');
ylabel('Produced water volume (m3)')
title('Produced water volume') 


%% Setup the fluid model
% Use the same two-phase fluid model for the matrix and the fractures
fluidMatrix = initSimpleADIFluid('phases','WO',    ... % Fluid phases: water and oil
                           'mu',  [muw, muo],      ... % Viscosities
                           'rho', [rhow,rhoo],     ... % Surface densities [kg/m^3]
                           'c',   [0, 0],          ... % Fluid compressibility[Cm, Cm],              
                           'pRef', 0,              ... % Reference pressure 
                           'cR',  0                ... % Rock compressibility
                           );
fluidFractures = fluidMatrix;                       

% Matrix relative permeabilities
fluidMatrix.krW  = coreyPhaseRelpermAD(nw, Swr, krw, Swr + Sor);
fluidMatrix.krO = coreyPhaseRelpermAD(no, Sor, kro, Swr + Sor);                                              

% Fractures relative permeabilities
fluidFractures.krW  = coreyPhaseRelpermAD(nw, Swr, krw, Swr + Sor);
fluidFractures.krO = coreyPhaseRelpermAD(no, Sor, kro, Swr + Sor);     

% Set the residual saturations which can be used in the transfer function                            
fluidMatrix.swr = Swr;
fluidMatrix.snr = Sor;
fluidFractures.swr = Swr;
fluidFractures.snr = Sor;
                       
% Initialize the 2 phase dual porosity - dual permeability model
model = TwoPhaseOilWaterModelDPDP(G, ...
                                {rock_fracture, rock_matrix}, ...  % Rock properties
                                {fluidFractures, fluidMatrix}, ... % Fluid properties  
                                'gravity', [0 0] ...
                                );
                           
% Calculate blocks' dimensions for using in EclipseTwoPhaseTransferFunction 
lx = (xmax-xmin) / Nx;
ly = (ymax-ymin) / Ny;
block_dimension = repmat([lx,ly,1],G.cells.num,1);

% Initialize the shape factors using the values from fine-scale simulations
model.transfer_model_object = myEclipseTwoPhaseTransferFunction('VariableShapeFactor', ...
   [rock.smin block_dimension]);

% Matrix and fracture capillary pressure
model.fluid_matrix.pcOW  = grauePcAD(Swr, Swr + Sor, model);
model.fluid.pcOW  = grauePcAD(Swr, Swr + Sor, model);


%% Initial conditions
clear state
state.pressure = ones(G.cells.num,1) * p0;  % Fracture pressure
state.pressure_matrix = ones(G.cells.num,1) * p0;       % Matrix pressure
state.s = repmat([Swr 1-Swr],G.cells.num,1);
state.sm = repmat([Swr 1-Swr],G.cells.num,1);

%% Boundary conditions

% Get the volumetric boundary rate
L = ymax-ymin;  
Qf = - qf * L / rhow;   % Water injection ("-" refers to DuMuX convention)

% Neumann BC at the left, Dirichlet BC at the right, no-flow conditions otherwise
% The initial pressure is interpreted as the right BC pressure 
bc = pside([], G, 'East', p0, 'sat', [1 0]);   
bc = fluxside(bc, G, 'West', Qf, 'sat', [1 0]);

% Set up the BC for fracture and matrix continua separately under the
% constraint that the total injection rate in the both continua is Qf.

% For each flux boundary cell, split the fraction of Qf to go to fracture and 
% matrix continua using the normal (x-) fracture permeability in this cell.
isSF = reshape(strcmpi(bc.type, 'flux'), [], 1);
N = G.faces.neighbors(bc.face(isSF),:);
BCcells = sum(N, 2);
Kbc = rock_fracture.perm(BCcells, 1);
Kfac = (Kbc - min(Kbc)) / (max(Kbc) - min(Kbc));
Qc = Qf / numel(Kbc);
bc.value(isSF) = Qc * Kfac;
bc.value_matrix = bc.value;
bc.value_matrix(isSF) = Qc * (1 - Kfac);

%% Schedule
% Use the same time instants as specified in the fine scale results
%dt = diff([0 t']);
dt = diff([0:max(t)/1000:max(t)]);
schedule = struct(...
    'control', struct('W',{[]}), ...    % Apparently need a fictious well
    'step',    struct('control', ones(numel(dt),1), 'val', dt));
schedule.control.bc = bc;

%% Run the model

% Format the simulation step output
tm = [0 ; reshape(cumsum(dt), [], 1)];
nSteps  = numel(tm) - 1;   
nDigits = floor(log10(nSteps)) + 1;
nChar = numel(formatTimeRange(tm(end), 2));
step_header = @(i) ...
  fprintf('Solving timestep %0*d/%0*d: %-*s -> %s\n', ...
          nDigits, i, nDigits, nSteps, nChar, ...
          formatTimeRange(tm(i + 0), 2), ...
          formatTimeRange(tm(i + 1), 2));

% Initialize the nonlinear solver with BackslashSolverAD as a linear solver
solver = NonLinearSolver('maxTimestepCuts', 20);

% Check if model and state are self-consistent and set up for current BC type
ctrl = schedule.control(schedule.step.control(1));
[forces, fstruct] = model.getDrivingForces(ctrl);
model = model.validateModel(fstruct);

% Get the indices of the inflow boundary
bwest = boundaryFaceIndices(model.G, 'West');

% The initial WIP and OIP
wipf0 = sum(model.operators.pv .* state.s(:, 1));
oipf0 = sum(model.operators.pv .* state.s(:, 2));
wipm0 = sum(model.operators_matrix.pv .* state.sm(:, 1));
oipm0 = sum(model.operators_matrix.pv .* state.sm(:, 2));

% The time loop
prevControl = nan;
for i = 1:nSteps
    step_header(i);
    state0 = state;
    
    currControl = schedule.step.control(i);
    if prevControl ~= currControl
        [forces, fstruct] = model.getDrivingForces(schedule.control(currControl));
        [model, state0]= model.updateForChangedControls(state, fstruct);
        prevControl = currControl;
    end
    
    extraArg = {}; % initialGuess
    
    % Solve the time step
    [state, report] = solver.solveTimestep(state0, dt(i), model,...
                                    forces{:}, 'controlId', currControl, extraArg{:});    
    
    % Check for convergence
    if ~report.Converged
        warning('NonLinear:Failure', ...
               ['Nonlinear solver aborted, ', ...
                'returning incomplete results']);
        failure = true;
        break;
    end
    
    % Update the secondary oil saturation for the matrix (for the
    % fracture this is done in updateState() of ReservoirModel
    state.sm(:, 2) = 1 - state.sm(:, 1);
    
    % WIP and OIP
    wipf = sum(model.operators.pv .* state.s(:, 1));
    oipf = sum(model.operators.pv .* state.s(:, 2));
    wipm = sum(model.operators_matrix.pv .* state.sm(:, 1));
    oipm = sum(model.operators_matrix.pv .* state.sm(:, 2));
    
    % Influxes
    qwf_w = sum(state.flux(bwest, 1));
    qof_w = sum(state.flux(bwest, 2));
    qwm_w = sum(state.flux_matrix(bwest, 1));
    qom_w = sum(state.flux_matrix(bwest, 2));    
        
    % Get the fracture & matrix & total phase fluxes at the outflow boundary
    % from mass conservation
    Twm = sum(state.Twm);
    Tom = sum(state.Tom);
    qwf_o = qwf_w - Twm - 1/dt(i)*(wipf - wipf0);
    qof_o = qof_w - Tom - 1/dt(i)*(oipf - oipf0);    
    qwm_o = qwm_w + Twm - 1/dt(i)*(wipm - wipm0);
    qom_o = qom_w + Tom - 1/dt(i)*(oipm - oipm0);  
    
    qw_o = qwf_o + qwm_o;
    qo_o = qof_o + qom_o;
    
    % Updating the previous WIP and OIP
    wipf0 = wipf;
    oipf0 = oipf;
    wipm0 = wipm;
    oipm0 = oipm;   
    
    % Get the cumulative fluxes
    if (i == 1)
        time_prev = state.time;        
        Qw_o = qw_o * dt(i);            
        Qo_o = qo_o * dt(i);            
        Qw_o_prev = Qw_o;   
        Qo_o_prev = Qo_o; 
    else
        Qw_o = state.Qw_o + qw_o * dt(i);            
        Qo_o = state.Qo_o + qo_o * dt(i); 
    end    
   
    state.Qw_o = Qw_o;
    state.Qo_o = Qo_o;

    % Plot the DPDP results together with the reference DFM solutions
    subplot(2, 1, 1)
    hold on
    plot([time_prev state.time]/timescale, [Qo_o_prev Qo_o], 'r')
    if i == 1, legend('Fine-scale simulation', 'MRST'); end
    
    subplot(2, 1, 2)
    hold on
    plot([time_prev state.time]/timescale, [Qw_o_prev Qw_o], 'r')  
    if i == 1, legend('Fine-scale simulation', 'MRST'); end
    
    drawnow
    
    time_prev = state.time;
    Qw_o_prev = Qw_o;   
    Qo_o_prev = Qo_o; 
    
end
    
%{
Copyright 2022 Geological Survey of Denmark and Greenland (GEUS).

Author: Nikolai Andrianov, nia@geus.dk.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

