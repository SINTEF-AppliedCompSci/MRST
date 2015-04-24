require vertical-equil-lab
require mex/libgeometry
require coarsegrid

clear;
close all;

%% Constants, etc.
%data_directory = 'modules/bitbucket/co2lab-matlab/vertical-equil-lab/examples/igems/data/surfaces';
mydir=mfilename('fullpath');i=find(mydir=='/',1,'last');mydir=mydir(1:i);
data_directory = fullfile(mydir,'../../examples/igems/data/surfaces');
surface_name   = 'OSSNP2';
surface_number = 1;
coarsening = [2 2]; % [4 4]
%coarsening = [4 4]; % [4 4]
perm = 500 * milli * darcy;
poro = 0.2;
z_levels = 10;  % 5
z_height = 100;
gravity on;

%keyboard;
% %{
%% Loading file, and creating associated rock structure
disp('loading file');
G = readIGEMSSurface(data_directory, ...
                     surface_name, ...
                     surface_number, ...
                     z_levels, ...
                     z_height, ...
                     coarsening);
G.cartDims
disp('finished');
rock = setupRockStructure(G, perm, poro);

%% Creating 2D version of the grid and the rock structure
% %}
Gt     = topSurfaceGrid(G);
rock2D = averageRock(rock, Gt);

%% @@ VISUALIZE/EXPLAIN TOP SURFACE
%-----------------------------------
Gt_lifted = Gt;
Gt_lifted.nodes.z = Gt.nodes.z - 100;
plotGrid(G, 'edgeAlpha', 0.05); 
plotGrid(Gt_lifted, 'faceColor', 'blue','edgeAlpha', 0.1); view([-20 10]);

disp('Hit a key to continue');
pause;

%% Run the trapping structure analysis
disp('running trapping structure analysis');
tic;
[z_spill_loc, g_trap, trap_level, z_spill_level, z_spill_loc_level, Gnew] = ...
    findTrappingStructure_dfs(Gt);
toc;
disp('finished');

% Below, the connection between traps is determined, as well as gathering
% together traps in groups having common highest levels.
tic;
[cell_lines, traps] = findCellLines(Gt, z_spill_loc);
toc;
%keyboard;
% Print basic information about the traps found
num_shown_traps = 10
table = printTrapInfoTable(Gt, z_spill_loc, trap_level, traps, -4, num_shown_traps, rock2D);


%% @@ ADD PIE CHART OF TRAPS
% - use colors to refer to entries in TrapInfoTable
%--------------------------
n_largest_vols = table(1:min(size(table,1), num_shown_traps), 4);
sum_rest = sum(table(:, 4)) - sum(n_largest_vols);
trap_labels = cellstr(num2str([1:num_shown_traps]')); 
trap_labels{num_shown_traps+1} = 'REST';
figure(1);
pie([n_largest_vols;sum_rest], trap_labels);

plot_trap_colors = zeros(size(traps));
% assigning colors oto explicitly labeled traps
for i = 1:num_shown_traps
    plot_trap_colors(find(traps==table(i, 1))) = i;
end
% assigning a common color to the remaining traps
plot_trap_colors(find(plot_trap_colors == 0 & traps~=0)) = num_shown_traps+1;

figure(2);
plotGrid(Gt, 'FaceAlpha', 0, 'EdgeAlpha', 0.1);
plotGrid(Gt, [cell_lines{:}], 'FaceColor', 'k'); % plotting 'rivers'
plotCellData(Gt, plot_trap_colors, find(traps~=0)); % plotting traps
view(-13, 54);



%% Visualizing individual traps
num_traps = max(traps);
while true
    fprintf('Locate trap with index (1 to %d).\n', num_traps)
    choice = input('Outside bounds means end: ');
    
    % end if choice is outside range
    if isempty(choice) || ~ismember(choice, 1:num_traps) break; end; 
                                                  
    % otherwise, plot result
    clf; plotCellData(Gt, double(traps~=0) + double(traps==choice));
end


%% @@ SPILL POINT ANALYSIS
%--------------------------


%% @@ Running simulation

% time parameters
T          = 750 * year();
stopInject = 150 * year();
dT         =   1 * year()/12;
dTplot     =   2 * dT();

% fluid parameters
% fluidVE_h = initVEFluidHForm(Gt, ...
%                              'mu' , [0.056641 0.30860] .* centi * poise, ...
%                              'rho', [686.54 975.86] .* kilogram/meter^3, ...
%                              'sr' , 0.2, ...
%                              'sw' , 0.1, ...
%                              'kwm', [0.2142 0.85]);
fluidVE = initSimpleVEFluid_s('mu',  [0.056641 0.30860] .* centi * poise, ...
                              'rho', [686.54 975.86] .* kilogram/meter^3, ...
                              'height', Gt.cells.H, ...
                              'sr', [0.2, 0.1]) % [co2, water]

% well 
wellIx = [G.cartDims(1:2)/5, G.cartDims([3 3])];
rate = 2.8e4*meter^3/day;
% WVE_old = convertwellsVE(verticalWell([], G, rock, wellIx(1), wellIx(2), ...
%                                       wellIx(3):wellIx(4), 'Type', 'rate', 'Val', ...
%                                       rate, 'Radius', 0.1, 'comp_i', [1 0], 'name', ...
%                                       'I'), ...
%                          G, Gt, rock2D);

[i2, j2] = ind2sub(Gt.cartDims, Gt.cells.indexMap);
wellCellIx = find(i2 == Gt.cartDims(1)/5 & j2 == Gt.cartDims(2)/5);
WVE = addWell([], Gt, rock2D, wellCellIx, 'Type','rate', 'Val', rate, ...
              'Radius', 0.1, 'comp_i', [1 0], 'InnerProduct','ip_simple', ...
              'name', 'I');
WVE.h = Gt.cells.H(WVE.cells); % VE-specific parameter


% boundary conditions
i = any(Gt.faces.neighbors==0, 2);  % find all outer faces
bcIxVE = Gt.cells.faces(i(Gt.cells.faces(:,1)), 1);
[~, rho] = fluidVE.properties();
bcVE  = addBC([], bcIxVE, 'pressure', Gt.faces.z(bcIxVE) * rho(2) * ...
              norm(gravity));
bcVE  = rmfield(bcVE, 'sat'); % no inflow faces, so not relevant
bcVE.h = zeros(size(bcVE.face)); % needed in 'computePressureRHSVE'

% inner product matrices
%SVE = computeMimeticIPVE(Gt, rock2D, 'Innerproduct', 'ip_simple');
%preComp = initTransportVE(Gt, rock2D);
SVE = computeTrans(Gt, rock2D); 
% multiply by height to get real transmisibilities for 2D case
SVE = SVE .* Gt.cells.H(rldecode(1:Gt.cells.num, diff(Gt.cells.facePos), 2)');


% solution structure
sol = initResSolVE(Gt, 0, 0, 'use_s_form', true);
sol.wellSol = initWellSol(WVE, 300*barsa());  % @@ need to use 3D well here?
% sol.s = height2Sat(sol, Gt, fluidVE); 

% setup plotting @@ Add plotting of percent physically trapped
[~, ~, res_sat] = fluidVE.properties();
opts = {'slice', wellIx, 'Saxis', [0 1-res_sat(2)], ...
        'maxH', 200, 'Wadd', 1000, 'wireH', true};

% main loop
t = 0;
trap_cells = find(z_spill_loc); %% identify cells that physically trap fluid
z_spill_depth = zeros(size(z_spill_loc, 1),1);
z_spill_depth(trap_cells) = z_spill_loc(trap_cells) - Gt.cells.z(trap_cells);

while t < T
    % advance solution: compute pressure, then transport
    % sol = solveIncompFlowVE(sol, Gt, SVE,rock,fluidVE,'bc',bcVE,'wells', WVE);
    % sol = explicitTransportVE(sol,Gt,dT,rock,fluidVE,'bc',bcVE,'wells',WVE, ...
    %                           'preComp', preComp,'intVert', false);
    sol = incompTPFA(sol, Gt, SVE, fluidVE, 'bc', bcVE,'wells', WVE);
    sol = implicitTransport(sol, Gt, dT, rock2D, fluidVE, 'wells', WVE);
    
    
    % reconstruct 'saturation' @@ PERHAPS WE NEED TO CONSIDER DEPTH OF 3D
    %sol.s = height2Sat(sol, Gt, fluidVE);
    %assert(max(sol.s(:,1))<1+eps && min(sol.s(:,1))>-eps);
    
    % reconstruct height
    sol.h = fluidVE.sat2height(sol);
    sol.h_max = max(sol.h, sol.h_max);
    
    % check if injection has ended
    if t>=stopInject
        WVE = []; bcVE = []; dT=5*year(); dTplot = dT;
    end
    
    % compute trapped and free volumes of CO2 @@ MODIFY THIS ONE
   fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d years', convertTo(t,year));

   res_trappedVol = sum((sol.h_max-sol.h).*rock2D.poro.*Gt.cells.volumes)* res_sat(1);
   % str_trappedVol = (rock2D.poro(trap_cells) .* Gt.cells.volumes(trap_cells))' * ...
   %                   min(sol.h(trap_cells), z_spill_depth(trap_cells)) * (1-res_sat(2));    

   str_trappedVol = (rock2D.poro(trap_cells) .* Gt.cells.volumes(trap_cells))' * ...
                     min(sol.h(trap_cells), z_spill_depth(trap_cells)) * (1-res_sat(2));    
      
   freeVol = sum(sol.h.*rock2D.poro.*Gt.cells.volumes)*(1-res_sat(2)) ...
             - str_trappedVol;

   
% $$$    trappedVol = sum((sol.h_max-sol.h).*rock2D.poro.*Gt.cells.volumes)* ...
% $$$                 fluidVE.sr;
% $$$    totVol = trappedVol + freeVol;    
   
   % @@ PHYSICAL TRAPPING OVERRIDE
% $$$    trappedVol = (rock2D.poro(trap_cells) .* Gt.cells.volumes(trap_cells))' * ...
% $$$                  min(sol.h(trap_cells), z_height);
   
% $$$    freeVol = totVol - trappedVol;
   
   
   
    % plotting
    if mod(t,dTplot)==0 || t >=T
        %plotPanelVEWithStructuralTrapping(G, Gt, WVE, sol,
        %t,[freeVol,trappedVol,totVol], opts{:});
        
        % for plotting, we need the full 3D profile of the saturation, so we generate
        % it here.  @@ Find a more elegant solution!
        sol_dummy = sol;
        sol_dummy.s = height2Sat(sol, Gt, fluidVE);
        
        plotPanelVEWithStructuralTrapping(G, Gt, WVE, sol_dummy, t, ...
                                          [freeVol,res_trappedVol, str_trappedVol], ...
                                          opts{:});
        
        drawnow;
    end
    
    t = t+dT;
end
