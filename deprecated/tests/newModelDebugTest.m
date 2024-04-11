mrstModule add ad-core ad-props co2lab
gravity reset on;

%% fluid parameters
g       = gravity;
rhow    = 1000;                                 % water density (kg/m^3)
co2     = CO2props();                           % CO2 property functions
p_ref   = 30 *mega*Pascal;                      % reference pressure
t_ref   = 94+273.15;                            % reference temperature
co2_rho = co2.rho(p_ref, t_ref);                % CO2 density
co2_c   = co2.rhoDP(p_ref, t_ref) / co2_rho;    % CO2 compressibility
wat_c   = 0;                                    % water compressibility
c_rock  = 4.35e-5 / barsa;                      % rock compressibility
srw     = 0.1;                                 % residual water
src     = 0.2;                                 % residual CO2
pe      = 5 * kilo * Pascal;                    % capillary entry pressure
muw     = 8e-4 * Pascal * second;               % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity

%% grid and rock
xres = 200;
zres = 200;

xlen = 10 * kilo * meter;
zlen = 100 * meter;

G = cartGrid([xres, 1, zres], [xlen, 1, zlen]);

% make grid sloping
slope = 0.05;
G.nodes.coords(:,3) = 2000 + G.nodes.coords(:,3) - G.nodes.coords(:,1) * slope;
G = computeGeometry(G);
[Gt, G] = topSurfaceGrid(G);

rock.poro = 0.2 * ones(G.cells.num, 1);
rock.perm = 300 * milli * darcy * ones(G.cells.num, 1);
rock2D = averageRock(rock, Gt);

%% initial state

initState.pressure = rhow * g(3) * Gt.cells.z;
initState.s = repmat([1,0], Gt.cells.num, 1);
initState.sGmax = initState.s(:,2);
initState.rs = zeros(Gt.cells.num, 1);

%% make fluid model
invPc3D = @(pc) (1-srw) .* (pe./max(pc, pe)).^2 + srw;
kr3D    = @(s) max((s-src)./(1-src), 0).^2; % uses CO2 saturation
fluid   = makeVEFluid(Gt, rock2D, 'si_rugosity'     , ...
               'co2_mu_ref'  , muco2, ...
               'wat_mu_ref'  , muw, ...
               'co2_rho_ref' , co2_rho                , ...
               'wat_rho_ref' , rhow                   , ...
               'co2_rho_pvt' , [co2_c, p_ref]         , ...
               'wat_rho_pvt' , [wat_c, p_ref]         , ...
               'residual'    , [srw, src]             , ...
               'pvMult_p_ref', p_ref                  , ...
               'pvMult_fac'  , c_rock                 , ...
               'invPc3D'     , invPc3D                , ...
               'kr3D'        , kr3D                   , ...
               'dissolution' , true, ...
               'dis_rate'    , 0, ...
               'dis_max'     , 0.03);
               

%% Setup simulation schedule

% hydrostatic pressure at xmin and xmax, and closed everywhere else
bc = pside([], Gt, 'xmax', Gt.cells.z(end) * rhow * g(3));
%bc = pside(bc, Gt, 'xmin', Gt.cells.z(1) * rhow *g(3));
bc.sat = repmat([1 0], 2, 1);

% well
wcell = ceil(xres/4);
%inj_rate = 0.045 * mega * 1e3 / year/ co2_rho;
inj_rate = 0.02 * mega * 1e3 / year/ co2_rho;
W2D = addWell([], Gt, rock2D, wcell, 'type', 'rate', ...
              'val', inj_rate, 'comp_i', [0 1]);
W2D.WI = W2D.WI .* Gt.cells.H(wcell);
W2D.h = Gt.cells.H(wcell);
W2D.dZ = 0;

schedule.control    = struct('W', W2D, 'bc', bc);
schedule.control(2) = struct('W', W2D, 'bc', bc);
schedule.control(2).W.val = 0;
%schedule.control(2).W.status = false;
schedule.step.val = [repmat(day, 50, 1); 
                     repmat(200 * day, 50, 1)];
schedule.step.control = [ones(50, 1); 2 * ones(50, 1)];

%% create model and run simulation
% model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
% [wellSol, states] = simulateScheduleAD(initState, model, schedule);

model2 = CO2VEBlackOilTypeModelNew(Gt, rock2D, fluid);

model2.fluid.dis_rate = 1e-8;
model2.fluid.dis_max = 0.03;
[wellSol2, states2] = simulateScheduleAD(initState, model2, schedule);


% plot saturation, residual saturation and dissolved CO2
figure
for i = 1:numel(states2)
    clf; hold on;
    [h, hmax] = upscaledSat2height(states2{i}.s(:,2), states2{i}.sGmax, Gt, ...
                                   'resSat', [srw, src]);
    % @@@ investigate possibly negative h
    
    plot(h/zlen);
    plot(hmax/zlen);
    plot(states2{i}.rs*10);
    legend('sG', 'sGmax', 'rs');
    pause(0.4);
end






% plot residual saturation
disgas = isfield(model.fluid, 'dis_max');
for i = 1:numel(states)
    [h, hmax] = upscaledSat2height(states{i}.s(:,2), states{i}.sGmax, Gt, ...
                                   'resSat', [srw, src]);
    s = height2Sat(h, hmax, Gt, srw, src);
    
    if disgas
        rs = repmat(states{i}.rs(:), zres, 1) * 5;
        s = s + rs;
        s(s>1) = 1;
    end
    
    
    
    [h2, hmax2] = upscaledSat2height(states2{i}.s(:,2), states2{i}.sGmax, Gt, ...
                                     'resSat', [srw, src]);
    s2 = height2Sat(h2, hmax2, Gt, srw, src);

    if disgas
        rs = repmat(states2{i}.rs(:), zres, 1) * 5;
        s2 = s2 + rs;
        s2(s2>1) = 1;
    end
    
    clf;
    subplot(1,2,1); plotCellData(G, s, 'edgealpha', 0.1); view(0,0); 
    subplot(1,2,2); plotCellData(G, s2, 'edgealpha', 0.1); view(0,0); 
    %clf; plotCellData(G, s2, 'edgealpha', 0.1); view(0,0); 
    pause(0.05)
end


[h2, hmax2] = upscaledSat2height(states2{50}.s(:,2), states2{50}.sGmax, Gt, ...
                                 'resSat', [srw, src]);
s2 = height2Sat(h2, hmax2, Gt, srw, src);
subplot(1,2,1); plotCellData(G, s2, 'edgealpha', 0.1); view(0,0); 

[h2, hmax2] = upscaledSat2height(states2{51}.s(:,2), states2{51}.sGmax, Gt, ...
                                 'resSat', [srw, src]);
s2 = height2Sat(h2, hmax2, Gt, srw, src);
subplot(1,2,2); plotCellData(G, s2, 'edgealpha', 0.1); view(0,0); 