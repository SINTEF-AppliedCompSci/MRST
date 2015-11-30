function compareFaultedModelsJohansen( varargin )
% Compare un-faulted to faulted model of Johansen


moduleCheck('ad-core','mrst-gui')

opt.transMult = []; % 0 if sealing, between 0 and 1 if conducting
% Note: as transmissibility of fault-face is modified only, and not the
% surrounding rock properties, there is not much gained in increasing the
% fault transmissibility by 10 or 10000 times the original value. The
% impact on CO2 migration across the fault is noticable when comparing
% 0.00001 and 0.01 times the original transmissibility value.
opt = merge_options(opt, varargin{:});


%% Load sector model of Johansen
[ grdecl_sect, G_sect, Gt_sect, rock2D_sect ] = makeJohansenSector();

% Get the position of the well (data given from 'Sector5_Well.txt');
wi    = find(G_sect.cells.indexMap==sub2ind(G_sect.cartDims, 48, 48, 1));
wcent = G_sect.cells.centroids(wi,:);
d = sqrt(sum(bsxfun(@minus, Gt_sect.cells.centroids, wcent(1:2)).^2, 2));
[~,wi_sect] = min(d);

% Well on west side of fault line
wi_sect = getCellIndex(Gt_sect, 5.322e5, 6.707e6);


%% Fault location is specified
% Specify cells on opposite sides of fault:
cells2D_1 = find(Gt_sect.cells.ij(:,1) == 43 & Gt_sect.cells.ij(:,2) <= 44);
cells2D_2 = find(Gt_sect.cells.ij(:,1) == 44 & Gt_sect.cells.ij(:,2) <= 44);

% Find the faces at the fault:
facesMat = sparse(double(Gt_sect.faces.neighbors(:,1))+1, ...
   double(Gt_sect.faces.neighbors(:,2))+1, 1:Gt_sect.faces.num);
facesMat = facesMat + facesMat';
faultFaces2D = diag( facesMat(cells2D_1+1, cells2D_2+1) );

% Adjust elevations of cells on opposite sides of fault:
Gt_sect.cells.z([cells2D_1; cells2D_2]) = ...
   G_sect.faces.centroids(Gt_sect.cells.map3DFace([cells2D_1; cells2D_2]), 3);


%% Adjust fault-face transmissibilities

T_all       = getVerticallyAveragedTrans(Gt_sect, rock2D_sect);
T_all_orig  = T_all;
% 
%     Gt = Gt_sect;
%     rock          = rock2D_sect; 
%     rock_tmp.perm = rock.perm .* Gt.cells.H; 
%     T             = computeTrans(Gt, rock_tmp); 
%     cf            = Gt.cells.faces(:, 1); 
%     nf            = Gt.faces.num; 
%     T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);
%     T_all = T;

T_all( faultFaces2D )   = T_all( faultFaces2D ).*opt.transMult;

figure
plotGrid(G_sect, 'EdgeAlpha',0.1,'FaceColor','none')
plotGrid(Gt_sect, 'EdgeAlpha',0.1,'FaceColor','none','EdgeColor','b')

transData1 = zeros(Gt_sect.cells.num,1);
assert(numel(cells2D_1) == numel(faultFaces2D));

transData1( cells2D_1 ) = T_all_orig( faultFaces2D );
avgtransData1 = mean(transData1( cells2D_1 ));
plotCellData(Gt_sect, transData1, cells2D_1, 'EdgeAlpha',0.1)

transData2 = zeros(Gt_sect.cells.num,1);
assert(numel(cells2D_2) == numel(faultFaces2D));

transData2( cells2D_2 ) = T_all( faultFaces2D );
avgtransData2 = mean(transData2( cells2D_2 ));
plotCellData(Gt_sect, transData2, cells2D_2, 'EdgeAlpha',0.1)
colorbar

plotFaces(Gt_sect, faultFaces2D, 'EdgeColor','r','LineWidth',3,'FaceColor','none')
legend('3D grid','2D grid',...
    ['original trans along fault, avg=',num2str(avgtransData1)],...
    ['modified trans along fault, avg=',num2str(avgtransData2)],...
    'fault')
view(3)


%% Simulate

[ initState, fluid, rhow, rhoc ] = simpleSetUp( Gt_sect, rock2D_sect );

% passing in keyword 'trans' will prevent computation of transmissibilities
% based on rock2D, but rather assigns s.T_all = T_all
model = CO2VEBlackOilTypeModel_trans(Gt_sect, rock2D_sect, fluid, 'trans', T_all); 


% Scheduling: (total injected volume, inj/mig time, bdrys, ...)
itime   = 50*year;
istep   = 25;
qtot    = 1*1e9*kilogram*convertTo(itime,year)/(rhoc*kilogram/meter^3); % total m3 over entire injection time
mtime   = 1000*year;
mstep   = 10;

schedule = setSchedule(Gt_sect, rock2D_sect, wi_sect, ...
    qtot, istep, itime, mstep, mtime, true);    % @@ bug if mtime = 1 but mstep = 2;

bdryFaces   = find( Gt_sect.faces.neighbors(:,1).*Gt_sect.faces.neighbors(:,2) == 0 );
bdryVal     = Gt_sect.faces.z(bdryFaces) * rhow*kilogram/meter^3 * norm(gravity);
bc          = addBC( [], bdryFaces, 'pressure', bdryVal, 'sat', [1 0] );
for i = 1:numel(schedule.control)
    schedule.control(i).bc = bc;
end


% Call to solver
[wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule);


%% Assess results
figure; set(gcf,'name',['trans multiplier = ',num2str(opt.transMult)])
plotToolbar(Gt_sect, states)


% % plot co2sat at select reservoir times
% plotCO2sat( Gt_sect, states, sim_report, 1 );
% plotCO2sat( Gt_sect, states, sim_report, 25 );
% plotCO2sat( Gt_sect, states, sim_report, 50 );
% plotCO2sat( Gt_sect, states, sim_report, 1000 );

pause


end

% -------------------------------------------------------------------------

function T = getVerticallyAveragedTrans(Gt, rock)
% code taken from CO2VEBlackOilTypeModel()

    rock_tmp      = rock; 
    rock_tmp.perm = rock.perm .* Gt.cells.H; 
    T             = computeTrans(Gt, rock_tmp); 
    cf            = Gt.cells.faces(:, 1); 
    nf            = Gt.faces.num; 
    T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);
    
end

% -------------------------------------------------------------------------

function [ initState, fluid, rhow, rhoc ] = simpleSetUp( Gt, rock2D )
% Set up initial conditions and model for CO2 injection/migration scenario

    % Fluid data at p = 300 bar
    muw = 0.30860;  rhow = 975.86; sw    = 0.1;
    muc = 0.056641; rhoc = 686.54; srco2 = 0.2;
    kwm = [0.2142 0.85];

    % initial state
    initState.pressure  = Gt.cells.z * norm(gravity) * rhow*kilogram/meter^3;  % hydrostatic pressure, in Pa=N/m^2
    initState.s         = repmat([1 0], Gt.cells.num, 1);                      % sat of water is 1, sat of CO2 is 0
    initState.sGmax     = initState.s(:,2);                                         % max sat of CO2 is initially 0
    initState.rs        = 0 * initState.sGmax;                                      % initially 0 

    ref_p               = mean(initState.pressure);
    info                = getSeaInfo('NorthSea',760);
    caprock_temperature = 273.15 + info.seafloor_temp + ...
            (Gt.cells.z - info.seafloor_depth) / 1e3 * info.temp_gradient; % Kelvin

    fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
              'fixedT'      , caprock_temperature, ...
              'wat_mu_ref'  , muw * centi*poise, ...
              'co2_mu_ref'  , muc * centi*poise, ...
              'wat_rho_ref' , rhow * kilogram/meter^3, ...
              'co2_rho_ref' , rhoc * kilogram/meter^3, ...
              'wat_rho_pvt' , [], ...
              'pvMult_p_ref', ref_p, ...
              'pvMult_fac'  , 0, ...
              'residual'    , [sw, srco2] , ...
              'dissolution' , false, 'dis_max', []);

end

% -------------------------------------------------------------------------

function plotCO2sat( Gt, states, sim_report, year2plot )
% plot co2 saturation in grid at specified year (since simulation start)

    %[inx] = find(sim_report.ReservoirTime == convertFrom(year2plot,year));
    plotTime_s = convertFrom( year2plot, year );
    diff = abs(bsxfun(@minus, sim_report.ReservoirTime, plotTime_s));
    minDiff = min( diff );
    inx = find(diff == minDiff, 1); % returns first index (earlier state)

    
    co2sat = states{inx}.s(:,2);
    press  = states{inx}.pressure;
    
    figure; set(gcf,'Position',[1 1 1500 500])
    
    subplot(1,2,1)
    title(['CO_2 saturation after ',num2str(year2plot),' year(s)'])
    plotGrid(Gt, 'EdgeAlpha',0.1,'FaceColor','none')
    plotCellData(Gt, co2sat, co2sat>0.01, 'EdgeAlpha',0.1)
    axis tight off; view(3)
    colorbar
    
    subplot(1,2,2)
    title(['Pressure after ',num2str(year2plot),' year(s)'])
    plotGrid(Gt, 'EdgeAlpha',0.1,'FaceColor','none')
    plotCellData(Gt, press, 'EdgeAlpha',0.1)
    axis tight off; view(3)
    colorbar

end

% -------------------------------------------------------------------------

function cellIndex = getCellIndex(Gt, Xcoord, Ycoord)
% Closest cell index of grid Gt corresponding to physical coordinate (X,Y)

    dv        = bsxfun(@minus, Gt.cells.centroids(:,1:2), [Xcoord, Ycoord]);
    [v, ind]  = min(sum(dv.^2, 2));
    cellIndex = ind; 
    % or Gt.cells.indexMap(i);
    
end

