%% test pressure-controlled wells


% determine the bhp of the wells based on the initial pressure in the
% formation (at the caprock), and the overpressure allowed (i.e., the
% difference between the bhp and the initial pressure). As formation
% pressure increases during the injection period, the overpressure will be
% decreasing. Thus the initial condition is used to set the bhp's of the wells.

[Gt, rock2D] = getFormationTopGrid('Stofm',5);
seainfo = getSeaInfo('Stofm',760);
trapCapacities = getTrappingInfo(Gt, rock2D, seainfo, 'plotsOn',false);

%wellinfo = getSnohvitInjectionInfo();
%wellinfo = rmfield(wellinfo,'inj_rate_MtperYr');
%wellinfo = rmfield(wellinfo,'wellCoords_prod');
% for i = 1:size(wellinfo.wellCoords_inj,1)
%     wcellinx(i) = getCellIndex(Gt, ...
%         wellinfo.wellCoords_inj(i,1), wellinfo.wellCoords_inj(i,2));
% end

wellinfo = getWellInfo(Gt, trapCapacities, 'prod',false, 'limits','none');
%maxOverP = 10 * mega * Pascal;
%wellinfo.BHP = trapCapacities.caprock_pressure(wellinfo.cinx_inj) + maxOverP;

% pass in wellinfo.initOverP and BHP is computed internally
wellinfo.initOverP = 1 * mega * Pascal;
[wellSols, states, sim_report, opt, var] = ...
    runGenericInjectionScenario(Gt, rock2D, seainfo, wellinfo, ...
                                        'wellType',     'bhp');
