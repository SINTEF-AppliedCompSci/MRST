%% Compare results obtained using different boundary conditions.
% The Sto formation is studied. An array of injection wells is placed to
% cover the entire formation. Well rates are computed based on the trapping
% capacity surrounding the well. (NB: only well coordinates need to be
% passed in by wellinfo, as well cell indices are computed internally)

gravity on

[Gt, rock2D]    = getFormationTopGrid('Stofm',5);
seainfo         = getSeaInfo('Stofm',760);
trapCapacities  = getTrappingInfo('Stofm',5, 'plotsOn',false);
wellinfo        = getWellInfo(Gt, trapCapacities, ...
                        'limits','none', ...
                        'prod',false, ...
                        'setInjRates',true);
                    
% A) Pressure boundaries.
[ wellSols_p, states_p, sim_report_p, opt_p, var_p ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType','pressure');


% B) No-flow boundaries.
[ wellSols_f, states_f, sim_report_f, opt_f, var_f ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType','flux');
                             
                             
% Compare the results:
compareTwoSims( 'pressure', Gt, rock2D, wellSols_p, states_p, sim_report_p, opt_p, var_p, ...
   'no flow', Gt, rock2D, wellSols_f, states_f, sim_report_f, opt_f, var_f);


%%% Instead of using an array of wells, we can use one injector
% We use the Snohvit injection details, however we increase the injection
% rate by a factor of 10 in order to make the impact of the two boundary
% conditions more noticable.
wellinfo = getSnohvitInjectionInfo();
wellinfo.inj_rate_MtperYr = wellinfo.inj_rate_MtperYr .* 10;


% A) Pressure boundaries.
[ wellSols_p, states_p, sim_report_p, opt_p, var_p ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType','pressure');


% B) No-flow boundaries.
[ wellSols_f, states_f, sim_report_f, opt_f, var_f ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType','flux');
                             
                             
% Compare the results:
compareTwoSims( 'pressure', Gt, rock2D, wellSols_p, states_p, sim_report_p, opt_p, var_p, ...
   'no flow', Gt, rock2D, wellSols_f, states_f, sim_report_f, opt_f, var_f);


%% Single injector with faulted formation
% Now we represent faults in the Sto formation by detecting fault faces,
% and passing in a transMult value, which then causes transmissibilities to
% be updated for these fault faces. Simulations are then run using the
% different boundary conditions.

% NB: grid, rock, and faultFace inx may change during call to detectFaults()
[ Gt, faultFaces2D, ~, ~, rock2D ] = detectFaults(Gt, rock2D, 'gradTol',100);

% grid and rock may have changed during call to detectFaults(), however the
% impact on wellinfo will only be felt if a port of the grid/rock has been
% cut or isolated by faults. Thus recompute wellinfo if using an array of
% wells and/or using trapCapacities. Recompute trapCapacities first.

% A) Pressure boundaries.
[ wellSols_p, states_p, sim_report_p, opt_p, var_p ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType',    'pressure', ...
                                 'faultFaces',  faultFaces2D, ...
                                 'transMult',   0.01);


% B) No-flow boundaries.
[ wellSols_f, states_f, sim_report_f, opt_f, var_f ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType',    'flux', ...
                                 'faultFaces',  faultFaces2D, ...
                                 'transMult',   0.01);
                             
% Compare the results:
compareTwoSims( 'pressure', Gt, rock2D, wellSols_p, states_p, sim_report_p, opt_p, var_p, ...
   'no flow', Gt, rock2D, wellSols_f, states_f, sim_report_f, opt_f, var_f);















