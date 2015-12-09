%% Compare results obtained using different boundary conditions.
% The Sto formation is studied. An array of injection wells is placed to
% cover the entire formation. Well rates are computed based on the trapping
% capacity surrounding the well. (NB: only well coordinates need to be
% passed in by wellinfo, as well cell indices are computed internally)

gravity on

[Gt, rock2D]    = getFormationTopGrid('Utsirafm',5);
seainfo         = getSeaInfo('Utsirafm',760);

trapCapacities  = getTrappingInfo(Gt, rock2D, seainfo, 'plotsOn',false);
wellinfo        = getWellInfo(Gt, trapCapacities, ...
                        'limits','none', ...
                        'prod',false, ...
                        'setInjRates',true, ...
                        'buffer', 5000, ...
                        'DX', 6*5000, ...
                        'DY', 6*5000 );
                    
% A) Pressure boundaries.
[ wellSols_p, states_p, sim_report_p, opt_p, var_p ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType','pressure');


% B) No-flow boundaries.
[ wellSols_f, states_f, sim_report_f, opt_f, var_f ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType','flux');
                             
                             
% Compare the results:
% compareTwoSims( 'pressure, array', Gt, rock2D, wellSols_p, states_p, sim_report_p, opt_p, var_p, ...
%    'no flow,, array', Gt, rock2D, wellSols_f, states_f, sim_report_f, opt_f, var_f);


%%% Instead of using an array of wells, we can use one injector
% We use the Snohvit injection details, however we increase the injection
% rate by a factor of 10 in order to make the impact of the two boundary
% conditions more noticable.
wellinfo = getSnohvitInjectionInfo();
wellinfo.inj_rate_MtperYr = wellinfo.inj_rate_MtperYr .* 10;


% A) Pressure boundaries.
[ wellSols_p2, states_p2, sim_report_p2, opt_p2, var_p2 ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType','pressure');


% B) No-flow boundaries.
[ wellSols_f2, states_f2, sim_report_f2, opt_f2, var_f2 ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo, ...
                                 'bdryType','flux');
                             
                             
% Compare the results:
% compareTwoSims( 'pressure, single', Gt, rock2D, wellSols_p2, states_p2, sim_report_p2, opt_p2, var_p2, ...
%    'no flow, single', Gt, rock2D, wellSols_f2, states_f2, sim_report_f2, opt_f2, var_f2);


%% Single injector with faulted formation
% Now we represent faults in the Sto formation by detecting fault faces,
% and passing in a transMult value, which then causes transmissibilities to
% be updated for these fault faces. Simulations are then run using the
% different boundary conditions.

% NB: grid, rock, and faultFace inx may change during call to detectFaults()
[ Gt_fault, faultFaces2D, ~, ~, rock2D_fault ] = detectFaults(Gt, rock2D, 'gradTol',100);

% grid and rock may have changed during call to detectFaults(), however the
% impact on wellinfo will only be felt if a port of the grid/rock has been
% cut or isolated by faults. Thus recompute wellinfo if using an array of
% wells and/or using trapCapacities. Recompute trapCapacities first.

% A) Pressure boundaries.
[ wellSols_p3, states_p3, sim_report_p3, opt_p3, var_p3 ] = ...
    runGenericInjectionScenario( Gt_fault, rock2D_fault, seainfo, wellinfo, ...
                                 'bdryType',    'pressure', ...
                                 'faultFaces',  faultFaces2D, ...
                                 'transMult',   0.01);


% B) No-flow boundaries.
[ wellSols_f3, states_f3, sim_report_f3, opt_f3, var_f3 ] = ...
    runGenericInjectionScenario( Gt_fault, rock2D_fault, seainfo, wellinfo, ...
                                 'bdryType',    'flux', ...
                                 'faultFaces',  faultFaces2D, ...
                                 'transMult',   0.01);
                             
% Compare the results:
% compareTwoSims( 'pressure, single, faulted', Gt, rock2D, wellSols_p, states_p, sim_report_p, opt_p, var_p, ...
%    'no flow, single, faulted', Gt, rock2D, wellSols_f, states_f, sim_report_f, opt_f, var_f);


% use array of wells again, for this faulted grid.
trapCapacities  = getTrappingInfo(Gt_fault, rock2D_fault, seainfo, 'plotsOn',false);
wellinfo        = getWellInfo(Gt_fault, trapCapacities, ...
                        'limits','none', ...
                        'prod',false, ...
                        'setInjRates',true, ...
                        'buffer', 5000, ...
                        'DX', 6*5000, ...
                        'DY', 6*5000 );
      

% A) Pressure boundaries.
[ wellSols_p4, states_p4, sim_report_p4, opt_p4, var_p4 ] = ...
    runGenericInjectionScenario( Gt_fault, rock2D_fault, seainfo, wellinfo, ...
                                 'bdryType',    'pressure', ...
                                 'faultFaces',  faultFaces2D, ...
                                 'transMult',   0.01);


% B) No-flow boundaries.
[ wellSols_f4, states_f4, sim_report_f4, opt_f4, var_f4 ] = ...
    runGenericInjectionScenario( Gt_fault, rock2D_fault, seainfo, wellinfo, ...
                                 'bdryType',    'flux', ...
                                 'faultFaces',  faultFaces2D, ...
                                 'transMult',   0.01);
                             
% Compare the results:
% compareTwoSims( 'pressure, array, faulted', Gt, rock2D, wellSols_p, states_p, sim_report_p, opt_p, var_p, ...
%    'no flow, array, faulted', Gt, rock2D, wellSols_f, states_f, sim_report_f, opt_f, var_f);


%%% Compare the results:
% 1.
compareTwoSims( 'pressure, array', Gt, rock2D, wellSols_p, states_p, sim_report_p, opt_p, var_p, ...
   'no flow,, array', Gt, rock2D, wellSols_f, states_f, sim_report_f, opt_f, var_f);

% 2.
compareTwoSims( 'pressure, single', Gt, rock2D, wellSols_p2, states_p2, sim_report_p2, opt_p2, var_p2, ...
   'no flow, single', Gt, rock2D, wellSols_f2, states_f2, sim_report_f2, opt_f2, var_f2);


% 3.
compareTwoSims( 'pressure, single, faulted', Gt_fault, rock2D_fault, wellSols_p3, states_p3, sim_report_p3, opt_p3, var_p3, ...
   'no flow, single, faulted', Gt_fault, rock2D_fault, wellSols_f3, states_f3, sim_report_f3, opt_f3, var_f3);

                             
% 4.
compareTwoSims( 'pressure, array, faulted', Gt_fault, rock2D_fault, wellSols_p4, states_p4, sim_report_p4, opt_p4, var_p4, ...
   'no flow, array, faulted', Gt_fault, rock2D_fault, wellSols_f4, states_f4, sim_report_f4, opt_f4, var_f4);









