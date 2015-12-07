function wellinfo = getSnohvitInjectionInfo()
% Returns well info corresponding to Snohvit injection project

% SEE ALSO:
%   runSnohvitInjectionScenario.m


    % Physical coordinate of Snohvit injection (approx, see chp 6 pg 135 in
    % Atlas):
    wellinfo.wellCoords_inj       = [9.225e5, 7.988e6];
    wellinfo.wellCoords_prod      = [];
    
    wellinfo.inj_rate_MtperYr     =  0.767; % Mt/yr; (for 30 years)


end