function surf = getSurfactant(surfNo, varargin)
%(Almost) uncodumented utility function in the co2-foam module. Used to
%retrieve property parameters from laboratory characterisation of
%CO2-brine surfactant systems.
%
% Set up fluid equilibrium properties for the surfactant system.
%
% Cases:
% 1: Test case
% 2-5: Reserved for cases in Vassenden and Holt (2000)
% 6-7: Eclipse field-scale simulations
% 10-15: Commercial surfactants tested in researcher project in 2018-2019.
%
% Defaults to 1.
% Partition coefficient (first column) and 
% maximum adsorption (second column)
    S= [ ...
        1.0  0.2; ... % 10. C14-16 AOS
        1.45 1.0; ... % 11. Tergitol 15-S-9
        0.87 0.2; ... % 12. Tergitol TMN 10
        0.10 0.5; ... % 13. Tergitol NP 10
        0.22 0.2; ... % 14. Igepal CO 720
        0.02 0.1];    % 15. Brij L23

    opt = struct('noPart',false);
    opt = merge_options(opt,varargin{:});

    if surfNo < 1 || surfNo > 15
        warning('Unknown fluid no., defaulting to 1');
        surfNo = 1;
    end

    switch surfNo
        case 1
            Cpart  = 1;
            adsMax = 1e-3;
            adsSat = 5e-4;
        case {2,3,4,5}
            Cpart  = 0;
            adsMax = 1e-3;
            adsSat = 5e-4;
        case {6,7}
            Cpart  = 0;
            adsMax = 1e-3;
            adsSat = 5e-4;
        case {10,11,12,13,14,15}
            Cpart  = S(surfNo-9,1);
            adsMax = S(surfNo-9,2)*1e-3;
            adsSat = 5e-4;
    end
    surf.Cpart = Cpart;
    surf.adsMax = adsMax;
    surf.adsSat = adsSat;
end
