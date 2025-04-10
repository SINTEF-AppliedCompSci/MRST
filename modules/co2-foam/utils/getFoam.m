function foam = getFoam(foamNo, varargin)
%Set up foam properties for the CO2 foam module 
%using examples from literature.
%
% SYNOPSIS:
%       foam = getFoam(n)
%
% PARAMETERS:
%   n  - Retrieve parameters for a SINTEF foam model.
%        Use parameter set number n.
%
% RETURNS:
%   foam - Structure with fields
%           sF - 
%           s1 -
%           s2 -
%           F0 -
%           v0 -

% 1:   Test case for CO2 Foam development
% 2-5: From Table 1 in Vassenden and Holt (2000)
% 6-7: Match Eclipse simulations
% 10-15: Parameters from testing of commercial CO2/brine foam surfactants
%        in 2018-2019. Used in 2019 NCCS deliverable.
% Defaults to 1 if input argument is <1 or >20.

opt = struct('noShear',false);
opt = merge_options(opt,varargin{:});

    if foamNo < 1 || foamNo > 20
        warning('Unknown foam no., switching to foam 1');
        foamNo = 1;
    end
    
    switch foamNo
        
        % From SINTEF petroleum foam example.
        case 1
            
            sF =   0.1  ;
            s1 =  50    ;
            s2 =   1    ;
            F0 =   1e-3 ;
            v0 =   1/day;
        
        % From table 1 in Vassenden and Holt (2000)
        case 2

            sF =        0;
            s1 =      500;
            s2 =       60;
            F0 =     2e-4;
            v0 = 0.25/day;

        
        case 3

            sF =     0.13;
            s1 =   237.5 ;
            s2 =     1   ;
            F0 =     9e-5;
            v0 = 0.26/day;

        case 4

            sF =    0.02;
            s1 = 1743   ;
            s2 =    1   ;
            F0 =    5e-6;
            v0 =   2/day;

        case 5
            
            sF =   0.005;
            s1 = 525    ;
            s2 =   1    ;
            F0 =   3e-5 ;
            v0 = 6.6/day;
            
        % Modified from case 1 (PB)
        case 6
            
            sF =   0.1;
            s1 =   50 ;
            s2 =   1  ;
            F0 =   1e-1 ;
            v0 = 6/day;

        % Modified from case 1 (AAG). Try to match Eclipse models
        % For use with a no shear model
        case 7
            sF =   0.1;
            s1 =   70 ;
            s2 =   0  ;
            F0 =   1e-2 ;
            v0 =   1;
 
		case 10
            % AOS anionic surfactant
            sF = 0.07;
            s1 = 250;
            s2 = 0;
            F0 = 1/4; % Originally 1/40
            v0 = 0.5/day; % m/day
        case 11
            % Tergitol 15-S-9
            sF = 0.12;
            s1 = 250;
            s2 = 0;
            F0 = 1/2;  % Originally 1/5
            v0 = 0.5/day; % m/day
        case 12
            % Tergitol TMN10
            sF = 0.142;
            s1 = 300;
            s2 = 0;
            F0 = 1/50;
            v0 = 0.5/day; % m/day
        case 13
            % Tergitol NP10
            sF = 0.14;
            s1 = 400;
            s2 = 0;
            F0 = 1/200;
            v0 = 0.5/day; % m/day
        case 14
            % Igepal CO 720
            sF = 0.115;
            s1 = 250;
            s2 = 0;
            F0 = 1/100;
            v0 = 0.5/day; % m/day
        case 15
            % Brij L23
            sF = 0.125;
            s1 = 300;
            s2 = 0;
            F0 = 1/10;  % Originally 1/30
            v0 = 0.5/day; % m/day
   end
 
	if opt.noShear
		v0 = 1;
	end
    foam = struct(        ...
           'sF'   , sF  , ... % Foam effective for sW >= sF
           's1'   , s1  , ... % Reduction factor slope for capillary regime
           's2'   , s2  , ... % Reduction factor slope for pres grad regime
           'F0'   , F0  , ... % Reduction factor for sW >> sF
           'v0'   , v0  , ... % Reference gas velocity (absoule value)
           'cfrac', 0.1);     % Foam effective for c >= cfrac*cmax
    % Quadratic permeability dependence. Reference permeability 500 mD.
    foam.PermDep = @(k) (k/(500*milli*darcy)).^2; 
end
        