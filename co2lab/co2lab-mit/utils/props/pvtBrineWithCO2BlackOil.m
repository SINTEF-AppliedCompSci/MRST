function [t, rho_co2_s, rho_brine_s, t2] = pvtBrineWithCO2BlackOil(T, P, S, ...
                                        saltVar, vapH2O, unsatVals, figs)
% SYNOPSIS:
% [t, rho_co2_s, rho_brine_s, t2] = pvtBrineWithCO2BlackOil(T, P, S, saltVar, vapH2O, unsatVals, figs)
%
% Lluís Saló-Salgado (lsalo@mit.edu), September 2019
% Updates:
% LS, July 2022:  adjust max p vals and update plots
% LS, April 2023: updates to consider fully ionized salts, no iteration in
%                 mixing rules, and update CO2 equilibrium constants 
%                 (kap0_co2) as well as molar volumes in the RK1949 EoS
% LS, July 2023:  update w_salt for rho_brine calculation (no CO2) and 
%                 output variables
%
% DESCRIPTION:
% Compute and save a table with R_s, P_aq, B_aq, mu_aq for a saline aqueous
% solution with dissolved CO2, and another table with (R_s), P_gas, B_gas 
% and mu_gas for a CO2-rich gaseous phase (with or without dissolved H2O).
%
% Intended for usage with a black-oil formulation, within the range tested
% by Spycher et al. (references below): T up to ~100 ºC and P up to ~600 bar.
% At higher pressures and temperatures, the assumptions used below are no
% longer valid, and the solubility model will start to lose accuracy,
% particularly in the composition of the gaseous phase.
%
% MAIN REFERENCES:
% Duan and Sun, Chem. Geol. (2003)
% Spycher et al., Geochim. Cosmochim. Acta (2003)
% Spycher & Pruess, Geochim. Cosmochim. Acta (2005)
% Hassanzadeh et al., IJGGC (2008)
% ECLIPSE reference manual by Schlumberger
%
% INPUT:
% T:    cell indicating reference temperature (typically temperature at
%       reservoir depth, since it is where it matters most. Note that the 
%       black-oil formulation is isothermal) and units. Accepted C, K, F.
%       Structured as a {'key',value} pair. 
%       Example: T = {'C', 50};
%
% P:    cell indicating pressure values at which Rs and Bw should be
%       computed, units (Accepted Pa, MPa, bar, psi), and type.
%       type: 'mMn'  --> [min p, max p, num vals]
%             'vals' --> array of size (1, n) containing all pressure
%                        values of interest.
%       Structured as {'units','type', array}. 
%       Example: pval = {'bar','mMn',[1, 300, 24]};
%
% S:    cell indicating reference salinity and units (accepted 'ppm' 
%       (mg salt/kg solution), 'g/kg' (g salt /kg solution), 'sal' (ppt in 
%       g salt / kg solution), 'ppml' (mg salt / l solution), and 'm' 
%       (molality)).
%       Structured as {'units','type', value}. 'type' indicates the type of
%       fluid for ion and salt breakdown assignment, based on the input
%       total salinity value. Accepted options for 'type' are 'NaCl',
%       'seawater' and 'brine' for nonmolal 'unit' and 'NaCl' or 'CaCl2'
%       for 'm' unit.
%       Example: S = {'ppm', 'brine', 10^5};
%
% saltVar: boolean indicating whether the change in salt mole fraction due
%          to CO2 dissolution into the brine should be accounted for.
%          Example: saltVar = 0; (0 for "no" or 1 for "yes")
%          Note that, if yes, we consider the salts to be fully ionized
%          (see comments in NOTES section below).
%
% vapH2O: boolean indicating whether vaporization of water in the gas phase
%         is desided in the output data. This does not affect the
%         computation in the solubility model, where we always consider
%         vaporization, but it does affect the reporting of quantities
%         related to the gas phase.
%         Example: vapH20 = 0;  (O for "no", in which case the ECLIPSE
%                                keyword PVDG will be used; 1 for "yes", in
%                                which case the keyword PVTG will be used).
%
% unsatVals: boolean indicating whether unsaturated values (undersaturated
%            oil and gas) should be provided at pressures below the maximum
%            pressure.
%
% figs:   double indicating whether figures illustrating the results
%         should or should not be plotted.
%         Example: figs = 0; (0 for "no" or 1 for "yes")
%
% Optional argument:
% directory:  string indicating the full path to the directory where output
%             table file should be saved.
%             Example: 'C:\Users\lsalo\matlab\mrst-2019a\myutils\'
%
% OUTPUT:
% t: Structure specifying temperature, pressures, brine salinity and
%    aqueous and gase phases values of interest. These are (see below for
%    naming conventions):
%    t.aq:  m_co2, x_co2, X_co2, X_salt, rho_b, rho_aq, Rs, B and mu
%    t.gas: Z_co2, phi_co2, gamma_co2, y_h2o, Y_h2o, rho_g, Rv, B and mu
%
% If directory is passed:
% Saves a table in .txt format in the specified directory. Units are 
% ECLIPSE 'METRIC' units. The table is saved as:
% - Aqueous phase: a table called fileName of size (n, 4), where n is the 
%   number of pressure values of interest. Rs is the 1st column, p is the 
%   2nd column and the corresponding B_aq and viscosity are the 3rd and 4th
%   columns. fileName = PVTO_Tval_PminPMax_StypeValUnit.txt
% - Gaseous phase: a table called fileName of size (n, 3/4) with the same
%   format as above for the aqueous phase. In case Rv = 0 (dry gas),
%   fileName = PVDG_... and size is (n, 3); in case Rv \neq 0 (wet gas) 
%   fileName = PVTG_... and size is (n, 4).
%
% NOTES:
% - Two additional values of pressure above maximum input pressure are
%   added in order to provide two values for undersaturated aqueous
%   (gaseous) phase(s) for the highest Rs (Rv). This is done because it is
%   the standard ECLIPSE format for .DATA files. For these two values, the
%   gas compressibility factor and fugacity coefficients are not computed
%   since they are not needed, and that is why they will appear as NaN in
%   the output table file (NaNs must be deleted when adding this to your
%   .DATA input file).
%
% - Iteration: In principle, the mixing rules and fugacity coefficients in
%   the gas phase depend on the mixture composition. But, in turn, the
%   fugacity coefficients are needed to determine the mixture composition,
%   which requires an iterative procedure. As suggested by Spycher et al.
%   (2003, 2005) and Hassanzadeh et al. (2008), we can assume infinite
%   dilution of water in the gas phase (y_h2o = 0) when applying the mixing
%   rules, such that the fugacity coefficients can be directly obtained.
%   This is reasonable given the small amounts of water vaporized in the
%   gas phase in the P, T range of interest here.
%   Because of this assumption, the molar volume of the gas phase
%   is also calculated assuming pure CO2, which affects the gas density.
%   Therefore, our accuracy in gaseous phase properties will decrease as 
%   more water vaporizes into the gas phase.
%   This function, however, is prepared for an iterative procedure if that
%   becomes of interest later. In this case, the water and CO2 
%   Redlich-Kwong molecular interaction params involving water (a_h20,
%   b_h20 and a_h2o-co2) should be re-calculated or obtained from literature
%   (based on experimental data) without assuming infinite dilution (values
%   used here are those by Spycher et al., 2003, Table 1).
%
% - if saltVar = 0, we assume that the salts are not fully ionized when
%   calculating salt mole fractions, to match Hassanzadeh et al. (2008).
%   if saltVar = 1, we currently assume that the salts present in the 
%   aqueous phase are fully ionized, independently of the salt composition. 
%   Strictly, this is true for strong electrolytes like NaCl and CaCl2, 
%   and, to a lesser degree, for KCl and MgCl2. Other species like CaSO4 or 
%   MgSO4 are weak electrolytes, so this assumption is not correct; however, 
%   given that we only consider solutions of NaCl, seawater or typical 
%   brines, strong electrolytes are much more abundant.
%
% VARIABLE NAME CONVENTIONS:
% T            temperature
% _K _C        Kelvin, Celsius
% P, p         pressure
% S            salinity
% rho          mass density
% rhox         molar density
% mu           viscosity
% Rs, Rv       solution CO2-brine ratio and vaporized H2O-CO2 ratio
% B            formation volume factor
% Z            gas compressibility factor
% phi          fugacity coefficient
% gamma        activity coefficient
% Mw_k         molar mass of k
% m_k          molality of species k in solution (mol k / kg solvent)
% w_k          mass fraction of species k in aqueous phase
% v_k          mass fraction of species k in gaseous phase
% x_k          mole fraction of species k in water phase
% X_k          mole fraction of species k in aqueous phase
% y_k, Y_k     mole fraction of species k in gaseous (CO2 rich) phase
% _b, _brine   brine (H2O + salt(s))
% _gas         gaseous phase (CO2-rich)
% _aq          aqueous saline solution with CO2 (H2O + salt(s) + CO2)
% _s           ECLIPSE standard conditions of T, P (15.56C, 1atm)
% _1, _2       solute and solvent, respectively
% _pt          partial (for partial molar volume)
% nu           stoichiometric number of ions contained in the dissolved salt
% 
%                       ----------------------
directory = fullfile(getDatasetPath('co2labmit'), 'fluid_props');
%% Correlations & EoS
% Pick co2, brine, aqueous (brine + CO2) and gas (CO2 + H2O) models for 
% density and viscosity. 
s_unit  = S{1};  s_type = S{2};   s_val  = S{3};

% Density
if strcmp(s_type, 'NaCl') || s_val == 0
    rho_b_fcn = @(T, P, S) pvtBrineRoweChou1970(T, P, S);             % [kg/m^3]
else
    rho_b_fcn = @(T, P, S) pvtBrineBW1992(T, P, S);                   % [ " ]
end
rho_co2_fcn = @(T, P, varargin) pvtCO2RK1949(T, P, varargin);             % [ " ]
rho_aq_fcn  = @(Mw_1, Mw_2, x_1, x_2, V_pt, rho_1) ...                    % Garcia (2001)
              (1 + (Mw_1/Mw_2)*(x_1/x_2)) / ((V_pt/Mw_2)*(x_1/x_2) + (1/rho_1));
rho_gas_fcn = @(Mw_gas, V) Mw_gas/V;
  
    
% Viscosity
mu_co2_fcn  = @(T, rho) viscCO2F1998(T, rho);                             % [Pa*s]
mu_aq_fcn   = @(T, P, m_S, w_co2) viscBrineCO2IC2012(T, P, m_S, w_co2);   % [ " ]
mu_gas_fcn  = @(x, Mw, mu) viscGasMixtD1993(x, Mw, mu);                    


%% Check and organize inputs
% ------------------------- Temperature -----------------------------------
T_unit = T{1}; T_val  = T{2};
if strcmp(T_unit, 'K'),     T_K = T_val;
elseif strcmp(T_unit, 'C'), T_K = T_val + 273.15;
elseif strcmp(T_unit, 'F'), T_K = (T_val + 459.67)*(5/9);
else, error("Temperature units must be 'K', 'C' or 'F'.")
end
T_C   = T_K - 273.15;                                                     
T_K_s = 273.15 + 15.56;                                                   

% ---------------------------- Pressure -----------------------------------
p_unit = P{1}; p_val  = P{3}; p_typ  = P{2};
assert(max(strcmp(p_unit, {'Pa','MPa','bar','psi'})))
switch p_typ
    case 'mMn', assert(numel(p_val)==3)
        p = logspace(log10(p_val(1)), log10(p_val(2)), p_val(3))';
    case 'vals'
        p = p_val';
    otherwise, error("p{2} ('type') must be 'mMn' or 'vals'.")
end
if strcmp(p_unit, 'Pa'),        p = p*10^(-5);
elseif strcmp(p_unit, 'MPa'),   p = p*10;
elseif strcmp(p_unit, 'psi'),   p = p*0.06894757;
end
P_bar_s = 1.013529;                                                       
pmax    = max(p);                                  
if unsatVals
    assert(numel(p) > 5, " Provide at least 4 pressure values for unSatVals")
    nv = round((numel(p))/5); 
    id_unsat = unique([round(linspace(1, numel(p), nv)) numel(p)]);     % id at which unSatVals will be provided
    p_unsat = cell(1,numel(id_unsat));
    n_unsat = 0;
    for n=1:numel(id_unsat)-1
        p_unsat{n} = logspace(log10(p(id_unsat(n))), log10(pmax), 4);       % 4 values total, log spaced
        n_unsat = n_unsat + numel(p_unsat{n});
    end
    p_unsat{end} = [p(end); 1.1*pmax; 1.25*pmax];
    n_tot = numel(p) + n_unsat + numel(p_unsat{end}) - numel(id_unsat);
else
    id_unsat = numel(p);
    p_unsat{1} = [p(end); 1.1*pmax; 1.25*pmax];                                  % always unSatVals for highest Rs
end

% ---------------------------- Salinity -----------------------------------
% Check input 
if numel(S)~=3
    error("number of elements in S cell must be 3")
end

% Initial variables
maxits      = 8;                                                          % max N iterations
saltSpecies = {'NaCl', 'KCl', 'CaSO4', 'CaCl2', 'MgSO4', 'MgCl2'};
nu          = [  2,     2,     2,       3,     2,      3];                % [-]
ions        = {'Na_+', 'K_+', 'Ca_2+', 'Mg_2+', 'Cl_-', 'SO4_2-'};
%         [ NaCl,     KCl,    CaSO4,  CaCl2,  MgSO4,  MgCl2]
Mw_salt = [58.44277 74.5513 136.1406 110.984 120.3676 95.211]/10^3;       % [kg/mol]
%         [ Na(+),   K(+),  Ca(2+), Mg(2+), Cl(-), SO4(2-)]
Mw_io   = [22.98977 39.0983 40.078 24.305 35.453 96.0626]/10^3;           % ["]

% Compute salt molalities
if ~strcmp(S{1}, 'm')
    switch s_unit
        case {'ppm','ppml'}
            kg_salt = s_val/10^6;                                         % [kg salt / kg solution]
        case {'g/kg', 'sal'}
            kg_salt   = s_val/10^3;                                       % [kg salt / kg solution]                 
        otherwise 
            error("Salinity units must be 'ppm(l)', 'g/l', 'sal', 'm'")
    end
    kg_wat  = 1-kg_salt;                                                  % [kg] water / kg solution
    if strcmp(s_unit,'ppml')                                              % correct for volume
        dif = 1;    tols = 1e-4;    it = 1;
        while dif > tols && it < maxits                                   % find kg water
            rho_bs  = rho_b_fcn(T_K_s, P_bar_s, kg_salt/(kg_wat+kg_salt));
            kg1_w = rho_bs/10^3 - kg_salt;                                % kg h2o / l sol
            dif = abs(kg1_w - kg_wat);
            kg_wat = kg1_w;
            if dif < tols
                break
            elseif it < maxits
                it = it + 1;
                %disp(dif)
            else, error('Iteration to correct for volume did not converge')
            end
        end
    end
    switch s_type
        case 'NaCl'
            m_NaCl  = kg_salt/(kg_wat*Mw_salt(1));                        % [m] 
            m_salt  = [m_NaCl 0 0 0 0 0];                                 % [m]
        case 'seawater'
            % https://www.britannica.com/science/seawater for data
            w_io = [10.68 0.40 0.41 1.28 19.16 2.68]/34.61;
            kg_io = w_io*kg_salt;
            m_io = (kg_io./(Mw_io*kg_wat));                               % [m]
            m_salt = [m_io(1)  m_io(2)  m_io(3)/2  m_io(3)/2 ...
                      m_io(6)-(m_io(3)/2)  m_io(4)-m_io(6)-(m_io(3)/2)]; 
        case 'brine'
            % Data from Hyeong & Capuano (2001);
            w_io = [23.5 0.181 0.676 0.184 36.1 0]/10^3;                  % [kg io / l sol]
            w_io = w_io/sum(w_io);                                        % [-]                     
            kg_io = w_io*kg_salt;   
            m_io = (kg_io./(Mw_io*kg_wat));                               % [m]
            m_salt = [m_io(5)-m_io(2)-2*(sum(m_io(3:4)))  m_io(2)  0  m_io(3) 0  m_io(4)];
        otherwise, error('This salt input is not supported.')
    end
else
    switch s_type
        case 'NaCl'
        m_salt  = [S{3}, 0, 0, 0, 0, 0];                                  % [m]
        kg_salt = S{3}*Mw_salt(1);
        case 'CaCl2'
        m_salt  = [0, 0, 0, S{3}, 0, 0];
        kg_salt = S{3}*Mw_salt(4);
        otherwise, error('This salt is not supported for molal inputs')
    end
    kg_wat = 1;
end

% ------------------- Ion molalities and salt mass fraction -------------------
% The total number of moles of ions in a solution (for fully ionized salts)
% will be a multiplier of the number of moles of a given salt species. This
% multiplier is given by the stoichiometric number (nu). Therefore, when
% using the molality of a fully dissociated salt, one can either
% use the molality of the salt species times the stoichiometric number, or
% the molality of the ions.                    
w_salt = kg_salt/(kg_wat+kg_salt);                                        % [-] 
if ~strcmp(s_type,'seawater') && ~strcmp(s_type,'brine')
    m_io   = [m_salt(1) m_salt(2) m_salt(3)+m_salt(4) m_salt(5)+m_salt(6) ...
              sum(m_salt(1:2))+2*(m_salt(4)+m_salt(6)) m_salt(3)+m_salt(5)];
end
   
% Check validity range
if pmax > 600 && pmax < 700
    warning('Input pressure contains values above 600 bar. Model may not be accurate.')
elseif pmax > 700
    error('Input pressure contains values significantly above 600 bar. Formulation possibly not accurate.')
end
if T_C > 100 && T_C <= 125
    warning('Input temperature above 100 C. Model may not be accurate.')
elseif T_C > 125
    error('Input temperature significantly above 100 C. This formulation is not accurate.')
end
if max(m_salt) > 4 && max(m_salt) < 6
    warning('Input salinity above 4 m. Model may not be accurate.')
elseif max(m_salt) > 6
    error('Input salinity above halite saturation. This formulation is not accurate.')
end

%% Preliminary variables 
R       = 83.1447;                                          % bar*cm^3/mol*K
Mw_h2o  = 18.01528/10^3;                                    % [kg/mol]
Mw_co2  = 44.0095/10^3;                                     %   "

% Compute mole fractions
m_h2o   = 1/Mw_h2o;                                         % mol H20/ kg H2O
vm = nu.*m_salt; 
if saltVar
    x_salt  = (vm) ./ (sum(vm) + m_h2o);    % [-] mole fractions (total), fully ionized salt
    x_salts = (vm) ./ sum(vm);              % [-] (of salts)
else
    x_salt  = m_salt ./ (sum(m_salt) + m_h2o);              % [-] Not fully ionized (Hassanzadeh et al., 2008)
    x_salts = m_salt./sum(m_salt);                          % [-] 
end
if all(isnan(x_salts))
    assert(S{3}==0)
    x_salts(:) = 0;
end
x_h2o   = 1 - sum(x_salt);                                  % [-]
Mw_salts = sum(Mw_salt.*x_salts);                           % [kg/mol] avg

% Compute mixture Mw (Mw_brine) in kg/mol
Mw_b = sum(Mw_salt.*x_salt) + Mw_h2o*x_h2o;                 % [kg/mol]

% Water and CO2 Redlich-Kwong molecular interaction params (Spycher et al., 2003) 
% The parameters involving H20 were calculated assuming infinite dilution
% of water in the gas phase, so they should be recalculated if this
% assumption is relaxed with this solubility model. In this case, a_h20
% would also be needed (not reported by Spycher et al., 2003, due to this
% assumption).
a_co2       = 7.54*10^7 - 4.13*10^4*T_K;                    % bar*cm^6*K^0.5/mol^2
a_h2o_co2   = 7.89*10^7;                                    % "
b_co2       = 27.8;                                         % cm^3/mol
b_h2o       = 18.18;                                        % "
a_h2o       = 0;                                            % not used here

% Average partial molar volumes of each pure condensed component, they are
% assumed constant in the p interval of interest (Spycher et al., 2003)
Vp_co2      = 32.6;                                         % cm^3/mol (both (g) and (l))
Vp_h2o      = 18.1;                                         % "

% True equilibrium ctnts. at p0=1bar (/kappa params; Spycher et al., 2003)
p0          = 1;                                            % bar
kap0_co2_g  = 10^(1.189 + 1.304*10^(-2)*T_C - 5.446*10^(-5)*T_C^2);
kap0_co2_l  = 10^(1.169 + 1.368*10^(-2)*T_C - 5.380*10^(-5)*T_C^2);
kap0_h2o    = 10^(-2.209 + 3.097*10^(-2)*T_C - 1.098*10^(-4)*T_C^2 ...
                  + 2.048*10^(-7)*T_C^3);

% Partial molar volume of CO2 (Garcia, 2001)
V_phi     = 37.51 - 9.585*10^(-2)*T_C + 8.740*10^(-4)*T_C^2 ...
            - 5.044*10^(-7)*T_C^3;                          % [cm^3/mol]
V_phi     = V_phi/10^6;                                     % [m^3/mol]

% CO2 and brine molar densities at standard conditions
[~, rhox_co2_s, rho_co2_s] = rho_co2_fcn(T_K_s, P_bar_s);   % [mol/m^3]
% With CoolProp library (EoS from Span & Wagner, 1996)
% rhox_co2_s   = py.CoolProp.CoolProp.PropsSI('DMOLAR',...
%                'T', T_K_s, 'P', P_Pa_s, 'CO2');           % [mol/m^3]
rho_brine_s   = rho_b_fcn(T_K_s, P_bar_s, w_salt);
rhox_brine_s  = rho_brine_s/Mw_b;                           % [mol/m^3]
rho_h2o_s     = pvtBrineRoweChou1970(T_K_s, P_bar_s, 0);    
rhox_h2o_s    = rho_h2o_s/Mw_h2o;                           % ["]


%% Calculate
% Initialize mole fractions in gaseous phase assuming infinite dilution of
% water
Y0_h2o      = 0;
Y_h2o_out   = 0;
Y0_co2      = 1;
Y_co2_out   = 1;

% Initialize result table
hdrs.gas    = {'Z_co2', 'phi_co2', 'gamma_co2', 'y_h2o', 'Y_h2o', ...
               'rho_co2_kgm3', 'rho_g_kgm3', 'Rv_Sm3Sm3', 'B_g_m3Sm3', ...
               'mu_co2_cP', 'mu_g_cP'};
hdrs.aq     = {'m_co2', 'x_co2', 'X_co2', 'X_salt', 'rho_b_kgm3', ...
               'rho_aq_kgm3', 'Rs_Sm3Sm3', 'B_b_m3Sm3', 'B_aq_m3Sm3', ...
               'mu_b_cP', 'mu_aq_cP'};
vartyp.gas  = cellstr(repmat('double',numel(hdrs.gas),1));
vartyp.aq   = cellstr(repmat('double',numel(hdrs.aq),1));
t.gas       = table('Size',[numel(p), numel(hdrs.gas)], ...
                    'VariableTypes', vartyp.gas, 'VariableNames', hdrs.gas);
t.aq        = table('Size',[numel(p), numel(hdrs.aq)], ...
                    'VariableTypes', vartyp.aq, 'VariableNames', hdrs.aq);
t.P_bar     = unique([p; p_unsat{end}]); 
t.T_Kelvin  = T_K; 
t.S_molal.salts = saltSpecies; t.S_molal.saltVals = m_salt;
t.S_molal.ions = ions;         t.S_molal.ionVals  = m_io;

if unsatVals
    hdrs2.gas = {'P_bar', 'Rv_Sm3Sm3', 'B_g_m3Sm3', 'mu_g_cP'};
    hdrs2.aq = {'Rs_Sm3Sm3', 'P_bar', 'B_aq_m3Sm3', 'mu_aq_cP'};
    vartyp.gas  = cellstr(repmat('double',numel(hdrs2.gas),1));
    vartyp.aq   = cellstr(repmat('double',numel(hdrs2.aq),1));
    t2.gas      = table('Size',[n_tot, numel(hdrs2.gas)], ...
                        'VariableTypes', vartyp.gas, 'VariableNames', hdrs2.gas);
    t2.aq       = table('Size',[n_tot, numel(hdrs2.aq)], ...
                    'VariableTypes', vartyp.aq, 'VariableNames', hdrs2.aq);
else
    t2 = [];
end

% Compute
maxits = 15;
iterate = 0;
idt2 = 1;
idp = 0;
for n=1:numel(p)
    it = 1;
    P  = p(n);
    while it <= maxits
        % Mixing rules (Prausnitz et al., 1986)
        a_m = Y0_h2o^2*a_h2o + 2*Y0_co2*Y0_h2o*a_h2o_co2 + ...
            Y0_co2^2*a_co2;
        b_m = Y0_co2*b_co2 + Y0_h2o*b_h2o;

        % 1. Gas molar volume (Redlich and Kwong (1949) EoS, = RK EoS)
        [V, ~, rho_co2] = rho_co2_fcn(T_K, P, a_m, b_m);              % [cm^3/mol, ~, kg/m^3]

        % 2. Gas compressibility factor
        Z = P*V/(R*T_K);                                              % [-]

        % 3. Fugacity coefficients (Spycher et al., 2003)
        %   3.1 CO2
        fa = log(V/(V-b_m));
        fb = b_co2/(V-b_m);
        fc = (2*(Y0_co2*a_co2 + Y0_h2o*a_h2o_co2)/(R*T_K^1.5*b_m)) * log((V+b_m)/V);
        fd = a_m*b_co2/(R*T_K^1.5*b_m^2) * (log((V+b_m)/V) - b_m/(V+b_m));
        fe = log(Z);
        phi_co2 = exp(fa + fb - fc + fd - fe);                        % [-]

        %   3.2 Water
        fbw = b_h2o/(V-b_m);
        fcw = (2*(Y0_co2*a_h2o_co2 + Y0_h2o*a_h2o)/(R*T_K^1.5*b_m)) * log((V+b_m)/V);
        fdw = a_m*b_h2o/(R*T_K^1.5*b_m^2) * ...
            (log((V+b_m)/V) - b_m/(V+b_m));                         % typo in Hassanzadeh et al., 2008
        phi_h2o = exp(fa + fbw - fcw + fdw - fe);                     % [-]

        % 4. CO2 molality in pure H2O at P, T conditions (m0_co2)
        if T_C < 31 && V < 94                                         % Spycher et al., 2003
            kap0_co2 = kap0_co2_l;
        else
            kap0_co2 = kap0_co2_g;
        end
        B = phi_co2*P/(m_h2o*kap0_co2) * exp(-(P-p0)*Vp_co2/(R*T_K));
        A = kap0_h2o/(phi_h2o*P) * exp((P-p0)*Vp_h2o/(R*T_K));
        % y_i = mole fraction of component i in the gaseous phase
        y_h2o = (1-B)/(1/A-B);                                        % [-] Spycher et al., 2003
        % x_i = mole fraction of component i in the aqueous phase
        x_co2 = B*(1-y_h2o);                                          % [-]
        % Check values within range
        if x_co2 < 0, warning(['x_co2 = ' num2str(x_co2) '. Set to 0'])
            x_co2 = 0;
        elseif x_co2 > 1, warning(['x_co2 = ' num2str(x_co2) '. Set to 1'])
            x_co2 = 1;
        end
        if y_h2o < 0, warning(['y_h2o = ' num2str(y_h2o) '. Set to 0'])
            y_h2o = 0;
        elseif y_h2o > 1, warning(['y_h2o = ' num2str(y_h2o) '. Set to 1'])
            y_h2o = 1;
        end
        x_h2o = 1 - x_co2;                                            % [-]
        % molality
        m0_co2 = m_h2o*x_co2/x_h2o;                                   % [m]

        % 5. CO2 molality in saline solution at P, T conditions (m_co2)
        %   5.1 CO2 activity coefficient (Duan and Sun, 2003, as in Spycher
        %   & Pruess, 2005; typos in Hassanzadeh et al., 2008)
        gamma_co2 = activityCO2DS2003(T_K, P, m_io);

        %   5.2 CO2 molality
        m_co2 = m0_co2/gamma_co2;

        % 6. H2O and CO2 mole fractions and equilibrium ratios in aqueous
        % saline solution with dissolved CO2.
        X_co2 = m_co2/(m_co2 + m_h2o + sum(vm));              % [-] mole fraction
        if saltVar
            X_salt = sum(vm)/(m_co2 + m_h2o + sum(vm));       % [-] fully ionized
        else
            X_salt = sum(x_salt);                                     % [-] Not fully ionized
        end
        X_solv    = 1 - X_co2;                                        % [-]
        X_h2o     = X_solv - X_salt;                                  % [-] fully ionized
        Y_h2o     = A*X_h2o;                                          % [-] fully ionized
        Y_co2     = 1-Y_h2o;                                          % [-]
        %K_co2     = Y_co2/X_co2;                                      % [-] equilibrium rat.
        %K_h2o     = Y_h2o/X_h2o;                                      % [-]

        if iterate == 1
            res = abs(Y_h2o - Y0_h2o);
            Y0_h2o = Y_h2o;
            Y0_co2 = Y_co2;
            if  res < 1e-6
                disp(['Iterations at p=' num2str(P) 'bar: ' num2str(it)]);
                break
            elseif it == maxits
                error('Iterative loop did not converge')
            end
            it = it + 1;
        else
            break
        end
    end

    if vapH2O
        Y_h2o_out = Y_h2o;
        Y_co2_out = Y_co2;
    end

    % 7. Compute Rs and Rv
    Rs = rhox_brine_s*X_co2/(rhox_co2_s*(1-X_co2));                   % [V CO2 Sm^3/ V br Sm^3] at sc.
    Rv = rhox_co2_s*Y_h2o_out/(rhox_h2o_s*Y_co2_out);                 % [V H2O Sm^3/ V CO2 Sm^3]
    
    pnum = 1;
    if ismember(P, p(id_unsat))
        idp = idp + 1;
        pnum = numel(p_unsat{idp});
    end
    for k=1:pnum
        if k > 1
            P = p_unsat{idp}(k);
            [V, ~, rho_co2] = rho_co2_fcn(T_K, P, a_m, b_m);                  % [cm^3/mol, ~, kg/m^3]
        end
        
        % 8. Molar masses and mass fractions
        Mw_aq   = Mw_salts*X_salt + Mw_h2o*X_h2o + Mw_co2*X_co2;              % [kg/mol] avg Molar mass
        Mw_solv = (Mw_salts*X_salt + Mw_h2o*X_h2o)/X_solv;
        Mw_gas  = Mw_co2*Y_co2_out + Mw_h2o*Y_h2o_out;
        %if saltVar == 1; w_salt = X_salt*(Mw_salts/Mw_aq); end
        w_co2 = X_co2*(Mw_co2/Mw_aq);                                         % [-] mass frac.
        v_h2o = Y_h2o_out*(Mw_h2o/Mw_gas);
           
        % 9. Aqueous phase (CO2 saturated; with dissolved salts) density
        %    (Garcia, 2001) and gaseous phase density
        rho_brine = rho_b_fcn(T_K, P, w_salt);                                % density of the brine (no CO2) [kg/m^3, 1/kPa]
        rho_aq    = rho_aq_fcn(Mw_co2, Mw_solv, X_co2, X_solv, V_phi, rho_brine);
        % We use the molar volume of the gas calculated above, which does not
        % yet account for the presence of water in the gas phase (mixing rules
        % assume Y_h2o = 0)
        rho_gas = rho_gas_fcn(Mw_gas, V/10^6); 
        
        % 10. FVF aqueous saline solution with CO2 (= B_aq) and gaseous
        % (CO2-rich) solution with H2O
        B_b   = rho_brine_s/rho_brine; 
        B_aq  = rho_brine_s/(rho_aq*(1 - w_co2));                             % [m^3/Sm^3]
        B_gas = rho_co2_s/(rho_co2*(1 - v_h2o));                              % [ " ] we use rho_co2 since diff is negligible
        
        % 11. Viscosity of aqueous and gaseous solutions
        mu_b   = mu_aq_fcn(T_K, P, sum(m_salt), 0)*1e3;                       % [cP]
        mu_aq  = mu_aq_fcn(T_K, P, sum(m_salt), w_co2)*1e3;                   % [ cP ]
        mu     = [mu_aq_fcn(T_K, P, 0, 0), mu_co2_fcn(T_K, rho_co2)];         % [ Pa*s ]
        mu_gas = mu_gas_fcn([Y_h2o_out, Y_co2_out], [Mw_h2o, Mw_co2], mu)*1e3; % [ cP ]

        if P < pmax && k == 1 || P >= pmax
            idt = n;
            if P>pmax, idt = n + (k-1); end
            % Assign variables to fields in result table
            [t.gas.Z_co2(idt), t.gas.phi_co2(idt), t.gas.gamma_co2(idt), ...
                t.gas.y_h2o(idt), t.gas.Y_h2o(idt), t.gas.rho_co2_kgm3(idt), t.gas.rho_g_kgm3(idt), ...
                t.gas.Rv_Sm3Sm3(idt), t.gas.B_g_m3Sm3(idt), t.gas.mu_co2_cP(idt), t.gas.mu_g_cP(idt)]  ...
                = deal(Z, phi_co2, gamma_co2, y_h2o, Y_h2o_out, rho_co2, rho_gas, Rv, B_gas, mu(2)*1e3, mu_gas);

            [t.aq.m_co2(idt), t.aq.x_co2(idt), t.aq.X_co2(idt), t.aq.X_salt(idt), ...
                t.aq.rho_b_kgm3(idt), t.aq.rho_aq_kgm3(idt), t.aq.Rs_Sm3Sm3(idt), ...
                t.aq.B_b_m3Sm3(idt), t.aq.B_aq_m3Sm3(idt), t.aq.mu_b_cP(idt), t.aq.mu_aq_cP(idt)]  ...
                = deal(m_co2, x_co2, X_co2, X_salt, rho_brine, rho_aq, Rs, B_b, B_aq, mu_b, mu_aq);
            
        end
        if unsatVals
            [t2.gas.P_bar(idt2), t2.gas.Rv_Sm3Sm3(idt2), t2.gas.B_g_m3Sm3(idt2), t2.gas.mu_g_cP(idt2)] = deal(P, Rv, B_gas, mu_gas);
            [t2.aq.Rs_Sm3Sm3(idt2), t2.aq.P_bar(idt2), t2.aq.B_aq_m3Sm3(idt2), t2.aq.mu_aq_cP(idt2)] = deal(Rs, P, B_aq, mu_aq);
            idt2 = idt2 + 1;
        end
    end
    
end


%% Plots
if figs == 1
    col_co2  = [180, 0, 0; ...
                100, 100, 100; ...
                0, 99, 67]./255;
    col_brine  = [0, 71, 148]./255;
    %axisarg    = {'FontSize', 11, 'TickLabelInterpreter','latex'};
    latx       = {'Interpreter','latex'};
    fontsz_tit = {'fontsize',14};
    
    h1 = figure(1); % Compare with Fig. 8 in Spycher et al. (2003)
    yyaxis left
    plot(t.P_bar,t.gas.Z_co2, 'color', col_co2(1,:),'linewidth', 1); grid on; ylim([0 1]);
    ylabel('$Z$ [-]','fontsize', 14, latx{:})
    yyaxis right
    plot(t.P_bar,t.gas.phi_co2, 'color', col_co2(1,:), 'LineStyle', '--', 'linewidth', 1); grid on;
    ylim([0 1]);
    ylabel('$\phi$ [-]','fontsize', 14, latx{:})
    xlabel('$p$ [bar]','fontsize', 14, latx{:})
    ax = gca;
    ax.YAxis(1).TickValues = 0:.1:1;
    ax.YAxis(2).TickValues = 0:.1:1;
    ax.YAxis(1).Color = col_co2(1,:); ax.YAxis(1).FontSize = 12;
    ax.YAxis(2).Color = col_co2(1,:); ax.YAxis(2).FontSize = 12;
    ax.XAxis.FontSize = 12;
    ax.XTick = [0 100 200 300 400 500 600];
    xlim([0 pmax]);
    legend('$Z$', '$\phi$', latx{:}, 'location', 'southwest', 'box', 'off', ...
           'fontsize', 14)
    title(['CO$_2$ coefs., $T$=' num2str(T_C) '$^\circ$C'], latx{:}, ...
        fontsz_tit{:})
    set(h1, 'Position', [600, 600, 325, 325])
    
    h2 = figure(2); % Compare with Fig. 2 in Spycher and Pruess (2005)
    yyaxis left
    plot(t.P_bar,t.aq.m_co2,'color', col_co2(1,:)); grid on; ylim([0 2]);
    ylabel('CO$_{2(aq)}$ [m]','fontsize', 14, latx{:})
    yyaxis right
    plot(t.P_bar,t.gas.Y_h2o*1000,'color', col_brine, 'LineStyle', '--'); grid on; ylim([4 14]);
    ylabel('$\Psi_{\mathrm{H}_2\mathrm{O}}\times10^3$ [-]','fontsize', 14, ...
        latx{:})
    xlabel('$p$ [bar]','fontsize', 14, latx{:})
    ax = gca;
    ax.YAxis(1).TickValues = [0, 0.5, 1.0, 1.5, 2.0];
    ax.YAxis(2).TickValues = 4:14;
    ax.YAxis(1).Color = col_co2(1,:); ax.YAxis(1).FontSize = 12;
    ax.YAxis(2).Color = col_brine; ax.YAxis(2).FontSize = 12;
    ax.XAxis.FontSize = 12;
    ax.XTick = [0 100 200 300 400 500 600];
    xlim([0 pmax]);
    legend('CO$_{2\mathrm{(aq)}}$', '$\Psi_{\mathrm{H}_2\mathrm{O}}$', ...
        latx{:}, 'location', 'northeast', 'box', 'off')
    title(['CO$_{2\mathrm{(aq)}}$ and H$_2$O$_{\mathrm{(g)}}$, $T$=' ...
        num2str(T_C) '$^\circ$C'], latx{:},  fontsz_tit{:})
    set(h2, 'Position', [600, 600, 300, 260])
    
    h3 = figure(3); % Compare with Fig. 4 in Hassanzadeh et al. (2008)
    yyaxis left
    plot(t.P_bar/10,t.aq.x_co2,'color', col_co2(1,:)); grid on; ylim([0.005 0.03]);
    ylabel('$\chi_{\mathrm{CO}_2}$ [-]','fontsize', 14, latx{:})
    yyaxis right
    plot(t.P_bar/10,t.gas.y_h2o,'color', col_brine, 'LineStyle', '--'); grid on; ylim([0 0.03]);
    ylabel('$\psi_{\mathrm{H}_2\mathrm{O}}\times10^3$ [-]','fontsize', 14, ...
        latx{:})
    xlabel('$p$ [MPa]','fontsize', 14, latx{:})
    ax = gca;
    ax.YAxis(1).TickValues = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03];
    ax.YAxis(2).TickValues = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03];
    ax.YAxis(1).Color = col_co2(1,:); ax.YAxis(1).FontSize = 12;
    ax.YAxis(2).Color = col_brine; ax.YAxis(2).FontSize = 12;
    ax.XAxis.FontSize = 12;
    ax.XTick = [0 10 20 30 40 50 60 70 80];
    xlim([0 pmax/10]);
    legend('$\chi_{\mathrm{CO}_2}$', '$\psi_{\mathrm{H}_2\mathrm{O}}$', ...
        latx{:}, 'location', 'northeast', 'box', 'off')
    title(['CO$_{2\mathrm{(aq)}}$ and H$_2$O$_{\mathrm{(g)}}$, $T$=' ...
        num2str(T_C) '$^\circ$C'], latx{:},  fontsz_tit{:})
    set(h3, 'Position', [600, 600, 300, 260])
    
    h4 = figure(4); % Compare with Fig. 5 in Hassanzadeh et al. (2008). Typo in
    % their caption, T should be 45C.
    yyaxis left
    plot(t.P_bar/10, t.aq.Rs_Sm3Sm3, 'color', col_co2(1,:),'linewidth', 1);
    %ylim([0 20])
    ylabel('R$_s$ [Sm$^3$ CO$_2$/Sm$^3$ brine]','fontsize', 14, latx{:})
    %set(gca,'Ytick',[0 5 10 15 20])
    ax = gca;
    ax.YAxis(1).Color = col_co2(1,:); ax.YAxis(1).FontSize = 12;
    ytl = get(gca, 'YTick');
    yyaxis right
    plot(t.P_bar/10, t.aq.B_aq_m3Sm3, 'color', col_brine, 'LineStyle', '--','linewidth',1); grid on;
    %ylim([1.012 1.032]); 
    %set(gca,'Ytick',[1.012 1.016 1.020 1.024 1.028 1.032])
    ylabel('B$_b$ [m$^3$/Sm$^3$]','fontsize', 14, latx{:})
    ax = gca;
    ax.YAxis(2).Color = col_brine; ax.YAxis(2).FontSize = 12;
    ytr = get(gca, 'YTick');
    ytrv = linspace(min(ytr), max(ytr), numel(ytl));
    ytrc = compose('%.3f', ytrv);
    set(gca, 'YTick', ytrv, 'YTickLabel', ytrc)
    
    xlabel('$p$ [MPa]','fontsize', 14, latx{:})
    xlim([0 50]);
    ax.XAxis.FontSize = 12;
    ax.XTick = [0 10 20 30 40 50];
    grid on; 
    legend('R$_s$', 'B$_b$', latx{:}, 'location', 'southeast', ...
           'box', 'off', 'fontsize', 14)
    title(['$T$=' num2str(T_C) 'C, $S$=' num2str(S{end}) S{1}], ...
        latx{:},  fontsz_tit{:})
    set(h4, 'Position', [600, 600, 400, 325])
    
    f5 = figure(5); % compare with Fig. 5b, d in Hassanzadeh et al. (2008)
    % subplot 1, density
    subplot(1,2,1)
    hold on
    p1 = plot(t.P_bar, t.aq.rho_b_kgm3, '--', 'color', col_brine, ...
        'linewidth', 1, 'DisplayName', '$\rho_\mathrm{b}$');
    p2 = plot(t.P_bar, t.aq.rho_aq_kgm3, '-', 'color', col_brine, ...
        'linewidth', 1, 'DisplayName', '$\rho_\mathrm{b}$, CO$_2$ sat.');
    p3 = plot(t.P_bar, t.gas.rho_g_kgm3, '-', 'color', col_co2(1,:), ...
        'linewidth', 1, 'DisplayName', '$\rho_\mathrm{g}$');
    hold off
    grid on
    xlabel('$p$ [bar]', latx{:}, 'fontsize', 14)
    ylabel('$\rho_\alpha$ [kg/m$^3$]', latx{:}, 'fontsize', 14)
    xlim([0 t.P_bar(end-2)])
    ylim([0 1200])
    legend([p1 p2, p3], latx{:}, 'fontsize', 14, 'location', 'southeast')
    % subplot 2, viscosity
    subplot(1,2,2)
    hold on
    p1 = plot(t.P_bar, t.aq.mu_aq_cP, '-', 'color', col_brine, ...
        'linewidth', 1, 'DisplayName', '$\mu_\mathrm{b}$, CO$_2$ sat.');
    p2 = plot(t.P_bar, t.gas.mu_g_cP, '-', 'color', col_co2(1,:), ...
        'linewidth', 1, 'DisplayName', '$\mu_\mathrm{g}$');
    hold off
    grid on
    xlabel('$p$ [bar]', latx{:}, 'fontsize', 14)
    ylabel('$\mu_\alpha$ [cP]', latx{:}, 'fontsize', 14)
    legend([p1 p2], latx{:}, 'fontsize', 14, 'location', 'southeast')
    set(gca,'yscale','log')
    ylim([1e-2 2])
    xlim([0 t.P_bar(end-2)])
    set(f5, 'Position', [200, 200, 600, 325])
end


%% Save PVT table
% Write to file
if nargin > 7  
    % Gas phase
    if max(t.gas.Rv_Sm3Sm3) > 0
        fileName_gas = ['PVTG_T' num2str(T_K) '_P' num2str(p(1)) 'P' ...
                        num2str(pmax) '_S' s_type  num2str(s_val) s_unit];
        if unsatVals
            writetable(t2.gas, [directory fileName_gas '.txt'], 'Delimiter', 'tab');
        else
            tgas = table(t.P_bar, t.gas.Rv_Sm3Sm3, t.gas.B_g_m3Sm3, t.gas.mu_g_cP);
            tgas.Properties.VariableNames = {'P_bar', 'Rv_Sm3Sm3', 'B_g_m3Sm3', 'mu_g_cP'};
            writetable(tgas, [directory fileName_gas '.txt'], 'Delimiter', 'tab');
        end
        disp('-----------------IMPORTANT NOTE (PVTG, Gas) ----------------------------')
        disp('The last 2 rows in PVTG must be substituted by a single row.')
        disp('This single row should be for last (max) p value used as input and should be:')
        disp('Rv = 0, and corresponding Bgas and viscosity.')
        disp('Refer to ECLIPSE Reference Manual for Keyword format.')
        disp('To obtain, run this function with vapH2O = 0 (PVDG) and get corresponding row.')
        disp('-----------------------------------------------------------')
    else
        fileName_gas = ['PVDG_T' num2str(T_K) '_P' num2str(p(1)) 'P' ...
                        num2str(pmax) '_S' s_type  num2str(s_val) s_unit];
        tgas = table(t.P_bar, t.gas.B_g_m3Sm3, t.gas.mu_g_cP);
        tgas.Properties.VariableNames = {'P_bar', 'B_g_m3Sm3', 'mu_g_cP'};
        writetable(tgas, [directory fileName_gas '.txt'], 'Delimiter', 'tab');
    end
    
    % Oleic phase (brine)
    fileName_aq  = ['PVTO_T' num2str(T_K) '_P' num2str(p(1)) 'P' ...
                    num2str(pmax) '_S' s_type  num2str(s_val) s_unit];
    if unsatVals
        writetable(t2.aq, [directory fileName_aq '.txt'], 'Delimiter', 'tab');
    else
        taq = table(t.aq.Rs_Sm3Sm3, t.P_bar, t.aq.B_aq_m3Sm3, t.aq.mu_aq_cP);
        taq.Properties.VariableNames = {'Rs_Sm3Sm3' 'P_bar' 'B_aq_m3Sm3' 'mu_b_cP'};
        writetable(taq, [directory fileName_aq '.txt'], 'Delimiter', 'tab');
    end

end