function parameter_struct = LET_imbibition_pc_handle(model, x, row_names)
%
% DESCRIPTION: builds the parameter input for the saturation function
%              builder based on the excel table given in the input for the
%              history matching
%
% SYNOPSIS:
%   parameter_struct = LET_imbibition_pc_handle(model, x, row_names)
%
% PARAMETERS:
%   - model: struct used for the history matching 
%   - x: the parameters used for the history match iterations
%   - row_names: the name corresponding to the parameters in x read from
%   the input excel
%
% RETURNS:
%   parameter_struct: struct used to build the saturation function for
%   history matching iterations
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
pc_input_struct = model.experiment.satfun.pc;
kr_input_struct = model.experiment.satfun.kr;

% handle swc
swc_pc_mask = strcmpi(row_names, 'swc_pc');
swc_kr_mask = strcmpi(row_names, 'swc_kr');
if any(swc_pc_mask) && any(swc_kr_mask)
    parameter_struct.swc_pc = x(swc_pc_mask);
    parameter_struct.swc_kr = x(swc_kr_mask);
elseif any(swc_pc_mask) && not(any(swc_kr_mask))
    parameter_struct.swc_pc = x(swc_pc_mask);
    parameter_struct.swc_kr = x(swc_pc_mask);
elseif not(any(swc_pc_mask)) && any(swc_kr_mask)
    parameter_struct.swc_pc = x(swc_kr_mask);
    parameter_struct.swc_kr = x(swc_kr_mask);
elseif not(any(swc_pc_mask)) && not(any(swc_kr_mask))
    if isfield(pc_input_struct, 'Swc')
        parameter_struct.swc_pc = pc_input_struct.Swc;
    else
        parameter_struct.swc_pc = 0;
    end
    if isfield(kr_input_struct, 'Swc')
        parameter_struct.swc_kr = kr_input_struct.Swc;
    else
        parameter_struct.swc_kr = 0;
    end
end

% handle sor
sor_pc_mask = strcmpi(row_names, 'sor_pc');
sor_kr_mask = strcmpi(row_names, 'sor_kr');
if any(sor_pc_mask) && any(sor_kr_mask)
    parameter_struct.sor_pc = x(sor_pc_mask);
    parameter_struct.sor_kr = x(sor_kr_mask);
elseif any(sor_pc_mask) && not(any(sor_kr_mask))
    parameter_struct.sor_pc = x(sor_pc_mask);
    parameter_struct.sor_kr = x(sor_pc_mask);
elseif not(any(sor_pc_mask)) && any(sor_kr_mask)
    parameter_struct.sor_pc = x(sor_kr_mask);
    parameter_struct.sor_kr = x(sor_kr_mask);
elseif not(any(sor_pc_mask)) && not(any(sor_kr_mask))
    if isfield(pc_input_struct, 'Sor')
        parameter_struct.sor_pc = pc_input_struct.Sor;
    else
        parameter_struct.sor_pc = 0;
    end
    if isfield(kr_input_struct, 'Sor')
        parameter_struct.sor_kr = kr_input_struct.Sor;
    else
        parameter_struct.sor_kr = 0;
    end
end

%handle spontaneous_multiplier
spontaneous_multiplier_mask = strcmpi(row_names, 'spontaneous_multiplier');
if any(spontaneous_multiplier_mask)
    parameter_struct.spontaneous_multiplier = x(spontaneous_multiplier_mask);
else
    if isfield(pc_input_struct, 'spontaneous_multiplier')
        parameter_struct.spontaneous_multiplier = pc_input_struct.spontaneous_multiplier;
    else
        parameter_struct.spontaneous_multiplier = 1;
    end
end

%handle forced_multiplier
forced_multiplier_mask = strcmpi(row_names, 'forced_multiplier');
if any(forced_multiplier_mask)
    parameter_struct.forced_multiplier = x(forced_multiplier_mask);
else
    if isfield(pc_input_struct, 'forced_multiplier')
        parameter_struct.forced_multiplier = pc_input_struct.forced_multiplier;
    else
        parameter_struct.forced_multiplier = 1;
    end
end

%handle min_pc
min_pc_mask = strcmpi(row_names, 'min_pc');
if any(min_pc_mask)
    parameter_struct.min_pc = x(min_pc_mask);
else
    if isfield(pc_input_struct, 'min_pc')
        parameter_struct.min_pc = pc_input_struct.min_pc;
    else
        parameter_struct.min_pc = -2;
    end
end

%handle max_pc
max_pc_mask = strcmpi(row_names, 'max_pc');
if any(max_pc_mask)
    parameter_struct.max_pc = x(max_pc_mask);
else
    if isfield(pc_input_struct, 'max_pc')
        parameter_struct.max_pc = pc_input_struct.max_pc;
    else
        parameter_struct.max_pc = 2;
    end
end

%handle L_spont
L_spont_mask = strcmpi(row_names, 'L_spont');
if any(L_spont_mask)
    parameter_struct.L_spont = x(L_spont_mask);
else
    if isfield(pc_input_struct, 'L_spont')
        parameter_struct.L_spont = pc_input_struct.L_spont;
    else
        parameter_struct.L_spont = 1;
    end
end

%handle E_spont
E_spont_mask = strcmpi(row_names, 'E_spont');
if any(E_spont_mask)
    parameter_struct.E_spont = x(E_spont_mask);
else
    if isfield(pc_input_struct, 'E_spont')
        parameter_struct.E_spont = pc_input_struct.E_spont;
    else
        parameter_struct.E_spont = 60;
    end
end

%handle T_spont
T_spont_mask = strcmpi(row_names, 'T_spont');
if any(T_spont_mask)
    parameter_struct.T_spont = x(T_spont_mask);
else
    if isfield(pc_input_struct, 'T_spont')
        parameter_struct.T_spont = pc_input_struct.T_spont;
    else
        parameter_struct.T_spont = 1;
    end
end

%handle L_forced
L_forced_mask = strcmpi(row_names, 'L_forced');
if any(L_forced_mask)
    parameter_struct.L_forced = x(L_forced_mask);
else
    if isfield(pc_input_struct, 'L_forced')
        parameter_struct.L_forced = pc_input_struct.L_forced;
    else
        parameter_struct.L_forced = 0.5;
    end
end

%handle E_forced
E_forced_mask = strcmpi(row_names, 'E_forced');
if any(E_forced_mask)
    parameter_struct.E_forced = x(E_forced_mask);
else
    if isfield(pc_input_struct, 'E_forced')
        parameter_struct.E_forced = pc_input_struct.E_forced;
    else
        parameter_struct.E_forced = 30;
    end
end

%handle T_forced
T_forced_mask = strcmpi(row_names, 'T_forced');
if any(T_forced_mask)
    parameter_struct.T_forced = x(T_forced_mask);
else
    if isfield(pc_input_struct, 'T_forced')
        parameter_struct.T_forced = pc_input_struct.T_forced;
    else
        parameter_struct.T_forced = 0.5;
    end
end

%handle sw_pc0
sw_pc0_mask = strcmpi(row_names, 'sw_pc0');
if any(sw_pc0_mask)
    parameter_struct.sw_pc0 = x(sw_pc0_mask);
else
    if isfield(pc_input_struct, 'sw_pc0')
        parameter_struct.sw_pc0 = pc_input_struct.sw_pc0;
    else
        parameter_struct.sw_pc0 = 0.5;
    end
end