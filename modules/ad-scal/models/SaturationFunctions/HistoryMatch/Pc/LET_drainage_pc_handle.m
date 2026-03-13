function parameter_struct = LET_drainage_pc_handle(model, x, row_names)
%
% DESCRIPTION: builds the parameter input for the saturation function
%              builder based on the excel table given in the input for the
%              history matching
%
% SYNOPSIS:
%   parameter_struct = LET_drainage_pc_handle(model, x, row_names)
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

%handle entry_multiplier
entry_multiplier_mask = strcmpi(row_names, 'entry_multiplier');
if any(entry_multiplier_mask)
    parameter_struct.entry_multiplier = x(entry_multiplier_mask);
else
    if isfield(pc_input_struct, 'entry_multiplier')
        parameter_struct.entry_multiplier = pc_input_struct.entry_multiplier;
    else
        parameter_struct.entry_multiplier = 0;
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

%handle entry_pc
entry_pc_mask = strcmpi(row_names, 'entry_pc');
if any(entry_pc_mask)
    parameter_struct.entry_pc = x(entry_pc_mask);
else
    if isfield(pc_input_struct, 'entry_pc')
        parameter_struct.entry_pc = pc_input_struct.entry_pc;
    else
        parameter_struct.entry_pc = 0.005;
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

%handle L_entry
L_entry_mask = strcmpi(row_names, 'L_entry');
if any(L_entry_mask)
    parameter_struct.L_entry = x(L_entry_mask);
else
    if isfield(pc_input_struct, 'L_entry')
        parameter_struct.L_entry = pc_input_struct.L_entry;
    else
        parameter_struct.L_entry = 1;
    end
end

%handle E_entry
E_entry_mask = strcmpi(row_names, 'E_entry');
if any(E_entry_mask)
    parameter_struct.E_entry = x(E_entry_mask);
else
    if isfield(pc_input_struct, 'E_entry')
        parameter_struct.E_entry = pc_input_struct.E_entry;
    else
        parameter_struct.E_entry = 60;
    end
end

%handle T_entry
T_entry_mask = strcmpi(row_names, 'T_entry');
if any(T_entry_mask)
    parameter_struct.T_entry = x(T_entry_mask);
else
    if isfield(pc_input_struct, 'T_entry')
        parameter_struct.T_entry = pc_input_struct.T_entry;
    else
        parameter_struct.T_entry = 1;
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
