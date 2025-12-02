function parameter_struct = modified_skjaeveland_pc_handle(model, x, row_names)
%
% DESCRIPTION: builds the parameter input for the saturation function
%              builder based on the excel table given in the input for the
%              history matching
%
% SYNOPSIS:
%   parameter_struct = modified_skjaeveland_pc_handle(model, x, row_names)
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

% cwi, coi, ri, bi, swd, sod

%handle cwi
cwi_mask = strcmpi(row_names, 'cwi');
if any(cwi_mask)
    parameter_struct.cwi = x(cwi_mask);
else
    if isfield(pc_input_struct, 'cwi')
        parameter_struct.cwi = pc_input_struct.cwi;
    else
        parameter_struct.cwi = 0.005;
    end
end

%handle coi
coi_mask = strcmpi(row_names, 'coi');
if any(coi_mask)
    parameter_struct.coi = x(coi_mask);
else
    if isfield(pc_input_struct, 'coi')
        parameter_struct.coi = pc_input_struct.coi;
    else
        parameter_struct.coi = -0.003;
    end
end

%handle ri
ri_mask = strcmpi(row_names, 'ri');
if any(ri_mask)
    parameter_struct.ri = x(ri_mask);
else
    if isfield(pc_input_struct, 'ri')
        parameter_struct.ri = pc_input_struct.ri;
    else
        parameter_struct.ri = -0.24;
    end
end

%handle bi
bi_mask = strcmpi(row_names, 'bi');
if any(bi_mask)
    parameter_struct.bi = x(bi_mask);
else
    if isfield(pc_input_struct, 'bi')
        parameter_struct.bi = pc_input_struct.bi;
    else
        parameter_struct.bi = 0.08;
    end
end

%handle swd
swd_mask = strcmpi(row_names, 'swd');
if any(swd_mask)
    parameter_struct.swd = x(swd_mask);
else
    if isfield(pc_input_struct, 'swd')
        parameter_struct.bi = pc_input_struct.swd;
    else
        parameter_struct.swd = 0.4;
    end
end

%handle sod
sod_mask = strcmpi(row_names, 'sod');
if any(sod_mask)
    parameter_struct.sod = x(sod_mask);
else
    if isfield(pc_input_struct, 'sod')
        parameter_struct.sod = pc_input_struct.sod;
    else
        parameter_struct.sod = 0.65;
    end
end