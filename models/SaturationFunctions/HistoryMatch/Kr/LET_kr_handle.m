function parameter_struct = LET_kr_handle(model, x, row_names)
%
% DESCRIPTION: builds the parameter input for the saturation function
%              builder based on the excel table given in the input for the
%              history matching
%
% SYNOPSIS:
%   parameter_struct = LET_kr_handle(model, x, row_names)
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

% handle krwSor
krwSor_mask = or(strcmpi(row_names, 'krwSor'), strcmpi(row_names, 'krw@Sor'));
if any(krwSor_mask)
    parameter_struct.krwSor = x(krwSor_mask);
else
    if isfield(kr_input_struct, 'krwSor')
        parameter_struct.krwSor = kr_input_struct.krwSor;
    else
        parameter_struct.krwSor = 1;
    end
end

% handle kroSwc
kroSwc_mask = or(strcmpi(row_names, 'kroSwc'), strcmpi(row_names, 'kro@Swc'));
if any(kroSwc_mask)
    parameter_struct.kroSwc = x(kroSwc_mask);
else
    if isfield(kr_input_struct, 'kroSwc')
        parameter_struct.kroSwc = kr_input_struct.kroSwc;
    else
        parameter_struct.kroSwc = 1;
    end
end

% handle Lw
Lw_mask = strcmpi(row_names, 'Lw');
if any(Lw_mask)
    parameter_struct.Lw = x(Lw_mask);
else
    if isfield(kr_input_struct, 'Lw')
        parameter_struct.Lw = kr_input_struct.Lw;
    else
        parameter_struct.Lw = 2;
    end
end

% handle Ew
Ew_mask = strcmpi(row_names, 'Ew');
if any(Ew_mask)
    parameter_struct.Ew = x(Ew_mask);
else
    if isfield(kr_input_struct, 'Ew')
        parameter_struct.Ew = kr_input_struct.Ew;
    else
        parameter_struct.Ew = 2;
    end
end

% handle Tw
Tw_mask = strcmpi(row_names, 'Tw');
if any(Tw_mask)
    parameter_struct.Tw = x(Tw_mask);
else
    if isfield(kr_input_struct, 'Tw')
        parameter_struct.Tw = kr_input_struct.Tw;
    else
        parameter_struct.Tw = 2;
    end
end

% handle Lnw
Lnw_mask = strcmpi(row_names, 'Lnw');
if any(Lnw_mask)
    parameter_struct.Lnw = x(Lnw_mask);
else
    if isfield(kr_input_struct, 'Lnw')
        parameter_struct.Lnw = kr_input_struct.Lnw;
    else
        parameter_struct.Lnw = 2;
    end
end

% handle Enw
Enw_mask = strcmpi(row_names, 'Enw');
if any(Enw_mask)
    parameter_struct.Enw = x(Enw_mask);
else
    if isfield(kr_input_struct, 'Enw')
        parameter_struct.Enw = kr_input_struct.Enw;
    else
        parameter_struct.Enw = 2;
    end
end

% handle Tnw
Tnw_mask = strcmpi(row_names, 'Tnw');
if any(Tnw_mask)
    parameter_struct.Tnw = x(Tnw_mask);
else
    if isfield(kr_input_struct, 'Tnw')
        parameter_struct.Tnw = kr_input_struct.Tnw;
    else
        parameter_struct.Tnw = 2;
    end
end