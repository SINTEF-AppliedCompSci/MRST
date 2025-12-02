function parameter_struct = brooks_corey_kr_handle(model, x, row_names)
%
% DESCRIPTION: builds the parameter input for the saturation function
%              builder based on the excel table given in the input for the
%              history matching
%
% SYNOPSIS:
%   parameter_struct = brooks_corey_kr_handle(model, x, row_names)
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

% handle lambda
lambda_pc_mask = strcmpi(row_names, 'lambda_pc');
lambda_kr_mask = strcmpi(row_names, 'lambda_kr');
if any(lambda_pc_mask) && any(lambda_kr_mask)
    parameter_struct.lambda_pc = x(lambda_pc_mask);
    parameter_struct.lambda_kr = x(lambda_kr_mask);
elseif any(lambda_pc_mask) && not(any(lambda_kr_mask))
    parameter_struct.lambda_pc = x(lambda_pc_mask);
    parameter_struct.lambda_kr = x(lambda_pc_mask);
elseif not(any(lambda_pc_mask)) && any(lambda_kr_mask)
    parameter_struct.lambda_pc = x(lambda_kr_mask);
    parameter_struct.lambda_kr = x(lambda_kr_mask);
elseif not(any(lambda_pc_mask)) && not(any(lambda_kr_mask))
    if isfield(pc_input_struct, 'lambda')
        parameter_struct.lambda_pc = pc_input_struct.lambda;
    else
        parameter_struct.lambda_pc = 2;
    end
    if isfield(kr_input_struct, 'lambda')
        parameter_struct.lambda_kr = kr_input_struct.lambda;
    else
        parameter_struct.lambda_kr = 2;
    end
end