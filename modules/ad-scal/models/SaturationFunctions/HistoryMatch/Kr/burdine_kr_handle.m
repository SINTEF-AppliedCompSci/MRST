function parameter_struct = burdine_kr_handle(model, x, row_names)
%
% DESCRIPTION: builds the parameter input for the saturation function
%              builder based on the excel table given in the input for the
%              history matching
%
% SYNOPSIS:
%   parameter_struct = burdine_kr_handle(model, x, row_names)
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
kr_input_struct = model.experiment.satfun.kr;

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

%handle kroSwc
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