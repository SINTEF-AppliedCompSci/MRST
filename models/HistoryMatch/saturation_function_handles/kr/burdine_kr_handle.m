function parameter_struct = burdine_kr_handle(model, x, row_names)

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