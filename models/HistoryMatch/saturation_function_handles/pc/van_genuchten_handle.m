function parameter_struct = van_genuchten_handle(x, row_names)

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

%handle alpha
alpha_mask = strcmpi(row_names, 'alpha');
if any(alpha_mask)
    parameter_struct.alpha = x(alpha_mask);
else
    if isfield(pc_input_struct, 'alpha')
        parameter_struct.alpha = pc_input_struct.alpha;
    else
        parameter_struct.alpha = 10;
    end
end

% handle n
n_pc_mask = strcmpi(row_names, 'n_pc');
n_kr_mask = strcmpi(row_names, 'n_kr');
if any(n_pc_mask) && any(n_kr_mask)
    parameter_struct.n_pc = x(n_pc_mask);
    parameter_struct.n_kr = x(n_kr_mask);
elseif any(n_pc_mask) && not(any(n_kr_mask))
    parameter_struct.n_pc = x(n_pc_mask);
    parameter_struct.n_kr = x(n_pc_mask);
elseif not(any(n_pc_mask)) && any(n_kr_mask)
    parameter_struct.n_pc = x(n_kr_mask);
    parameter_struct.n_kr = x(n_kr_mask);
elseif not(any(n_pc_mask)) && not(any(n_kr_mask))
    if isfield(pc_input_struct, 'n')
        parameter_struct.n_pc = pc_input_struct.n;
    else
        parameter_struct.n_pc = 2;
    end
    if isfield(kr_input_struct, 'n')
        parameter_struct.n_kr = kr_input_struct.n;
    else
        parameter_struct.n_kr = 2;
    end
end