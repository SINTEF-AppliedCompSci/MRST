function alpha = get_alpha(x, model)

x = x(:);
row_names = model.history_match.RowNames;
rock = model.experiment.rock;

% handle alpha
alpha_mask = strcmpi(row_names, 'alpha');
if any(alpha_mask)
    alpha = x(alpha_mask);
else
    if isfield(rock, 'alpha')
        alpha = rock.alpha.value;
    else
        alpha = 1;
    end
end
