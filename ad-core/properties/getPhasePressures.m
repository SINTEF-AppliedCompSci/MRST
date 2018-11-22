function p_alpha = getPhasePressures(p, pc)
    nph = numel(pc);
    p_alpha = cell(1, nph);
    for i = 1:nph
        if isempty(pc{i})
            p_alpha{i} = p;
        else
            p_alpha{i} = p + pc{i};
        end
    end
end