function c = colorizeWell(type, index, D)
    switch lower(type)
        case 'prod'
            cmap = jet(numel(D.prod));
            c = cmap(index,:);
        case 'inj'
            cmap = gray(numel(D.inj));
            c = cmap(index,:);
        case 'global'
            if any(D.prod == index)
                c = colorizeWell('prod', find(D.prod == index), D);
            else
                c = colorizeWell('inj', find(D.inj == index), D);
            end
    end
end
