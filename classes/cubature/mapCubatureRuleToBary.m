function [lambda, w] = mapCubatureRuleToBary(rule, type)

    if strcmpi(type, 'triangle')
        t = '%f%f';
    elseif strcmpi(type, 'tetrahedron')
        t = '%f%f%f';
    end

    fid = fopen([rule, '_x.txt']);
    x = textscan(fid, t, 'HeaderLines', 0, 'CollectOutput', 1);
    x = x{:};
    fclose(fid);
    
    fid = fopen([rule, '_w.txt']);
    w = textscan(fid, '%f', 'HeaderLines', 0, 'CollectOutput', 1);
    w = w{:};
    fclose(fid);
    
    fid = fopen([rule, '_r.txt']);
    v = textscan(fid, t, 'HeaderLines', 0, 'CollectOutput', 1);
    v = v{:};
    fclose(fid);
    
    n = numel(w);
    
    x = [x, ones(size(x,1),1)];
    R = [v, ones(size(v,1),1)];
    lambda = x/R;

end