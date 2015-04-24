function plotCurveTube(curve, radius, varargin)

    opt = struct('data',    [], ...
                 'interpMethod', 'linear', ...
                 'interpDataMethod', 'linear', ...
                 'EdgeColor', {'none'}, ...
                 'color',   {[1 0 0]}, ...
                 'refine',  1);
    opt = merge_options(opt, varargin{:});
    
    [x, y, z] = cylinder(1);
    C = [z(1, :)', y(1, :)', x(1, :)'];
    
    nel = size(curve, 1);

    
    if size(radius, 1) == 1 && nel ~= 1
        radius = repmat(radius, nel, 1);
    end
    
    data = opt.data;
    if opt.refine ~= 1 && nel > 1
        d = cumsum([0; sqrt(sum(diff(curve).^2, 2))]);
        len = d(end);
        d_new = (0:len/(nel*opt.refine):len)';
        % Include old data points
        d_new = unique([d_new; d]);
        d_new = sort(d_new);

        curve  = interp1(d, curve,  d_new, opt.interpMethod);
        radius = interp1(d, radius, d_new, opt.interpMethod);
        if ~isempty(data)
            data = interp1(d, data, d_new, opt.interpDataMethod);
        end
        nel = size(curve, 1);
    end
    
    if size(curve, 2) == 2
        curve = [curve, zeros(nel, 1)];
        radius = [radius, mean(radius, 2)];
    end
    
    hasData = ~isempty(data);
    
    if nel == 1
        % We were only given a single data point - extend it to two almost
        % identical data points
        vectors = [0 0 1; 0 0 1];
        curve = [curve; curve];
        curve(:, 3) = curve(:, 3) +  sqrt(eps)*[-1; 1];
        if hasData
            data = [data; data];
        end
        radius = [radius; radius];
        nel = 2;
    else
        % Tangential vectors along the curve using first order
        % approximation
        vectors = diff(curve, 1, 1);
        % First element repeated to ensure that we have a vector for each
        % data point
        vectors = [vectors(1, :); vectors];
        vectors = normalize(vectors);
    end
    [x, y, z] = deal(zeros(nel, size(C, 1)));
    
    for i = 1:nel
        % Loop through all the segments and rotate the reference circles,
        % adding them to the list of cylinders as we go
        v = vectors(i, :);
        % Rotate
        pts = rotate(C, v);
        pts = bsxfun(@mtimes, pts, radius(i, :));
        % Translate
        pts = bsxfun(@plus, pts, curve(i, :));

        sub = i;
        x(sub, :) = pts(:, 1)';
        y(sub, :) = pts(:, 2)';
        z(sub, :) = pts(:, 3)';
    end
    
    if hasData
        d = repmat(data, 1, size(z, 2));
        surf(x, y, z, d, 'EdgeColor', opt.EdgeColor)
    else
        surf(x, y, z, 'EdgeColor', opt.EdgeColor, 'FaceColor', opt.color, 'CData', [])
    end
    
    % Plot circles at the endpoints
    for i = [1, size(x, 1)];
        if hasData
            patch(x(i, :), y(i, :), z(i, :), d(i, :))
        else
            patch(x(i, :), y(i, :), z(i, :), nan, 'FaceColor', opt.color)
        end
    end
end
        
function pts = rotate(pts, newvec)
    angle = acos(newvec);
    M = rotateX(angle(1))*rotateY(angle(2));
    for i = 1:size(pts, 1)
        pts(i, :) = pts(i, :)*M;
    end
end

function M = rotateX(o)
    M = [1,      0,        0; ...
         0, cos(o),  -sin(o); ...
         0, sin(o),   cos(o);];
end

function M = rotateY(o)
    M = [cos(o),  0, sin(o); ...
         0,     1,     0; ...
         -sin(o), 0, cos(o);];
end

function M = rotateZ(o)
    M = [cos(o),  -sin(o), 0; ...
         sin(o), cos(o), 0; ...
         0, 0, 1;];
end

function vectors = normalize(vectors)
    vectors = bsxfun(@rdivide, vectors, sqrt(sum(vectors.^2, 2)));
end
