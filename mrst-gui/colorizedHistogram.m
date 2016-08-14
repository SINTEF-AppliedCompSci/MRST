function colorizedHistogram(data, n_sample)
    if nargin < 2
        n_sample = 10;
    end;
    assert(size(data, 1) == 1 || size(data, 2) == 1, ...
        'Matrix input not supported. Please provide a vector as input.')
    data = double(data);
    data = data(all(isfinite(data), 2), :);
    % Find data
    [n, X] = hist(data, n_sample);

    h = X(2) - X(1);
    for i = 1:n_sample
        plotRectangle(X(i) - h, n(i), h, X(i));
    end
end


function h = plotRectangle(x, height, width, value, varargin)
    h = patch([x; x+width; x+width; x], [0; 0; height; height], value, varargin{:});
end