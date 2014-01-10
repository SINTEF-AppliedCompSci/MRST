function N = getGridNeighbors(N)
persistent x
if isempty(x)
    if nargin >= 1
        x = N;
    else
        x = false;
    end
else
    if nargin >= 1
        warning('Input ignored, clear function to reset output')
    end
end
N = x;
end
