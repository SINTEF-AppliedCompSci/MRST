function rectangular = findRectangularFractures(fractures)
rectangular = [];
for i = 1:numel(fractures)
    if size(fractures(i).points,1) == 4
        p = [fractures(i).points;fractures(i).points(1,:)];
        d = sqrt(sum(diff(p,1).^2,2));
        if d(1) == d(3) && d(2) == d(4) && d(1) + d(2) == d(3) + d(4)
            rectangular = [rectangular;i]; %#ok
        end
    end
end
            