function points=removeexcess(points,trigger)
% Removes excess preallocated rows. These rows have values of trigger.

xcoord=points(:,1);

points=points(xcoord>trigger,:);


end

