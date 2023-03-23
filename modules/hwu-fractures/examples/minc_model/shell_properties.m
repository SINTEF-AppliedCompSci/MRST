function [ dist, area ] = shell_properties( volume_fractions,block_dimensions )
%shell_properties Summary of this function goes here
%   Detailed explanation goes here

x = get_shell_boundary_pos_from_volume_fractions( volume_fractions(1:end-1), block_dimensions );

dist = zeros(length(volume_fractions),1);
area = zeros(length(volume_fractions),1);

for i=1:length(volume_fractions)-1
    dist(i) = (x(i+1)-x(i))/2;
    area(i) = minc_proximity_function_dev(x(i),block_dimensions)*prod(block_dimensions);
end

dist(end) = x(end)-x(end-1);
area(end) = minc_proximity_function_dev(x(end-1),block_dimensions)*prod(block_dimensions);

end

