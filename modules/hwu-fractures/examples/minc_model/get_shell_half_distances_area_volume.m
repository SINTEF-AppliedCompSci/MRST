function [ dist, area, volume ] = get_shell_half_distances_area_volume( volume_fractions,L )
%GET_SHELL_HALF_DISTANCES_AND_AREA Summary of this function goes here
%   Detailed explanation goes here

x = get_shell_boundary_pos_from_volume_fractions( volume_fractions, L );

dist = zeros(length(volume_fractions)+1,1);
area = zeros(length(volume_fractions)+1,1);
volume = zeros(length(volume_fractions)+1,1);

for i=1:length(volume_fractions)
    dist(i) = (x(i+1)-x(i))/2;
    area(i) = minc_proximity_function_dev(x(i),L);
    volume(i) = (minc_proximity_function(x(i+1),L)...
        - minc_proximity_function(x(i),L))*prod(L);
end

dist(end) = x(end)-x(end-1);
area(end) = minc_proximity_function_dev(x(end-1),L);
volume(end) = (minc_proximity_function(x(end),L)...
    - minc_proximity_function(x(end-1),L))*prod(L);
end

