function [ x ] = get_shell_boundary_pos_from_volume_fractions( volume_fractions, L )
%GET_SHELL_BOUNDARY_POS Summary of this function goes here
%   Detailed explanation goes here
%   sum(volume_fractions) <1
%   1..n-1 entries volume fractions of shells
%   the last shell volume fraction is computed as 1- sum(volume_fractions)
% volume_fractions(end+1) = 1-sum(volume_fractions);

x = zeros(length(volume_fractions)+1,1);
x(1) = 0;
x(end+1) = min(L)/2.0;

p = @(x) minc_proximity_function(x,L);

accum_vol = 0;
for i = 1:length(volume_fractions)
    accum_vol = accum_vol + volume_fractions(i);
    find_x = @(x) p(x)-accum_vol;
    x0 = 1.0;
    x(i+1) = fzero(find_x,x0);
    while isnan(x(i+1))
        x(i+1) = fzero(find_x,x0*rand);
    end
end

