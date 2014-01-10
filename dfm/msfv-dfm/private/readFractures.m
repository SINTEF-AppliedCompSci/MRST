function [vertices,fractures] = readFractures(input,box,precision)
% read fractures from a multiple of inputs, merge them
% and output a list of vertices and fractures
%
% SYNOPSIS
% [vertices,fractures] = readFractures(input,box,precision)
%
% PARAMETERS:
%       input - a cell of structures where the fracture data is stored
%       box   - bounding box of the domain
%       precision - the relative precision of the grid
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.


vertices = [];
fractures = [];
for i = 1 : numel(input)

    % build the fractures using the input data
    disp('Get fine-scale fractures')
    [vert_new, frac_new] = build_fractures_mod (input{i});

    % all boundaries in the fine-scale is removed.
    frac_new = frac_new( frac_new(:,3)>=0,:);

    if isfield(input{i},'scale')
        vert_new = vert_new .* input{i}.scale;
    end
    % we may want to manage the tags our self
    if isfield(input{i},'tag')
        frac_new(frac_new(:,3)>0 ,3) = input{i}.tag;
    end
    % the 4th column is used for the boundary tags
    frac_new = [frac_new zeros(size(frac_new,1),1)];

    % make sure the vertices are within the box
    vert_new = snap2box(vert_new,box);

    % merge the new fractures with the old ones
    [vertices,fractures] = mergeFractures(vertices,fractures,vert_new,frac_new,box,precision);

    % remove duplicated fractures. Keep the old ones
    [vertices,fractures] = removeDuplicates(vertices,fractures);
    clear frac_new vert_new
end

end

%% helpers
function coords = snap2box(coords,box)
% snap the points to the bounding box
coords(:,1) = max(coords(:,1),box(1,1));
coords(:,2) = max(coords(:,2),box(1,2));
coords(:,1) = min(coords(:,1),box(2,1));
coords(:,2) = min(coords(:,2),box(2,2));
end