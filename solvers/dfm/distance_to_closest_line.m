function varargout = distance_to_closest_line(varargin)
% Copyright 2011-2012 University of Bergen
%
% This file is licensed under the GNU General Public License v3.0.

buildmex distance_to_closest_line.cpp
[varargout{1:nargout}] = distance_to_closest_line(varargin{:});