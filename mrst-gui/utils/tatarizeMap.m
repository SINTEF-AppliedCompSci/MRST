function map = tatarizeMap(n)
%TATARIZEMAP Independent colors maximally dissimilar from all previous colors
%   TATARIZEMAP returns an 64-by-3 matrix containing a colormap
%   The colors have been computed by a brute force approach that generates
%   all possible colors and then picks the most distinct.
%
%   TATARIZEMAP(M) returns an M-by-3 matrix which either consists of the
%   first M colors in the original 64-map if M<=64, or is obtained by
%   linear interpolation if M>64
%
%   EXAMPLE
%
%   This example shows how to reset the colormap of the current figure.
%
%       colormap(tatarizeMap)
%
%   See also AUTUMN, BONE, COLORCUBE, COOL, COPPER, FLAG, GRAY, HOT, HSV,
%   JET, LINES, PARULA, PINK, PRISM, SPRING, SUMMER, WHITE, WINTER,
%   COLORMAP, RGBPLOT.

values = [
0 0 0
1 0 103
213 255 0
255 0 86
158 0 142
14 76 161
255 229 2
0 95 57
0 255 0
149 0 58
255 147 126
164 36 0
0 21 68
145 208 203
98 14 0
107 104 130
0 0 255
0 125 181
106 130 108
0 174 126
194 140 159
190 153 112
0 143 156
95 173 78
255 0 0
255 0 246
255 2 157
104 61 59
255 116 163
150 138 232
152 255 82
167 87 64
1 255 254
255 238 232
254 137 0
189 198 255
1 208 255
187 136 0
117 68 177
165 255 210
255 166 254
119 77 0
122 71 130
38 52 0
0 71 84
67 0 44
181 0 255
255 177 103
255 219 102
144 251 146
126 45 210
189 211 147
229 111 254
222 255 116
0 255 120
0 155 255
0 100 1
0 118 255
133 169 0
0 185 23
120 130 49
0 255 198
255 110 65
232 94 190
];

values = values/255;
if nargin < 1
   map = values;
elseif n<=size(values,1)
   map = values(1:n,:);
else
   P = size(values,1);
   map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');
end