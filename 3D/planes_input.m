function fractures = planes_input()
fractures = struct;

% fractures(1).points =  [15 10 3;
%                         60 75 7;
%                         60 75 10;
%                         15 10 7];
%                     %{
%                     [20 2 0;
%                         20 75 0;
%                         20 75 5;
%                         20 2 5];
%                         %}
% fractures(2).points = [70 10 10;
%                         70 70 10;
%                         90 70 10;
%                         85 10 10];  
% 
% fractures(1).points =  [20 20 25;
%                         20 20 75;
%                         20 80 75;
%                         20 80 25];
% fractures(2).points = [50 20 50;
%                         80 20 50;
%                         80 80 50;
%                         50 80 50]; 
% 
% fractures(1).points = [5 5 2;
%                        5 5 4;
%                        10 20 4;
%                        10 20 2];
% fractures(2).points = [5 10 2;
%                        5 10 4;
%                        20 25 4;
%                        20 25 2];
                   
% fractures(3).points = [10 5 1;
%                        10 5 5;
%                        25 20 5;
%                        25 20 1]; 
% fractures(3).aperture = 1/25;
% fractures(1).aperture = 1/25;
% fractures(2).aperture = 1/25;

% 
% fractures(1).points = [7 5  2;
%                        7 25 2;
%                        5 25 4;
%                        5 5  4];
% fractures(2).points = [11 5  2;
%                         11 25 2;
%                         13 25 4;
%                         13 5  4];
% fractures(3).points = [19 5  2;
%                         19 25 2;
%                         17 25 4;
%                         17 5  4];
% fractures(4).points = [23 5  2;
%                         23 25 2;
%                         25 25 4;
%                         25 5  4];
% fractures(1).aperture = 1/25;
% fractures(2).aperture = 1/25;
% fractures(3).aperture = 1/25;
% fractures(4).aperture = 1/25;


fractures(1).points = [15 5  1;
                       15 25 1;
                       15 25 5;
                       15 5  5];

fractures(1).aperture = 1/25;

return