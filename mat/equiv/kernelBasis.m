clc; clear; close all;

q1 = [ 1 -1  1 -1  1 -1  1 -1];
q2 = [ 1 -1 -1  1 -1  1  1 -1];
q3 = [ 1  1 -1 -1 -1 -1  1  1];
q4 = [-1  1 -1  1  1 -1  1 -1];

% h = 2;
% alpha = sqrt(9/8);
% beta = sqrt(27/8);
% 
% h = 1;
% alpha = sqrt(144);
% beta = sqrt(1728);

% Q = [alpha*q1', alpha*q2', alpha*q3', beta*q4'];

Q = [q1', q2', q3', q4'];

% S = 1/diag([h^2, h^2, h^2, h^3],0);

S = 1/(2*sqrt(2))*eye(4);

S'*Q'*Q*S
