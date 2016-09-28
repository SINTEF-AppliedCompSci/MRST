clc; clear; close all;

fprintf('--------------------------------------------------------------\n');
fprintf('Test 1: Consistency. Errors should be close to machine epsilon\n');
fprintf('--------------------------------------------------------------\n');

run unitTest.m

fprintf('\n------------------------------------------------------------------ \n');
fprintf('Test 2: Comparison with MFD for soruce/pressure bc/outflow bc problem\n');
fprintf('---------------------------------------------------------------------\n');

run unitTest2.m

fprintf('\n---------------------------------------------------------\n');
fprintf('Test 3: Comparison with MPFA for pressure drop problem bc/ \n');
fprintf('-----------------------------------------------------------\n');

run unitTest3.m