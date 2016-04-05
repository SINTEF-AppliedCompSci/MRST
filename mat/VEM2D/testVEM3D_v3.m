clc; clear all; close all;

addpath('../');
run('../startup.m');

n = 10;

G = cartGrid([n,n], [1,1]);
f = @(X) zeros(size(X,1),1);
k = 2;
G = computeVEM2DGeometry(G,f,k,1);

bc = {};

sol = VEM2D_v3(G,f,k,bc)