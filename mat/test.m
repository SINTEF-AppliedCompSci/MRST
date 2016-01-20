clc; clear all; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

G = cartGrid([2,2,2],[1,1,1]);
% G.nodes.coords(1:2,:) = G.nodes.coords(1:2,:) -0.5;
G = computeGeometry(G);
G = mrstGridWithFullMappings(G);
G = computeVEMGeometry(G);

VEM3D_loc_v2(G,1)

% I = faceInt2(G) 

%VEM3D_loc(G,1)
% 
% nF = G.faces.num;
% 
% figure;
% axis([-1 2 -1 2 -1 2]);
% view(3);
% hold on;
% for i = 1:nF
%     edgeNum = G.faces.edgePos(i):G.faces.edgePos(i+1)-1;
%     edges = G.faces.edges(edgeNum);
%     nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
%     nodes = G.edges.nodes(nodeNum);
%     plot3(G.nodes.coords(nodes,1), G.nodes.coords(nodes,2), G.nodes.coords(nodes,3), '-ob');
%     nE = numel(edges);
%     edgeNormals = G.faces.edgeNormals(edgeNum,:);
%     for j = 1:nE
%         plot3([G.edges.centroids(edges(j),1), ...
%                G.edges.centroids(edges(j),1) + edgeNormals(j,1)], ...
%               [G.edges.centroids(edges(j),2), ...
%                G.edges.centroids(edges(j),2) + edgeNormals(j,2)],...
%               [G.edges.centroids(edges(j),3), ...
%                G.edges.centroids(edges(j),3) + edgeNormals(j,3)],'--r');
% %            pause
%     end
% %     pause
% end
