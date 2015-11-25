clc; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

f =  @(X)  4*pi^2*sin(X(:,1)*pi).*cos(X(:,2)*pi);
gD = @(X)  2*sin(X(:,1)*pi).*cos(X(:,2)*pi) - log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
gN = @(X) -2*pi*cos(X(:,1)*pi).*cos(X(:,2)*pi) - 2*(X(:,1)+0.1)./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2);

meshSize = [8, 16, 32, 64];
cellVec = zeros(4,1);
hVec = zeros(4,1);
errVec = zeros(4,1);
for i = 1:4
    n = meshSize(i);
    nx = n; ny = n;
    G = unitSquare(nx, ny);
    %G = cartGrid([nx, ny], [1,1]);
    G = sortEdges(G);
    G = mrstGridWithFullMappings(G);
    G = computeGeometry(G);

    tol = 1e-6;
    boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
    bNeu = abs(G.faces.centroids(boundaryEdges,1)) < tol;
    bc = struct('bcFunc', {{gN, gD}}, 'bcFaces', {{boundaryEdges(bNeu), boundaryEdges(~bNeu)}}, 'bcType', {{'neu', 'dir'}});
    %bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});
    
    U = VEM2D(G,f,bc);

    X = [G.nodes.coords ; G.faces.centroids];
    Uexact = gD(X);

    Nc = G.cells.num;
    hM = 0;
    for c = 1:Nc
        nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
        nodes = G.cells.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);
        hM = hM + cellDiameter(X);
    end
    hM = hM/Nc;
    cellVec(i) = Nc;
    hVec(i) = hM;
    
    errVec(i) = hM*norm(U(1:end-G.cells.num) - Uexact, inf);
    
    clf;
    
    plotGrid(G, 'facecolor', 'w');
    xlabel('x');
    ylabel('y');
    axis equal off;
    fileName = strcat('../tex/project/Figures/2D_grid_',num2str(i),'.eps');
    saveas(gcf, fileName, 'eps');
    
    clf;
    
    plotVEM(G, U, '');
    view(-45,38);
    axis([0, 1, 0, 1, -4, 1.5]);
    set(gcf,'Units','normalized')
    fileName = strcat('../tex/project/Figures/2D_sol_',num2str(i),'.eps');
    saveas(gcf, fileName, 'eps');
    
    
end

slope = polyfit(log(hVec), log(errVec), 1);
slope = slope(1);
loglog(hVec, errVec, 'k');
str = strcat('Slope = ', num2str(slope));
fontSize = 18;
text(0.15,1e-5,str, 'FontSize', fontSize);
xlabel('log h', 'FontSize', fontSize);
ylabel('log e_h', 'FontSize', fontSize);
axis equal;
fileName = '../tex/project/Figures/2D_conv.eps';
saveas(gcf, fileName, 'eps');
