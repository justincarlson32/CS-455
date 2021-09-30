clc,clear
close all

d = 15;
k = 1.2;
r = k*d; 
rPrime = .25*k*r;
epsilon = 0.1;
numNodes = 100;
n = 2;

graphSize = 50;

nodes = graphSize.*rand(numNodes, n) + graphSize.*repmat([0 1], numNodes, 1);
pNodes = zeros(numNodes,n);

deltaT = 0.05;
t = 0:deltaT:5;

oldNodes = nodes;

qMean = zeros(size(t,2), n); %position of nodes
pMean = zeros(size(t,2), n); %velocity of nodes


connectivity = []; %connectivity instantiation


qNodesCell = cell(size(t,2), numNodes); %cells of positions of nodes
pNodesCell = cell(size(t,2), numNodes); %cells of velocityies of nodes

nFrames = 20;
mov(1:nFrames) = struct('cdata', [],'colormap', []);

for i = 1:length(t)
    [nAgent, A] = getNeighbors(nodes, r);
    [Ui] = algorithm1(numNodes, nodes, nAgent, n, epsilon, r, d, pNodes); 

    pNodes = (nodes - oldNodes)/deltaT; %velocity of nodes
    
    pNodesCell{i} = pNodes; %updating cell of velocity of nodes
    
    oldNodes = nodes;
    nodes = oldNodes + pNodes*deltaT + 0.5 * Ui* deltaT * deltaT;
    qMean(i,:) = mean(nodes);

    
    pMean(i,:) = mean(pNodes); %average velocity for graphing
    
    qNodesCell{i} = nodes; %cell array for nodes
    
    
    connectivity(i)= (1 / numNodes) * rank(A); %connectivity for graphing

    
    %plotting points
    plot(nodes(:,1),nodes(:,2), '.')
    hold on
    
    
    %connecting nodes
    plot(nodes(:,1),nodes(:,2), 'k>','LineWidth',.2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5) 
    hold off
    
    for node = 1:numNodes
        tmp=nodes(nAgent{node},:);
        for j = 1:size(nodes(nAgent{node},1))
            line([nodes(node,1),tmp(j,1)],[nodes(node,2),tmp(j,2)])
        end
    end
    
    mov(i) = getframe;
    hold off

end