clc,clear
close all

%ALGORITHM SELECT
% 0 = algorithm 1
% 1 = algorithm 2
% 2 = algorithm 3
% 3 = algorithm 4

algType = 1;


%constants
d = 15;
k = 1.2;
r = k*d; 
rPrime = .22*k*r;
epsilon = 0.1;

numNodes = 30;
n = 2;

graphSize = 150;

%target
q1 = [];
p1 = [];

%-------------------Algorithm Type-specific Setup--------------------------


if (algType == 1 || algType == 0)
    graphSize = 50;
end

%target creating
if (algType ~= 0) 
    q1 = [150 150];
	p1 = [0 0]; 
end

%-------------------end----------------------------------------------------

nodes = graphSize.*rand(numNodes, n) + graphSize.*repmat([0 1], numNodes, 1);
pNodes = zeros(numNodes,n);

deltaT = 0.01;
t = 0:deltaT:5;

oldNodes = nodes;

qMean = zeros(size(t,2), n); %position of nodes
pMean = zeros(size(t,2), n); %velocity of nodes


connectivity = []; %connectivity instantiation


qNodesCell = cell(size(t,2), numNodes); %cells of positions of nodes
pNodesCell = cell(size(t,2), numNodes); %cells of velocities of nodes

%screen updating
nFrames = 20;
mov(1:nFrames) = struct('cdata', [],'colormap', []);



for i = 1:length(t)
    
    subplot(2,2,1) % main plot
    
    %--------trajectories creations--------------
    if (algType ~= 0)
        
        qX1 = 0;
        qY1 = 0;
                
        if (algType == 2) %sine-wave
           qX1 = 50 + 50*t(i);
           qY1 = 295 - 50*sin(t(i));
        end
        
        if (algType == 3) %circle 
            qX1 = 310 - 160*cos(t(i)); 
            qY1 = 255 - 160*sin(t(i)); 
        end

        if (algType == 2 || algType == 3)
            q1(i,:) = [qX1, qY1];
            if (i > 1)
                p1(i,:) = (q1(i,:) - q1(i-1,:)) / deltaT;
            else
                continue 
            end 
        end
        
        plot(q1(:,1),q1(:,2), 'ro', 'LineWidth', 2, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4)
        hold on
    end
    %----------------------------------------------
    
    %---------------Algorithm-type Update----------
    
    if algType == 0 %flocking
        [nAgent, A] = getNeighbors(nodes, r);
        [Ui] = algorithm1(numNodes, nodes, nAgent, n, epsilon, r, d, pNodes); 
    elseif algType == 1 % static target
        [nAgent, A] = getNeighbors(nodes, r);
        [Ui] = algorithm2(algType, numNodes, nodes, nAgent, n, epsilon, r, d, q1, p1, pNodes); 
    elseif algType == 2 || algType == 3 %dynamic targeting
        [nAgent, A] = getNeighbors(nodes, r);
        [Ui] = algorithm2(algType, numNodes, nodes, nAgent, n, epsilon, r, d, q1(i,:), p1(i,:), pNodes); 
    end
    
    %-----------------------------------------------
    
    %---------------Plot updating------------------
    
    pNodes = pNodes + Ui * deltaT;%(nodes - oldNodes) / deltaT; % compute velocities of nodes
    pNodesCell{i} = pNodes;
    oldNodes = nodes;
    nodes = oldNodes + (pNodes * deltaT) + (((deltaT^2)/2) * Ui);
    qMean(i,:) = mean(nodes); % compute position of nodes
    
    if (algType ~= 0)
        plot(qMean(:,1),qMean(:,2),'ro', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
        hold on
    end
    
    pMean(i,:) = mean(pNodes); %Compute velocity of COM of MSN
    qNodesCell{i} = nodes; 
    Connectivity(i)= (1 / numNodes) * rank(A);
    
    for node = 1:numNodes
        tmp = nodes(nAgent{node},:);
        for j = 1:size(nodes(nAgent{node},1))
            line([nodes(node,1),tmp(j,1)],[nodes(node,2),tmp(j,2)])
        end
    end
   
    plot(nodes(:,1),nodes(:,2), '.')
    hold on
    plot(nodes(:,1),nodes(:,2), 'k>','LineWidth',.2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5) 
    hold off
    
     
    mov(i) = getframe;
    hold off
    
    title('MSN Graph')
%----------------------Plotting velocity------------------
    subplot(2,2,2) % velocity plot
    title('Velocity')
    pNodei = [];
    for i2 = 2:i %iterates through the timesteps for the history cell matrix
      tmp = pNodesCell{i2};
      for j = 1:numNodes
        if j == 1 %Plot velociy of sensor node 1; you can change this number to plot for other nodes 
            pNodei(i2) = norm(tmp(j,:));
            plot(pNodei, 'b')
            hold on
        end
      end
    end
%----------------Plotting Connectivity----------------------   
   subplot(2,2,3) %connectivity plot
   plot(Connectivity)
   title('Connectivity')
   
%----------------Plotting Trajcetory------------------------
    subplot(2,2,4) %traj plot
    for i2 = 2:i
     tmp = qNodesCell{i2};
     plot(tmp(:,1), tmp(:,2), 'k.') 
     hold on
    end
    hold on
    plot(nodes(:,1), nodes(:,2), 'm>', 'LineWidth', .2, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 5)
    title('Trajectory')
%-----------------------------------------------
end



