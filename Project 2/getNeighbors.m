function [nAgent, A] = getNeighbors(nodes, r)

numNodes = size(nodes,1);
dif = cell(numNodes,1); 

distanceAlpha = zeros(numNodes,numNodes);
   
nAgent = cell(numNodes,1); 
                                
for i = 1:numNodes
    dif{i} = repmat(nodes(i,:),numNodes,1) - nodes;
    tmp = dif{i};
    
    for j = 1:numNodes
        dTemp(j,:) = norm(tmp(j,:));
    end
    
    distanceAlpha(i,:)= dTemp; 
end

for k = 1:numNodes
    nAgent{k} = find(distanceAlpha(:,k) < r & distanceAlpha(:,k) ~= 0);
end

A = zeros(numNodes,numNodes);

for i = 1:numNodes
    for j = 1:numNodes
        if i ~= j
            dist_2nodes = norm(nodes(j,:) - nodes(i,:));
            if dist_2nodes < r && dist_2nodes ~= 0
                A(i,j) = 1; 
            end
        end
    end 
end