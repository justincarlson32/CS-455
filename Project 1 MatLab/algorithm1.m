function [Ui] = algorithm1(numNodes, nodes, nAgent, n, epsilon, r, d, pNodes)

    c1 = 30;
    c2 = 2*sqrt(c1);
    
    a = 5;
    b = 5;
    c = abs(a - b) / sqrt(4*a*b);
    
    rSigma = sigmaNorm(r);
    dSigma = sigmaNorm(d);

    nIJ = zeros(numNodes,numNodes,n);
    for i = 1:numNodes
        for j = 1:numNodes
            q = norm(nodes(j,:) - nodes(i,:));
            sigmaGradient = (nodes(j,:) - nodes(i,:)) / (1 + epsilon * sigmaNorm(nodes(j,:) - nodes(i,:)));
            if q < r && q ~= 0
                nIJ(i,j,:) = sigmaGradient;
            end   
        end 
    end
   
    U = zeros(numNodes, n);
    
    gradient = 0;
    conscensus = 0;
    adjacencyIJ = zeros(numNodes,numNodes);
    
    for i = 1:numNodes
        for j = 1:size(nAgent{i})
            nVal = nAgent{i}(j);
            if(i ~= nVal) 
                z = sigmaNorm(nodes(nVal,:) - nodes(i,:));
                zPhi = z - dSigma;
                rhoH = bumpFunction(z / rSigma);
                sigmoidal = (zPhi + c) / sqrt(1 + (zPhi + c)^2);
                phi = .5 * ((a + b) * sigmoidal + (a - b));
                
                phiAlpha = rhoH * phi;
                
                adjacencyIJ(i,nVal) = rhoH;
               
                gradient = phiAlpha * [nIJ(i,nVal,1) nIJ(i,nVal,2)];

                conscensus = adjacencyIJ(i,nVal) * (pNodes(nVal,:) - pNodes(i,:));
            end
        end
        
        U(i,:) = (c1 * gradient) + (c2 * conscensus);
    end
    
    Ui = U;
end