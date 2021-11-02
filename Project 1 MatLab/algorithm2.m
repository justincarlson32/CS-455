function [Ui] = algorithm2(algType, numNodes, nodes, nAgent, n, epsilon, r, d, q1, p1, pNodes) 
    

    c1 = 30;
    c2 = 2*sqrt(c1);
    
    cm1 = 5.1;
    cm2 = 2*sqrt(cm1);
    
    a = 5;
    b = 5;
    c = abs(a - b) / sqrt(4*a*b);
    U = zeros(numNodes, n);

    for i = 1:numNodes
        gradient = 0;
        conscensus = 0;
        for j = 1:numNodes
            dist = norm(nodes(j,:) - nodes(i,:));
            if(dist < r && i ~= j) 
                z = nodes(j,:) - nodes(i,:);          
                normZ = norm(z);
                phiAlphaVal = phiAlpha(sigmaNorm(normZ),a,b,c,d,r);
                gradient = phiAlphaVal * getNIJ(i,j,nodes, epsilon);
                conscensus = getAIJ(i,j,nodes,r) * (pNodes(j,:) - pNodes(i,:));
            end
        end
        
    
        fG = cm1 * (nodes(i,:) - q1);
        fA = (c1 * gradient) + (c2 * conscensus);
        
        U(i,:) = fA - fG;
        
        if (algType ~= 1)
            U(i,:) = fA - fG - cm2 * (pNodes(i,:) - p1);
        end
    end
    
    
    Ui = U;

end


function returnVal = sigma1(z)

    returnVal = 1 + z^2;
    returnVal = sqrt(returnVal);
    returnVal = z/returnVal;
    
end

function phi = phi(z,c,a,b)

    val1 = a + b;
    val2 = sigma1(z + c);
    val3 = a - b;
    
    val = val1 * val2 + val3;
    val = val/2;

    phi = val;
    
end

function phiAlpha = phiAlpha(z,a,b,c,d,r)

    rSigm = z / sigmaNorm(r);
    dSigm = z - sigmaNorm(d);
    val1 = bumpFunction(rSigm);
    val2 = phi(dSigm,a,b,c);
    phiAlpha = val1 * val2;

end

function getAIJ = getAIJ(i,j, nodes, r)

    val1 = nodes(j,:) - nodes(i,:);
    normal = norm(val1);
    val2 = sigmaNorm(normal)/sigmaNorm(r);
    getAIJ = bumpFunction(val2);

end

function getNIJ = getNIJ(i,j,nodes, epsilon)

    val1 = nodes(j,:) - nodes(i,:);
    normal = norm(val1);
    val2 = 1 + epsilon * normal^2;
    getNIJ = val1/sqrt(val2);

end