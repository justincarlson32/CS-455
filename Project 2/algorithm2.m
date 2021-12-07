function [Ui] = algorithm2(nodes, nAgent, numNodes, r, d, epsilon, pNodes, actions) 
    c1a = 6;
    c2a = 2*sqrt(c1a);
    c1mt = 1.1;
    Ui = zeros(numNodes, 2);
    
    a = 5;
    b = 5;
    c = abs(a-b) / sqrt(4*a*b);
    
    for i = 1:numNodes
        
        gradient = 0;
        consensus = 0;
        feedback = 0;
        
        q_mt = pointFromAction(actions(i)); % get target point to swap out with original static
        
        for j = 1:length(nAgent{i})
            z = nodes(nAgent{i}(j),:) - nodes(i,:);          
            normZ = sigmaNorm(z);
            phiAlphaVal = phiAlpha(normZ,a,b,c,d,r);
            gradient = gradient + phiAlphaVal * getNIJ(nodes(i,:), nodes(nAgent{i}(j),:), epsilon);
            consensus = consensus + getAIJ(nodes(i,:), nodes(nAgent{i}(j),:), r) * (pNodes(nAgent{i}(j),:) - pNodes(i,:));
        end
        feedback = nodes(i,:) - q_mt;
        
        Ui(i,:) = ((c1a * gradient) + (c2a * consensus) - (c1mt * feedback)); 
    end

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

function getAIJ = getAIJ(i,j, r)

    val1 = j - i;
    normal = sigmaNorm(val1);
    val2 = normal/sigmaNorm(r);
    getAIJ = bumpFunction(val2);

end

function getNIJ = getNIJ(i,j,epsilon)

    val1 = j - i;
    normal = sigmaNorm(val1);
    val2 = 1 + epsilon * normal;
    getNIJ = val1/sqrt(val2);

end