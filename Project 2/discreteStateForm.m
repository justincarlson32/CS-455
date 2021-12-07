function state = discreteStateForm(nodes, curNodeIndex, neighbors)
    
    if isempty(neighbors)
        state = 1;
        return;
    end
      
    minDist = 1000000; % really large initialized value
    
    minIndex = 1; % lowest possible index
    for i = 1:length(neighbors)
        dist = norm(nodes(curNodeIndex, :) - nodes(neighbors(i),:));
        if dist < minDist
            minDist = dist;
            minIndex = neighbors(i);
        end
    end
    
    state = (minIndex + 1) * 2;
end