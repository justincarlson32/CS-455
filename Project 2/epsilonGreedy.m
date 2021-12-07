function epsilonGreed = epsilonGreedy(qTable, state, explorationRate)
    
    epsilonGreed = rand();

    if (epsilonGreed < explorationRate)
        epsilonGreed = randi([1 4]);
    else
        [reward, action] = max(qTable(state,:));
        epsilonGreed = action(1);
    end

    return;
end