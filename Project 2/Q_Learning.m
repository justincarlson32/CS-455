function [Q_update, Connectivity, R_nodes, R_sum_all, A_sum_cooQ, mean_Delta_Q, q_nodes_all]  = Q_Learning(Q_update, statelist, actionlist, nstates, nactions, num_nodes, n, nodes, epsilon_learning, t, delta_t, safePoints, learningType)
    
    d = 30;% distance between nodes
    k = 1.1;% scaling factor
    r = d*k;% interaction range
    epsilon = 0.1;% sigma norm
    alpha = 0.2;
    gamma = 0.9;
    
    w = 0.5;% fudge it, 50/50 neighbors sounds good
    
    
    R_nodes = zeros(length(t), num_nodes);% store rewards
    A_sum_cooQ = 1:length(t);% action value of nodes
    q_nodes_all = cell(length(t),1);% position of nodes 
    mean_Delta_Q = 1:length(t);% mean change in Q values
    p_nodes = zeros(num_nodes,n);%velocity of nodes
    nodes_old = nodes; % old nodes
    Connectivity = 1:length(t); % connectivity
    R_sum_all = 1:length(t);% sum of rewards


    center_of_mass = zeros(size(t,2),n);% center of mass
    
    curStates = [];% cur states
    curActions = [];% cur actions
    nextStates = [];% next states
    nextActions = [];% next actions
    
    [Nei_agent, A] = getNeighbors(nodes, r);% generate neighbors and adj matrix
    
    for i = 1:num_nodes
        curStates(i) = discreteStateForm(nodes, i, Nei_agent{i});
        if (1/(num_nodes))*rank(A) == 1 % if fully connected 
            curStates(:) = 2;
        end
        curActions(i) = epsilonGreedy(Q_update{i}, curStates(i), epsilon_learning);
    end
    
    Q_prev = cell(1, num_nodes);% previous iterators of qtables for graphing
    
    %episode iterators
    for iterator = 1:length(t)
        
        % plotting safe spots
        plot(safePoints(:,1),safePoints(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.5)
        hold on
        
        A_sum_cooQ(iterator) = 0;% action sum = 0

        [Ui] = algorithm2(nodes, Nei_agent, num_nodes, r, d, epsilon, p_nodes, curActions);
        p_nodes = (nodes - nodes_old)/delta_t; % set velocities
        nodes_old = nodes; % set old nodes
        nodes = nodes_old + p_nodes*delta_t  + Ui*delta_t*delta_t /2; % get new pos of nodes
        center_of_mass(iterator,:) = mean(nodes); % COM of nodes
        
        %Plot nodes
        plot(center_of_mass(:,1),center_of_mass(:,2),'ro','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k','MarkerSize',4.2)
        hold on
        q_nodes_all{iterator} = nodes; %q_nodes_all{iterator}(:,:) = nodes;
        
        Connectivity(iterator)= (1/(num_nodes))*rank(A);
        [Nei_agent, A] = getNeighbors(nodes, r);   %Determine neighbors for each node after taking action
        
        R_sum_all(iterator) = 0;   %Initialize sum of rewards
        sum_Delta_Q_row = 1:num_nodes; %Initialize sum of mean delta Q for a row of Q table
        
        %Observe new state and reward after taking action
        for i = 1:num_nodes
            nextStates(i) = discreteStateForm(nodes, i, Nei_agent{i});    %Update next state
            
            %Compute reward
            if length(Nei_agent{i}) <= 6 %&& Connectivity(iterator) ~= 1
                reward = length(Nei_agent{i});
            else
                reward = 6;
            end
            
            %Save reward values
            R_nodes(iterator, i) = reward;
            R_sum_all(iterator) = R_sum_all(iterator) + reward;
            %
            
            %determine next action
            nextActions(i) = epsilonGreedy(Q_update{i}, nextStates(i), epsilon_learning);
            
            %Save next action
            A_sum_cooQ(iterator) = A_sum_cooQ(iterator) + nextActions(i);  
            
            %store previous Q value
            Q_prev{i} = Q_update{i};
                
            %Compute Q_update value
            qUpdate(i, w, curStates, curActions, nextActions, Q_update, reward, alpha, gamma, Nei_agent, nstates, nactions, learningType);            
            % 
            
            %Set next action as current action
            curActions(i) = nextActions(i);
            
            %Set next state as current state
            curStates(i) = nextStates(i);
                    
            %Determine Delta Q
            Delta_Q = Q_update{i} - Q_prev{i};
            sum_Delta_Q_col = sum(abs(Delta_Q));
            sum_Delta_Q_row(i) = sum(sum_Delta_Q_col);
        end
        
        
        %get average DeltaQ value
        mean_Delta_Q(iterator) = sum(sum_Delta_Q_row(:))/num_nodes;
        
        % Plot nodes with links
        plot(nodes(:,1),nodes(:,2), '.')
        hold on
        plot(nodes(:,1),nodes(:,2), 'k>','LineWidth',.25,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 4.5)
        hold off
        for node_i = 1:num_nodes
            buffer = nodes(Nei_agent{node_i},:);
            for j = 1:size(nodes(Nei_agent{node_i},1))
                line([nodes(node_i,1),buffer(j,1)],[nodes(node_i,2),buffer(j,2)]) 
            end
        end
        drawnow;
        hold off
        %
    end
end


function value = qUpdate(i, w, curStates, curActions, nextStates, Q_update, reward, alpha, gamma, Nei_agent, nstates, nactions, learningType)

    if (learningType == 0) % ind
        max1 = max(Q_update{i}(nextStates(i),:));% get max of rewards
        Q_update{i}(curStates(i),curActions(i)) = Q_update{i}(curStates(i),curActions(i)) + alpha * (reward + gamma * max1 - Q_update{i}(curStates(i),curActions(i)));% update Q table
    end         
            
    if (learningType == 1) % coop
            max1 = max(Q_update{i}(nextStates(i),:));% get max of rewards
            Qval = Q_update{i}(curStates(i),curActions(i)) + alpha * (reward + gamma*max1 - Q_update{i}(curStates(i),curActions(i)));
            Qi = Q_update{i};
            Qi(curStates(i),curActions(i)) = Qval;

            Q_neigh_sum = zeros(nstates, nactions);
            for j = 1:length(Nei_agent{i}) % update neighbors Q table
                Q_neigh_sum = Q_neigh_sum + Q_update{Nei_agent{i}(j)};
            end
            Q_update{i} = 0.01 * ((w * Qi) + ((1-w)*Q_neigh_sum));  %Update actual Q table
    end   

end