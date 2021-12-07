clc,clear
close all

%======LEARNING SELECT======
% 0 = individual learning  =
% 1 = cooperative learning =
%                          =
    learningType = 0;
%                          =
%END SELECT                =

maxepisodes = 6; %num episodes
num_nodes = 15;  %Randomly generate nodes
n = 2; %number of dimensions
nodes = 90.*rand(num_nodes,n)+90.*repmat([0 1],num_nodes,1);
statelist   = BuildStateList(num_nodes);  % the list of states
actionlist  = BuildActionList(); % the list of actions
nstates     = size(statelist,1);
nactions    = 4; % dont know why this broke but eh

safePoints = generateSafePoints(nactions);

epsilon_learning     = 0.0001;   % probability of a random action selection for e-greedy policy
delta_t = 0.03; 
t = 0:delta_t:3.1;
Q_initial = load('Qcell_4actions2.mat'); %FOR 4 ACTIONS CASE
Q_update = Q_initial.Q;

%SAVE DATA FOR EVALUATION
Connectivity_episodes = cell(1, maxepisodes );
Connectivity_episodes_learning = cell(1, maxepisodes );
R_ind_episodes = cell(1, maxepisodes);
R_all_episodes = cell(1, maxepisodes );
A_sum_cooQ_episodes = cell(1, maxepisodes );
Topo_eva_all_epi = cell(1, maxepisodes );
mean_Delta_Q_epi = cell(1, maxepisodes );
q_nodes_epi = cell(1, maxepisodes);

for i=1:maxepisodes
    nodes = 90.*rand(num_nodes,n)+90.*repmat([0 1],num_nodes,1);
    %Training
    [Q_update, Connectivity, R_nodes, R_sum_all, A_sum_cooQ, mean_Delta_Q, q_nodes_all]  = Q_Learning(Q_update, statelist, actionlist, nstates, nactions, num_nodes, n, nodes, epsilon_learning, t, delta_t, safePoints, learningType);
    %Save data
    Connectivity_episodes{i} = Connectivity;
    R_ind_episodes{i} = R_nodes;
    R_all_episodes{i} = R_sum_all;
    A_sum_cooQ_episodes{i} = A_sum_cooQ;
    q_nodes_epi{i} = q_nodes_all;
    mean_Delta_Q_epi{i} = mean_Delta_Q; 
    disp(['Espisode: ',int2str(i),]) 
    %decrease probability of a random action selection (for e-greedy selection)
    epsilon_learning = epsilon_learning * 0.99;
end


% part 2 (trajectory plotting) 1,2,last

% first 
q_nodes_mat = cell2mat(q_nodes_epi{1});
figure(2),plot(q_nodes_mat(:,1),q_nodes_mat(:,2),'k.')
hold on
plot(q_nodes_epi{1}{length(t)}(:,1),q_nodes_epi{1}{length(t)}(:,2), 'm>','LineWidth',.25,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize', 4.5)
title('Node trajectory in first episode')

% second
q_nodes_mat = cell2mat(q_nodes_epi{2});
figure(3),plot(q_nodes_mat(:,1),q_nodes_mat(:,2),'k.')
hold on
plot(q_nodes_epi{2}{length(t)}(:,1),q_nodes_epi{2}{length(t)}(:,2), 'm>','LineWidth',.25,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize', 4.5)
title('Node trajectory in second episode')

% last
q_nodes_mat = cell2mat(q_nodes_epi{maxepisodes});
figure(4),plot(q_nodes_mat(:,1),q_nodes_mat(:,2),'k.')
hold on
plot(q_nodes_epi{maxepisodes}{length(t)}(:,1),q_nodes_epi{maxepisodes}{length(t)}(:,2), 'm>','LineWidth',.25,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize', 4.5)
title('Node trajectory in last episode')

% part 3 (individual reward)
R_each_node = zeros(length(t)*maxepisodes, num_nodes);
for ep = 1:maxepisodes
    temp = R_ind_episodes{ep};
   for i = 1:length(t)
      for j = 1:num_nodes
          R_each_node((ep-1)*length(t) + i, j) = temp(i, j);
      end
   end
end
figure(5), plot(R_each_node);
title('Individual Reward over learning episodes')
grid on

% part 4 (sum of rewards)
R_all_epi_mat = cell2mat(R_all_episodes);
[R_all_diff0, index_R]= find(R_all_epi_mat>0);
figure(6), plot(R_all_epi_mat(index_R));
title('Total reward over learning episodes')
grid on

% part 5 (action selection)
A_cooQ_matrix = cell2mat(A_sum_cooQ_episodes);
[A_diff0, index_A_cooQ] = find(A_cooQ_matrix>0);
figure(7), plot(A_cooQ_matrix(index_A_cooQ))
title('Action Selection over learning episodes')
grid on 

% part 6 (average delta_Q)
mean_Delta_Q_mat = cell2mat(mean_Delta_Q_epi);
[delta_Q_diff, index_delta_Q] = find(mean_Delta_Q_mat > 0);
figure(8), plot(mean_Delta_Q_mat(index_delta_Q))
title('Average delta_Q over learning episodes')
grid on

%=========================== END PLOTTING================

function states = BuildStateList(n)
    states = zeros(n*2 +2, 2);
    for row = 1:n*2 +2
       for col = 1:2
           if col == 1
               if mod(row, 2) == 0 % if even
                   states(row,col) = row/2 - 1;
               else % if odd
                   states(row,col) = floor(row/2);
               end
           else
              states(row,col) = mod(row-1, 2);
           end
       end
    end
end


function actions = BuildActionList()
    actions = 1:4; % safe places
end

function safePoints = generateSafePoints(nactions)
    safePoints = [];
    for actions = 1:nactions
        safePoints(actions,:) = pointFromAction(actions);
    end
end
