
% Set Q-values (taken from Harry, start of session 5)
Q = [0.0312, 0.7120, 0.7113;
      0.7787, 0.0000, 0.7041;
      0.7533, 0.7167, 0.0000];

% Set chain of actions (arbitrary)
actions = [1 3 1 2 1 3 1 3 2 1 3];

% Set rewards for actions (arbitrary)
rewards = [0 1 1 1 1 0 0 0 1 0];

% Infer state
states = actions(1:end-1);
actions = actions(2:end);

% Other variables (taken from Harry)
params.alpha = 0.0111;
params.discount = 0.6805;
params.rpe_recency = 0.0142;
params.recency = 7.0509;
sa_pairs = reshape([1:9],3,3);
for i = 1:9
    samples(i).samples = zeros(length(actions), 6);
    samples(i).dn = 1;
end

% Calculate RPE for each trial
for iTrial = 1:length(actions)
    
    a = actions(iTrial);
    s = states(iTrial);
    r = rewards(iTrial);
    s_prime = a;
    sa = sa_pairs(s,a);
    
    % Q-values for possible next state-action pairs
    available_sa_prime = sa_pairs(s_prime,:);
    Qvals = Q(available_sa_prime);
        
    % On-policy most likely next action
     [maxval, ~] = max(Qvals);
     [bestA] = find(Qvals==maxval);
    if length(bestA) > 1
        a_prime_star = bestA(1);
    else
        a_prime_star = bestA;
    end
    sa_prime_star = sa_pairs(s_prime,a_prime_star);
        
    % Reward prediction error for reward on this trial
    rpe(iTrial) = r + params.discount*Q(sa_prime_star) - Q(sa);
    
    % Update Q-value for state-action pair
    Q(sa) = (1-params.alpha)*Q(sa) + params.alpha*(r + params.discount*Q(sa_prime_star));
    
end

% Add trials to memory buffer
for iTrial = 1:length(actions)
    
    a = actions(iTrial);
    s = states(iTrial);
    r = rewards(iTrial);
    s_prime = a;
    sa = sa_pairs(s,a);
    RPE = rpe(iTrial);
    
    % Add current trial to bank of samples
    samples(sa).samples(samples(sa).dn,:) = [s, a, r, s_prime, RPE, iTrial];
    
    % Record how many times this state-action pair has been experienced
    samples(sa).dn = samples(sa).dn + 1;
    
end

% Calculate probability of replay for each state-action pair
u_sa = find([samples.dn] > 1);

% Random replay: all state-action pairs in memory buffer equal
sa_replay_prob.random = zeros(1, 9);
sa_replay_prob.random(u_sa) = 1/length(u_sa);

% Reward-biased replay: in proportion to Q-value
sa_replay_prob.rewardbiased = zeros(1, 9);
sa_replay_prob.rewardbiased(u_sa) = Q(u_sa)/sum(Q(u_sa));

% RPE-prioritised replay: weight state-action pairs in memory buffer by
% their RPE and replay only the top one
sa_replay_prob.RPEprioritised = zeros(1, 9);
rpe_weighting = NaN(1, 9);
for i = 1:length(u_sa)
    % Samples from this SA pair
    sample_list = samples(u_sa(i)).samples(find(samples(u_sa(i)).samples(:,1)),:);
    % Take recency-weighted mean of RPE
    rpe_weighting(u_sa(i)) = mean(abs(sample_list(:,5)).*flip((params.rpe_recency)).^(1./(1:(size(sample_list,1)))'));
end 
% Choose transition with greatest recent RPE to replay
[~, max_rpe] = max(abs(rpe_weighting));
sa_replay_prob.RPEprioritised(max_rpe) = 1;

% RPE-proportional replay: weight state-action pairs in memory buffer by
% their RPE and take probabilities in proportion to this
sa_replay_prob.RPEproportional = zeros(1, 9);
rpe_weighting = zeros(1, 9);
for i = 1:length(u_sa)
    % Samples from this SA pair
    sample_list = samples(u_sa(i)).samples(find(samples(u_sa(i)).samples(:,1)),:);
    % Take recency-weighted mean of RPE
    rpe_weighting(u_sa(i)) = mean(abs(sample_list(:,5)).*flip((params.rpe_recency)).^(1./(1:(size(sample_list,1)))'));
end 
% Choose transition with greatest recent RPE to replay
sa_replay_prob.RPEproportional = rpe_weighting / sum(rpe_weighting);

% Calculate probability of replay for each trial
trial_replay_prob.random = zeros(1, length(actions));
trial_replay_prob.rewardbiased = zeros(1, length(actions));
trial_replay_prob.RPEprioritised = zeros(1, length(actions));
trial_replay_prob.RPEproportional = zeros(1, length(actions));
for sa = 1:9
    
    % Find trials that belong to the state-action pair
    sample_list = samples(sa).samples(find(samples(sa).samples(:,1)),:);
    
    % Weight them by recency
    pdf = sample_list(:, 6) .^ params.recency / sum(sample_list(:, 6) .^ params.recency);
    
    trial_replay_prob.random(sample_list(:, 6)) = sa_replay_prob.random(sa) * pdf;
    trial_replay_prob.rewardbiased(sample_list(:, 6)) = sa_replay_prob.rewardbiased(sa) * pdf;
    trial_replay_prob.RPEprioritised(sample_list(:, 6)) = sa_replay_prob.RPEprioritised(sa) * pdf;
    trial_replay_prob.RPEproportional(sample_list(:, 6)) = sa_replay_prob.RPEproportional(sa) * pdf;

end

% Print table showing results for every trial
behaviour = table([1:length(actions)]', states', actions', rewards', rpe', 'VariableNames', {'Trial', 'State', 'Action', 'Reward', 'RPE'})

% Print table showing probability of replay for every state-action pair
sa_state = {'High', 'High', 'High', 'Medium', 'Medium', 'Medium', 'Low', 'Low', 'Low'}';
sa_action = {'High', 'Medium', 'Low', 'High', 'Medium', 'Low', 'High', 'Medium', 'Low'}';
replay_prob = table(sa_state, sa_action, [samples(:).dn]'-1, Q(:), rpe_weighting',...
    sa_replay_prob.random', sa_replay_prob.rewardbiased', sa_replay_prob.RPEprioritised', sa_replay_prob.RPEproportional',...
    'VariableNames', {'State', 'Action', 'Number_trials', 'Q_value', 'Recency_weighted_mean_absolute_RPE',...
    'Random_replay', 'Reward_biased_replay', 'RPE_prioritised_replay', 'RPE_proportional_replay'})

% Plot probability of replay of each state-action pair
colour_scheme;
axes_properties;
figure('Position', [680 558 680 420]); hold on
b = bar([sa_replay_prob.random' sa_replay_prob.rewardbiased' sa_replay_prob.RPEprioritised' sa_replay_prob.RPEproportional'], 'FaceColor', colourscheme.random);
b(1).FaceColor = colourscheme.random;
b(2).FaceColor = colourscheme.expectedReward;
b(3).FaceColor = colourscheme.RPE;
b(4).FaceColor = colourscheme.RPEProb;
for i = 1:4
    b(i).EdgeColor = 'none';
end

% Formatting
xlim([0.5 9.5])
xlabel('State-action pair')
ylabel('Probability of replay')
l = legend({'Random', 'Reward-biased', 'RPE-prioritised', 'RPE-proportional'});
l.Location = 'northwest';
l.Title.String = 'Replay policy';
l.Box = 'off';
set(gca, 'YTIck', [0 1]); ylim([0 1])

% Xtick labels
arms = {'high', 'mid', 'low'};
ticklabels = {};
for i = 1:3
    for j = 1:3
        ticklabels{(i-1)*3+j} = sprintf([arms{j} '\\newline' arms{i}]);
    end
end
set(gca, 'XTick', 1:9)
set(gca, 'XTickLabel', ticklabels)

% Plot probability of replay of each trial
colour_scheme;
axes_properties;
figure('Position', [680 558 680 420]); hold on
b = bar([trial_replay_prob.random' trial_replay_prob.rewardbiased' trial_replay_prob.RPEprioritised' trial_replay_prob.RPEproportional'], 'FaceColor', colourscheme.random);
b(1).FaceColor = colourscheme.random;
b(2).FaceColor = colourscheme.expectedReward;
b(3).FaceColor = colourscheme.RPE;
b(4).FaceColor = colourscheme.RPEProb;
for i = 1:4
    b(i).EdgeColor = 'none';
end

% Formatting
xlim([0.5 length(actions)+0.5])
xlabel('Trial number')
ylabel('Probability of replay')
l = legend({'Random', 'Reward-biased', 'RPE-prioritised', 'RPE-proportional'});
l.Location = 'northwest';
l.Title.String = 'Replay policy';
l.Box = 'off';
set(gca, 'YTIck', [0 1]); ylim([0 1])
