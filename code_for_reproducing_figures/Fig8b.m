
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
policies = {'noReplay', 'random', 'expectedReward', 'RPE', 'RPEProb'};
shuffled = true;
n_replays = 15;

% Plot
[h] = plot_reliability_error_over_sessions(rat_names, policies, shuffled, n_replays);

% Remove legend
hLegend = findobj(gcf, 'Type', 'Legend');
hLegend.Visible = 'off';