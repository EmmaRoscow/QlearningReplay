
colour_scheme;

rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
policies = {'random', 'expectedReward', 'RPE', 'RPEProb'};
shuffled = false;

% Create plot
[h, pval, ksamples] = plot_reliability_error_over_replays(rat_names, policies, shuffled);

% Adjust size
h.Position(3) = h.Position(3)/3*2;
y = [0.25 1.5]; ylim(y);
set(gca, 'YTick', [y(1) 1 y(2)])

% Plot where significantly below 1
for iPol = 1:length(policies)
    pvalues = pval(iPol, :);
    sig = pvalues < 0.05 & ksamples > 0;
    sig_replays = ksamples(sig);
    evalc(['c = colourscheme.' policies{iPol}]);
    for iSig = 1:length(sig_replays)
        text(sig_replays(iSig), y(2)-0.05*iPol, '*', 'FontSize', 24, 'HorizontalAlignment', 'center', 'Color', c)
    end
end

% Adjust legend position
hLegend = findobj(gcf, 'Type', 'Legend');
hLegend.Location = 'southwest';
hLegend.Position(2) = hLegend.Position(2) - 0.02;