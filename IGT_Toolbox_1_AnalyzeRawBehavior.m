
%%%%%%%%%%%%%%%% IGTtoolbox Results - Behavioral Measures %%%%%%%%%%%%%%%%%
% This script outputs all the behavioral measures which are independent of
% computational modeling (e.g net scores, directed exploration etc.).
% It computes automatically the statistics according to the type of
% analysis performed:
% - takes into account the independent vs dependent nature of your observations
% - takes into account the number of groups/conditions to decide for ex. if
% t-tests or ANOVAs should be run
% - tests whether your data is normally distributed using a lillietest and
% output both parametric and nonparametric statistics in any case.
% It also plots the different measurements (plot2svg can be used to obtain
% vectorial images).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Romain Ligneul / romain.ligneul@gmail.com / Dec 2018.
% v1.2

% clean the environment
clear all; close all;clc;
% get useful function
addpath(genpath('Tools/OTHERS'));
% load relevant dataset (matching its folder name in IGTdata/)
model_dir = [uigetdir('IGTdata') '/'];

load([model_dir 'IGTdata.mat']);

% set color map
colormap_custom = 'brewer1';

%% Loop over subjects and compute behavioral measures
% in order to avoid polluting the workspace with useless variable, all work
% is saved in structure

for s = 1:length(data)
    
    disp(['subject n°' num2str(s)])
    
    % subject id and condition id
    subj_id(s,1) = s;
    cond_id(s,1) = data{s}.cond;
    try
        cond_label{s,1} = data{s}.cond_label;
    catch
        cond_label{s,1} = num2str(cond_id(s,1));
    end
    
    % compute net score
    netscore_time(s,:) = cumsum(ismember(data{s}.deck, [3 4]))-cumsum(ismember(data{s}.deck, [1 2]));
    netscore_total(s,:) = netscore_time(s,end);
    
    for d=1:4
        times_chosen(d,:) = cumsum(data{s}.deck==d);
        dumpayoff = cumsum(double(data{s}.deck==d).*(data{s}.win-data{s}.lose))';
        avg_payoff(d,:) = dumpayoff./times_chosen(d,:);
    end
    avg_payoff = [[NaN;NaN;NaN;NaN],avg_payoff(:,1:end-1)];
    times_chosen = [[0;0;0;0],times_chosen(:,1:end-1)];

    % compute win-stay lose-shift metrics
    stay = data{s}.deck(2:end)==data{s}.deck(1:end-1);
    prvloss = abs(data{s}.lose(1:end-1))>0;
    WSLS(s,:) = [mean(stay(~prvloss)) mean(~stay(prvloss))];
    clear stay prvloss
    
    % evaluate repeat after bigpun
    bigpun_ind = find(data{s}.lose<=-1000);
    bigpun_number(s,1) = numel(bigpun_ind);
    if ~isempty(bigpun_ind) && bigpun_ind(end)==100
        bigpun_ind(end)=[];
    end
    if ~isempty(bigpun_ind)
        repeat_after_bigpun(s,1) = mean(double(data{s}.deck(bigpun_ind)==data{s}.deck(bigpun_ind+1)));
    else
        repeat_after_bigpun(s,1) = NaN;
    end
    clear bigpun_ind
    
    % compute choice entropy and mutual information between successive
    % choices
    H_total(s,1) = h(data{s}.deck);
    MI2choices_total(s,1) = mi(data{s}.deck(1:end-1),data{s}.deck(2:end));
    
    % compute directed exploration indexes
    for t = 1:98
        directed_exploration3(s,t) = double(numel(unique(data{s}.deck(t:t+2)))==3);
        if t<98
            directed_exploration4(s,t) = double(numel(unique(data{s}.deck(t:t+3)))==4);
        end
    end
    tt=0;
    for t = 1:4:97
        tt=tt+1;
        directed_exploration4_fixed(s,t) = double(numel(unique(data{s}.deck(t:t+3)))==4);
    end
    % permute to obtain empirical chance levels on exploration indexes
    for r = 1:round(5000/length(data))
        perm_deck = data{s}.deck(randperm(length(data{s}.deck)));
        for t = 1:length(perm_deck)-2
            perm_DE3(r,t) = double(numel(unique(perm_deck(t:t+2)))==3);
        end
        for t = 1:length(perm_deck)-3
            perm_DE4(r,t) = double(numel(unique(perm_deck(t:t+3)))==4);
        end
        perm_deck = data{s}.deck(randperm(length(data{s}.deck)));        tt=0;
        for t = 1:4:length(perm_deck)-2
            tt=tt+1;
            perm_DE3_fixed(r,tt) = double(numel(unique(perm_deck(t:t+2)))==3);
        end
        tt=0;
        for t = 1:length(perm_deck)-3
            tt=tt+1;
            perm_DE4_fixed(r,tt) = double(numel(unique(perm_deck(t:t+3)))==4);
        end
    end
    chance_DE3(s) = mean2(perm_DE3);
    chance_DE4(s) = mean2(perm_DE4);
    chance_DE3_fixed(s) = mean2(perm_DE3_fixed);
    chance_DE4_fixed(s) = mean2(perm_DE4_fixed);
end
nsubj = s;
ntrial = 100;

%% statistics
% test effects on parameters
[condition_list, ia] = unique(cond_label);
condition_number = cond_id(ia);
% ask user to select conditions to use for statistical inference
condition_selection = listdlg('PromptString', 'Select the conditions you want to include in the analysis:', ...
    'ListString', condition_list, 'ListSize', [300 200]);
condition_ids = cellstr(char(condition_list(condition_selection)));
condition_number = condition_number(condition_selection);
% ask user to define test type
if numel(condition_ids)>1
    dumstr = questdlg('Are the conditions reflecting repeated measures or independent observation units?', 'Test type', 'Repeated', 'Independent', 'Independent');
    test_type = strmatch(dumstr, {'Repeated', 'Independent'});
    % ask user to chose what to plot
    dumstr = questdlg('Plot all or just selected condtions/groups?', 'Plot type', 'All', 'Selected', 'Selected');
    plot_type = strmatch(dumstr, {'All', 'Selected'});
else
    test_type=0;
    plot_type=1;
end

%%% get means, std and normality tests of relevant variables
gmat = [];
netscore_ymat = [];
WS_ymat = []; LS_ymat = [];
MI_successive_choices_ymat = [];
H_choices_ymat = [];
DE3_ymat = [];
DE4_ymat = [];
gmat = [];
% loop over conditions
for c = 1:length(condition_selection);
    gmat = [gmat; condition_selection(c)*ones(sum(condition_number(c)==cond_id),1)];
    % net score
    stats_all.netscore.mat{c} = netscore_total(condition_number(c)==cond_id);
    stats_all.netscore.normality(c) = 1-lillietest(stats_all.netscore.mat{c});
    netscore_ymat = [netscore_ymat; stats_all.netscore.mat{c}];
    % win stay
    stats_all.WS.mat{c} = WSLS(condition_number(c)==cond_id,1);
    stats_all.WS.normality(c) = 1-lillietest(stats_all.WS.mat{c});
    WS_ymat = [WS_ymat; stats_all.WS.mat{c}];
    % lose shift
    stats_all.LS.mat{c} = WSLS(condition_number(c)==cond_id,2);
    stats_all.LS.normality(c) = 1-lillietest(stats_all.LS.mat{c});
    LS_ymat = [LS_ymat; stats_all.LS.mat{c}];
    % MI successive choices
    stats_all.MI_successive_choices.mat{c} = MI2choices_total(condition_number(c)==cond_id);
    stats_all.MI_successive_choices.normality(c) = 1-lillietest(stats_all.MI_successive_choices.mat{c});
    MI_successive_choices_ymat = [MI_successive_choices_ymat; stats_all.MI_successive_choices.mat{c}];
    % choice entropy total
    stats_all.H_choices.mat{c} = H_total(condition_number(c)==cond_id);
    stats_all.H_choices.normality(c) = 1-lillietest(stats_all.H_choices.mat{c});
    H_choices_ymat = [H_choices_ymat; stats_all.H_choices.mat{c}];
    % directed exploration 3
    stats_all.DE3.mat{c} = nanmean(directed_exploration3(condition_number(c)==cond_id,:),2);
    stats_all.DE4.normality(c) = 1-lillietest(stats_all.DE3.mat{c});
    DE3_ymat = [DE3_ymat; stats_all.DE3.mat{c}];
    % directed exploration 4
    stats_all.DE4.mat{c} = nanmean(directed_exploration4(condition_number(c)==cond_id,:),2);
    stats_all.DE4.normality(c) = 1-lillietest(stats_all.DE4.mat{c});
    DE4_ymat = [DE4_ymat; stats_all.DE4.mat{c}];
    
end
stats_all.gmat = gmat;
stats_all.netscore.means = meanbycond(netscore_ymat, gmat,[]);
stats_all.netscore.stds = stdbycond(netscore_ymat, gmat,[]);
stats_all.WS.means = meanbycond(WS_ymat, gmat,[]);
stats_all.WS.stds = stdbycond(WS_ymat, gmat,[]);
stats_all.LS.means = meanbycond(LS_ymat, gmat,[]);
stats_all.LS.stds = stdbycond(LS_ymat, gmat,[]);
stats_all.MI_successive_choices.means = meanbycond(MI_successive_choices_ymat, gmat,[]);
stats_all.MI_successive_choices.stds = stdbycond(MI_successive_choices_ymat, gmat,[]);
stats_all.H_choices.means = meanbycond(H_choices_ymat, gmat,[]);
stats_all.H_choices.stds = stdbycond(H_choices_ymat, gmat,[]);
stats_all.DE3.means = meanbycond(DE3_ymat, gmat,[]);
stats_all.DE3.stds = stdbycond(DE3_ymat, gmat,[]);
stats_all.DE4.means = meanbycond(DE4_ymat, gmat,[]);
stats_all.DE4.stds = stdbycond(DE4_ymat, gmat,[]);

% report the type of tests used
if test_type==2
    stats_all.nonparam_tests = 'Wilcoxon unpaired (2 groups) or Kruskall-Wallis (>2 groups)';
    stats_all.param_tests = 't-test unpaired (2 groups) or anova1 (>2 groups)';
else
    stats_all.netscore.nonparam_tests = 'Wilcoxon paired (2 groups) or Friedmann (>2 groups)';
    stats_all.netscore.param_tests = 't-test paired (2 groups) or repeated measure ANOVA (>2 groups)';
end

%%% compute statistics on all measures of interest
if test_type==2
    if length(condition_number)==2
        % netscore
        [stats_all.netscore.nonparam_p, ~, stats_all.netscore.nonparam_stats] = ranksum(netscore_ymat(gmat==1),netscore_ymat(gmat==2));
        [~, stats_all.netscore.param_p,  ~, stats_all.netscore.param_stats] = ttest2(netscore_ymat(gmat==1),netscore_ymat(gmat==2));
        % WS
        [stats_all.WS.nonparam_p, ~, stats_all.WS.nonparam_stats] = ranksum(WS_ymat(gmat==1),WS_ymat(gmat==2));
        [~, stats_all.WS.param_p, stats_all.WS.param_stats] = ttest2(WS_ymat(gmat==1),WS_ymat(gmat==2));
        % LS
        [stats_all.LS.nonparam_p, ~, stats_all.LS.nonparam_stats] = ranksum(LS_ymat(gmat==1),LS_ymat(gmat==2));
        [~, stats_all.LS.param_p,  ~, stats_all.LS.param_stats] = ttest2(LS_ymat(gmat==1),LS_ymat(gmat==2));
        % MI_successive_choices
        [stats_all.MI_successive_choices.nonparam_p, ~, stats_all.MI_successive_choices.nonparam_stats] = ranksum(MI_successive_choices_ymat(gmat==1),MI_successive_choices_ymat(gmat==2));
        [~, stats_all.MI_successive_choices.param_p, ~,  stats_all.MI_successive_choices.param_stats] = ttest2(MI_successive_choices_ymat(gmat==1),MI_successive_choices_ymat(gmat==2));
        % H_choices
        [stats_all.H_choices.nonparam_p, ~, stats_all.H_choices.nonparam_stats] = ranksum(H_choices_ymat(gmat==1),H_choices_ymat(gmat==2));
        [~, stats_all.H_choices.param_p, ~,  stats_all.H_choices.param_stats] = ttest2(H_choices_ymat(gmat==1),H_choices_ymat(gmat==2));
        % DE3
        [stats_all.DE3.nonparam_p, ~, stats_all.DE3.nonparam_stats] = ranksum(DE3_ymat(gmat==1),DE3_ymat(gmat==2));
        [~, stats_all.DE3.param_p, stats_all.DE3.param_stats] = ttest2(DE3_ymat(gmat==1),DE3_ymat(gmat==2));
        % DE4
        [stats_all.DE4.nonparam_p, ~, stats_all.DE4.nonparam_stats] = ranksum(DE4_ymat(gmat==1),DE4_ymat(gmat==2));
        [~, stats_all.DE4.param_p, ~, stats_all.DE4.param_stats] = ttest2(DE4_ymat(gmat==1),DE4_ymat(gmat==2));
    elseif length(condition_number)>2
        % netscore
        [stats_all.netscore.param_p, stats_all.netscore.param_table, stats_all.netscore.stats] = anova1(netscore_ymat, gmat, 'off');
        [stats_all.netscore.nonparam_p,  stats_all.netscore.nonparam_table, stats_all.netscore.nonparam_stats] = kruskalwallis(netscore_ymat, gmat, 'off');
        % WS
        [stats_all.WS.param_p, stats_all.WS.param_table, stats_all.WS.stats] = anova1(WS_ymat, gmat, 'off');
        [stats_all.WS.nonparam_p,  stats_all.WS.nonparam_table, stats_all.WS.nonparam_stats] = kruskalwallis(WS_ymat, gmat, 'off');
        % LS
        [stats_all.LS.param_p, stats_all.LS.param_table, stats_all.LS.stats] = anova1(LS_ymat, gmat, 'off');
        [stats_all.LS.nonparam_p,  stats_all.LS.nonparam_table, stats_all.LS.nonparam_stats] = kruskalwallis(LS_ymat, gmat, 'off');
        % MI_successive_choices
        [stats_all.MI_successive_choices.param_p, stats_all.MI_successive_choices.param_table, stats_all.MI_successive_choices.stats] = anova1(MI_successive_choices_ymat, gmat, 'off');
        [stats_all.MI_successive_choices.nonparam_p,  stats_all.MI_successive_choices.nonparam_table, stats_all.MI_successive_choices.nonparam_stats] = kruskalwallis(MI_successive_choices_ymat, gmat, 'off');
        % H_choices
        [stats_all.H_choices.param_p, stats_all.H_choices.param_table, stats_all.H_choices.stats] = anova1(H_choices_ymat, gmat, 'off');
        [stats_all.H_choices.nonparam_p,  stats_all.H_choices.nonparam_table, stats_all.H_choices.nonparam_stats] = kruskalwallis(H_choices_ymat, gmat, 'off');
        % DE3
        [stats_all.DE3.param_p, stats_all.DE3.param_table, stats_all.DE3.stats] = anova1(DE3_ymat, gmat, 'off');
        [stats_all.DE3.nonparam_p,  stats_all.DE3.nonparam_table, stats_all.DE3.nonparam_stats] = kruskalwallis(DE3_ymat, gmat, 'off');
        % DE4
        [stats_all.DE4.param_p, stats_all.DE4.param_table, stats_all.DE4.stats] = anova1(DE4_ymat, gmat, 'off');
        [stats_all.DE4.nonparam_p,  stats_all.DE4.nonparam_table, stats_all.DE4.nonparam_stats] = kruskalwallis(DE4_ymat, gmat, 'off');
    else
        warning('No statistics computed: only one group or condition was defined!');
    end
else
    if length(condition_number)==2
        % netscore
        [stats_all.netscore.nonparam_p, ~, stats_all.netscore.nonparam_stats] = signrank(netscore_ymat(gmat==1),netscore_ymat(gmat==2));
        [~, stats_all.netscore.param_p, stats_all.netscore.param_stats] = ttest(netscore_ymat(gmat==1),netscore_ymat(gmat==2));
        % WS
        [stats_all.WS.nonparam_p, ~, stats_all.WS.nonparam_stats] = signrank(WS_ymat(gmat==1),WS_ymat(gmat==2));
        [~, stats_all.WS.param_p, stats_all.WS.param_stats] = ttest(WS_ymat(gmat==1),WS_ymat(gmat==2));
        % LS
        [stats_all.LS.nonparam_p, ~, stats_all.LS.nonparam_stats] = signrank(LS_ymat(gmat==1),LS_ymat(gmat==2));
        [~, stats_all.LS.param_p, stats_all.LS.param_stats] = ttest(LS_ymat(gmat==1),LS_ymat(gmat==2));
        % MI_successive_choices
        [stats_all.MI_successive_choices.nonparam_p, ~, stats_all.MI_successive_choices.nonparam_stats] = signrank(MI_successive_choices_ymat(gmat==1),MI_successive_choices_ymat(gmat==2));
        [~, stats_all.MI_successive_choices.param_p, stats_all.MI_successive_choices.param_stats] = ttest(MI_successive_choices_ymat(gmat==1),MI_successive_choices_ymat(gmat==2));
        % H_choices
        [stats_all.H_choices.nonparam_p, ~, stats_all.H_choices.nonparam_stats] = signrank(H_choices_ymat(gmat==1),H_choices_ymat(gmat==2));
        [~, stats_all.H_choices.param_p, stats_all.H_choices.param_stats] = ttest(H_choices_ymat(gmat==1), H_choices_ymat(gmat==2));
        % DE3
        [stats_all.DE3.nonparam_p, ~, stats_all.DE3.nonparam_stats] = signrank(DE3_ymat(gmat==1),DE3_ymat(gmat==2));
        [~, stats_all.DE3.param_p, stats_all.DE3.param_stats] = ttest(DE3_ymat(gmat==1), DE3_ymat(gmat==2));
        % DE4
        [stats_all.DE4.nonparam_p, ~, stats_all.DE4.nonparam_stats] = signrank(DE4_ymat(gmat==1),DE4_ymat(gmat==2));
        [~, stats_all.DE4.param_p, stats_all.DE4.param_stats] = ttest(DE4_ymat(gmat==1), DE4_ymat(gmat==2));
    elseif length(condition_number)>2
        % netscore
        [stats_all.netscore.param_p, stats_all.netscore.param_table] = anova_rm(cell2mat(stats_all.netscore.mat), 'off');
        [stats_all.netscore.nonparam_p, stats_all.netscore.nonparam_table, stats_all.netscore.nonparam_stats] = friedman(cell2mat(stats_all.netscore.mat{c}));
        % WS
        [stats_all.WS.param_p, stats_all.WS.param_table] = anova_rm(cell2mat(stats_all.WS.mat), 'off');
        [stats_all.WS.nonparam_p, stats_all.WS.nonparam_table, stats_all.WS.nonparam_stats] = friedman(cell2mat(stats_all.WS.mat{c}));
        % LS
        [stats_all.LS.param_p, stats_all.LS.param_table] = anova_rm(cell2mat(stats_all.LS.mat), 'off');
        [stats_all.LS.nonparam_p, stats_all.LS.nonparam_table, stats_all.LS.nonparam_stats] = friedman(cell2mat(stats_all.LS.mat{c}));
        % MI_successive_choices
        [stats_all.MI_successive_choices.param_p, stats_all.MI_successive_choices.param_table] = anova_rm(cell2mat(stats_all.MI_successive_choices.mat), 'off');
        [stats_all.MI_successive_choices.nonparam_p, stats_all.MI_successive_choices.nonparam_table, stats_all.MI_successive_choices.nonparam_stats] = friedman(cell2mat(stats_all.MI_successive_choices.mat{c}));
        % H_choices
        [stats_all.H_choices.param_p, stats_all.H_choices.param_table] = anova_rm(cell2mat(stats_all.H_choices.mat{c}), 'off');
        [stats_all.H_choices.nonparam_p, stats_all.H_choices.nonparam_table, stats_all.H_choices.nonparam_stats] = friedman(cell2mat(stats_all.H_choices.mat));
        % DE3
        [stats_all.DE3.param_p, stats_all.DE3.param_table] = anova_rm(cell2mat(stats_all.DE3.mat{c}), 'off');
        [stats_all.DE3.nonparam_p, stats_all.DE3.nonparam_table, stats_all.DE3.nonparam_stats] = friedman(cell2mat(stats_all.DE3.mat));
        % DE4
        [stats_all.DE4.param_p, stats_all.DE4.param_table] = anova_rm(cell2mat(stats_all.DE4.mat{c}), 'off');
        [stats_all.DE4.nonparam_p, stats_all.DE4.nonparam_table, stats_all.DE4.nonparam_stats] = friedman(cell2mat(stats_all.DE4.mat));
    else
        warning('No statistics computed: only one group or condition was defined!');
    end
end

% log test type
if test_type==2
    stats_all.nonparam_tests = 'Wilcoxon unpaired (2 groups) or Kruskall-Wallis (>2 groups)';
    stats_all.param_tests = 't-test unpaired (2 groups) or anova1 (>2 groups)';
else
    stats_all.netscore.nonparam_tests = 'Wilcoxon paired (2 groups) or Friedmann (>2 groups)';
    stats_all.netscore.param_tests = 't-test paired (2 groups) or repeated measure ANOVA (>2 groups)';
end

%% plot and perform statistical test on NET scores
% in time
figure('Name', 'Net Score','units','normalized','position',[0 0 .45 .45]);
clear g;
x=repmat(1:ntrial,nsubj,1);
y=netscore_time;
color=repmat(cond_label,1,ntrial);
if plot_type==1
    g(1,1) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(1,1) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(1,1).stat_summary('type', 'sem', 'geom', 'area', 'setylim', true);
g(1,1).set_title('NET score in time');
g(1,1).set_names('x', 'trial', 'y', 'NET score');
g(1,1).set_color_options('map', colormap_custom);
% average
y=netscore_total;
color=cond_label;
x = ones(size(color,1), size(color,2));
if plot_type==1
    g(2,1) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(2,1) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(2,1).stat_summary('type', 'sem', 'geom', 'bar', 'dodge', 0.7, 'width', 0.6);
g(2,1).stat_summary('type', 'sem', 'geom', 'errorbar', 'dodge', 0.7, 'width', 0.6, 'setylim', true);
g(2,1).axe_property('Xlim', [0 2], 'xtick', [], 'xticklabel', []); ;
g(2,1).set_title('Final NET score');
g(2,1).set_names('x', 'group', 'y', 'NET score');
g(2,1).set_color_options('map', colormap_custom);
% display panels
g.draw();

%% plot and perform statistical test on other behavioral measures
% in time
figure('Name', 'Other behavioral measures','units','normalized','position',[0.55 0 .45 .45]);
clear g;
% average Win Stay
y=WSLS(:,1);
color=cond_label;
x = ones(size(color,1), size(color,2));
if plot_type==1
    g(1,1) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(1,1) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(1,1).stat_summary('type', 'sem', 'geom', 'bar', 'dodge', 0.7, 'width', 0.6);
g(1,1).stat_summary('type', 'sem', 'geom', 'errorbar', 'dodge', 0.7, 'width', 0.6, 'setylim', true);
g(1,1).axe_property('Xlim', [0 2], 'xtick', [], 'xticklabel', []); ;
g(1,1).set_title('Win Stay');
g(1,1).set_names('x', 'group', 'y', 'frequency');
g(1,1).set_color_options('map', colormap_custom);
% average Win Stay
y=WSLS(:,2);
color=cond_label;
x = ones(size(color,1), size(color,2));
if plot_type==1
    g(1,2) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(1,2) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(1,2).stat_summary('type', 'sem', 'geom', 'bar', 'dodge', 0.7, 'width', 0.6);
g(1,2).stat_summary('type', 'sem', 'geom', 'errorbar', 'dodge', 0.7, 'width', 0.6, 'setylim', true);
g(1,2).axe_property('Xlim', [0 2], 'xtick', [], 'xticklabel', []); ;
g(1,2).set_title('Lose Shift');
g(1,2).set_names('x', 'group', 'y', 'frequency');
g(1,2).set_color_options('map', colormap_custom);
% choice entropy
y=H_total;
color=cond_label;
x = ones(size(color,1), size(color,2));
if plot_type==1
    g(2,1) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(2,1) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(2,1).stat_summary('type', 'sem', 'geom', 'bar', 'dodge', 0.7, 'width', 0.6);
g(2,1).stat_summary('type', 'sem', 'geom', 'errorbar', 'dodge', 0.7, 'width', 0.6, 'setylim', true);
g(2,1).axe_property('Xlim', [0 2], 'xtick', [], 'xticklabel', []); ;
g(2,1).set_title('H(choice)');
g(2,1).set_names('x', 'group', 'y', 'bits');
g(2,1).set_color_options('map', colormap_custom);
% mutual information of successive choices (perseveration measure)
y=MI2choices_total;
color=cond_label;
x = ones(size(color,1), size(color,2));
if plot_type==1
    g(2,2) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(2,2) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(2,2).stat_summary('type', 'sem', 'geom', 'bar', 'dodge', 0.7, 'width', 0.6);
g(2,2).stat_summary('type', 'sem', 'geom', 'errorbar', 'dodge', 0.7, 'width', 0.6, 'setylim', true);
g(2,2).axe_property('Xlim', [0 2], 'xtick', [], 'xticklabel', []); ;
g(2,2).set_title('MI(t,t+1)');
g(2,2).set_names('x', 'group', 'y', 'bits');
g(2,2).set_color_options('map', colormap_custom);
% draw
g.draw();

%% plot and perform statistical test on directed exploration
% in time
figure('Name', 'Directed exploration','units','normalized','position',[0 0.45 .45 .45]);
clear g;
x=repmat(1:ntrial-2,nsubj,1);
y=directed_exploration3;
color=repmat(cond_label,1,ntrial-2);
if plot_type==1
    g(1,1) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(1,1) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(1,1).stat_summary('type', 'sem', 'geom', 'area', 'setylim', true);
g(1,1).set_title('DE3 in time');
g(1,1).set_names('x', 'trial', 'y', 'pattern frequency');
g(1,1).set_color_options('map', colormap_custom);
g(1,1).geom_hline('yintercept',mean(chance_DE3));
% average
y=nanmean(directed_exploration3,2);
color=cond_label;
x=ones(size(cond_label));
if plot_type==1
    g(1,2) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(1,2) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(1,2).stat_summary('type', 'sem', 'geom', 'bar', 'dodge', 0.7, 'width', 0.6);
g(1,2).stat_summary('type', 'sem', 'geom', 'errorbar', 'dodge', 0.7, 'width', 0.6, 'setylim', true);
g(1,2).axe_property('Xlim', [0 2], 'xtick', [], 'xticklabel', []); ;
g(1,2).set_title('DE3 average');
g(1,2).set_names('x', 'group', 'y', 'frequency');
g(1,2).set_color_options('map', colormap_custom);
g(1,2).geom_hline('yintercept',mean(chance_DE3));
% DE4 in time
x=repmat(1:ntrial-3,nsubj,1);
y=directed_exploration4;
color=repmat(cond_label,1,ntrial-3);
if plot_type==1
    g(2,1) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(2,1) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(2,1).stat_smooth();%'type', 'sem', 'geom', 'area', 'setylim', true);
g(2,1).set_title('DE4 in time');
g(2,1).set_names('x', 'trial', 'y', 'pattern frequency');
g(2,1).set_color_options('map', colormap_custom);
g(2,1).geom_hline('yintercept',mean(chance_DE4));
g(2,1).axe_property('ylim',[0 0.4]);
% DE4 average
y=nanmean(directed_exploration4,2);
color=cond_label;
x = ones(size(cond_label));
if plot_type==1
    g(2,2) = gramm('x', x(:), 'y', y(:), 'color', color(:));
else
    subset = repmat(ismember(cond_label, condition_list(condition_selection)),1,size(y,2));
    g(2,2) = gramm('x', x(:), 'y', y(:), 'color', color(:), 'subset', subset(:));
end
g(2,2).stat_summary('type', 'sem', 'geom', 'bar', 'dodge', 0.7, 'width', 0.6);
g(2,2).stat_summary('type', 'sem', 'geom', 'errorbar', 'dodge', 0.7, 'width', 0.6, 'setylim', true);
g(2,2).axe_property('Xlim', [0 2], 'xtick', [], 'xticklabel', []); ;
g(2,2).set_title('DE4 average');
g(2,2).set_names('x', 'group', 'y', 'frequency');
g(2,2).set_color_options('map', colormap_custom);
g(2,2).geom_hline('yintercept',mean(chance_DE4));
% draw
g.draw();

