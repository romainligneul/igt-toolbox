%%%%%%%%%%%%%%%%%%%%%%% IGTtoolbox Results - Modeling %%%%%%%%%%%%%%%%%%%%%
% This script outputs all the behavioral measures which are dependent of
% computational modeling (e.g best fitting parameters, recovery analyses,
% model comparison, etc.).
% It computes automatically the statistics according to the type of
% analysis performed:
% - takes into account the independent vs dependent nature of your observations
% - takes into account the number of groups/conditions to decide for ex. if
% t-tests or ANOVAs should be run
% - tests whether your data is normally distributed using a lillietest and
% output both parametric and nonparametric statistics in any case.
% It also plots everything (plot2svg can be used to obtain vectorial images).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Romain Ligneul / romain.ligneul@gmail.com / August 2018.
% v1.1

% clean the environment
clear all;close all;clc;

% ask user for directory (what should be analyzed)
model_dir = [uigetdir('IGTmodelrecovery', 'Select folder containing the models of interest') '/'];

% get all models
dummy_list = dir(model_dir)
for d = 1:length(dummy_list); model_list{d,1} = dummy_list(d).name; end; model_list(1:2)=[];

% select models to use for the comparison
model_selection = listdlg('PromptString', 'Select the models you want to compare:', ...
    'ListString', model_list, 'ListSize', [300 200], 'InitialValue', 1:length(model_list));
model_list = model_list(model_selection);

% load relevant dataset (matching its folder name in IGTdata/)
data_dir = [uigetdir('IGTdata', 'Select folder corresponding to the modeling analysis') '/'];
load([data_dir 'IGTdata.mat']);

% start looping over models list and logging info
for m = 1:size(model_list,1)
    
    disp(['model ' num2str(m) ' = ' model_list(m)]);

    % log subject id and condition id
    if m == 1
        for s = 1:length(data)
            subj_id(s,1) = s;
            cond_id(s,1) = data{s}.cond;
            try
                cond_label{s,1} = data{s}.cond_label;
            catch
                cond_label{s,1} = num2str(cond_id(s,1));
            end
        end
    end
    
    % list dir inside the model dir
    submodellist = dir([model_dir model_list{m}]);% '/fitted_model.mat']
    submodellist(1:2)=[];
    
    for mm=1:size(model_list,1)
        
        % load model
        load([model_dir model_list{m} '/' submodellist(mm).name '/fitted_model.mat'],'R');
        
        % get Goodness of fit metrics
        AIC{m}(mm, :) = R.GoF(:,3);
        F{m}(mm, :) = R.GoF(:,1);
        BIC{m}(mm, :) = R.GoF(:,2);
        LL{m}(mm, :) = R.GoF(:,4);
        
        % make short name
        underscore_ind = strfind(model_list{m}, '_');
        model_short_name{mm} = model_list{mm}(underscore_ind(1)+1:underscore_ind(2)-1);
        
        
    end
    
%     % make group comparison analysis on BIC
%     % recovery of the good model based on the simulation analysis
%     options.figName = model_list{m};
%     options.modelNames = model_short_name;
%     VBA_groupBMC(LL{m},options);
    
    subplot(floor(size(model_list,1)/2),ceil(size(model_list,1)/2),m)
%     figure('Name', ['Goodness of fit ' model_list{m}]);
    diff_metrics = [-2*sum(BIC{m},2)-(-2*sum(BIC{m}(m,:),2)), -2*sum(AIC{m},2)-(-2*sum(AIC{m}(m,:),2)) -2*sum(F{m},2)-(-2*sum(F{m}(m,:),2))];
    bar(diff_metrics);
    set(gca,'xticklabels', model_short_name)
    set(gcf, 'color', 'w');
    ylabel('Diff GoF');
    ylim([-6000 8000])
    title(model_short_name{m})
    
%     close all
    
end

saveas(gcf,'recovery', 'svg')
