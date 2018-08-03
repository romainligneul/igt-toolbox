clear all;
files = dir('*xlsx');


[xchoices] = xlsread([files(1).name]);
[xlosses] = xlsread([files(2).name]);
[xwins] = xlsread([files(3).name]);
[study] = importdata(['index_100.txt']);
study = study(2:end);

% recompile study labels (and age condition for Wood et al.)
wood_ind = 0;
for f = 1:length(study)
    dumind = strfind(study{f,1}, '"');
    study{f,1} = study{f,1}(dumind(3)+1:dumind(4)-1);
    if strcmp(study{f,1},'Wood');
        wood_ind = wood_ind+1;
        if wood_ind < 91
            study{f,1} = 'Wood_Young';
        else
            study{f,1} = 'Wood_Old';
        end
    end
end
    
for f = 1:size(xchoices,1)
    
    data{f}.id = f;
    data{f}.cond = 0;
    data{f}.cond_label =  study{f,1};
    if strcmp(study{f,1}, 'Wood_Old')
        data{f}.cond = 2;
    elseif strcmp(study{f,1}, 'Wood_Young')
        data{f}.cond = 1;
    end
    
    data{f}.trial(:,1) = 1:100;
%     for t = 1:size(xnum,1)
        data{f}.deck(:,1) = xchoices(f,:); %strmatch(xtxt(t,2),{'A', 'B', 'C', 'D'})
%     end
    data{f}.win = xwins(f,:)'; % xnum(:,3);
    data{f}.lose = xlosses(f,:)';%xnum(:,4);
    data{f}.score = nan(size(xchoices,2), 1); %sxnum(:,5);
    data{f}.rt = nan(size(xchoices,2), 1);
       
end

save('IGTdata.mat', 'data')