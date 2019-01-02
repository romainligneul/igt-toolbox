clear all;
files = dir('IGTdata/*txt');


[raw{1}] = importdata(['IGTdata_amphetamine.txt']);
[raw{2}] = importdata(['IGTdata_heroin.txt']);
[raw{3}] = importdata(['IGTdata_healthy_control.txt']);
cond_labels = {'amphetamine', 'heroin', 'controls'};
ss=0
for f = 1:3;%size(xchoices,1);
    raw{f}.data;
    subjects = unique(raw{f}.data(:,end));
    for s = subjects'
        if length(raw{f}.data(raw{f}.data(:,end)==s,1))==100
        ss = ss+1;
        data{ss}.id = s;
        data{ss}.cond = f;
        data{ss}.cond_label = cond_labels{data{ss}.cond};
        
        data{ss}.trial = raw{f}.data(raw{f}.data(:,end)==s,1);
        %     for t = 1:size(xnum,1)
        data{ss}.deck(:,1) = raw{f}.data(raw{f}.data(:,end)==s,2); %strmatch(xtxt(t,2),{'A', 'B', 'C', 'D'})
        %     end
        data{ss}.win = raw{f}.data(raw{f}.data(:,end)==s,3); %xwins(f,:)'; % xnum(:,3);
        data{ss}.lose = raw{f}.data(raw{f}.data(:,end)==s,4);
%         data{ss}.score = nan(size(xchoices,2), 1); %sxnum(:,5);
%         data{ss}.rt = nan(size(xchoices,2), 1);  
        end
    end
end

save('IGTdata.mat', 'data')