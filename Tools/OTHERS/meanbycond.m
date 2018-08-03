function meanbycond = meanbycond(data, cond0, filter)
% require unidimensionnal data and cond vectors of same length
% and eventually a filter rule:
% - filter alone must be a logical 0/1 array. It will keep only the
% data where filter == 1.
% remark, meanbycond are ordered by ascending order of cond numbers/letters
% if sum(isnan(cond0))>0
%     cond0(isnan(cond0))= pi;
%     display('NaN in the cond data');
% end

condtype = unique(cond0);

if ~isempty(filter)
    data = data(filter);
    cond0 = cond0(filter);
end

condtype = sort(condtype);


for c = 1:length(condtype)
    meanbycond(1,c) = nanmean(data(cond0 == condtype(c)));
end;

end

