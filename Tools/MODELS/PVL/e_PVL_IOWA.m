function  fx = e_PVL_IOWA(x,P,u,in)
%%%% PVL model / Evolution function
% P(1) corresponds to value sensitivity, with more linear representations of
% values being associated with a parameter closer to 1. [0-1]
% P(2) corresponds the loss aversion parameter. Higher values indicate 
% more intense reactivity to losses. [0 5]
% P(3) corresponds the decay parameter, that is, the speed at which values 
% are forgotten. A parameter closer to 1 indicate less forgetting. [0 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For trial 0, the evolution function uses the starting values provided in 
% in.hs.initval

%% parameter transformation 
for pp = 1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end;

%% write current on next
fx = x;

%% update

% deck selected
d = u(2);

% compute utility of the outcome harvested
r = u(3)-abs(u(4));
if r>=0
    v = r^P(1);
else
    v = -P(2)*(abs(r)^P(1));
end

% update relevant state (except trial 0 in which we just map the initial values)
if u(1)==0
else
    for de = 1:4
        if d~=de
            fx(de)= x(de)*P(3);
        else
            fx(de)= x(de)*P(3)+v;
        end
    end
end;

end