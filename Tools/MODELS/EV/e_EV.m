function  [fx] = e_EV_IOWA(x,P,u,in)
%%%% EV model / Evolution function
% P(1) corresponds to weight placed on losses relative to gains. A value
% closer to one means that learning is more driven by losses. [0-1]
% P(2) corresponds to the learning rate of the model. [0-1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For trial 0, the evolution function uses the starting values provided in 
% in.hs.initval

%% parameter transformation
% Parameters are transformed in order to vary within their specific range.
for pp = 1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end;

%% write current on next
fx = x;

% deck selected
d = u(2);
% compute utility of the outcome harvested
v = (1-P(1))*u(3)- P(1)*abs(u(4));

% update relevant value
if u(1)==0
    fx= x;
else
    fx(d)= x(d)+P(2)*(v-x(d));
end;

end