function  [fx] = e_VPP_IOWA(x,P,u,in)
%%%% VPP model / Evolution function
% P(1) corresponds to value sensitivity, with more linear representations of
% values being associated with a parameter closer to 1. [0-1]
% P(2) corresponds the loss aversion parameter. Higher values indicate 
% more intense reactivity to losses. [0 5]
% P(3) corresponds to the learning rate of the model. [0-1]
% P(4) corresponds to the decay parameter of the perseveration values. [0 1]
% P(5) corresponds to the tendency to persevere after a gain [-1 1]
% P(6) corresponds to the tendency to persevere after a loss [-1 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For trial 0, the evolution function uses the starting values provided in 
% in.hs.initval

%% Specific comments:
% in first version of the evolution function:
% - one single learning rate for SAS & SS
% - one learning rate for controllability

%% parameter transformation / should always be performed.
% P(1) = Shape of the value function [0 1]
% P(2) = Loss aversion [0 5];
% P(3) = updating [0 1];
% raw parameters correspond to the x=x transformation.
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
    v = r.^P(1);
else
    v = -P(2)*(abs(r)^P(1));
end

% update relevant state (except trial 0 in which we just map the initial values)
if u(1)==0
else
    fx(d)= x(d)+P(3)*(v-x(d));
end;

% update perseveration value
if u(1)==0
else
    for de = 1:4
        if d==de
             if r>=0
                fx(de+4)= P(4)*x(de+4) + P(5);
             else
                fx(de+4)= P(4)*x(de+4) + P(6);
             end
         else
            fx(de+4)= P(4)*x(de+4);
        end
    end
end
            
end