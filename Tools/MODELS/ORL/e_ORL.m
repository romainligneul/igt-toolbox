function  [fx] = e_EXPLORE(x,P,u,in)
%%%% EXPLORE model / Evolution function
% see run_EXPLORE for a description of the model
% see documentation of VBA_toolbox for a description of inputs / outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R.Ligneul 06/17

%% parameter transformation
for pp = 1:length(P)
    P(pp) = in.param_transform{pp}(P(pp));
end

%% write current on next
fx = x;

%% update
% deck selected
d = u(2);
outcome = u(3)-abs(u(4));
pers_param = (3^P(3))-1;

if u(1)==0 % initialization
    fx=[0 0 0 0 0.5 0.5 0.5 0.5 1 1 1 1];
    
else       % update
    for de = 1:4
        % selected
        if d==de
            if outcome>=0
                fx(de)= x(de)+P(1)*(outcome-x(de));
                fx(de+4)=x(de+4)+P(1)*(sign(outcome)-x(de+4));
                fx(de+8)=1/(1+pers_param);
            else
                fx(de)= x(de)+P(2)*(outcome-x(de));
                fx(de+4)=x(de+4)+P(2)*(sign(outcome)-x(de+4));
            end
        else % counterfactual learning only on outcome frequencies
            if outcome>=0
                fx(de+4)=x(de+4)+P(2)*((-sign(outcome)/3)-x(de+4));
            else
                fx(de+4)=x(de+4)+P(1)*((-sign(outcome)/3)-x(de+4));
            end
            fx(de+8)=x(de+8)/(1+pers_param);
        end
    end
end

end