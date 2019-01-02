function  [gx] = o_PVL_IOWA(x,P,u,in)
%%%% PVL model / Observation function
% P(1) corresponds to the consistency of the model. It is transformed in a 
% way which may choices more deterministic as a function of time. [-2 2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parameter transformation
for pp =1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end;

%% compute probability of each choice given hidden states
decks_ind = in.hs.decks_ind;
consistency = 3^P(1)-1;
gx = sub_softmax(x(decks_ind), consistency);

%% softmax subfunction
function p = sub_softmax(collapsed, consistency)
    for i = 1:length(collapsed)
      p(i,1) = exp(collapsed(i)*consistency)/(exp(collapsed(1)*consistency)+exp(collapsed(2)*consistency)+exp(collapsed(3)*consistency)+exp(collapsed(4)*consistency));
    end
end   

end