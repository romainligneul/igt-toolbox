function  [gx] = o_evolution_template(x,P,u,in)
%%%% VPP model / Observation function
% P(1) corresponds to the consistency of the model. [0 5]
% P(2) corresponds to the 'expectancy weight', that it, to the tendency to
% rely on learned values rather than perseverative tendency when making
% decision. [0 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parameter transformation
for pp =1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end;

%% compute probability of each choice given hidden states
decks_ind = in.hs.decks_ind;
for d=1:4
weighting(d) = x(d)*P(1) + (1-P(1))*x(d+4);
end
consistency = 3^P(2)-1;
gx = sub_softmax(weighting, consistency);

%% softmax subfunction
function p = sub_softmax(collapsed, consistency)
    for i = 1:length(collapsed)
      p(i,1) = exp(collapsed(i)*consistency)/(exp(collapsed(1)*consistency)+exp(collapsed(2)*consistency)+exp(collapsed(3)*consistency)+exp(collapsed(4)*consistency));
    end
end

end   
