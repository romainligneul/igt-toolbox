function  gx = o_EXPLORE(x,P,u,in)
%%%% EXPLORE model / Observation function
% see run_EXPLORE for a description of the model
% see documentation of VBA_toolbox for a description of inputs / outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R.Ligneul 06/17

%% parameter transformation / should always be performed.
for pp =1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end

%% compute probability of each choice given hidden states
gx = sub_softmax(x(1:4)+P(1)*x(5:8)+P(2)*x(9:12), 1);

%% softmax subfunction
function p = sub_softmax(collapsed, consistency)
    for i = 1:length(collapsed)
      p(i,1) = exp(collapsed(i)*consistency)/(exp(collapsed(1)*consistency)+exp(collapsed(2)*consistency)+exp(collapsed(3)*consistency)+exp(collapsed(4)*consistency));
    end
end

end