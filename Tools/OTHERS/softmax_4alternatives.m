function p = softmax_4alternatives(x, consistency)
   % temp = 3^consistency-1;
   temp = consistency;
    for i = 1:length(x)
      p(i,1) = exp(x(i)*temp)/(exp(x(1)*temp)+exp(x(2)*temp)+exp(x(3)*temp)+exp(x(4)*temp));
    end
end   