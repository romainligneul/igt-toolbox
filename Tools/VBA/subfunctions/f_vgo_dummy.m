function [fx] = f_vgo_dummy(x,P,u,in)

fx = x + sig(P(1))*(u-x) + sig(P(2))*x;