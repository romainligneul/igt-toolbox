function [fb] = h_truefalse(deck_chosen,t,in)
% compares the entry yt with a stored reference answer u0(t)
% if t == 1
%     fb = [0 0];
% else
ind = randi(length(in{find(deck_chosen)}));
fb = in{find(deck_chosen)}(ind,:);
% end
