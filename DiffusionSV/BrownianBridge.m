function [BB]=BrownianBridge(t_1,t_2,b_1,b_2,m)

% Imputes a Brownian bridge with m points, excluding the initial
% and final points, between the times t_1 and t_2 with 
% corresponding values b_1 and b_2.

B=BrownianMotion(t_1,t_2,m+1); % B is a Browian motion with initial value B(t_1)=0
t=(0:m+1)./(m+1);
BB=b_1*ones(1,size(B,2))+B+(b_2-B(m+2)-b_1)*t;
