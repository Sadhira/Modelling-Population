%% Q1

t = 1:0.1:100;
alpha = 1;
beta = 1;
gamma = 0.3;

[T,X] = ode45(@xdot_m4, [0 10], [1 2]');  
plot(T,X)



% @ is a function handle that tells matlab to use the following inputs
%as inputs for that function
% [0 10] is my time span
% 1 and 2 are my arbitrary "positive initial condition"s for x and y



% Q2

% nullclines: 
% x = 1/gamma * alpha/(1+y^4)
% y = 1/gamma * beta/(1+x^4)

   



%Q3

% previous is plotting x(x,y,t) and y(x,y,t) against t
% now plot dxdt(x,y) and dydt(x,y) against x and y 


%x = [0 1]

  
% for 
%     dxdt = (alpha /(1 + x(2)^4)) - (gamma * x(1));
% 
%     dydt = (beta/(1 + x(1)^4) - (gamma * x(2)));
    
    
