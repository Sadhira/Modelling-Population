%% need to do 2b and 3

%% Q2 - analytical solution, using x

k=3/5;
t=1:11;
x=zeros(1:11);
x(1) = 4;

%  dx/dt = -kx

figure(1)
x = x(1)*exp(-k*(t-1));
plot(0:10, x); 


%% Q1 - numerical solution, using dx/dt

%  M(t,y)*y' = f(t,y)
%  dx/dt = -kx


%"expeq" is defined in another script as 
% 
% function xp = expeq(t, x);
% k = 3/5;
% xp = -k * x;

figure(2)
%[t,x]=ode45(’F’,[t0,tf],[x10,x20]);
%in this case, there is only x, not x1 and x2
%so for the last parameter, just use x0 which is 4
[t,x] = ode45('expeq', [0 10], 4);
plot(t, x);

%% Q2b - mean squared error of both solutions

% [x,t] = house_dataset;
% net = feedforwardnet(10);
% net.performFcn = 'mse';  % Redundant, MSE is default
% net.performParam.regularization = 0.01;
% net = train(net,x,t);
% y = net(x);
% perf  = perform(net,t,y);

% figure(1)
% perf = mse(x,t,y);
% plot(t, mse);


%MSE = mean((desired - mean).^2);

%1 - analytical
k=3/5;
t=1:11;
x1=zeros(1:11);
x1(1) = 4;

figure(3)
x1 = x1(1)*exp(-k*(t-1));
hold on;
plot(0:10, x1); 

%2 - numerical
hold on;
figure(4)
[t,x2] = ode45('expeq', 0:10, 4);
plot(t, x2);



%reshape(mean(mean((x1 - x2).^2, 2),1)[1,3]);

%%error = mse(x2, t, x1)
% not sure if this is correct?
x1 = x1';
msevals = zeros(11:1);
%msevals = mean((x2(1:11)-x1(1:11)).^2);
msevals = (x2-x1).^2;
mseval = mean(msevals)
mseval2 = immse(x1, x2)

% figure(4)
% plot((0:10), msevals(1:11));

%[T numerical_solution_values] = ode45(@(t,x) expel(t,x,k), boundaries, initial_value);
%plot(T, 

% %% --
% figure;
% plot(t, x, r)

%% Q3 

% [T numerical_solution_values] = ode45(@(t,x) training1numerical(t,x,k), boundaries, initial_value);



%% Q3 final answer explanation

% ode45 takes smaller steps at the beginning 
% because small changes in t result in large changes 
% in x. as the line flattens, larger steps can be taken
% because there are smaller differences in x. the 
% function then takes smaller steps again to reach a 
% more a accurate value at the end
% accurate 
%