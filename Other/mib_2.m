%% Q1 - numerical, Euler

%%x(t + h) = x(t) + h [?k x(t)],  t ? [0, 10]       equation

k = 0.25;
h = 0.01;

t = [0:h:10];
%x = [0:0.01:10];
xout(1) = 5;
 
% for i = 1:10
%     x(i+h) = x(i) + -h(k*x(i));       false attempt 1
% end

for n=1:length(t)-1
        xout(n+1)= eulermib2(xout(n), h, k);
end
  

% for  t = 0:h:10
%     x = eulermib2(x(i), h, k);        false attempt 2
%     i = i + 1;
% end


%% Q1b - analytical

%could have used dsolve for mib_1, too
    
[SOLUTIONS]= dsolve('Dx=-k*x','x(0)=5');
Solution=eval(SOLUTIONS); 
T=t;

figure(1)
plot(T, Solution);


%alternative solution from answers:
%xA1=xE1(1)*exp(-k*t1);


%% Q2 h = 0.001

figure(2)

h = 0.01;
t1 = [0:h:10];
xout1(1)=5;

for n=1:length(t1)-1
        xout1(n+1)= eulermib2(xout1(n), h, k);
end

plot(xout1, t1, 'r');


%dsolve seems to only use t as its time input, so it is probably safe
%to define them specifically using 't' before using the function
%(not doing so did cause an error where t=1:0.01:10 was used

t=0:0.01:10; 
[SOLUTIONS]= dsolve('Dx=-k*x','x(0)=5');
Solution1=eval(SOLUTIONS); 
T1=t1;

plot(T1, Solution1);


hold on

h = 0.001;
t2 =0:h:10;
xout2(1)=5;

for n=1:length(t2)-1
        xout2(n+1)= eulermib2(xout2(n), h, k);
end

plot(t2, xout2, 'c');
   
t=0:0.001:10;
[SOLUTIONS]= dsolve('Dx=-k*x','x(0)=5');
Solution2=eval(SOLUTIONS); 
hold on
plot(t, Solution2, 'g');      

%need to zoom in a lot to clearly see the differences between the two lines


mseval1 = immse(Solution1, xout1)
mseval2 = immse(Solution2, xout2)

%there is a smaller error using smaller values of k because 

