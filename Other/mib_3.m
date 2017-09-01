%% Q1

% dx = ?kxdt + ?dW
% x(t + h) = x(t) + h [?k x(t)] + ? root-h * randn


k = 3/16;
h = 0.01;
xout(1) = 6;
xanalyt(1) = 6;
t = [0:h:10];
sigma = 0.2;

xsum = 1:zeros(length(t+1));
xsum(1) = xout(1);

% a) superimpose 20 runs of the function

figure(1)

for i = 1:20
    for n=1:length(t)-1
            xout(n+1)= eulermib2(xout(n), h, k) + (sigma * h^0.5 * randn);
            xanalyt(n+1)= eulermib2(xout(n), h, k);
            
            %xsum(n+1)= xsum(n+1)+xout(n+1);            
            
    end
    
    %xsum = xsum + xout; apparently doesn't work
    xsums(i,:)=xout(:)
    
    hold on
    plot (t, xout)    
end
    
    

% for i=2:20
%     xsum(n) = xsum(n) / 20;
% end

%xmean = xsum /20;

xmean = mean(xsums, 1); %averages along columns
                        %mean(xsums, 2) would average along rows

hold on
    plot (t, xmean, '-.or');

mseval = immse(xmean, xanalyt)


%% Q2 - use greater t range


k = 3/16;
h = 0.01;
t2 = [0:h:200];

x2out = [1:length(t2)];   %it's best to always define the size
x3out = [1:length(t2)];
xvals = [2:(length(t2)+1)]; %don't know why,length doesn't match without +1

x2out(1) = 0;
x3out(1) = 0;
%xanalyt(1) = 0;

sigmaa = 0.1;
sigmab = 5;


for n=1:length(t2)-1
    x2out(n+1)= eulermib2(x2out(n), h, k) + (sigmaa * h^0.5 * randn);
    x3out(n+1)= eulermib2(x3out(n), h, k) + (sigmab * h^0.5 * randn);
end

figure(2)
subplot(1,2,1)
hist(x2out);

subplot(1,2,2)
%figure(3)
%hold on
hist(x3out);


xvals(1,:) = x2out;
xvals(2,:) = x3out;
xmean2 = mean(xvals,1); 

figure(4)
hist(xmean2);   %not necessary to plot

meanstd2 = std(x2out)
meanstd3 = std(x3out)


%the increase in width is caused by the increase in noise aka stochastic
%variation, and so there is also a higher standard deviation





