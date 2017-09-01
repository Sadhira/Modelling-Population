%% Q1

N = 276;
p = 0.05;

A = zeros(N);

for i=1:N
	for j=1:N
    	pval = rand;
    	if (pval <= p)
        	A(i,j) = 1;
    	end
	end
end

spy(A);



[x,y]=GenerateRingCoordinates(N);

T = [x',y'];

gplot(A,T);



%% Q2a

SNO = load('SocialNetwork.mat');  %loading struct from file
SN = SNO.A;                   	%getting array from struct

figure(1)
spy(SN);

figure(2)
gplot(SN, T);




%Q2b i

for i=1:N  
    	%dist = (dijkstra(SN,i)), i)
    	mydist(:,i) = dijkstra(SN,i);
    	 
end


%Q2b ii

diam = 0;

for i = 1:276
    for j = i:276                   %ignores values below the diagonal
        diam = diam + mydist(i,j);
    end
end

netdiam = diam/(276^2);   
                      %denom is squared because of the 276x276 matrix size
                      
                       %need to plot histogram
                   %%   
%Q2b iii
                       
RN = zeros(276);

tL = 276*(276-1)

p = L / tL;

%p = 0.85;

% 


for i = 1:276
    for j = 1:276
        
        pval = rand;
    	if (pval <= p)
        	RN(i,j) = 1;
        end
    end
end

figure(8)
spy(RN);

distr = zeros(276);

for i=1:N  
    	%dist = (dijkstra(SN,i)), i)
    	distr(:,i) = dijkstra(RN,i);
    	 
end

diamr = 0;

for i = 1:276
    for j = i:276                   %ignores values below the diagonal
        diamr = diamr + distr(i,j);
    end
end

netdiamr = diamr/(276*276); 
%netdiamr = diamr/100000000; 


%%

%Q2c

%d = diameter
d0 = mean(diam); %averages array in one dimension to give row vector
d = mean(d0); %averages row vector to give final single mean


F = sum(SN, 2); %number of friends per person

L = sum(F, 1); % number of links

figure(3)
hist(F);

%NEED TO DO THIRD PART
%I do have a code for the clustering coefficient in part 3


%% test

TSN = [0 1 1 1 0 0 1 1;    % TSN = test social network
      1 0 1 0 0 1 0 0;
      1 1 0 0 0 0 0 0;
      1 0 0 0 0 0 0 0;
      0 0 0 0 0 1 0 1;
      0 1 0 0 1 0 1 0;
      1 0 0 0 0 1 0 0;
      1 0 0 0 1 0 0 0;];
  
%% test

N = zeros(8,1);
M = zeros(8,8);

for i = 1:8

    for j = 1:8

        if TSN(i,j) > 0

            for k = 1:8

                if TSN(i,k) == TSN(j,k)
                    
                    if TSN(i,k)>0

                        N(i,1) = N(i,1) + 1;
                        M(i,k) = 1;
                    else M(i,k) = 0;
                    end
                end
            end
        end

    end

    N(i,1) = N(i,1) * 0.5;

end

%% test

%Ci = 2ni/k(k-1)   ----> ??

%note nA = 1 in lecture slide

% for row = 1:156
% 	for col = 1:256
%    	 
%     	if SN(row, col)>0
%        	 
%        	 
% 	C(i)
% C = F


N = zeros(276,1);
M = zeros(276,276);

for i = 1:276

    for j = 1:276

        if SN(i,j) > 0

            for k = 1:276

                %if SN(i,k) == SN(j,k)
                if SN(i,k)>0 && SN(j,k)>0   %some values are 2 instead of 1
                    
                    if SN(i,k)>0

                        N(i,1) = N(i,1) + 1;
                        M(i,k) = 1;
                    else M(i,k) = 0;
                    end
                end
            end
        end

    end

    N(i,1) = N(i,1) * 0.5;

end
    
 

%Ci = 2ni/k(k-1)

%for i = 1:276
    
%    CC(i,1) = 2 * 


%% Q3a

%DN = SN;  %disease network seperate to not alter original

%add extra dimension to indicated noninfected (0) vs infected (1)

%infect one person (at random?)

%loop ten times: go to each infected node and infect its neighbours


DN = zeros(276,11) ;  % disease network
                      % consists of initial states and 10 cycles

DN(27,1) = 1;  %infects one initial node identified by row number

p = 1; %proportion of populatation infected

n = 1; %number of cycles including start initialisation
K = 0; %row vector of number of people infected per cycle
k = 1;                %number of infected people for last cycle
t = round(p * 276);   %target of infected people



while t>k
    
    for i = 1:276

        if DN(i,n) == 1

            for j = 1:276

                if SN(i,j) > 0

                    if DN(j,n+1) == 0

                        DN(j,n+1) = 1;
                    end
                end
            end
        end
    end
    
    K = sum(DN);
    k = K(1,n+1);
    n = n+1;
    
end

n = n-1  %infection cycles without counting the initialisation

P = K./276.*100;

figure(1)
bar(0:10, P);
%plot(0:10, P);
ylabel('infected population %');
xlabel('number of cycles');

% a.i) it takes 5 steps to infect the whole network from node 1

% a.ii) change p from 1 to to 0.5
%       it takes 2 steps to infect half the network from node 1

% a.iii) it does matter slightly who is chosen initially
%        eg it takes 4 steps to infect the whole network from node 27
%        but still takes 2 steps to infect half the network from node 27
        

%% Q3b

DN = zeros(276,11) ;  % disease network

DN(27,1) = 1;  %infects one initial node

p = 1; %proportion of populatation infected

n = 1; %number of cycles including start initialisation
K = 0; %row vector of number of people infected per cycle
k = 1;                %number of infected people for last cycle
t = round(p * 276);   %target of infected people


figure(7)

for d = 1:276
    
    DN = zeros(276,11);
    DN(d,1) = 1;
    %DN(276,1) = 1;
    n=1;
    
    while n<11

        for i = 1:276

            if DN(i,n) == 1

                for j = 1:276

                    if SN(i,j) > 0

                        if DN(j,n+1) == 0

                            DN(j,n+1) = 1;
                        end
                    end
                end
            end
        end

        n = n+1;

    end
    
    K = sum(DN);
        %k = K(1,n+1);
    P = K./276.*100;
        
    hold on;
        
    plot(0:10, P);
    xlabel('number of cycles');
    ylabel('infected % of population');

    
end

% the entire network will always be infected by the 7th cycle

% for some reason, it seems that there are infected population values that 
% jump to and stay at 0 for every cycle. I currently have no idea why, but 
% that's what the downward sloping lines in the graph are.



%% Q4












