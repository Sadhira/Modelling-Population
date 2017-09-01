%% mib coursework

%% q2a

BCM = load('fMRIBrainConnectivityMatrix.mat'); % brain connectivity matrix

CD = BCM.Coord;                 	% coordinates original

WCO = BCM.WeightedConnectivity;  	% weighted connectivity original

WC  = WCO;      	

N = size(WC,1);    %N = 638;
                	

for i = 1:N    	
	WC(i,i) = 0;        	% sets diagonals to 0 so nodes don't connect to themselves
                	
	for j = 1:N      	
%     	if WC(i,j) > 0   	% sets values to only be 0 or 1
%         	WC(i,j) = 1;  
%     	end
    	
    	WC(j,i) = WC(i,j);  % makes values symmetric along diagonal
        	
	end
end



%% q2b

spy(WC);
title('WCacency Matrix of Weighted Connectivity');

%% number of nodes
N = size(WC,1);

nodeConnected = zeros(N,1);
connections = 0;

for i = 1:N
   for j = 1:N
    	if WC(i,j) > 0
        	connections = connections + 1;
    	end
   end
  
   if connections > 0
   	nodeConnected(i,1) = 1;
   end
  
   connections = 0;
end


nodesActive = 0;
nodesActive = sum(nodeConnected);

%%
nodesDeleted = 0;
%nodesDeleted corrects for the offset in the index caused by deletion 

for i = 1:N
    if nodeConnected(i,1) == 0
        WC(i - nodesDeleted,:) = [];  
        WC(:,i - nodesDeleted) = [];
        nodesDeleted = nodesDeleted + 1;
 
    end
end
 
% figure(2)
% spy(WC);
% title('WCacency Matrix of Weighted Connectivity');

%% diameter

N = size(WC,1);

dist = zeros(N);

for i=1:N 
	 %dist = (dijkstra(WC,i)), i)
	 dist(:,i) = dijkstra(WC,i);
	 
end

diam = 0;

for i = 1:N
	for j = i:N                  	%ignores values below the diagonal
        if dist(i,j) ~=inf
            diam = diam + dist(i,j);
        end
	end
end

netDiam = diam/(N^2);
%netdiamr = diamr/100000000;

%%
% Ratty's code 
% 
% shortest_distance=Inf(N); %%all nodes are initially assumed to be unconnected
% for i=1:1:N                     %%This loop uses Dijkstra's algorithm to find the shortest
%     shortest_distance(i,:)=dijkstra(WC,i);  %%distance of a selected node to all other nodes
% end
% total_length=0;  %%variable for the sum of all the shortest lengths
% for i=1:1:N       %%loop for calculating the sum for all shortest lengths
%     for j=i:1:N
%         if shortest_distance(i,j)~=Inf & i~=j
%            total_length=total_length+shortest_distance(i,j);
%         end
%     end
% end
% b=1:1:287;
% b=sum(b,2)-13*287; %%-13*287 to exclude the unconnected nodes.
% diameter=total_length/b; %%Here the variable b represents total number of connections 
% 

%% number of links

BC = WC; % BC = basic connections of 1s and 0s, where they are unweighted

for i = 1:N
    for j = 1:N
        if BC(i,j) > 0
            BC(i,j) = 1;
        end
    end
end

links2 = sum(BC,2);         % column vector of the sum of each BC row
links1 = sum(links2,1);     % sum of the above column vector
links = links1/2            % sum divided by 2

%%

nodeD = sum(BC,2);      % column vector of the sum of each BC row

figure(1)
hist(nodeD,15)


[freq, histEdges0] = histcounts(nodeD,15);
 
histEdges = histEdges0;
histEdges(:,1) = [];       %removes starting 0 value 
 
figure(2)
scatter(freq, histEdges)
set(gca,'xscale','log')
set(gca,'yscale','log')




%% q2c

%plot3(WC,BCM,WC)

%seperation of data into x, y, and z node coordinates

nX = CD(:,1);
nY = CD(:,2);
nZ = CD(:,3);

%plot3(nX, nY, nZ, '.')      %plots nodes


%
subplot(1,4,1)
plot3(nX, nY, nZ, '.')
subplot(1,4,2)
plot3(nX, nY, nZ, '.')
subplot(1,4,3)
plot3(nX, nY, nZ, '.')
subplot(1,4,4)
plot3(nX, nY, nZ, '.')
title('Visualisation in Anatomical Space');


% visualisation of only nodes





subplot(2,2,1)
plot3(nX, nY, nZ, '.')
title('X-Y-Z view')

subplot(2,2,2)
plot3(nX, nY, nZ, '.')
title('X-Y view')

subplot(2,2,3)
plot3(nX, nY, nZ, '.')
title('X-Z view')

subplot(2,2,4)
plot3(nX, nY, nZ, '.')
title('Y-Z view')

for s = 1:4                     
    subplot(2,2,s)
    
    plot3(nX, nY, nZ, '.')        % plots nodes
    
    switch s 
        case 1
            title('X-Y-Z view');
        case 2
            title('X-Y view');
        case 3
            title('X-Z view');
        case 4
            title('Y-Z view');
    end
end

% visualisation of nodes with edges
subplot(2,2,1)
plot3(nX, nY, nZ, '.')
line(nX, nY, nZ)
title('X-Y-Z view')

subplot(2,2,2)
plot3(nX, nY, nZ, '.')
line(nX, nY, nZ)
title('X-Y view')

subplot(2,2,3)
plot3(nX, nY, nZ, '.')
line(nX, nY, nZ)
title('X-Z view')

subplot(2,2,4)
line(nX, nY, nZ)
plot3(nX, nY, nZ, '.')
line(nX, nY, nZ)
title('Y-Z view')


figure(8)
for s = 1:4                     
    subplot(2,2,s)
    
    plot3(nX, nY, nZ, '.')        % plots nodes
    line(nX, nY, nZ)              % plots edges
    
    switch s 
        case 1
            title('X-Y-Z view');
        case 2
            title('X-Y view');
        case 3
            title('X-Z view');
        case 4
            title('Y-Z view');
    end
end
    
%
figure(2)
plot3(nX, nY, nZ)



%% q2c physical diameter


DC = zeros(N);          % DC = distance connectivity

for i = 1:N
    for j = i:N
        
        if BC(i,j) > 0 
                     
            xdist = (CD(i,1) + CD(j,1))^2;
            ydist = (CD(i,2) + CD(j,2))^2;
            zdist = (CD(i,3) + CD(j,3))^2;
            
            DC(i,j) = sqrt(xdist + ydist + zdist);      % calculates distance between each node pair
            
            DC(j,i) = DC(i,j);                          % mirrors along diagonal
            
        end
    end
end



pDist = zeros(N);                   %pD = physical dist

for i=1:N 
	 %dist = (dijkstra(WC,i)), i)
	 pDist(:,i) = dijkstra(WC,i);
	 
end

pDiam = 0;

for i = 1:N
	for j = i:N                  	%ignores values below the diagonal
        if dist(i,j) ~=inf
            pDiam = pDiam + dist(i,j);
        end
	end
end

pNetDiam = pDiam/(N);



%% take 2 ^

%% q2c physical diameter

N = 638;
DC = zeros(N);          % DC = distance connectivity

for i = 1:N
    for j = i:N
        
        if WCO(i,j) > 0 
                     
            xdist = (CD(i,1) - CD(j,1))^2;
            ydist = (CD(i,2) - CD(j,2))^2;
            zdist = (CD(i,3) - CD(j,3))^2;
            
            DC(i,j) = sqrt(xdist + ydist + zdist);      % calculates distance between each node pair
            
            DC(j,i) = DC(i,j);                          % mirrors along diagonal
            
        end
    end
end



pDist = zeros(N);                   %pDist = physical dist

for i=1:N 
	 %dist = (dijkstra(WC,i)), i)
	 pDist(:,i) = dijkstra(DC,i);
	 
end

pDiam = 0;
count = 0;

for i = 1:N
	for j = i:N                  	%ignores values below the diagonal
        if pDist(i,j) ~= Inf||0
            pDiam = pDiam + pDist(i,j);
            count = count + 1;
        end
	end
end

pNetDiam = pDiam/(count);



%%
% determine spatial embedded diameter--------------------------------------
 
% calculate the distance between the shortest path
 
distance_matrix = zeros(638,638);%store distance here
 
for i =1:638
    for j=1:638
        if WCO(i,j)>0
            distance_matrix(i,j)= sqrt(((CD(i,1)-CD(j,1))^2) + ((CD(i,2)-CD(j,2))^2) + ((CD(i,3)-CD(j,3))^2) );
        end
    end
end
 
% take the shortest distance
for i = 1:638;
    shortest_distance(i,:) = dijkstra(distance_matrix ,i);
end
 
%eliminate the 0s and any inf and take the average.
sum1 = 0; 
count1 = 0;
for i=1:638; %writting for loop like this only takes the upper triangular matrix. no need to create u.t matrix for shortest distance
    for j =i:638;
        if shortest_distance(i,j) ~= Inf||0 ;
            sum1 = sum1 + shortest_distance(i,j);
            count1=count1+1;
        end
    end
end
mean_dist_mills = (sum1/count1)


%%


nX = CD(:,1);
nY = CD(:,2);
nZ = CD(:,3);

%figure(10)
%plot3(nX, nY, nZ, '.')  

%s loop is for plotting 4 graphs to manually display different views

% for s = 1:4                     
%     subplot(2,2,s)
%     plot3(nX, nY, nZ, '.')
% 
%     for i = 1:N
% 
%         edge(1,:) = CD(i,:);
%         edge(2,:) = [0,0,0];
% 
%         for j = 1:N
%             if WC(i,j) > 0
% 
%                 edge(2,:) = CD(j,:);
%                 line(edge(:,1), edge(:,2), edge(:,3));     
%             end
%         end
%     end
%     
%      switch s 
%         case 1
%             title('X-Y-Z view');
%         case 2
%             title('X-Y view');
%         case 3
%             title('X-Z view');
%         case 4
%             title('Y-Z view');
%      end
% end





%% q2di

% average degree and degree distribution:
% Global topological characteristics of networks are the Average node degree and more
% informative - the node degree distribution. The latter is represented by a histogram,
% showing what proportion of nodes over the total number of nodes have a specific degree

% The degree of a node is the number of edges connected to the node.
% ---
 
% figure(2)
% hist(links2,15)
                             % w = weighted
wlinks2 = sum(WC,2);         % column vector of the sum of each WC row
wlinks1 = sum(links2,1);     % sum of the above column vector
wlinks = links1/2;           % sum divided by 2

figure(1)
hist(wlinks2,15)


%scatter(links2, links2)
%scatter(1:10,1:10)
%set(gca,'xscale','log')
%[nwfreq, nwhistEdges0] = histcounts(links2,15);

% [freq, histEdges0] = histcounts(wlinks2,15);
% 
% histEdges = histEdges0;
% histEdges(:,1) = [];


[freq, histMids] = hist(wlinks2,15);


figure(2)
scatter(freq, histMids)
set(gca,'xscale','log')
set(gca,'yscale','log')


hold on
% myfit = polyfit(freq, histMids);
% plot(myfit);



%% q2dii


LC = zeros(N,1);          % LC = linked connectivity

for i = 1:N
    
    for j = 1:N
        
        if BC(i,j) > 0
            
            for s = 1:N
                
                if BC(i,s)>0 && BC(j,s)>0
                    
                    LC(i,1) = LC(i,1) + 0.5; 
                    
                end
            end
        end
    end
end
         

% CC = local clustering coefficient   = C
% LC = linked connectivity (links between neighbours)   = n
% nodeD = node degree   = k

% for each node, C = 2n/k(k-1)
% therefore CC = 2*LC/nodeD(nodeD-1)


CC = zeros(N:1);
CC = 2.*LC./(nodeD.*(nodeD-1));

%scatter(nodeD,CC);
plot(nodeD,CC, 'o');
xlabel('node degree k');
ylabel('local clustering coefficient C');


CC2 = CC;
nodeD2 = nodeD;

nodesDeleted = 0;

CC2(isnan(CC2)) = -100

for i = 1:N
    
    
    if CC2(i - nodesDeleted, 1) == -100
        
        CC2(i - nodesDeleted, :) = [];
        nodeD2(i - nodesDeleted, :) = [];
        
        nodesDeleted = nodesDeleted + 1;
    end
end
        
    



% scatter(nodeD2,CC2);
% 
% hold on
% myfit = polyfit(nodeD2, CC2);
% plot(myfit);




figure(5)
scatter(nodeD2,CC2);
hold on ;
my_poly=polyfit(nodeD2,CC2,1);
X2= 1:200; % X data range 
Y2=polyval(my_poly,X2);                 % my_poly(1,1) = 1.7015e-04

plot(X2,Y2,'-');



% 
% figure(6)
% scatter(nodeD,CC);
% hold on ;
% my_poly=polyfit(nodeD,CC,1);
% X2= 1:200; % X data range 
% Y2=polyval(my_poly,X2);
% 
% plot(X2,Y2,'-');

%%scatter(LC, nodeD);
%%xlabel('links between neighbours n');
%%ylabel('node degree k');


%%

[lCC, gCC] = ClusteringCoefficient(WC);

lCC = lCC';


scatter(nodeD,lCC);
hold on ;
my_poly=polyfit(nodeD,lCC,1);
X2= 1:200; % X data range 
Y2=polyval(my_poly,X2);                 % my_poly(1,1) = 2.6389e-04

plot(X2,Y2,'-');


%%

[lCC, gCC] = ClusteringCoefficient(WCO);

lCC = lCC';

nodeDO = sum(WCO, 2);

figure(7)
scatter(nodeDO,lCC);                     
hold on ;
my_poly=polyfit(nodeDO,lCC,1);
X2= 1:200; % X data range 
Y2=polyval(my_poly,X2);                 % my_poly(1,1) = 5.2551e-04

plot(X2,Y2,'-');


%% q2diii


% T = total number of possible links
% L = total number of actual links
% p = probabilty needed to get actual number of links
% p = L/T


%linkCapacity = combntns(1:N,2)
linkCapacity = nchoosek(N,2);

prob = links/linkCapacity;

RC = zeros(N);              % RC = random connectivity


for i = 1:N
    for j = i:N
    
        randVal = rand;
        
        if randVal <= prob
            
            RC(i,j) = 1;
            RC(j,i) = 1;
            
        end
    end
end
            
%% cont ^


rLC = zeros(N,1);          % rLC = random linked connectivity

for i = 1:N
    
    for j = 1:N
        
        if RC(i,j) > 0
            
            for s = 1:N
                
                if RC(i,s)>0 && RC(j,s)>0
                    
                    rLC(i,1) = rLC(i,1) + 0.5; 
                    
                end
            end
        end
    end
end
         

% CC = local clustering coefficient   = C
% LC = linked connectivity (links between neighbours)   = n
% nodeD = node degree   = k

% for each node, C = 2n/k(k-1)
% therefore CC = 2*LC/nodeD(nodeD-1)


rNodeD = sum(BC,2);      % column vector of the sum of each RC row

rCC = zeros(N:1);
rCC = 2.*rLC./(rNodeD.*(rNodeD-1));





% for i = 1:N
%     if CC(N,1) = inf  
%         CC(N-i             %remove infinities?




rLCC = mean(rCC);
LCC = mean(CC);

if LCC > rLCC
    smallNet = 1
end

%%
% BW = rand([50 50 50])>.25; % your 3d matrix
%  i1 = [3, 2, 20];           % coordinates of initial point
%  i2 = [6, 10, 25];          % coordinates of end point
%  n = max(abs(i2-i1))+1;     % number of steps
%  i = arrayfun(@(a,b)round(linspace(a,b,n)),i1,i2,'uni',0);
%  idx = sub2ind(size(BW),i{:});
%  sumBW = nnz(BW(idx));
%  disp(cell2mat(i'));        % display trajectory
%  disp(sumBW);               % display number of 1's in path
%  
 %%
 
 
% BW = rand([50 50 50])>.25; % your 3d matrix
%  i1 = [3, 2, 20];           % coordinates of initial point
%  i2 = [6, 10, 25];          % coordinates of end point
%  n = max(abs(i2-i1))+1;     % number of steps
%  i = arrayfun(@(a,b)round(linspace(a,b,n)),i1,i2,'uni',0);
%  idx = sub2ind(size(BW),i{:});
%  sumBW = nnz(BW(idx));
%  disp(cell2mat(i'));        % display trajectory
%  disp(sumBW);               % display number of 1's in path
