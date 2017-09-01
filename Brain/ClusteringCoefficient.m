 function [Clocal,Cglobal] = ClusteringCoefficient(A) 
 % Note, that A is assumed to have no self-connections (diagonal is zero) 
 % This code also works for directed networks. 
 A = A-diag(diag(A)) % Remove any self-connections (just in case) 
 
 N = size(A,1); 
 Clocal = zeros(1,N); 
 
 % loop over all edges 
 for v=1:N 
     [Nv] = find(A(v,:) + A(:,v)'); 
     numFriends = length(Nv); 
     if (numFriends>1) 
         % Local clustering coefficient: Number of friends that are friends,
         % over the possible number of friendships between numFriends friends
         Clocal(v) = sum(sum(A(Nv,Nv)))./(numFriends ^2-numFriends); 
     end; 
     if (numFriends ==1) 
         Clocal(v) = 0; 
     end; 
 end;
 Cglobal = mean(Clocal); 