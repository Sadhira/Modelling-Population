function [x,y] = GenerateRingCoordinates(numNodes)


for ll=0:(numNodes-1)
x(ll+1) = cos(2*pi*ll/numNodes);
y(ll+1) = sin(2*pi*ll/numNodes);
end

