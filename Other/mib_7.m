% snakes and ladders
% 3x3 board (9 squares)
% biased coin with probability p for head
% heads advances player 1 place
% tails advances player 2 places


% Assignment:
% 1. Model the Game as a Markov process.
% 2. How many turns on average does it take to complete the game?
% Use the Fundamental matrix to solve the problem (in Matlab).
% Use Monte-Carlo simulation of the game. 


%Q1

p = 0.4;
q = 1 - p;


    % 1 2 3 4 5 6 7 8 9
    
T = [ 0 0 0 q p 1 0 0 0 ;   %1
      0 0 0 0 0 0 0 0 0 ;   %2
      0 0 0 0 0 0 0 0 0 ;   %3
      0 0 0 0 0 0 p 1 0 ;   %4
      q 0 1 p 0 0 0 0 0 ;   %5
      0 0 0 0 0 0 0 0 0 ;   %6
      p 1 0 0 q 0 0 0 0 ;   %7
      0 0 0 0 0 0 0 0 0 ;   %8
      0 0 0 0 0 0 q 0 1 ];  %9
  
  
%Q2  

W = inv((eye(9) - T'))   %and then...?