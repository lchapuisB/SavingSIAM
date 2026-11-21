function [Z,PI] = discretize_normal(N,mu,rho,sigma)
% Code to approximate AR(1) process using the Rouwenhorst method 

%Purpose:    Finds a Markov chain whose sample paths approximate those of
%            the AR(1) process
%                z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
%            where eps are normal with stddev sigma
%
%
%Input:      N       scalar, number of nodes for Z
%            mu      scalar, unconditional mean of process
%            rho     scalar
%            sigma   scalar, std. dev. of epsilons
%
%Output:     Z       N*1 vector, nodes for Z
%            PI      N*N matrix, transition probabilities
%


sigmaz = sigma / sqrt(1-rho^2);

p  = (1+rho)/2;
PI = [p 1-p; 1-p p];

for n = 3:N
    PI = p*[PI zeros(n-1,1); zeros(1,n)] + ...
         (1-p)*[zeros(n-1,1) PI; zeros(1,n)] + ...
         (1-p)*[zeros(1,n); PI zeros(n-1,1)] + ...
         p*[zeros(1,n); zeros(n-1,1) PI];
    PI(2:end-1,:) = PI(2:end-1,:)/2;
end

fi = sqrt(N-1)*sigmaz;
Z  = linspace(-fi,fi,N)';
Z  = Z + mu;
if(N==1)
    PI=1;
end

