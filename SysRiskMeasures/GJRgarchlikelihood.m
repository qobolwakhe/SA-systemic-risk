% PURPOSE:
%     Likelihood for GJR-garch estimation
% 
% USAGE:
%     [LLF, h, likelihoods] = GJRgarchlikelihood(parameters , data , q , p, stdEstimate, T)
% 
% INPUTS:
%     parameters	=  A vector of GARCH process params + (error param) of the form [constant, arch, gamma, garch]
%     data          =  A set of zero mean residuals
%     q             =  The lag order length for ARCH
%     p             =  The lag order length for GARCH
%     stdEstimate	=  The conditional std deviation of the data
%     T             =  Length of data
% 
% OUTPUTS:
%     LLF           =  Minus 1 times the log likelihood
%     h             =  Time series of conditional volatilities
%     likelihoods   =  Time series of likelihoods
% 
% COMMENTS:
%
% Author: Sylvain Benoit    Date Revision: 17/09/2014
% Inspiration from: Kevin Sheppard 

function [LLF, h, likelihoods] = GJRgarchlikelihood(parameters , data , q , p, stdEstimate, T)

alpha=parameters(1:q);            %récupère l'initialisation de la composante ARCH;
gamma=parameters(q+1:2*q);        %récupère l'initialisation de la composante ARCH quand Epsilont < 0;
beta=parameters(2*q+1:2*q+p);     %récupère l'initialisation de la composante GARCH;
 
m = max(p,q);                       %récupère le nb de retard maximum

h=zeros(size(data));                %initialisation de Ht avec que des 0 [T by 1];

h(1:m,1)=stdEstimate^2;             %initialisation des m premières valeurs par la variance non conditionnel;

for t = (1 + m):T

    datalag = data(t-q:t-1);
    
    h(t) = (1-(sum(alpha)+0.5*sum(gamma)+sum(beta)))*var(data) + (datalag.*datalag)'*alpha + (datalag.*datalag.*(datalag<0))'*gamma + h(t-p:t-1)'*beta; %si erreur ici penser au flip left right

end;

t = (1 + m):T;

LLF  =  sum(log(h(t))) + sum((data(t).^2)./h(t));
LLF  =  0.5 * (LLF  +  (T - m)*log(2*pi));

if nargout > 2
    likelihoods = 0.5 * ((log(h(t))) + ((data(t).^2)./h(t)) + log(2*pi));
    likelihoods = -likelihoods;    
end

h=h(t);

end