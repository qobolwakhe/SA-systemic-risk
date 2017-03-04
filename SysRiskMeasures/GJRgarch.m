% PURPOSE:
%     GJR-garch(P,Q) parameter estimation with the Normal Error Distribution
% 
% USAGE:
%     [parameters, likelihood, stderrors, robustSE, ht, scores] = fattailed_garch(data , q , p , startingvals, options)
% 
% INPUTS:
%     data          = A single column of zero mean random data, normal or not for quasi likelihood
%     Q             = Non-negative, scalar integer representing a model order of the ARCH process
%     P             = Positive, scalar integer representing a model order of the GARCH process: P is the number of lags of the lagged conditional 
%                           variances included, Can be empty([]) for ARCH process
%     startingvals	= A (1+q+q+p) (plus 1 if STUDENTST or GED is selected for the nu parameter) vector of starting vals or [] => default
%     options       = Default options are use, otherwise you can provide an options vector.  See HELP OPTIMSET
% 
% OUTPUTS:
%     parameters    = a [1+q+q+p X 1] column of parameters with omega, alpha1, alpha2, ..., alpha(p), beta1, beta2, ... beta(q)
%     likelihood    = the loglikelihood evaluated at he parameters
%     robustSE      = QuasiLikelihood std errors which are robust to some forms of misspecification(see White 94)
%     stderrors     = the inverse analytical hessian, not for quasi maximum liklihood
%     ht            = the estimated time varying VARIANCES
%     scores        = The numberical scores(# fo params by t) for M testing   
% 
% 
% COMMENTS:
%   GJR-GARCH(P,Q) constraints
%     (1) Omega > 0
%     (2) Alpha(i) >= 0 for i = 1,2,...Q
%     (3) Gamma(i) >= 0 for i = 1,2,...Q    
%     (4) Beta(j)  >= 0 for j = 1,2,...P
%     (5) sum(Alpha(i) + 0.5*Gamma(i) + Beta(j)) < 1 for i = 1,2,...Q and j = 1,2,...P
%
%   The time-conditional variance, H(t), of a GJR-GARCH(P,Q) process is modeled as follows:
%
%     H(t) = Omega + Alpha(1)*r_{t-1}^2 + Alpha(2)*r_{t-2}^2 +...+ Alpha(Q)*r_{t-q}^2+...
%                   + Gamma(1)*(r_{t-1}^2)*I-_{t-1} + Gamma(2)*(r_{t-2}^2)*I-_{t-2} +...+ Gamma(Q)*(r_{t-q}^2)*I-_{t-q}+...
%                      Beta(1)*H(t-1)+ Beta(2)*H(t-2)+...+ Beta(0)*H(t-p)
%
% Modifications: Sylvain Benoit    Date Revision: 17/09/2014
% Initial codes Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 4/1/2004

function [parameters, likelihood, stderrors, robustSE, ht, scores] = GJRgarch(data , q , p , startingvals, options)

t = size(data,1);

m = max(p,q);   

if isempty(startingvals)
   alpha  =  0.05*ones(q,1)/q;                  %l'initialisation pourrait être plus propre!!
   gamma  =  0.05*ones(q,1)/q;
   beta   =  0.75*ones(p,1)/p;
else
    alpha=startingvals(1:q);
    gamma=startingvals(q+1:2*q);
    beta=startingvals(2*q+1:2*q+p);
end
                        
LB = zeros(1,q+q+p)+2*options.TolCon;     
UB = [];     
sumA =  [-eye(q+q+p) ; ones(1,q) 0.5*ones(1,q) ones(1,p)]; %constraint sumA*X <= sumB
sumB =  [zeros(q+q+p,1)+ 2*options.TolCon ; 1 - 2*options.TolCon];

startingvals = [alpha ; gamma ; beta];

% Estimate the parameters
stdEstimate = std(data,1);                      %calcul de la standard deviation
data        = [stdEstimate(ones(m,1)) ; data];  %initialisation des m premiers éléments
T           = size(data,1);

[parameters, LLF, EXITFLAG, OUTPUT, LAMBDA, GRAD] =  fmincon('GJRgarchlikelihood', startingvals ,sumA  , sumB ,[] , [] , LB , UB, [], options, data, q , p, stdEstimate, T);

hess = hessian_2sided('GJRgarchlikelihood',parameters,data,q,p,stdEstimate,T);

[likelihood, ht] = GJRgarchlikelihood(parameters,data,q,p,stdEstimate,T);
likelihood = - likelihood;

stderrors = inv(hess);                          %ou hess^(-1)

%% Compute robust STD by finite difference

if nargout > 5
   h=max(abs(parameters/2),1e-2)*eps^(1/3);
   hplus=parameters+h;
   hminus=parameters-h;
   likelihoodsplus=zeros(t,length(parameters));
   likelihoodsminus=zeros(t,length(parameters));
   for i=1:length(parameters)
      hparameters=parameters;
      hparameters(i)=hplus(i);
      [HOLDER, HOLDER1, indivlike] = GJRgarchlikelihood(hparameters, data, p, q, stdEstimate, T);
      likelihoodsplus(:,i)=indivlike;
   end
   for i=1:length(parameters)
      hparameters=parameters;
      hparameters(i)=hminus(i);
      [HOLDER, HOLDER1, indivlike] = GJRgarchlikelihood(hparameters, data, p, q, stdEstimate, T);
      likelihoodsminus(:,i)=indivlike;
   end
   scores=(likelihoodsplus-likelihoodsminus)./(2*repmat(h',t,1));
   scores=scores-repmat(mean(scores),t,1);
   B=scores'*scores;
   robustSE=stderrors*B*stderrors;
end