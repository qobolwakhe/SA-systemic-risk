% PURPOSE:
%        Full likelihood for use in the DCC_MVGARCH estimation and
%        returns the likelihood of the QMLE estimates of the DCC parameters
% 
% USAGE:
%        [logL, Rt, likelihoods, Qt, Zt]=dcc_mvgarch_full_likelihood(parameters, data, archQ, garchP, dccQ, dccP)
% 
% INPUTS:
%    parameters  = A k+sum(archQ)+sum(garchP)+dccQ+dccP vector of parameters of the form
%                  [GarchParams(1) GarchParams(2) ... GarchParams(k) DCCParams]
%                  where the garch parameters from each estimation are of the form
%                  [omega(i) alpha(i1) alpha(i2) ... alpha(ip(i)) beta(i1) beta(i2) ... beta(iq(i))]
%                  and DCCparams are [DCCa DCCb]
%    data        = A t by k matrix of zero mean residuals
%    archQ       = A vector of arch innovation lag lengths
%    garchP      = A vector of Garch AR lag lengths
%    dccQ        = A scalar of the DCC innovations lag length
%    dccP        = A scalar of the DCC AR lag lengths
% 
% OUTPUTS:
%    logL        = The estimated log likelihood
%    Rt          = The estimates covariances
%    likelihoods = The likelihoods (t by 1)
%    Qt          = A k by k by t array of Qt elements
% 
% COMMENTS: pour que les estimateurs soient fully efficient
%  
% Modifications: Sylvain Benoit    Date Revision: 17/09/2014
% Initial codes Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 4/1/2004


function [logL, Rt, likelihoods, Qt]=dcc_mvgarch_full_likelihood(parameters, data, archQ,garchP,dccQ,dccP)

% First we need to make the T by K matrix of variances.
[t,k]=size(data);
H=zeros(size(data));            %initialisation de Ht avec que des 0 [T by 1];

index=1;

for i=1:k
    
    alpha = parameters(index:index+archQ(i)-1);
    gamma = parameters(index+archQ(i):index+archQ(i)+archQ(i)-1);
    beta  = parameters(index+archQ(i)+archQ(i):index+archQ(i)+archQ(i)+garchP(i)-1);
                      
    m = max(archQ(i),garchP(i));         

    h=zeros(size(data,1),1);                   %initialisation de Ht avec que des 0 [T by 1];
    h(1:m,1)=std(data(:,i))^2;                 %initialisation des m premières valeurs par la variance non conditionnel;
    
%écriture de ma variance conditionnelle

for T = (1 + m):t

    datalag = data(T-archQ(i):T-1,i);
    
    h(T,1) = (1-(sum(alpha)+0.5*sum(gamma)+sum(beta)))*var(data(:,1)) + (datalag.*datalag)'*alpha + (datalag.*datalag.*(datalag<0))'*gamma + h(T-garchP(i):T-1,1)'*beta; %si erreur ici penser au flip left right

end;
     
    H(:,i)=h;
    
    index=index+archQ(i)+archQ(i)+garchP(i); 
        
end

stdresid=data./sqrt(H);

%First compute Qbar, the unconditional Correlation Matrix

Qbar=cov(stdresid);

stdresid=[ones(max(dccQ,dccP),k);stdresid]; %on rajoute

a=parameters(index:index+dccQ-1);
b=parameters(index+dccQ:index+dccQ+dccP-1);
sumA=sum(a);
sumB=sum(b);

% Next compute Qt
m=max(dccQ,dccP);

Qt=zeros(k,k,t+m);
Qt(:,:,1:m)=repmat(Qbar,[1 1 m]);
Rt=zeros(k,k,t+m);
Rt(:,:,1:m)=repmat(Qbar,[1 1 m]);

logL=0;
likelihoods=zeros(t+m,1);
H=[zeros(m,k);H];
Q=dccQ;
P=dccP;

for j=(m+1):t+m
      Qt(:,:,j)=Qbar*(1-sumA-sumB);
   for i=1:Q
      Qt(:,:,j)=Qt(:,:,j) + a(i)*(stdresid(j-i,:)'*stdresid(j-i,:));
   end
   for i=1:P
      Qt(:,:,j)=Qt(:,:,j) + b(i)*Qt(:,:,j-i);
   end
   Rt(:,:,j)=Qt(:,:,j)./(sqrt(diag(Qt(:,:,j)))*sqrt(diag(Qt(:,:,j)))');
   likelihoods(j) = k*log(2*pi) + sum(log(H(j,:))) + log(det(Rt(:,:,j))) + stdresid(j,:)*inv(Rt(:,:,j))*stdresid(j,:)'; %vraisemblance complète!
   logL=logL + likelihoods(j);
end;

Qt=Qt(:,:,(m+1:t+m));           %output de Qt matrix [k by k by t]
Rt=Rt(:,:,(m+1:t+m));
logL=(1/2)*logL;
likelihoods=(1/2)*likelihoods(m+1:t+m);