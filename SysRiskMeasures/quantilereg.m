%PURPOSE, compute Delta CoVaR (= < VaR)
%INPUT :    Y       = A t by 1 matrix of dependant variable
%           X       = A t by 1 matrix of independant variable
%           theta   = the sample quantile [0,1]
%           Ht      = A k by k by t array of conditional variances covariances
%           ht_m    = A t by 1 market conditional volatility matrix
%           rho     = A t by 1 correlation conditional volatility matrix
%                       between maket and asset
%
%OUTPUT :   beta    = Estimated Coefficients
%           tstats  = T-Students of the coefficients
%
% Author: Sylvain Benoit    Date Revision: 17/09/2014

function [beta tstats]=quantilereg(y,x,theta)

ry=length(y);       %number of rows of dependant variable (ry=rx)
[rx cx]=size(x);    %number of rows and columns of independant variable
x=[ones(rx,1) x];   %add a constant to the independant variable
cx=cx+1;            %number of independant variables

% Finding first estimates by solving the system

itrat=0;
xstar=x;
diff=1;
beta=ones(cx,1);
while itrat<1000 && diff>1e-6
    itrat=itrat+1;
    beta0=beta;
    beta=inv(xstar'*x)*xstar'*y;
    resid=y-x*beta;
    resid(abs(resid)<.000001)=.000001;
    resid(resid<0)=theta*resid(resid<0);
    resid(resid>0)=(1-theta)*resid(resid>0);
    resid=abs(resid);
    z=[];
for i=1:cx 
    z0 = x(:,i)./resid;
    z=[z z0];
end

xstar=z;
beta1=beta;
diff=max(abs(beta1-beta0));

end

e=y-x*beta;

%Estimating variances based on Green 2008(quantile regression)

iqre=iqr(e);
if theta==0.5
  h=0.9*std(e)/(ry^0.2);
else
  h=0.9*min(std(e),iqre/1.34)/(ry^0.2);
end
u=(e/h);
fhat0=(1/(ry*h))*(sum(exp(-u)./((1+exp(-u)).^2)));
D(ry,ry)=0;
DIAGON=diag(D);
DIAGON(e>0)=(theta/fhat0)^2;
DIAGON(e<=0)=((1-theta)/fhat0)^2;
D=diag(DIAGON);
VCQ=(x'*x)^(-1)*x'*D*x*(x'*x)^(-1);


%Standard errors and t-stats

tstats=beta./diag(VCQ).^0.5;
stderrors=diag(VCQ).^0.5;
PValues=2*(1-tcdf(abs(tstats),ry-cx));

