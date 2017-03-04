%PURPOSE, compute MES
%INPUT :    data    = A t by 2 return with zero mean matrix
%           c       = A scalar which determine the cut
%           ht_m    = A t by 1 market conditional volatility matrix
%           ht_i    = A t by 1 asset conditional volatility matrix
%           rho     = A t by 1 correlation conditional volatility matrix
%                       between maket and asset
%
%OUTPUT :   MES = A t by 1 conditional MES matrix
%
% Author: Sylvain Benoit, Christophe Hurlin,    Date Revision: 17/09/2014

function [MES] = fct_MES(data,c,ht_m,ht_i,rho)

    em=(data(:,1))./ht_m; %market first column

    xi=((data(:,2)./ht_i)-rho.*em)./ sqrt(1-rho.^2); %asset second column

    bwd=1*(size(data,1)^(-0.2)); %Scaillet's bwd p21, I put 1 instead of the standard deviation because our shocks are iid with unit variance

    K1=sum(em.*(normcdf(((c./ht_m)-em)./bwd)))./(sum(normcdf(((c./ht_m)-em)./bwd)));

    K2=sum(xi.*(normcdf(((c./ht_m)-em)./bwd)))./(sum(normcdf(((c./ht_m)-em)./bwd)));

MES = (ht_i.*rho.*K1) + (ht_i.*sqrt(1-rho.^2).*K2);