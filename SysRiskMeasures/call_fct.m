% PURPOSE: Main Script
%       Run DCC with GJR Garch estimation (Maximum Likelihood in 2 step) and compute:
%           - MES, Delta CoVaR dcc, Delta CoVaR quantile and Delta CoVaR ols
% 
% USERS INPUTS:
%       index   = A t by 1 vector with index returns
%       asset   = A t by 1 vector with asset returns
%       alpha   = A scalar between [0,1], risk level of our VaR   
%
% TREATMENT:
%       data_center     = A zero mean t by 2 matrix of returns, per convention the first column is the market returns and the second is the asset
%
% GJR-GARCH and DCC:
%       dcc_mvgarch: call function to estimate DCC-GJRGARCH (Kevin Sheppard)
%
%   USERS OUTPUTS:
%       res.ht_i            = A t by 1 Asset conditional volatility matrix
%       res.ht_m            = A t by 1 Market conditional volatility matrix
%       res.rho             = A t by 1 conditional correlation matrix
%       res.Asset_VaR       = A t by 1 Asset_VaR matrix with empirical quantile at alpha
%       res.Market_VaR      = A t by 1 Market_VaR matrix with empirical quantile at alpha
%       res.beta            = A t by 1 conditional beta matrix
%       res.loglikelihood   = The log likelihood evaluated at the optimum
%       res.parameters      = A vector of parameters estimated for our data of the form
%                          [GarchParams(market) GarchParams(asset) DCCParams]
%                               where the garch parameters from each estimation are of the form
%                          [omega(i) alpha(i1) beta(i1))]
%       res.recap_stderror  = A vector of standard errors for each parameters
%       res.recap_t_stat    = A vector of T-Statistic for each parameters
%
% MES:
%       fct_MES: call function to estimate MES 
%       c             = A scalar which is used to define the systemic risk event, Market HS VaR here 
%                           it's a constant(-2% in MES paper) dividing by market conditional standard deviation
%   USERS OUTPUTS:
%       res.MES           = A t by 1 MES matrix
%
% SRISK:
%   USERS OUTPUTS:
%       SRISK         = A t by 1 SRISK matrix
%
% CoVaR:
%   USERS OUTPUTS:
%       res.Delta_CoVaR_dcc     = A t by 1 Delta CoVaR dcc matrix
%       res.CoVaR_quant         = A t by 1 CoVaR quantile matrix
%       res.Delta_CoVaR_quant   = A t by 1 Delta CoVaR quantile matrix
%       res.Delta_CoVaR_ols     = A t by 1 Delta CoVaR ols matrix
%       res.gam_quant           = Slope parameter estimated by Quantile regression
%       res.gam_ols             = Slope parameter estimated by OLS
%
% Author: Sylvain Benoit    Date Revision: 29/09/2015


function [res] = call_fct(index,asset,LTQ,MV,k,alpha)

warning off

%% Data treatment

data = [index asset];

data = data(~isnan(data(:,2).*data(:,1)),:); % select only elements which are not NaN

data_center = data - ones(size(data,1),1)*mean(data); %demeaned returns

%% GJR-GARCH and DCC

[parameters, loglikelihood, Ht, Rt, Qt, stdresid, likelihoods, stderrors, A,B, jointscores, H]=dcc_mvgarch(data_center,1,1,1,1);

ht_m=sqrt(H(:,1)); %market conditional volatility

ht_i=sqrt(H(:,2)); %asset conditional volatility

rho=squeeze(Rt(1,2,:)); %conditional correlation

stderrors2=(stderrors.*eye(size(parameters,1),size(parameters,1)));

recap_stderror=stderrors2(stderrors2~=0);

recap_t_stats(:,1)=parameters./(stderrors2(stderrors2>0));

%Conditional Asset VaR from dcc with Empirical Quantile at alpha
    Asset_VaR = ht_i*quantile(data_center(:,2)./ht_i,alpha);
        
%Conditional Market VaR from dcc with Empirical Quantile at alpha
    Market_VaR = ht_m*quantile(data_center(:,1)./ht_m,alpha);

%Conditional beta
    Beta = rho.*ht_i./ht_m;

    
%% Descriptives Statistics & Graph

%mean, std deviation, Skewness, Kurtosis, Max, Min,

%% Tails NonParametric estimation and MES (expected losses (in %) of the asset return when market return goes down of more than c%); Eq 9

disp(' ')
disp(' === Third Step : Estimation of the MES ===')
disp(' ')

%attention au cut!!!
c = quantile(data_center(:,1),alpha); % HS VaR (nonparametric), it's our systemic event, it's a scalar here

MES = - fct_MES(data_center,c,ht_m,ht_i,rho);

disp(sprintf('   Average MES = %3.4f',mean(MES))), disp(' ')

%% Delta CoVaR DCC when asset return is equal to his VaR; Eq 24

disp(' ')
disp(' === Fourth Step : SRISK computation ===')
disp(' ')

    LRMES = (1-exp(-18.*MES)); %without simulation
    
disp(sprintf('   Average LRMES = %3.4f',mean(LRMES))), disp(' ')

    SRISK = k.*LTQ - (1-k).*(1-LRMES).*MV;

disp(sprintf('   Average SRISK = %3.4f',mean(SRISK))), disp(' ')

%% Delta CoVaR DCC when asset return is equal to his VaR; Eq 24

disp(' ')
disp(' === Fifth Step : Estimation of CoVaR and DCoVaR ===')
disp(' ')

    Delta_CoVaR_dcc = - ( rho.*ht_m./ht_i.*Asset_VaR...
        - rho.*ht_m./ht_i.*repmat(median(data_center(:,2)),size(data_center,1),1) );

    disp(sprintf('   Average DCoVaR (DCC) = %3.4f',mean(Delta_CoVaR_dcc))), disp(' ')

%% Delta CoVaR Quantile when asset return is equal to his VaR; Eq 23

[beta tstats]=quantilereg(data_center(:,1),data_center(:,2),alpha); %Quant with constant

gam_quant = beta(2);

CoVaR_quant = - ( beta(1) + beta(2).*Asset_VaR );

    Delta_CoVaR_quant = - ( beta(2).*Asset_VaR ...
        - beta(2).*repmat(median(data_center(:,2)),size(data_center,1),1) );

        disp(sprintf('   Average DCoVaR (quantile regression) = %3.4f',mean(Delta_CoVaR_quant))), disp(' ')
        
%% Delta CoVaR OLS when asset return is equal to his VaR; Eq 25

    gam_ols = [data_center(:,2)] \ data_center(:,1); %OLS without constant (r_mt = gam_ols *r_it) cause demeaned returns
        
    Delta_CoVaR_ols = - ( gam_ols.*Asset_VaR ...
        - gam_ols.*repmat(median(data_center(:,2)),size(data_center,1),1) );
    
    disp(sprintf('   Average DCoVaR (OLS regression) = %3.4f',mean(Delta_CoVaR_ols))), disp(' ')
    
%% Clear and treatments on OTHERS OUTPUTS

clear Ht Rt Qt stdresid likelihoods stderrors A B jointscores H stderrors2 c beta tstats

Market_VaR = - Market_VaR;

Asset_VaR = - Asset_VaR;

%% Affectations 
res.ht_i=ht_i;
res.ht_m=ht_m;
res.rho=rho;
res.Asset_VaR=Asset_VaR;
res.Market_VaR=Market_VaR;
res.Beta=Beta;
res.loglikelihood=loglikelihood;
res.parameters=parameters;
res.recap_stderror=recap_stderror;
res.recap_t_stats=recap_t_stats;
res.MES=MES;
res.LRMES=LRMES;
res.SRISK=SRISK;
res.Delta_CoVaR_dcc=Delta_CoVaR_dcc;
res.CoVaR_quant=CoVaR_quant;
res.Delta_CoVaR_quant=Delta_CoVaR_quant;
res.Delta_CoVaR_ols=Delta_CoVaR_ols;
res.gam_quant=gam_quant;
res.gam_ols=gam_ols;


%% Graphs

x=linspace(1,size(asset,1),size(asset,1))'; %x axis!
    
close all
%Asset Conditional Volatility
f1=figure;

subplot(3,1,1)
   plot(ht_i);
            title('Conditional volatility of the asset returns','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[0 (max(ht_i)+0.1*max(ht_i))],...
            'ytick',linspace(0,(max(ht_i)+0.1*max(ht_i)),5));

%Market Conditional Volatility
subplot(3,1,2)
    plot(ht_m);
            title('Conditional volatility of the market returns','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[0 (max(ht_m)+0.1*max(ht_m))],...
            'ytick',linspace(0,(max(ht_m)+0.1*max(ht_m)),5));
            
%Conditional Correlation
subplot(3,1,3)
    plot(rho);
            title('Correlation market/asset returns','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[0 (max(rho)+0.1*max(rho))],...
            'ytick',linspace(0,(max(rho)+0.1*max(rho)),5));

%Asset VaR & Demeaned Returns
f2=figure; 
subplot(2,1,1)
    plot(x,data_center(:,2),'b',x,-Asset_VaR,'r');
            title('VaR of asset returns','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('VaR & Demeaned Returns','fontsize',8,'fontweight','b','color','k');
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[(min(data_center(:,2))+0.1*min(data_center(:,2))) (max(data_center(:,2))+0.1*max(data_center(:,2)))],...
            'ytick',linspace((min(data_center(:,2))+0.1*min(data_center(:,2))),(max(data_center(:,2))+0.1*max(data_center(:,2))),5));         
        legend('Returns','Value-at-Risk','location','northwest')

%Market VaR & Demeaned Returns
subplot(2,1,2)
plot(x,data_center(:,1),'b',x,-Market_VaR,'r');
            title('VaR of market returns','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('VaR & Demeaned Returns','fontsize',8,'fontweight','b','color','k');
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[(min(data_center(:,1))+0.1*min(data_center(:,1))) (max(data_center(:,1))+0.1*max(data_center(:,1)))],...
            'ytick',linspace((min(data_center(:,1))+0.1*min(data_center(:,1))),(max(data_center(:,1))+0.1*max(data_center(:,1))),5));             
        legend('Returns','Value-at-Risk','location','northwest')           
            
%MES
f3=figure;
    plot(x,MES);
            title('MES','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('%','fontsize',8,'fontweight','b','color','k');
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[0 (max(MES)+0.1*max(MES))],...
            'ytick',linspace(0,(max(MES)+0.1*max(MES)),5));
            
%SRISK
f4=figure;
    plot(x,SRISK);
            title('SRISK','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('million','fontsize',8,'fontweight','b','color','k');
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[0 (max(SRISK)+0.1*max(SRISK))],...
            'ytick',linspace(0,(max(SRISK)+0.1*max(SRISK)),5));
        
%Delta_CoVaR_dcc
f5=figure;
subplot(2,1,1)
    plot(x,Delta_CoVaR_dcc);
            title('\DeltaCoVaR (DCC)','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('%','fontsize',8,'fontweight','b','color','k');           
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[0 (max(Delta_CoVaR_dcc)+0.1*max(Delta_CoVaR_dcc))],...
            'ytick',linspace(0,(max(Delta_CoVaR_dcc)+0.1*max(Delta_CoVaR_dcc)),5));
            
%Delta_CoVaR_quantile
subplot(2,1,2)
    plot(x,Delta_CoVaR_quant);
            title('\DeltaCoVaR (quantile regression)','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('%','fontsize',8,'fontweight','b','color','k');
            set(gca,'xlim',[min(x) max(x)],'xtick',linspace(min(x),max(x),4),...
            'ylim',[0 (max(Delta_CoVaR_quant)+0.1*max(Delta_CoVaR_quant))],...
            'ytick',linspace(0,(max(Delta_CoVaR_quant)+0.1*max(Delta_CoVaR_quant)),5));

        
%Conditional MES VS Delta_CoVaR_dcc
f6=figure;         
        [A,H1A,H2A] = plotyy(x,MES,x,Delta_CoVaR_dcc);
            title('Comparison MES and \DeltaCoVaR (DCC)','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            set(get(A(1),'Ylabel'),'String','MES','fontsize',8,'fontweight','b','color','b')
            set(get(A(2),'Ylabel'),'String','DCC-\DeltaCoVaR','fontsize',8,'fontweight','b','color','r')
            set(A(1),'ylim',[0 (max(MES)+0.1*max(MES))],'ytick',linspace(0,(max(MES)+0.1*max(MES)),5),...
                'xlim',[min(x) max(x)], 'xtick',linspace(min(x),max(x),4))
            set(A(2),'ylim',[0 (max(Delta_CoVaR_dcc)+0.1*max(Delta_CoVaR_dcc))],'ytick',linspace(0,(max(Delta_CoVaR_dcc)+0.1*max(Delta_CoVaR_dcc)),5),...
                'xlim',[min(x) max(x)], 'xtick',[])                    
            legend('MES','\DeltaCoVaR (DCC)','location','northwest')