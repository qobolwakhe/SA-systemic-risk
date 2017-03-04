% PURPOSE: Main Script
%       Run DCC with GJR Garch estimation (Maximum Likelihood in 2 step) and compute:
%           - MES, Delta CoVaR dcc, Delta CoVaR quantile and Delta CoVaR ols
% 
% USERS INPUTS:
%       index   = A t by 1 vector with index returns
%       asset   = A t by 1 vector with firm's equity returns
%       LTQ     = A t by 1 vector with the total amount of liabilities of the firm
%       MV      = A t by 1 vector with the market capitalisation of the firm
%       alpha   = A scalar between [0,1], risk level of our VaR
%       k       = A scalar between [0,1], Prudential Capital Required (k*LTQ), usually sets at 4% or 8%
%
% TREATMENT:
%       data_center     = A zero mean t by 2 matrix of returns, per convention the first column is the market returns and the second is the asset
%
% GJR-GARCH and DCC:
%       dcc_mvgarch: call function to estimate DCC-GJRGARCH 
%   USERS OUTPUTS:
%       ht_i            = A t by 1 Asset conditional volatility matrix
%       ht_m            = A t by 1 Market conditional volatility matrix
%       rho             = A t by 1 conditional correlation matrix
%       Asset_VaR       = A t by 1 Asset_VaR matrix with empirical quantile at alpha
%       Market_VaR      = A t by 1 Market_VaR matrix with empirical quantile at alpha
%       beta            = A t by 1 conditional beta matrix
%       loglikelihood   = The log likelihood evaluated at the optimum
%       parameters      = A vector of parameters estimated for our data of the form
%                          [GarchParams(market) GarchParams(asset) DCCParams]
%                               where the garch parameters from each estimation are of the form
%                          [omega(i) alpha(i1) beta(i1))]
%       recap_stderror  = A vector of standard errors for each parameters
%       recap_t_stat    = A vector of T-Statistic for each parameters
%
% MES:
%       fct_MES: call function to estimate MES 
%       c             = A scalar which is used to define the systemic risk event, Market HS VaR here 
%                           it's a constant(-2% in MES paper) dividing by market conditional standard deviation
%   USERS OUTPUTS:
%       MES           = A t by 1 MES matrix
%
% SRISK:
%   USERS OUTPUTS:
%       LRMES         = A t by 1 LRMES (without simulation) matrix
%       SRISK         = A t by 1 SRISK matrix
%
% CoVaR:
%   USERS OUTPUTS:
%       Delta_CoVaR_dcc     = A t by 1 Delta CoVaR dcc matrix
%       CoVaR_quant         = A t by 1 CoVaR quantile matrix
%       Delta_CoVaR_quant   = A t by 1 Delta CoVaR quantile matrix
%       Delta_CoVaR_ols     = A t by 1 Delta CoVaR ols matrix
%       gam_quant           = Slope parameter estimated by Quantile regression
%       gam_ols             = Slope parameter estimated by OLS
%
% Author: Sylvain Benoit,    Date Revision: 29/09/2015

clear
warning off all
clc
close all

[data]=xlsread('Data_RMC.xls');

index = data(:,1);            % Returns of the market (system)

asset = data(:,2);            % Returns of the firm's equity

LTQ = data(:,3);              % Total amount of liabilities

MV = data(:,4);               % Market Capitalisation

alpha = 0.05;                 % Risk level of our VaR

k = 0.08;                     % Prudential Capital Required (k*LTQ)

res = call_fct(index,asset,LTQ,MV,k,alpha);
