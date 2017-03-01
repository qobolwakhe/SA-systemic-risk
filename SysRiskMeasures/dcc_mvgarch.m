% PURPOSE:
%        Estimates a multivariate GARCH model using the DCC estimator of Engle and Sheppard
% 
% USAGE:
%        [parameters, loglikelihood, Ht, Rt, Qt, likelihoods, stdresid, stderrors, A,B, jointscores]...
%                =dcc_mvgarch(data,dccQ,dccP,archQ,garchP)
% 
% INPUTS:
%      data          = A zero mean t by k (asset number) vector of residuals from some filtration [t by k]
%      dccQ          = The lag length of the innovation term in the DCC estimator (a scalar)
%      dccP          = The lag length of the lagged correlation matrices in the DCC estimator (a scalar)
%      archQ         = One of two things:     A scalar, q     in which case a p innovation model is estimated for each series
%                                             A k by 1 vector in which case the ith series has innovation terms p=archP(i)
%      garchP        = One of two things:     A scalar, p     in which case a q GARCH lags is used in estimation for each series
%                                             A k by 1 vector in which case the ith series has lagged variance terms q=archQ(i)
%
% OUTPUTS:
%      parameters    = A vector of parameters estimated form the model of the form
%                          [GarchParams(1) GarchParams(2) ... GarchParams(k) DCCParams]
%                          where the Garch Parameters from each estimation are of the form
%                          [omega(i) alpha(i1) alpha(i2) ... alpha(ip(i)) beta(i1) beta(i2) ... beta(iq(i))]
%      loglikelihood = The log likelihood evaluated at the optimum
%      Ht            = A k by k by t array of conditional variances
%      Rt            = A k by k by t array of Rt elements
%      Qt            = A k by k by t array of Qt elements
%      stdresid      = A [t by k] matrix of standardized residuals
%      likelihoods   = the estimated likelihoods t by 1
%      stderrors     = A length(parameters)^2 matrix of estimated correct standard errors
%      A             = The estimated A from the robust standard errors
%      B             = The estimated B from the standard errors
%      jointscores   = The estimated scores of the likelihood t by length(parameters)
%      H             = Conditional Volatility univariate
% 
% COMMENTS:
%
% Modifications: Sylvain Benoit    Date Revision: 17/09/2014
% Initial codes Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 4/1/2004

function [parameters, loglikelihood, Ht, Rt, Qt, stdresid, likelihoods, stderrors, A,B, jointscores, H]=dcc_mvgarch(data,dccQ,dccP,archQ,garchP)

[t,k]=size(data);

if length(archQ)==1
    archQ=ones(1,k)*archQ;
end

if length(garchP)==1
    garchP=ones(1,k)*garchP;
end

disp(' ')
disp(' === First Step : Estimation of GJR-GARCH ===')

%% Now lest do the univariate garching using fattailed_GJRgarch

options  =  optimset('fmincon');
options  =  optimset(options , 'TolFun'      , 1e-006);
options  =  optimset(options , 'Display'     , 'on');
options  =  optimset(options , 'Diagnostics' , 'off');
options  =  optimset(options , 'LargeScale'  , 'off');
options  =  optimset(options , 'MaxFunEvals' , 1000*max(archQ+garchP+1));
options  =  optimset(options , 'MaxIter'     , 1000*max(archQ+garchP+1));  
options  =  optimset(options , 'MaxSQPIter'  , 1000);

for i=1:k
    disp(' ')
    fprintf(1,' Estimating GJR-GARCH model for Series %d\n',i)
    
    [univariate{i}.parameters, univariate{i}.likelihood, univariate{i}.stderrors, univariate{i}.robustSE, univariate{i}.ht, univariate{i}.scores] ... 
        = GJRgarch( data(:,i) , archQ(i) , garchP(i) , [], options);
    
    stdresid(:,i)=data(:,i)./sqrt(univariate{i}.ht); %calcul des résidus standardisés!

end


%% Now lest do the multivariate garching using dcc_mvgarch_likelihood


options_base  =  optimset('fmincon');
options  =  optimset(options_base,'Algorithm','sqp',...
    'LargeScale', 'off', 'Display', 'on', 'Diagnostics', 'off');

dccstarting=[ones(1,dccQ)*.01/dccQ ones(1,dccP)*.97/dccP];

disp(' ')
disp(' === Second Step : Estimation of the DCC ===')
disp(' ')

[dccparameters,dccllf,EXITFLAG,OUTPUT,LAMBDA,GRAD]=fmincon('dcc_mvgarch_likelihood',dccstarting, ones(size(dccstarting)),[1-2*options.TolCon],[],[],...
    zeros(size(dccstarting))+2*options.TolCon,[],[],options,stdresid,dccQ,dccP);


%% We now have all of the estimated parameters

parameters=[];
H=zeros(t,k);
for i=1:k
    parameters=[parameters;univariate{i}.parameters];
    H(:,i)=univariate{i}.ht;
end
parameters=[parameters;dccparameters'];


%% We now have Ht and the likelihood but we want the FULL LIKELIHOOD

% options  =  optimset('fminsearch');
% options  =  optimset(options , 'Display'            , 'iter');
% options  =  optimset(options , 'Diagnostics'        , 'on');
% options  =  optimset(options , 'LevenbergMarquardt' , 'on');
% options  =  optimset(options , 'LargeScale'         , 'off');
% options  =  optimset(options , 'MaxIter'            , 1);
% % à utiliser avec parcimonie car on ne maîtrise plus les contraintes!
% [parameters,loglikelihood,exitflag]=fminsearch('dcc_mvgarch_full_likelihood', parameters, options, data, archQ, garchP, dccQ, dccP);
% [loglikelihood, Rt, likelihoods, Qt]=dcc_mvgarch_full_likelihood(parameters, data, archQ, garchP, dccQ, dccP);

[loglikelihood, Rt, likelihoods, Qt]=dcc_mvgarch_full_likelihood(parameters, data, archQ, garchP, dccQ, dccP);

likelihoods=-likelihoods;                                                   %output;
loglikelihood=-loglikelihood;                                               %output;
Ht=zeros(k,k,t);
Hstd=H.^(0.5);
for i=1:t
    Ht(:,:,i)=diag(Hstd(i,:))*Rt(:,:,i)*diag(Hstd(i,:));
    %stdresid(i,:)=data(i,:)*Ht(:,:,i)^(-0.5);
end
save tempHt Ht
clear Ht


if nargout >=7
    % How are we going to get STD errors?  Partitioned invers probably.  Well, we need to get the scores form the dcc model, the joint likelihood.
    % We then need to get A12 and A22 so we can have it all.  We also need to get A11 in the correct form.
    A=zeros(length(parameters),length(parameters));
    index=1;
    for i=1:k
        workingsize=size(univariate{i}.stderrors);
        A(index:index+workingsize-1,index:index+workingsize-1)=univariate{i}.stderrors^(-1); %début de la construction de la matrice A
        index=index+workingsize;
    end
   
    % Ok so much for a A11 and A12 and A22, as we have them all between whats above
%    fprintf(1,'\n\nCalculating Standard Errors, this can take a while\n');
    otherA=dcc_hessian('dcc_mvgarch_full_likelihood',parameters, dccQ + dccP, data, archQ, garchP, dccQ, dccP); %calcul la matrice hessienne pour le DCC
    A(length(parameters) - dccQ - dccP + 1:length(parameters),:)=otherA; %met la matrice A en bas à droite de la matrice A existante
    
    % We now need to get the scores for the DCC estimator so we can finish B
    jointscores=zeros(t,length(parameters));
    index=1;
    for i=1:k
        workingsize=size(univariate{i}.scores,2);
        jointscores(:,index:index+workingsize-1)=univariate{i}.scores;
        index=index+workingsize;
    end
    
    % Now all we need to do is calculate the scores form the dcc estimator and we have everything
    h=max(abs(parameters/2),1e-2)*eps^(1/3);
    hplus=parameters+h;
    hminus=parameters-h;
    likelihoodsplus=zeros(t,length(parameters));
    likelihoodsminus=zeros(t,length(parameters));
    for i=length(parameters)-dccQ-dccP+1:length(parameters)
        hparameters=parameters;
        hparameters(i)=hplus(i);
        [HOLDER, HOLDER1, indivlike] = dcc_mvgarch_full_likelihood(hparameters, data, archQ, garchP, dccQ, dccP);
        likelihoodsplus(:,i)=indivlike;
    end
    for i=length(parameters)-dccQ-dccP+1:length(parameters)
        hparameters=parameters;
        hparameters(i)=hminus(i);
        [HOLDER, HOLDER1, indivlike] = dcc_mvgarch_full_likelihood(hparameters, data, archQ, garchP, dccQ, dccP);
        likelihoodsminus(:,i)=indivlike;
    end
    DCCscores=(likelihoodsplus(:,length(parameters)-dccQ-dccP+1:length(parameters))-likelihoodsminus(:,length(parameters)-dccQ-dccP+1:length(parameters)))...
        ./(2*repmat(h(length(parameters)-dccQ-dccP+1:length(parameters))',t,1));
    jointscores(:,length(parameters)-dccQ-dccP+1:length(parameters))=DCCscores;
    B=cov(jointscores);
    A=A/t;
    stderrors=A^(-1)*B*A'^(-1)*t^(-1);
end
% Done!
load tempHt