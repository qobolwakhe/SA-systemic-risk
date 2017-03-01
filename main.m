% PURPOSE: Main Script
%       Run DCC with GJR Garch estimation and compute:
%           - MES and SRISK
%         Required functions available for download at:
%         http://www.runmycode.org/companion/view/175
%           Author: Sylvain Benoit, Christophe Hurlin
%         https://www.kevinsheppard.com/MFE_Toolbox
%           Author: Kevin Sheppard
% REQUIRED INPUTS:
%       Dta     = A t by 2 vector with index returns and firm's equity returns
%       TOTL    = A t by 1 vector with the total amount of liabilities of the firm
%       MCAP    = A t by 1 vector with the market capitalisation of the firm
%       alpha   = A scalar between [0,1], risk level of our VaR
%       k       = A scalar between [0,1], Prudential Capital Required (k*LTQ), usually sets at 4% or 8%
%       M       = Order of symmetric innovations in DCC model
%       L       = Order of asymmetric innovations in ADCC model
%       N       = Order of lagged correlation in DCC model
%       P       = Positive, scalar integer representing the number of symmetric innovations in the
%                 univariate volatility models
%       O       = Non-negative, scalar integer representing the number of asymmetric innovations in the
%                 univariate volatility models
%       Q       = Non-negative, scalar integer representing the number of conditional covariance lags in
%                 the univariate volatility models
%       GJRTYPE = Either 1 (TARCH/AVGARCH) or 2 (GJR-GARCH/GARCH/ARCH)


%       Data available from Thompson Reuters Eikon. Spreadsheet attached is
%       linked to the application and will automatically update when
%       connected. User is also to update the fields and obtain alternate
%       time series.
%
% OUTPUT
%       PARAMETERS   - Estimated parameters.  Output depends on METHOD.
%                    3-stage: [VOL(1) ... VOL(K) corr_vech(R)' vech(N)' alpha gamma beta]
%                    2-stage: [VOL(1) ... VOL(K) corr_vech(R)' alpha gamma beta]
%                    where VOL(j) is a (1+P(i)+O(i)+Q(i)) vector containing the parameters from
%                    volatility model i.
%       LL           - The log likelihood at the optimum
%       HT           - A [2 2 t] dimension matrix of conditional covariances
%       VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)/T)
%       SCORES       - A t by numParams matrix of individual scores
%       DIAGNOSTICS  - A structure containing further outputs.
%                    gScale - the scale used in estimation of asymmetric models
%       MES          -A t by 1 conditional MES matrix
%       SRISK        -A t by 1 SRISK matrix
% 
%  Script also outputs a csv file in the directory OUT\ with the dates, MES
%  and SRISK for each individual firm



%   The dynamics of a the correlations in a DCC model are:
%     3-stage:
%     Q(t) = R*(1-sum(a)-sum(b))-sum(g)*N + a(1)*e(t-1)'*e(t-1) + ... + a(m)*e(t-m)'*e(t-m)
%     + g(1)*v(t-1)'*v(t-1) + ... + g(l)*v(t-l)*v(t-l) + b(1)*Q(t-1) + ... + b(n)*Q(t-1)
%

%Author: Qobolwakhe Dube

%% Add functions to path:
%       Point directory to location of extracted files

addpath(genpath('C:\Users\Lwakhe\Dropbox\Project 3\RunMyCode_SysRiskMeasures_R2015b'))
addpath(genpath('C:\Users\Lwakhe\Dropbox\Project 3\MFE'))


mkdir('OUT')                                                    %Output will be saved here

%% Reading in data:
%   NB!! First check for duplicate rows wrt dates; after updating
%   spreadsheet, delete rows marked as duplicates in first column
%       this is a technical problem that will sometimes happen when
%       updating reuters data

data.Price = xlsread('MES_data.xlsx','Share prices','C5:R4320');
data.Index = xlsread('MES_data.xlsx','Index','C3:C4437');
data.Liabilities = xlsread('MES_data.xlsx','Liabilities','D3:S96');
data.MarkCap = xlsread('MES_data.xlsx','Market Cap','C3:R4330');
data.Equity = xlsread('MES_data.xlsx','Equity','C3:R95');
[~,~,Series] = xlsread('MES_data.xlsx','Share prices','C3:R4');                           % Cell containing tickers and names of respective equities
[~, ~, raw_dates, dates] = xlsread('MES_data.xlsx','Share prices','B5:B4320','',@convertSpreadsheetExcelDates);
dates = dates(:,1);
dates = datetime([dates{:,1}].', 'ConvertFrom', 'Excel');
n = length(Series);

    %cleaning data
[~,~,duplicates.Index] = xlsread('MES_data.xlsx','Index','A3:B4437');
df = strcmp('Duplicate',duplicates.Index);
data.Index(df(:,1))=[];
duplicates.Index(df(:,1),:)=[];


[~,~,duplicates.Price] = xlsread('MES_data.xlsx','Share prices','A5:B4320');
df = strcmp('Duplicate',duplicates.Price);
data.Price(df(:,1),:)=[];
dates(df(:,1))=[];
duplicates.Price(df(:,1),:)=[];


[~,~,duplicates.MarkCap] = xlsread('MES_data.xlsx','Market Cap','A3:B4330');
df = strcmp('Duplicate',duplicates.MarkCap);
data.MarkCap(df(:,1),:)=[];
duplicates.MarkCap(df(:,1),:)=[];


[f,g,h]=intersect(duplicates.Price(:,2),duplicates.MarkCap(:,2));
data.Price=data.Price(g,:);
dates=dates(g,:);
duplicates.Price=duplicates.Price(g,:);
data.MarkCap=data.MarkCap(h,:);

[f,g,h]=intersect(duplicates.Price(:,2),duplicates.Index(:,2));
data.Price=data.Price(g,:);
dates=dates(g,:);
data.MarkCap=data.MarkCap(g,:);
data.Index=data.Index(h,:);

% Assumption wrt liabilities is that value of debt at the last
% reported date 


[~, ~, tempraw_dates, tempdates] = xlsread('MES_data.xlsx','Liabilities','C3:C96','',@convertSpreadsheetExcelDates);
tempdates = tempdates(:,1);
tempdates = datetime([tempdates{:,1}].', 'ConvertFrom', 'Excel');
temp=[];

for i=1:n
    df=find(~isnan(data.Liabilities(:,i)));
    in=length(df);
    
    
    if dates(1)<tempdates(df(1))
        temp1=nan(sum(dates<tempdates(df(1))),1);
    else
        temp1=repmat(data.Liabilities(df(1),i),sum(dates<tempdates(df(1))),1);
    end

    for j=2:in
        temp1a=repmat(data.Liabilities(df(j-1),i),sum(dates<tempdates(df(j))& dates>=tempdates(df(j-1))),1);
        temp1=[temp1;temp1a];
    end
   temp1a=repmat(data.Liabilities(df(j),i),sum(dates>=tempdates(df(j))),1);
   temp1=[temp1;temp1a];
   temp=[temp temp1];
end

LTL=temp;

[~, ~, tempraw_dates, tempdates] = xlsread('MES_data.xlsx','Equity','B3:B95','',@convertSpreadsheetExcelDates);
tempdates = tempdates(:,1);
tempdates = datetime([tempdates{:,1}].', 'ConvertFrom', 'Excel');
temp=[];

for i=1:n
    df=find(~isnan(data.Equity(:,i)));
    in=length(df);
    
    
    if dates(1)<tempdates(df(1))
        temp1=nan(sum(dates<tempdates(df(1))),1);
    else
        temp1=repmat(data.Equity(df(1),i),sum(dates<tempdates(df(1))),1);
    end

    for j=2:in
        temp1a=repmat(data.Equity(df(j-1),i),sum(dates<tempdates(df(j))& dates>=tempdates(df(j-1))),1);
        temp1=[temp1;temp1a];
    end
   temp1a=repmat(data.Equity(df(j),i),sum(dates>=tempdates(df(j))),1);
   temp1=[temp1;temp1a];
   temp=[temp temp1];
end

EQ = temp ;


%% Computing log returns

Returns.Firm = data.Price(2:end,:)./data.Price(1:end-1,:);
Log_ret.Firm = log(Returns.Firm);

Returns.Index = data.Index(2:end,:)./data.Index(1:end-1,:);
Log_ret.Index = log(Returns.Index);

%% Computing each firm's MES and SRISK

M = 1;
L = 0;
N = 1;
P = 1;
O = 1;
Q = 1;
GJRTYPE = 2;

alpha= 0.05;
k = 0.08;

for i= 1:n
    
    TOTL =  LTL(2:end,i);
    MCAP = data.MarkCap(2:end,i);
    Dta = [Log_ret.Index Log_ret.Firm(:,i) ];
    
    R = [Dta TOTL MCAP];
    
    %deleting rows with NaN entries
    b = dates(2:end,1);                     
    b(any(isnan(R),2),:) = [];  
    eval(strcat('x.',Series{2,i},'=b;')); %dates of remaining observations for series i
    
    R(any(isnan(R),2),:) = []; 
    
    Dta = R(:,1:2);
    TOTL = R(:,3);
    MCAP = R(:,4);
             
    
    
    % covariance and volatility estimation 
    [PARAMETERS,LL,HT,VCV,SCORES,DIAGNOSTICS] = dcc(Dta,[],M,L,N,P,O,Q,GJRTYPE,'3-stage');
    
    
    %volatility model parameters
    eval(strcat('params.',Series{2,i},'.omega = PARAMETERS(5);'));
    eval(strcat('params.',Series{2,i},'.beta = PARAMETERS(8);'));
    eval(strcat('params.',Series{2,i},'.alpha = PARAMETERS(6);'));
    eval(strcat('params.',Series{2,i},'.lambda = PARAMETERS(7);'));
    
    
    rho=squeeze(HT(1,2,:))./sqrt(squeeze(HT(1,1,:)).*squeeze(HT(2,2,:)));  % conditional correlation matrix
    ht_i = sqrt(squeeze(HT(1,1,:)));                                       % firm conditional volatility matrix
    ht_m = sqrt(squeeze(HT(2,2,:)));                                       % market conditional volatility matrix
    
    
    c = quantile(Dta(:,1),alpha);
    
    %     Compute MES
    
    eval(strcat('MES.',Series{2,i},' = -fct_MES(Dta,c,ht_m,ht_i,rho);'));
    
    %     Compute SRISK
    
    eval(strcat('LRMES = (1-exp(-18.*MES.',Series{2,i},'));')) %without simulation assumption from Acharya et al 2012
    eval(strcat('SRISK.',Series{2,i},' = k.*TOTL - (1-k).*(1-LRMES).*MCAP;'));
    
end

%SRISK contribution
z= zeros(length(dates),n);

for i=1:n
    eval(strcat('g = x.',Series{2,i},';'));
    [E, ia, ib] = intersect(dates,g);
    
    eval(strcat('h = SRISK.',Series{2,i},';'));
    z(ia,i)= h;
end
z(z<0)=0;
SRISK_tot=sum(z,2);
D=z./repmat(SRISK_tot,1,n);

% Aligning dates
for i=1:n
    eval(strcat('g = x.',Series{2,i},';'));
    [E, ia, ib] = intersect(dates,g);
    
    w = D(:,i);
    w=w(ia);
    eval(strcat('SRISK_per.',Series{2,i},'=w;'))
end


%% Systemic risk ranking 


T = length(dates);
ranking=[];

% Daily rankings by SRISK contribution; 1=Most systematically important
% Institutions with no contribution to systemic risk (i.e. negative or Nan) on a given day assigned a rank of 0 

for i=1:T
    ind = zeros(1,19);
    [s, index]=sort(D(i,:),'descend');
    ind(index)=1:19;
    ind(D(i,:)==0)=0;
    ind(isnan(D(i,:)))=0;
    ranking = [ranking; ind];
end

% Aligning dates
for i=1:n
    eval(strcat('g = x.',Series{2,i},';'));
    [E, ia, ib] = intersect(dates,g);
    
    w = ranking(:,i);
    w=w(ia);
    eval(strcat('SRISK_rank.',Series{2,i},'=w;'))
end

%% Writing output to csv file

for i=1:n
    time_ = eval(strcat('cellstr(datestr(x.',Series{2,i},'))'));
    MESts = eval(strcat('MES.',Series{2,i}));
    MESts = num2cell(MESts);
    SRISKts = eval(strcat('SRISK.',Series{2,i}));
    SRISKts = num2cell(SRISKts);
    SRISKperts = eval(strcat('SRISK_per.',Series{2,i}));
    SRISKperts = num2cell(SRISKperts);
    SRISKrnkts = eval(strcat('SRISK_rank.',Series{2,i}));
    SRISKrnkts = num2cell(SRISKrnkts);
    
    TEMP = [time_ MESts SRISKts SRISKperts SRISKrnkts];

    
    
    fileID = fopen(eval(ID),'w');
    fprintf(fileID,'%s, %s, %s, %s, %s\n','Date','MES','SRISK', 'SRISK contribution', 'SRISK ranking');
    
    for g=1:length(SRISKts)
        fprintf(fileID,'%s, %f, %f, %f, %f\n' ,TEMP{g,:});
    end
   fclose('all');
end

%% Graphs
%   To switch between firms, change the index j for the respective series

j=1;
%   MES
figure
subplot(3,1,1)
    eval(strcat('plot(x.',Series{2,j},',MES.',Series{2,j},')'));
            title('MES','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('%','fontsize',8,'fontweight','b','color','k');
    eval(strcat('set(gca,''ylim'',[0 (max(MES.',Series{2,j},')+0.1*max(MES.',Series{2,j},'))]',...
            ',''ytick'',linspace(0,(max(MES.',Series{2,j},')+0.1*max(MES.',Series{2,j},')),5),''YTickLabel''',...
            ',100.*linspace(0,(max(MES.',Series{2,j},')+0.1*max(MES.',Series{2,j},')),5))'));

%   SRISK
subplot(3,1,2)
sc=1000000000;
    eval(strcat('plot(x.',Series{2,j},',SRISK.',Series{2,j},')'));
            title('SRISK','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('$Billions','fontsize',8,'fontweight','b','color','k');
    eval(strcat('set(gca,''ylim'',[0 (max(SRISK.',Series{2,j},')+0.1*max(SRISK.',Series{2,j},'))]',...
            ',''ytick'',linspace(0,(max(SRISK.',Series{2,j},')+0.1*max(SRISK.',Series{2,j},')),5),''YTickLabel''',...
            ',linspace(0,(max(SRISK.',Series{2,j},')+0.1*max(SRISK.',Series{2,j},'))./sc,5))'));
        

%   SRISK contribution
subplot(3,1,3)
    eval(strcat('plot(x.',Series{2,j},',SRISK_per.',Series{2,j},')'));
            title('SRISK contribution','FontSize',10,'fontweight','b','color','k');
            xlabel('','fontsize',8,'fontweight','b','color','k');
            ylabel('%','fontsize',8,'fontweight','b','color','k');
            set(gca,'ylim',[0 1],'ytick', linspace(0,1,5),'YTickLabel',100*linspace(0,1,5))
            
% identifying largest institutions
% NB! graph legend has to be changed manually
cap = data.MarkCap(end,:);
capm = cap;
capm(isnan(cap))=[];
capm = sort(capm,'descend');
cap5 = capm(1:5);
[E, ia, ib] = intersect(cap,cap5);
ia = sort(ia);

figure
hold on
for j=ia
eval(strcat('plot(x.',Series{2,j},',tsmovavg(100*SRISK_per.',Series{2,j},',''s'',30,1))'))
end
title('SRISK Contribution of 5 Biggest Institutions by Market Cap - 30 Day MA','FontSize',10,'fontweight','b','color','k')
xlabel('','fontsize',8,'fontweight','b','color','k');
ylabel('%','fontsize',8,'fontweight','b','color','k');
legend('Industrial & Commercial Bank of China','China Construction Bank corp',...
    'Agricultural Bank of China ltd','Bank of Communications ltd','China Merchants Bank')

hold off
saveas(gcf,'OUT\SRISKper.png')



% Look through data to identify largest institutiions within regions
figure
hold on
for j=[1 9 11 14 16]
eval(strcat('plot(x.',Series{2,j},',tsmovavg(100*SRISK_per.',Series{2,j},',''s'',30,1))'))
end
title('SRISK Contribution of 5 Biggest Institutions by Region - 30 Day MA','FontSize',10,'fontweight','b','color','k')
xlabel('','fontsize',8,'fontweight','b','color','k');
ylabel('%','fontsize',8,'fontweight','b','color','k');
legend('Industrial and commercial bank of china','OAO Sberbank of Russia',...
    'Itau Unibanco holdings SA','Barclays Africa','State bank of India')

hold off
saveas(gcf,'OUT\SRISKreg.png')
