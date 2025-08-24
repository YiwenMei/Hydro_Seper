% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut, created 2016
% Last updated on 8/24/2025

%% Functionality
% This code performs recession analysis by fitting a recession model on streamflow
% time series to
%  1) locate turning points of the streamflow series;
%  2) construct baseflow time series by connecting all turning points with straight
%     lines; and
%  3) calculate the recession coefficient and maximum baseflow index of the basin;

%% Inputs
%  Q : streamflow time series for a basin (mm/h);
% sc : number of time step in a day (e.g., if Q is hourly, sc is 24; if Q is
%       daily, sc is 1);
%  A : basin size (km^2);
% ofn: full name of the .mat file to store the flow time series characteristics
%       (i.e. change rate of flow dQ/dt, recession coefficient k, logarithmic
%       change rate of recesssion coefficient dk);

% pflg: parallel flag (false/true - squential/parallel, default is false);
% Rnc : threshold to define no changes in recession coefficient (defaults are
%        1e-4 for hourly Q, 1e-7 for minute Q, 1e-2 for daily Q);
% r2M : coefficient of determination threshold to define acceptable model fitting
%        peformance;
% BFIi: initial guess of baseflow index;
%  n  : parameter determines the order of the recession model fitted to the streamflow
%        time series (i.e., dQ/dt=-kQ^n);

%% Outputs
%  pt : potential time steps used to define start (RiP) and end (ReP) points
%        of flow events;
%  Qb : baseflow time series formed by connecting points in pt with straight
%        lines;
%  K  : recession coefficient of the basin;
% BFIm: maximum baseflow index of the basin.

function [pt,Qb,BFIm,K,r,LSP]=RCK(Q,sc,A,ofn,varargin)
%% Check the inputs
narginchk(4,9);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Q',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Q'));
addRequired(ips,'sc',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'sc'));
addRequired(ips,'A',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'A'));
addRequired(ips,'ofn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ofn'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg')); 
if sc==1
  addOptional(ips,'Rnc',1e-2,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Rnc'));
elseif sc>1 && sc<=24
  addOptional(ips,'Rnc',1e-4,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Rnc'));
elseif sc>24 && sc<=24*60
  addOptional(ips,'Rnc',1e-7,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Rnc'));
end
addOptional(ips,'r2M',.8,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'r2M'));
addOptional(ips,'BFIi',NaN,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'BFIi'));
addOptional(ips,'n',1,@(x) validateattributes(x,{'double'},{'scalar','integer'},mfilename,'n'));

parse(ips,Q,sc,A,ofn,varargin{:});
pflg=ips.Results.pflg;
Rnc=ips.Results.Rnc;
r2M=ips.Results.r2M;
BFIi=ips.Results.BFIi;
n=ips.Results.n;
clear ips varargin

%% Streamflow characteristics
% Determine the length of regression window
LRW_i=5; % Minimum size of regression window
LSP=.827*sc*A^.2; % Length of recession; .827 convert (km^2)^.2 to (mi^2)^.2
LRW=floor(LSP/5+1); % Maintain 80% of point
LRW=max(LRW_i,LRW);

% Fitting recession model to find dQ/dt, k, r2, and dk/dt 
Qs=[nan(floor((sc-1)/2),1);movmean(Q,[floor((sc-1)/2) ceil((sc-1)/2)],'Endpoints','discard');...
    nan(ceil((sc-1)/2),1)]; % Moving average of a daily window
Q(Q==0)=NaN;
Qs(Qs==0 | isnan(Q))=NaN;

if exist(ofn,'file')~=2
  dQdt=nan(size(Qs));
  k=nan(size(Qs));
  r2=nan(size(Qs));
  dk=nan(size(Qs));

  T=(LRW-1)/2+1:length(Qs)-(LRW-1)/2;
  switch pflg
    case true
      parfor t=1:length(T)
        [b1,b2,r]=RCK_sub(T,t,LRW,Qs,n);
        dQdt(t)=b1; % Change rate of flow dQ/dt (L/T^2)
        k(t)=b2; % Recession coefficient k (T^(n-2)/L^(n-1))
        r2(t)=r; % Coefficient of determination
      end

      parfor t=1:length(T)
        b=RCK_sub1(T,t,LRW,k);
        dk(t)=b; % change rate of k (T^(n-3)/L^(n-1))
      end

    case false
      for t=1:length(T)
        [b1,b2,r]=RCK_sub(T,t,LRW,Qs,n);
        dQdt(t)=b1; % Change rate of flow dQ/dt (L/T^2)
        k(t)=b2; % Recession coefficient k (T^(n-2)/L^(n-1))
        r2(t)=r; % Coefficient of determination
      end

      for t=1:length(T)
        b=RCK_sub1(T,t,LRW,k);
        dk(t)=b; % change rate of k (T^(n-3)/L^(n-1))
      end
  end

% Determine the response ratio
  benb=[nan(LRW,1);movmean(Q,[LRW 0],'omitnan','Endpoints','discard')];
  bena=[movmean(Q,[0 LRW],'omitnan','Endpoints','discard');nan(LRW,1)];
  Rr=bena./benb;

  save(ofn,'dQdt','k','r2','dk','Rr');
else
  load(ofn,'dQdt','r2','dk','Rr');
end

%% Find RiPs and RePs
if LSP<LRW_i
  sw=2*LRW_i;
else
  sw=round(2*LSP);
end
ben=[nan(sw,1);movmean(Q,[sw sw],'omitnan','Endpoints','discard');nan(sw,1)]; % Envelope of baseflow
if isnan(BFIi)
  Qbi=min([ben Q],[],2);
  k=isnan(Q) | isnan(Qbi);
  Qbi(k)=NaN;
  Q(k)=NaN;
  BFIi=sum(Qbi,'omitnan')/sum(Q,'omitnan');
end

T=(1:length(Q))';
RiP=T(dQdt>1e-10 & Q<=ben & Rr>1/BFIi);
ReP=T(dk<log(Rnc) & r2>r2M & Q<=ben);
% plot(Q);hold on;plot(pt.ReP,Q(pt.ReP),'*');hold on;plot(pt.RiP,Q(pt.RiP),'*');
% plot(Q);hold on;plot(Qbi);hold on;plot(ben)
pt.RiP=RiP;
pt.ReP=ReP;
clear RiP ReP bena benb Rr k

%% Baseflow time series
Qb=interp1(union(pt.RiP,pt.ReP),Q(union(pt.RiP,pt.ReP)),T);
Qb(Qb<0)=0;
Qb(Qb>Q)=Q(Qb>Q);
Qb(isnan(Q))=NaN;
BFIm=sum(Qb,'omitnan')/sum(Q,'omitnan'); % Baseflow index

% Recession coefficient
try
  id1=dk<=log(Rnc) & ~isnan(r2) & dQdt<-1e-8 & r2>r2M;
  [b,~,~,~,r]=regress(dQdt(id1),[ones(length(Qs(id1)),1) Qs(id1)]);
catch
  id1=~isnan(dk) & ~isnan(r2) & dQdt<-1e-8 & r2>r2M;
  [b,~,~,~,r]=regress(dQdt(id1),[ones(length(Qs(id1)),1) Qs(id1)]);
end
K=-b(2);

% tbl=array2table([log(Qs(id1)) log(-dQdt(id1)) dk(id1)]);
% scatter(tbl,'Var1','Var2','filled','ColorVariable','Var3');
% colorbar;
% colormap('spring');
% fprintf('%i\n',sum(id1));
end

function [b1,b2,r]=RCK_sub(T,t,LRW,Qs,n)
ti=(T(t)-(LRW-1)/2:T(t)+(LRW-1)/2)';

y=Qs(ti);
if length(find(~isnan(y)))>.5*LRW
  X=[ones(length(ti),1) ti];
  [b,~,~,~,~]=regress(y,X);
  b1=b(2); % Change rate of flow dQ/dt (L/T^2)
else
  b1=NaN;
end

if n==1
  y=log(Qs(ti)); % dQ/Q
else
  y=Qs.^(1-n)/(1-n); % dQ/Q^n
end
if length(find(~isnan(y)))>.5*LRW
  X=[ones(length(ti),1) ti];
  [b,~,~,~,r]=regress(y,X);
  b2=-b(2); % Recession coefficient k (T^(n-2)/L^(n-1))
  r=r(1); % Coefficient of determination
else
  b2=NaN;
  r=NaN;
end
end

function b=RCK_sub1(T,t,LRW,k)
ti=(T(t)-(LRW-1)/2:T(t)+(LRW-1)/2)';

y=k(ti);
if length(find(~isnan(y)))>.5*LRW
  X=[ones(length(ti),1) ti];
  [b,~,~,~,~]=regress(y,X);
  b=log(abs(b(2))); % change rate of k (T^(n-3)/L^(n-1))
else
  b=NaN;
end
end
   