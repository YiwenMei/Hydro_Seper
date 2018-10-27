% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut
% Last updated on 10/20/2017

%% Functionality:
% This code performs recession analysis by fitting a recession model on streamflow
% time series to
%  1) locate turning points of the streamflow series;
%  2) construct baseflow time series by connecting all turning points with straight
%     lines; and
%  3) calculate the recession coefficient and maximum baseflow index of the basin;

%% Inputs
%  Q : streamflow time series for a basin (mm/h);
%  A : basin size (km^2);
%  n : parameter determines the order of the recession model fitted to the streamflow
%      time series (i.e., dQ/dt=-kQ^n);
% Rnc: threshold to define no changes in recession coefficient;
% r2M: coefficient of determination threshold to define acceptable model fitting
%      peformance;
% LRW: length of regression window (h, if it is set to [], the code estimate
%      LRW based on the basin size);
% Obi: baseflow time series (if it is set to [], the moving average algorithm
%      is evoked);
% TC : full name of the .mat file to store the flow time series characteristics
%      (i.e. change rate of flow dQ/dt, recession coefficient k, logarithmic
%      change rate of recesssion coefficient dk).

%% Outputs:
% pt.RiP: potential time steps used to define start of flow events; 
% pt.ReP: potential time steps used to define end of flow events;
%   Qb  : baseflow time series formed by connecting points in pt with straight
%         lines;
%   K   : recession coefficient of the basin;
%  BFIm : maximum baseflow index of the basin.

function [pt,Qb,K,BFIm]=RCK(Q,A,n,Rnc,r2M,LRW,Qbi,TC)
Q(Q==0)=NaN;
T=1:length(Q);
LSP=.827*24*A^.2; % Empirical estimate of length of recession

% Determine the length of regression window
if isempty(LRW)
  idx=mod(LSP/5+1,2)<1;
  LRW=floor(LSP/5+1); % Maintain 80% of point
  LRW(idx)=LRW(idx)+1;
end

%% Fitting recession model
if exist(TC,'file')~=2
  dQdt=nan(size(Q));
  k=nan(size(Q));
  r2=nan(size(Q));
  dk=nan(size(Q));

  for t=(LRW-1)/2+1:length(Q)-(LRW-1)/2
    ti=(t-(LRW-1)/2:t+(LRW-1)/2)';
    y=Q(ti);
    if length(find(~isnan(y)))>.5*LRW
      X=[ones(length(ti),1) ti];

      [b,~,~,~,~]=regress(y,X);
      dQdt(t)=b(2); % Change rate of flow dQ/dt (L/T^2)
    end

    if n==1
      y=log(Q(ti)); % dQ/Q
    else
      y=Q.^(1-n)/(1-n); % dQ/Q^n
    end
    if length(find(~isnan(y)))>.5*LRW
      X=[ones(length(ti),1) ti];

      [b,~,~,~,r]=regress(y,X);
      k(t)=-b(2); % Recession coefficient k (T^(n-2)/L^(n-1))
      r2(t)=r(1); % Coefficient of determination
    end

    y=k(ti);
    if length(find(~isnan(y)))>.5*LRW
      X=[ones(length(ti),1) ti];

      [b,~,~,~,~]=regress(y,X);
      dk(t)=log(abs(b(2))); % change rate of k (T^(n-3)/L^(n-1))
    end
  end
  save(TC,'dQdt','k','r2','dk');

else
  load(TC);
end

%% Determine the turning points
% Envelope of baseflow
ben=[nan(2*round(LSP),1);movmean(Q,[2*round(LSP) 2*round(LSP)],'omitnan','Endpoints',...
    'discard');nan(2*round(LSP),1)];

% RiPs and RePs
if isempty(Qbi)
%   ben=[nan(2*round(LRE),1);movmean(Q,[2*round(LRE) 2*round(LRE)],'omitnan','Endpoints',...
%       'discard');nan(2*round(LRE),1)];
  Qbi=ben;
  Qbi(Qbi>Q)=Q(Qbi>Q);
% else
%   ben=Qbi;
end
BFI=nansum(Qbi)/nansum(Q);

benb=[nan(LRW,1);movmean(Q,[LRW 0],'omitnan','Endpoints','discard')];
bena=[movmean(Q,[0 LRW],'omitnan','Endpoints','discard');nan(LRW,1)];
Rr=bena./benb;

RiP=T(dQdt>1e-10 & Q<=ben & Rr>1/BFI);
ReP=T(dk<Rnc & r2>r2M & Q<=ben);
% plot(Q);hold on;plot(pt.ReP,Q(pt.ReP),'*');hold on;plot(pt.RiP,Q(pt.RiP),'*');
% plot(Q);hold on;plot(Qb);hold on;plot(ben)
pt.RiP=RiP;
pt.ReP=ReP;
clear RiP ReP

%% Baseflow time series
Qb=interp1(union(pt.RiP,pt.ReP),Q(union(pt.RiP,pt.ReP)),T');
Qb(Qb<0)=0;
Qb(Qb>Q)=Q(Qb>Q);
Qb(isnan(Q))=NaN;
BFIm=nansum(Qb)/nansum(Q); % Baseflow index

% Recession coefficient
id1=~isnan(dk) & ~isnan(r2) & dQdt<-1e-8;
[b,~,~,~,~]=regress(dQdt(id1),[ones(length(Q(id1)),1) Q(id1)]);
% loglog(Q(id1),dQdt(id1),'.')
K=-b(2);
end
