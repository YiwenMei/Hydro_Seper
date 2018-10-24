% Yiwen Mei
% University of Connecticut
% 10/04/2016

% Description:
% The code a) locate the turning points of the hydrograph; b) calculate the
% maximum baseflow index & recession coefficient and c) separate baseflow
% from streamflow by the RCK methods; 

% Inputs:
% flow: Streamflow time series in mm/h
%  p1 : Half-regression period (half length of the period to calculate the
%       recession constant)
%  p2 : Constant slope threshold (1e-7 min-2 is suggested in Blume et al.
%       2007)
%  p3 : Minimu acceptabe r2 (r2 statistics of the t vs. ln(Q) corellation).
%  p4 : Upper limit of baseflow rate
%  p5 : Low slope (-dQ/dt)

% Outputs:
%   TS   : Baseflow time series
% pt.RiPs: Location of points in rising limb
% pt.RePs: Location of points in recession limb
%   Sc   : Recession coefficient K
%  BFIm  : Maximum baseflow index

function [TS,pt,Sc,BFIm]=RCK(flow,p1,p2,p3,p4,p5)

k=NaN(length(flow),4);
for t=p1+1:length(flow)-p1
  x=(t-p1:t+p1)';
  X=[ones(length(x),1) x];
  y=flow(x);
  z=log(flow(x));
  if length(find(~isnan(y)))>4
    [b,~,~,~,~]=regress(y,X);
    k(t,1)=-b(2); % Slope for every time step

    [b1,~,~,~,r]=regress(z,X);
    k(t,2)=-b1(2); % Recession coefficient for every time step
    k(t,3)=r(1); % r2 stats for every time step
  end
end

k(1:size(k,1)-1,4)=diff(k(:,2)); % Change of slope in h-2
k(isnan(flow) | isnan(k(:,4)),:)=NaN;

c=find(abs(k(:,4))<=p2 & k(:,3)>=p3 & k(:,1)>p5 & ~isnan(flow)); % criteria p4 is not
[f,~,~,~,~]=regress(k(c,1),[ones(length(c),1) flow(c)]); % considered here
K=f(2);

Qck=flow;
Qck(abs(k(:,4))>p2 | k(:,3)<p3 | flow>p4 | k(:,1)<p5 | isnan(flow))=NaN; % Select the points
Qck([1:p1 length(flow)-p1+1:length(flow)])=NaN;

RePs=find(Qck>=0);

rp=find(flow<=p4 & k(:,3)<p3);
idx=0;
for t=2:length(rp)-1
  if flow(rp(t))<=flow(rp(t-1)) && flow(rp(t))<flow(rp(t+1))
    idx=idx+1;
    RiPs(idx,1)=rp(t);
  end
end

TPs=unique([1;RiPs;RePs;length(flow)]);

CK_sl=interp1(TPs,flow(TPs),(1:length(flow))'); % RCK Baseflow time series
CK_sl(CK_sl<0)=0;
CK_sl(CK_sl>flow)=flow(CK_sl>flow);
CK_sl(isnan(flow))=NaN;

pt.RiPs=RiPs;
pt.RePs=RePs;
Sc=K;
TS=CK_sl;
BFIm=nansum(CK_sl)/nansum(flow);
