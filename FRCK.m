% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut
% Last updated on 11/16/2017

%% Functionality:
% This code computes the filtered-baseflow hydrograph.

%% Input
% Q : flow time series (mm/h);
% BF: baseflow time series (mm/h) or baseflow index;
% K : Recession coefficient.

%% Output
%  Qb : Baseflow time series;
% BFIm: maximum baseflow index of the time peirod.

%% Additional note
% If BF is a time series, the code outputs the filtered baseflow time
%   series of BF;
% If BF is a baseflow index, the code outputs the filteres baseflow time
%   series of Q (in this case the code works as the recusive digital filter
%   documented by Eckhardt 2005);

function [Qb,BFIm]=FRCK(Q,BF,K)
a=exp(-K); % Recession constant

% Construct the baseflow time series
if numel(BF)>1 % if BF is a baseflow time series
  BFIm=nansum(BF)/nansum(Q);
  Qb=BF;

  for t=2:length(Q)
    if ~isnan(BF(t-1)) && ~isnan(Q(t))
      Qb(t)=((1-BFIm)*a*BF(t-1)+(1-a)*BFIm*Q(t))/(1-a*BFIm);
    end
  end

else % if BF is baseflow index
  Qb=Q;
  for t=2:length(Q)
    if ~isnan(Qb(t-1)) && ~isnan(Q(t))
      Qb(t)=((1-BF)*a*Qb(t-1)+(1-a)*BF*Q(t))/(1-a*BF);
    end
  end
end

Qb(Qb>Q)=Q(Qb>Q);
Qb(Qb<0)=0;
Qb(isnan(Q))=NaN;

BFIm=nansum(Qb)/nansum(Q);
end
