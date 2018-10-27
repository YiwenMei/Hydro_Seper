% Yiwen Mei
% University of Connecticut
% 11/17/2012

% Description:
% The code construct the baseflow time series by using the recursive
% digital filter (following Eckhardt, 2005)

% Inputs:
% flow: Streamflow time series
% BFIm: Maximum baseflow index
%  K  : Recession coefficient (note that a=exp(-Kt))

% Outputs:
% Qbr: Baseflow time series

function [Qb,BFI]=RDF(Q,BFIm,K)
a=exp(-K);

Q_flip=flipud(Q);
Qb_flip=Q_flip;
for t=2:length(Q)
  Qb_flip(t)=((1-BFIm)*a*Qb_flip(t-1)+(1-a)*BFIm*Q_flip(t))/(1-a*BFIm);
end

Qb=flipud(Qb_flip);
for t=2:length(Q)
  Qb(t)=((1-BFIm)*a*Qb(t-1)+(1-a)*BFIm*Q(t))/(1-a*BFIm);
end

for t=1:length(Q)
  if Qb(t)<0
    Qb(t)=0;
  elseif Qb(t)>Q(t)
    Qb(t)=Q(t);
  end
end

BFI=nansum(Qb)/nansum(Q);
end
