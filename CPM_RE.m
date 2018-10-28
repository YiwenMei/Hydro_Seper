% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut
% Last updated on 1/22/2012

%% Functionality:
% The code locate the beginning and end of all the rainfall event.

%% Input
%  P : Rainfall time series (mm/h);
% LMG: Maximum pause of rainfall (typically set to 0);
% r0 : Low-rain rate threshold (typically set to 0).

%% Output:
% RE: Beginnings, ends, duration, cumulative volume and centroids of raifnall
%     events.

function RE=CPM_RE(rain,LMG,r0)
if isempty(LMG)
  LMG=0;
end
if isempty(r0)
  r0=0;
end

rm=nanmean(rain);

P=rain;
fP=flipud(P);
for t=1:length(rain)-LMG
  if nansum(P(t:t+LMG))<=r0*(LMG+1)*rm
    P(t)=NaN;
  end

  if nansum(fP(t:t+LMG))<=r0*(LMG+1)*rm
    fP(t)=NaN;
  end
end
fP=flipud(fP);
P(isnan(fP))=NaN;

Pa=find(P>=0);
Pb=ones(length(Pa),1);
for x=2:length(Pa)
  if Pa(x)>Pa(x-1)+1
    Pb(x)=Pb(x-1)+1;
  else
    Pb(x)=Pb(x-1);
  end
end

Pc=zeros(max(Pb),5);
for x=1:max(Pb)
  Pc(x,1)=Pa(find(Pb==x,1,'first'));
  Pc(x,2)=Pa(find(Pb==x,1,'last'));
  Pc(x,3)=length(Pc(x,1):Pc(x,2));
  Pc(x,4)=nansum(rain(Pc(x,1):Pc(x,2)));
  Pc(x,5)=max(rain(Pc(x,1):Pc(x,2)));
end

tcr=zeros(size(Pc,1),1);
for x=1:size(Pc,1)
  t=(Pc(x,1):Pc(x,2))';
  tcr(x)=nansum(t.*rain(t))/nansum(rain(t));
end

RE=[Pc tcr];
end
