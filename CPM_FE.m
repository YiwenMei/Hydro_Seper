% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut
% Last updated on 11/17/2012

%% Functionality:
% The code identify potential hydrologic events.

%% Inputs
% peak: Location of peak flow;
%  pts: points in rising and recession limb;
%   Q : flow time series (mm/h).

%% Output
% FE: Locations (the n-th hour of the begin, peak and end) of flow events

function FE=CPM_FE(peak,pts,Q)

PE=nan(length(peak),3); % length(peak)-1
for x=1:length(peak) % length(peak)-1
  rip=find(peak(x)-pts.RiPs>0,1,'last');
  rep=find(peak(x)-pts.RePs<0,1,'first');
  if ~isempty(rip) && ~isempty(rep)
    PE(x,1)=pts.RiPs(rip);
    PE(x,2)=peak(x);
    PE(x,3)=pts.RePs(rep);
  end
end
PE(isnan(PE(:,1)),:)=[];

J1=zeros(size(PE,1),1);
J=zeros(size(PE,1)-1,1);
for x=1:size(PE,1)
  X=setdiff(1:size(PE,1),x);

  for y=1:length(X)
    if ~isempty(max(intersect(PE(x,1):PE(x,3),PE(X(y),1):PE(X(y),3))))
      J(y)=max(intersect(PE(x,1):PE(x,3),PE(X(y),1):PE(X(y),3)));
    else
      J(y)=0;
    end
  end
  if ~isempty(min(J(J~=0)))
    J1(x)=min(J(J~=0));
  else
    J1(x)=0;
  end
end
J=J1(J1~=0);

J=unique(J);
MPEs=zeros(length(J),1);
for y=1:length(J)
  MPEs(y,1)=min(PE(J1==J(y),1));
  MPEs(y,3)=max(PE(J1==J(y),3));

  pid=PE(J1==J(y),2);
  [~,i]=max(Q(pid));
  MPEs(y,2)=pid(i);
end

FE=[PE(J1==0,:);MPEs];
end
