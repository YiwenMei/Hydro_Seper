% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut
% Last updated on 5/17/2016

%% Functionality:
% The code find the location of peaks.

%% Input
%  Q : streamflow time series (mm/h);
%  A : basin size (km^2);
% pts: Location of turning points
% Qb: baseflow time series from FRCK (mm/h).

%% Output
% Qma: Location of the peaks

function Qma=CPM_peak(Q,A,Qb)
LSP=.827*24*A^.2;
rp=nansum(Q)/nansum(Qb); % Use the inverse of BFIm as response ratio

Qc=Q; % Initial and final Characteristic flow points
for t=2:length(Q)-1
  if Q(t)==Q(t-1) && Q(t)==Q(t+1)
    Qc(t)=NaN;
  end
end
Qc1(:,1:2)=[find(~isnan(Qc)) Q(~isnan(Qc))];

skm=0; % Maximum number of skiped time step
for x=1:size(Qc1,1)-1
  if Qc1(x,2)==Qc1(x+1,2)
    skm=skm+1;
  end
end

sk=0;
for x=1:size(Qc1,1)-skm-1
  if Qc1(x+sk,2)==Qc1(x+sk+1,2)
    Qc2(x,1)=floor((Qc1(x+sk,1)+Qc1(x+1+sk,1))/2); % Characteristic flow points
    Qc2(x,2)=Qc1(x+sk,2);
    sk=sk+1;
  else
    Qc2(x,:)=Qc1(x+sk,:);
  end
end

Qma=Qc2(:,1); % Maxima flow Characteristic points
for x=2:size(Qc2,1)-1
  if Qc2(x,2)>Qc2(x-1,2) && Qc2(x,2)>Qc2(x+1,2)
    Qma(x)=Qc2(x,1);
  else
    Qma(x)=NaN;
  end
end
Qma(isnan(Qma))=[];
Qma(length(Qma))=[];

% Qma(Q(Qma)./Qb(Qma)<=rp | Q(Qma)<=nanmean(Q) | Qma>max(pts.RePs)...
%     | Qma<min(pts.RiPs),:)=NaN;
Qma(Q(Qma)./Qb(Qma)<=rp | Q(Qma)<=nanmean(Q),:)=NaN;
Qma(isnan(Qma))=[];

Qma(:,2)=1; % Exclude the adjadcent Maxima
for x=2:length(Qma)
  if Qma(x,1)-Qma(x-1,1)>LSP
    Qma(x,2)=Qma(x-1,2)+1;
  else
    Qma(x,2)=Qma(x-1,2);
  end
end

Qma1=Qma;
for j=1:max(Qma(:,2))
  X=find(Qma(:,2)==j);
  if length(X)>2
    for x=2:length(X)-1
      if Q(Qma(X(x),1))<Q(Qma(X(x)-1,1)) && Q(Qma(X(x),1))<Q(Qma(X(x)+1,1))
        Qma1(X(x),:)=NaN;
      end
    end
  end
end
Qma=Qma1(:,1);
Qma(isnan(Qma))=[];

Qma(:,2)=1;
for x=2:length(Qma)
  if Qma(x,1)-Qma(x-1,1)>LSP
    Qma(x,2)=Qma(x-1,2)+1;
  else
    Qma(x,2)=Qma(x-1,2);
  end
end

Qma1=Qma;
sk=0;
for j=1:max(Qma(:,2))
  X=Qma(Qma(:,2)==j,1);
  x=find(Q(X)<max(Q(X)));
  Qma1(x+sk,:)=NaN;
  sk=sk+length(X);
  x=find(Qma1(:,2)==j);
  if length(x)>1
    x=x(x~=round(mean(x)));
    Qma1(x,:)=NaN;
  end
end
Qma=Qma1(:,1);
Qma(isnan(Qma))=[];
end
