% Yiwen Mei
% University of Connecticut
% 5/28/2013

% Description:
% The code calculate the baseflow time series based on UKIH method (Gutard
% 1992). This method is for daily flow series.

% Inputs:
% flow: flow time series;
%  p1 : Block periods
%  p2 : minima ratio

% Outputs:
% Qb: Baseflow time series

function [Qb,BFI]=UKIH(Q,A)

idx=mod(.827*24*A^.2,2)<1;
Lb=floor(.827*24*A^.2);
Lb(idx)=Lb(idx)+1;

Qmi=nan(length(Q),1);
Tmi=nan(length(Q),1);
for h=1:Lb:length(Q)
  if h+Lb-1<=length(Q)
    [q,i]=min(Q(h:h+Lb-1));
    Qmi(h)=q;
    Tmi(h)=h+i-1;
  end
end
Tmi(isnan(Qmi))=[];
Qmi(isnan(Qmi))=[];

% dQf=[diff(Qmi)./Qmi(2:length(Qmi));NaN]; % (Q_t-Q_t+1)/Q_t+1
dQf=[Qmi(1:length(Qmi)-1)./Qmi(2:length(Qmi));NaN]; % Q_t/Q_t+1
% Qmi=flipud(Qmi);
% dQmi=flipud([diff(Qmi)./Qmi(2:length(Qmi));NaN]); % (Q_t-Q_t-1)/Q_t-1
dQb=[NaN;Qmi(2:length(Qmi))./Qmi(1:length(Qmi)-1)]; % Q_t/Q_t-1
idx=~isnan(dQf) & ~isnan(dQb);
dQ=[dQf(idx) dQb(idx)];
[idx,Cg]=kmeans(dQ,2);

% plot(dQf(idx==1),dQb(idx==1),'.');hold on;plot(dQf(idx==2),dQb(idx==2),'.');

i=find(Cg(:,1)==min(Cg(:,1)) & Cg(:,2)==min(Cg(:,2)));
% i=find(Cg(:,1)==max(Cg(:,1)) & Cg(:,2)==max(Cg(:,2)));
idx=[NaN;idx;NaN];

% Qb=flipud(Qb);
Qmi=Qmi(idx==i);
Tmi=Tmi(idx==i);

Ti=1:length(Q);
Qb=interp1(Tmi,Qmi,Ti)';
Qb(Qb>Q)=Q(Qb>Q);

BFI=nansum(Qb)/nansum(Q);
end

% bm=floor(length(flow)/p1);
% 
% TP=ones(bm,1);
% for b=1:bm
%   [~,d]=min(flow(p1*b-(p1-1):p1*b));
%   TP(b)=p1*(b-1)+d;
% end
% 
% TP1=TP;
% for t=2:length(TP)-1
%   if flow(TP(t+1))/flow(TP(t))<=p2 && flow(TP(t-1))/flow(TP(t))<=p2
%     TP1(t)=NaN;
%   end
% end
% TP(isnan(TP1))=[];
% TP=(union([1 length(flow)],TP))';
% 
% Qbs=zeros(length(flow),1);
% for b=1:length(TP)-1
%   for t=TP(b):TP(b+1)
%     Qbs(t)=((flow(TP(b+1))-flow(TP(b)))/(TP(b+1)-TP(b)))*t+(TP(b+1)*flow(TP(b))-TP(b)*flow(TP(b+1)))/(TP(b+1)-TP(b));
%     
%     if Qbs(t)>flow(t)
%       Qbs(t)=flow(t);
%     end
%   end
% end
% 
% Qb=Qbs;
