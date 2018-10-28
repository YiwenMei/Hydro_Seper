addpath('XXX');
inpth='XXX';
load([inpth 'Demo_data'])

%% Recession analysis and baseflow separation
% Convert the units to mm/h
Q=3.6*Q/A;
Q(Q==0)=NaN;

n=1;
TCout='XXX.mat';

BFI=nan(6,1);
for y=1:6
% Find the rise and recession points  
  if y==1
    [pt, Qbf, K, ~]=RCK(Q, A, n, -11.5, .9, [], [], TCout);
  else
    [pt, Qbf, K, ~]=RCK(Q, A, n, -11.5, .9, [], Qbf, TCout);
  end

% Filter the baseflow time series 
  [Qbf, BFIm]=FRCK(Q, Qbf, K);

  BFI(y,:)=BFIm;
end

%% Event identification
% Locate the peak flow rate from the long-term hydrograph
Pks=CPM_peak(Q, A, Qbf);

% Formation of flow events
FE=CPM_FE(Pks, pt.RiP, pt.ReP, Q);

% Formation of rainfall episodes
RE=CPM_RE(P, [], []);
RE(RE(:,3)==1 | RE(:,4)./RE(:,3)<=.01, :)=[]; % Remove the minor episodes

% Event association
mlag=nan(4,1);
for y=1:4
  if y==1
    RFE=CPM(Q, P, Qbf, FE, RE, A, []);
  else
    RFE=CPM(Q, P, Qbf, FE, RE, A, mlag(y));
  end

  k=find(isnan(RFE.tE(:,1))); % Remove the non-associable events
  if ~isempty(k)
    RFE.D(k,:)=[];
    RFE.V(k,:)=[];
    RFE.tE(k,:)=[];
    RFE.R(k,:)=[];
  end

  lag=RFE.D(:,3);
  lag(RFE.D(:,3)<0 |  RFE.R(:,1)>1 | RFE.D(:,3)==max(RFE.D(:,3)))=[];
  mlag(y)=mean(lag); % update mlag for every iteration
end

%% Extra Filters
k=find(RFE.D(:,3)<0 | RFE.R(:,1)>1 | RFE.V(:,2)./RFE.D(:,2)<nanmean(Q));
RFE.D(k,:)=[];
RFE.V(k,:)=[];
RFE.tE(k,:)=[];
RFE.R(k,:)=[];

%% Save data
paTab=table(A,K,BFIm,LSP,mlag(length(mlag)),'VariableNames',{'Area','ReceCoef',...
    'meanBFI','SearchLen','meanTlag'},'RowNames',{'XXX'});
save('XXX','RFE','paTab','Qbf');
