function Demo()
addpath('C:\Fun_pool\Hysep_tool');
addpath('C:\Fun_pool\Drought_tool');

ipth='C:\Fun_pool\Hysep_tool';
fn=fullfile(ipth,'Demo_data.mat');
load(fn,'A','P','Q');

% Convert the units to mm/h
Q=3.6*Q/A;
Q(Q==0)=NaN;

%% Recession analysis and baseflow separation
n=1;
Rnc=10e-6;
r2x=.9;

ofn=fullfile(ipth,'TC.mat');

BFI=nan(6,1);
for y=1:6
% Find the rise and recession points
  if y==1
    [pt,Qbf,BFIx,K,r,LSP]=RCK(Q,24,A,ofn,false,Rnc,r2x,NaN,n);
  else
    [pt,Qbf,BFIx,K,r,LSP]=RCK(Q,24,A,ofn,false,Rnc,r2x,BFIx,n);
  end

% Filter the baseflow time series 
  [Qbf,BFIx]=RDF(Q,BFIx,K,Qbf);
  BFI(y,:)=BFIx;
end

%% Event identification
% Locate the peak flow rate from the long-term hydrograph
Pks=CPM_peak(Q,LSP);

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
k=find(RFE.D(:,3)<0 | RFE.R(:,1)>1 | RFE.V(:,2)./RFE.D(:,2)<mean(Q,'omitnan'));
RFE.D(k,:)=[];
RFE.V(k,:)=[];
RFE.tE(k,:)=[];
RFE.R(k,:)=[];

%% Save data
paTab=table(A,K,BFIx,LSP,mlag(length(mlag)),'VariableNames',{'Area','ResC',...
    'BFIx','LSP','Tl_m'},'RowNames',{'XXX'});
paTab.Properties.VariableUnits={'km^2','h^(n-2)/mm^(n-1)','-','h','h'};

ofn=fullfile(pth,'Demo_O.mat');
save(ofn,'RFE','paTab','Qbf');
end
