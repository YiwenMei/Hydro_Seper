path='C:\Users\Mei Yiwen\Desktop\Study\Y.Mei.2013\Usermanual\CPM\Demo_data\';
load([path 'Demo_data'])

%% Baseflow separation
% Convert the units to mm/h
Q=3.6*Q/A;
Q(Q==0)=NaN;

% Clustering of flow
q=Q(~isnan(Q));
idx=kmeans(q, 3);
mq=nan(3, 1);
for c=1:3
    mq(c)=max(q(idx==c));
end
ben=min(mq,[],1)';

% RCK Method
[Qbf_sl, pts, K, BFIm]=RCK(Q, 11, 5e-5, .9, ben, 1e-10);

% FRCK Method
BFIM=[BFIm;0;0;0;0;0;0]; % Check whether it converge
for y=1:6
    [Qbf, BFIm]=FRCK(Q, Qbf_sl, K);

    BFIM(y+1)=BFIm;
end

%% Event identification
% Locate the peak flow rate from the long-term hydrograph
LSP=.827*24*A^.2;
Pks=CPM_peak(Q, LSP, pts.RePs, Qbf);

% Formation of flow events
FE=CPM_FE(Pks, pts.RiPs, pts.RePs);

% Formation of rainfall episodes
RE=CPM_RE(P, 0, 0);

RE(RE(:, 3)==1 | RE(:, 4)./RE(:, 3)<=.01, :)=[]; % Remove the minor episodes

% Event association
mlag=[LSP;0;0;0;0]; % Initialize mlag by LSP
for y=1:4
    RFE=CPM(Q, P, Qbf, FE, RE, LSP, mlag(y));

    k=find(isnan(RFE.tE(:, 1))); % Remove the non-associable
    if ~isempty(k)               % flood events
        RFE.D(k, :)=[];
        RFE.V(k, :)=[];
        RFE.tE(k, :)=[];
        RFE.R(k, :)=[];
        RFE.GEs(k, :)=[];
    end

    lag=RFE.D(:, 3);
    lag(RFE.D(:, 3)<0 |  RFE.R(:, 1)>1 | ...
        RFE.D(:, 3)==max(RFE.D(:, 3)))=[];
    mlag(y+1)=mean(lag); % update mlag for every iteration
end

%% Extra Filters
k=find(RFE.D(:, 3)<0 | ...
    RFE.R(:, 1)>1 | ...
    RFE.V(:, 2)./RFE.D(:, 2)<nanmean(Q) | ...
    Q(RFE.tE(:, 7))./max(Q(RFE.tE(:, [2 4])), [], 2)<1/BFIm(length(BFIm)) | ...
    min(Qbf(RFE.tE(:, [2 4]))./Q(RFE.tE(:, [2 4])), [], 2)<BFIm(length(BFIm)));

RFE.D(k, :)=[];
RFE.V(k, :)=[];
RFE.tE(k, :)=[];
RFE.R(k, :)=[];

%% Concatenate the parameters
para=[A K BFIm ben LSP mlag(length(mlag))];

%% Save data
clear BFIm c FE idx k K lag LSP mlag mq Pks pts q Qbf_sl RE ben y
save([path 'Results'],'RFE','P','para','Q','Qbf');
