% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut
% Last updated on 10/6/2022

%% Functionality:
% This code calculate the baseflow time series using the Smooth Minima (or UKIH)
%  method (Gutard 1992). This function is adaptable for sub-daily streamflow.

%% Inputs
%  Q : streamflow time series for a basin (m3/time step or mm/time step);
% sc : number of time step in a day (e.g. if Q is hourly, sc is 24; if Q is daily,
%       sc is 1);
% N_d: number of days for each block;
% rt : response ratio to determine if a point is a turning point.

%% Outputs:
%  Qb : baseflow time series of the same resolution and unit as Q; 
% BFI : long-term baseflow index of the basin based on SMM;
% N_bp: number of time step with negative baseflow and baseflow larger than streamflow.

function [Qb,BFI,N_bp]=SMM(Q,sc,N_d,rt)
%% Check the inputs
narginchk(4,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Q',@(x) validateattributes(x,{'double'},{'vector','nonnegative'},mfilename,'Q'));
addRequired(ips,'sc',@(x) validateattributes(x,{'double'},{'scalar','integer','positive'},...
    mfilename,'sc'));
addRequired(ips,'N_d',@(x) validateattributes(x,{'double'},{'scalar','integer','positive'},...
    mfilename,'N_d'));
addRequired(ips,'rt',@(x) validateattributes(x,{'double'},{'scalar','positive'},mfilename,'rt'));

parse(ips,Q,sc,N_d,rt);
clear ips varargin

%% Determine the segments
TSi=~isnan(Q);
k=Run_Length(TSi,true,Q);
k=reshape(k',2,length(k)/2)';
k(k(:,1)==1,1)=1:sum(k(:,1));
TSi=Run_Length(reshape(k(:,1:2)',size(k,1)*2,1),false,[])';

%% Calculate baseflow time series
l=round(N_d*sc); % Block size

Qb=nan(size(Q));
N_bp=nan(max(TSi),1);
for i=1:max(TSi)
  q=Q(TSi==i);
  if length(q)>=3*l
    [qb,n_bp]=cal_BF(q,l,rt);
    Qb(TSi==i)=qb;

    N_bp(i)=n_bp;
  end
end
N_bp=sum(N_bp,'omitnan');

%% Calculate the baseflow index
k=any(isnan([Qb Q]),2);
BFI=sum(Qb(~k),'omitnan')/sum(Q(~k),'omitnan');
end

function [qb,n_bp]=cal_BF(q,l,rt)
%% Find minima of each block
N_blk=length(q)/l;
blk=(1:ceil(N_blk))'; % At least three blocks
id=repelem(blk,l,1);
id=id(1:length(q));

q_sts=nan(max(id),2);
for i=1:max(id)
  q_blk=q(id==i);
  if length(q_blk)==l
    [qm,j]=min(q_blk,[],'omitnan');
    qm_id=find(id==i,1)-1+j;
    q_sts(i,:)=[qm_id,qm];
  end
end

%% Filter the minima
k=rt*q_sts(2:end-1,2)<=min([q_sts(1:end-2,2) q_sts(3:end,2)],[],2); % Smooth minima criteria
k=[false;k;false];
q_sts(~k,:)=[];

%% Interpolate the filtered minima with straight lines
if size(q_sts,1)>=2
  Ti=(1:length(q))';
  qb=interp1(q_sts(:,1),q_sts(:,2),Ti);
  n_bp=sum(qb<0 | qb>q);
  qb(qb<0)=0; % baseflow must be greater than 0 and lower than total flow
  qb(qb>q)=q(qb>q);
else
  qb=nan(size(q));
  qb(q_sts(:,1))=q_sts(:,2);
  n_bp=NaN;
end
end
