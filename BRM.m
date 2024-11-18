% Yiwen Mei (yiwen.mei@uconn.edu)
% CIRCA, University of Connecticut
% Last updated on 11/15/2022

%% Functionality:
% This code implements the Bump and Rise Method (BRM, Stewart, 2015). It optimizes
%  the initial baseflow index through iterations.

%% Inputs
%  Q  : streamflow time series for a basin (m^3/s);
% f_bp: a constant fraction of the increase or decrease of streamflow (bump);
% k_rs: the slope of the dividing line (rise, in m^3/s/time step);

% thr_E: BFI difference threshold under which to terminate the iteration;
% Nitr : maximum number of iteration;
% BFI0 : initial baseflow index.

%% Outputs:
%  Qb : baseflow time series of the same resolution and unit as Q; 
% BFI : long-term baseflow index of the basin based on RDF;
% N_bp: number of time step with negative baseflow and baseflow larger than streamflow;
% MRE : final MRE of BFI after the iteration process.

function [Qb,BFI,N_bp,MRE]=BRM(Q,f_bp,k_rs,varargin)
%% Check the inputs
narginchk(3,6);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Q',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Q'));
addRequired(ips,'f_bp',@(x) validateattributes(x,{'double'},{'scalar','<',1,'>',0},mfilename,'f_bp'));
addRequired(ips,'k_rs',@(x) validateattributes(x,{'double'},{'scalar','positive'},mfilename,'k_rs'));

addOptional(ips,'thr_E',1e-6,@(x) validateattributes(x,{'double'},{'scalar','nonnegative'},...
    mfilename,'thr_E'));
addOptional(ips,'Nitr',5,@(x) validateattributes(x,{'double'},{'scalar','integer','positive'},...
    mfilename,'Nitr'));
addOptional(ips,'BFI0',.5,@(x) validateattributes(x,{'double'},{'scalar','positive'},mfilename,'BFI0'));

parse(ips,Q,f_bp,k_rs,varargin{:});
thr_E=ips.Results.thr_E;
Nitr=ips.Results.Nitr;
BFI0=ips.Results.BFI0;
clear ips varargin

%% Determine the segments
TSi=~isnan(Q);
k=Run_Length(TSi,true,Q);
k=reshape(k',2,length(k)/2)';
k(k(:,1)==1,1)=1:sum(k(:,1));
TSi=Run_Length(reshape(k(:,1:2)',size(k,1)*2,1),false,[])';

ct=1;
MRE=thr_E;
while abs(MRE)>=thr_E && ct<Nitr
%% Baseflow for every segment
  Qb=nan(size(Q));
  for i=1:max(TSi)
    q=Q(TSi==i);
    if length(q)>1

      qb=nan(size(q));
      qb(1)=BFI0*q(1);
      for j=2:length(q)
        if q(j)>qb(j-1)+k_rs*(j-(j-1))
          qb(j)=qb(j-1)+k_rs*(j-(j-1))+f_bp*(q(j)-q(j-1));
        else
          qb(j)=q(j);
        end
      end
      Qb(TSi==i)=qb;
    end
  end

  N_bp=sum(Qb<0 | Qb>Q);
  Qb(Qb<0)=0; % baseflow must be greater than 0 and lower than total flow
  Qb(Qb>Q)=Q(Qb>Q);

%% Calculate the baseflow index
  k=any(isnan([Qb Q]),2);
  BFI=sum(Qb(~k),'omitnan')/sum(Q(~k),'omitnan');
  MRE=(BFI0-BFI)/BFI0;
  BFI0=BFI;
  ct=ct+1;
end
end
