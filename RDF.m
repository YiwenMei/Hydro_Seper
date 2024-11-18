% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut
% Last updated on 10/6/2022

%% Functionality:
% This code implements the recursive digital filter (RDF, Eckhardt, 2005) and
%  the filtered UKIH method (FUKIH, Aksoy et al. 2009). In the second case, the
%  UKIH-based baseflow time series must be provided as an optional input.

%% Inputs
%  Q  : streamflow time series for a basin (m3/time step or mm/time step);
% BFIm: maximum baseflow index;
%  a  : recession constant (note that a=exp(-Kt))

% Qb0: initial baseflow time series of the same size as Q.

%% Outputs:
%  Qb : baseflow time series of the same resolution and unit as Q; 
% BFI : long-term baseflow index of the basin based on RDF;
% N_bp: number of time step with negative baseflow and baseflow larger than streamflow.

function [Qb,BFI,N_bp]=RDF(Q,BFIm,a,varargin)
%% Check the inputs
narginchk(3,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Q',@(x) validateattributes(x,{'double'},{'vector','nonnegative'},mfilename,'Q'));
addRequired(ips,'BFIm',@(x) validateattributes(x,{'double'},{'scalar','<',1,'>',0},mfilename,'BFIm'));
addRequired(ips,'a',@(x) validateattributes(x,{'double'},{'scalar','<',1,'>',0},mfilename,'a'));

addOptional(ips,'Qb0',[],@(x) validateattributes(x,{'double'},{'vector','nonnegative'},mfilename,'Qb0'));

parse(ips,Q,BFIm,a,varargin{:});
Qb0=ips.Results.Qb0;
clear ips varargin

%% Determine the segments
TSi=~isnan(Q);
k=Run_Length(TSi,true,Q);
k=reshape(k',2,length(k)/2)';
k(k(:,1)==1,1)=1:sum(k(:,1));
TSi=Run_Length(reshape(k(:,1:2)',size(k,1)*2,1),false,[])';

%% Baseflow for every segment
Qb=nan(size(Q));
for i=1:max(TSi)
  q=Q(TSi==i);
  if length(q)>1

% Without initial baseflow (RDF)
    qb=nan(size(q));
    if isempty(Qb0)
      qb(1)=BFIm*q(1);
      for j=2:length(q)
        qb(j)=((1-BFIm)*a*qb(j-1)+(1-a)*BFIm*q(j))/(1-a*BFIm);
      end

% With initial baseflow (FUKIH)
    else
      qb0=Qb0(TSi==i);
      for j=2:length(q)
        qb(j)=((1-BFIm)*a*qb0(j-1)+(1-a)*BFIm*q(j))/(1-a*BFIm);
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
end
