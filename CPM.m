% Yiwen Mei (yiwen.mei@uconn.edu)
% CEE, University of Connecticut
% Last updated on 7/25/2013

%% Functionality:
% The code associates the rainfall events with runoff events and calculate
% the correspoinding event properties.

%% Input
%  Q  : streamflow time series (mm/h);
%  P  : rainfall time series (mm/h);
%  Qb : baseflow time series (mm/h);
%  FE : timing of flow events (begin, end, and peak)
%  RE : timing of rainfall episodes (begin, end, and centroid)
%  A  : basin size (km^2);
% mlag: Mean event time lag (h).

%% Output
% Mat_RE.D : Durations of rainfall event, flow event and event time lag
% Mat_RE.V : Volume of event rainfall, flow and baseflow
% Mat_RE.tE: Timing of the matched rainfall and flow events (begins of rainfall &
%            flow, ends of rainfall & flow, centroids of rainfall & flow and max
%            flow peak)
% Mat_RE.R : Event-based runoff coefficient and baseflow index

function Mat_RE=CPM(Q,P,Qb,FE,RE,A,mlag)
LSP=.827*24*A^.2;

if ~isempty(FE(:,1))
  sh=sort(FE(:,1));
  eh=sort(FE(:,3));
  ph=sort(FE(:,2));
end

sr=RE(:,1); % starting points of hyetographs
er=RE(:,2); % ending points of hyetographs
cr=RE(:,6); % centroids of hyetographs

eh=[1;eh];
ch=nan(size(sh,1),1);
ccr=nan(size(sh,1),1);

pbE=cell(size(sh,1),1);
pcE=cell(size(sh,1),1);
precE=cell(size(sh,1),1);

GAEs=NaN*zeros(size(sh,1),2);
GBEs=NaN*zeros(size(sh,1),2);
for n=1:size(sh,1)
  th=sh(n):eh(n+1);
  ch(n)=nansum(Q(th)'.*th)/nansum(Q(th)); % Centroid of flow event

% Start of hydrograph after rainfall within a range
% Range is the mean time lag for wet period or the gap between two flow
% event, LSP is used for quality control
  pbE{n}=find(sh(n)-sr>=0 & sh(n)-sr<min([LSP,sh(n)-eh(n),mlag]));

% Centroid of hyetograph before hydrograph but after the last pbE{n} 
  if ~isempty(pbE{n})
    pcE{n}=find(cr>cr(pbE{n}(length(pbE{n}))) & ch(n)>=cr);
  else % Centroid of hyetograph before hydrograph but after start of hydrograph
    pcE{n}=find(cr>sh(n) & ch(n)>=cr);
  end
  precE{n}=union(pcE{n},pbE{n});

  if ~isempty(precE{n})
    tcr=sr(precE{n}(1)):er(precE{n}(length(precE{n})));
    ccr(n)=nansum(P(tcr)'.*tcr)/nansum(P(tcr));
  end

  if ~isempty(pbE{n})
    GAEs(n,:)=[sr(min(pbE{n})) er(max(pbE{n}))];
  end
  if ~isempty(pcE{n})
    GBEs(n,:)=[sr(min(pcE{n})) er(max(pcE{n}))];
  end
end

lag=ch-ccr;

postE=cell(size(sh,1),1);
aRE=cell(size(sh,1),1);

Dh=nan(size(sh,1),1);
Dcr=nan(size(sh,1),1);

Vcr=nan(size(sh,1),1);
Vt=nan(size(sh,1),1);
Vb=nan(size(sh,1),1);

tE=nan(size(sh,1),7);

GCEs=nan(size(sh,1),2);
for n=1:size(sh,1)
  postE{n}=find(ch(n)<cr & eh(n+1)>er & eh(n+1)-cr>lag(n));

  if ~isempty(postE{n})
    GCEs(n,:)=[sr(min(postE{n})) er(max(postE{n}))];
  end

  aRE{n}=union(precE{n},postE{n});

  if ~isempty(aRE{n})
    tcr=sr(aRE{n}(1)):er(aRE{n}(length(aRE{n}))); % Period of rainfall event
    th=sh(n):eh(n+1); % Period of flow event

    Dh(n)=length(~isnan(th)); % Duration of flow event
    Dcr(n)=length(tcr); % Duration of rainfall event

    Vcr(n)=nansum(P(tcr)); % Event rainfall volume
    Vt(n)=nansum(Q(th)); % Event flow volume
    Vb(n)=nansum(Qb(th)); % Event baseflow volume by FRCK

    ccr(n)=nansum(P(tcr)'.*tcr)/Vcr(n); % Centroid of rainfall event
    tE(n,:)=[min(tcr) min(th) max(tcr) max(th) ccr(n) ch(n) ph(n)]; % Event timing
  end
end

lag=ch-ccr; % Event time lag
BFI=Vb./Vt; % Event-based baseflow index by FRCK
RC=(Vt-Vb)./Vcr; % Event-based runoff coefficient by FRCK

Mat_RE.D=[Dcr Dh lag];
Mat_RE.V=[Vcr Vt Vb];
Mat_RE.tE=tE;
Mat_RE.R=[RC BFI];
end
