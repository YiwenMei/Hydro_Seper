% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last updated on 11/18/2017

%% Functionality:
% The code calculates the rainfall-runoff event properties.

%% Input
% Q: streamflow time series (mm/h);
% P: rainfall time series (mm/h);
% T: event timing stored in a n-by-4 matrix (n is the event number and the
%    four columns stand for 1)the begin of rainfall, 2) begin of flow, 3)
%    end of rainfall and 4) end of flow

%% Output
% Event_pro: event properties stored in the order of
%            1) Volume of rainfall,     2) volume of flow,
%            3) centroid of rainfall,   4) centroid of flow,
%            5) time of peak flow,      5) event time lag,
%            6) dispersion of rainfall, 6) dispersion of flow.

function Event_pro=hydro_pro(Q,P,T)
Vr=nan(size(T,1),1);
Vf=nan(size(T,1),1);

Cr=nan(size(T,1),1);
Cf=nan(size(T,1),1);
Pf=nan(size(T,1),1);

Dr=nan(size(T,1),1);
Df=nan(size(T,1),1);
for i=1:size(T,1)
  tr=T(i,1):T(i,3); % Period of rainfall event
  tf=T(i,2):T(i,4); % Period of flow event

  Vr(i)=nansum(P(tr)); % Event rainfall volume
  Vf(i)=nansum(Q(tf)); % Event flow volume

  Cr(i)=nansum(P(tr)'.*tr)/Vr(i); % Centroid of rainfall event
  Cf(i)=nansum(Q(tf)'.*tf)/Vf(i); % Centroid of flow event
  [~,j]=max(Q(tf));
  Pf(i)=min(tf(j)); % Peak of flow event

  Dr(i)=sqrt(nansum(P(tr)'.*(tr-Cr(i)).^2)/Vr(i));
  Df(i)=sqrt(nansum(Q(tf)'.*(tf-Cf(i)).^2)/Vf(i));
end

Event_pro=[Vr Vf Cr Cf Pf Cf-Cr Dr Df]; % Event properties
end
