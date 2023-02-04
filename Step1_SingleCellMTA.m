%--------------------------------------------------------------------------
% Generate single-cell MTA
%
% Written by K. Kunida & Y. Sakumura
% Feb. 5th, 2023
% Ver. 1.0    
%
% Kunida et al. Cell Reports, 2023, https://doi.org/10.1016/j.celrep.2023.112071
%--------------------------------------------------------------------------

clear variables;
close all;

%--------------------------------------------------------------------------
% Users should edit their own data files and parameters
%
% Variable names should be specified below:
% - 'V_ref'     : Reference velocities (RefID x time).
% - 'VelHeatmap': spatio-temporal heatmap of veloccity (space x time)
% - 'ActHeatmap': spatio-temporal heatmap of activity (space x time)
%
%  This code standardizes the activity heatmap.
%--------------------------------------------------------------------------

%---
% Specify the file of reference velocities
% Rows are time series, columns are reference velocity IDs
% The variable 'V_ref' (RefID x time) is stored in this file.

RefVelFile = 'SampleData/ReferenceVelocity.mat'; % variable 'V_ref' is stored

%---
% Specify the data dirctory name containing,
% - 'VelHeatmap.mat' (variable 'VelHeatmap' is stored)
% - 'ActHeatmap.mat' (variable 'ActHeatmap' is stored)
% 
% The output MTA file 'MTA.mat' is laso saved here.

%HeatmapDir    = 'SampleData/rac-08/';   % Figure 1 of Kunida et al. Cell Rep 2023.
%HeatmapDir    = 'SampleData/rac-01/';
HeatmapDir    = 'SampleData/rac-03/';

%---
% Specify Reference Vel indices with a small sample of histogram in Figure 1
% e.g. If you use all ref V, NotUsedRefVidx = [];

NotUsedRefVidx = [3, 6, 9, 12, 16];

%---
% Tolerance for classification as reference velocity

Tol = 10; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following need not be edited
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ActFile = sprintf('%s/ActHeatmap.mat', HeatmapDir);
VelFile = sprintf('%s/VelHeatmap.mat', HeatmapDir);
load(ActFile, 'ActHeatmap'); % Load ActHeatmap
load(VelFile, 'VelHeatmap'); % Load VelHeatmap
load(RefVelFile); % Load the reference velocity

% Data standardization (mean 0, std 1)
ActHeatmapNorm = (ActHeatmap - mean(ActHeatmap(:))) / std(ActHeatmap(:));

%--------------------------------------------------------------------------
% PARAMETERS
%--------------------------------------------------------------------------
refVelSize = size(V_ref,2); % Reference velocity length
refVelNum = size(V_ref,1); % Reference velocity number
refVelNum1 = refVelNum + 1; % One additional velocity class that was not classified in any reV (Others)

% Extract and save the MTA results excluding high velocity ID
usedRefVidx = 1:refVelNum1;
i = ismember(usedRefVidx,NotUsedRefVidx);
usedRefVidx = usedRefVidx(i==0);

%--------------------------------------------------------------------------
% DATA RESHAPE
%--------------------------------------------------------------------------
% Reconstruction of spatio-temporal data into time-series data.
HeatmapWidth = size(ActHeatmapNorm, 2);
Ncycle = (HeatmapWidth-(rem(HeatmapWidth,refVelSize)))./refVelSize;

Vel = [];
Act = [];
AR = [];
for c=1:Ncycle
  frange=[refVelSize*(c-1)+1:refVelSize*c];
  Vel = [Vel;VelHeatmap(:,frange)];      % velocity
  Act = [Act;ActHeatmapNorm(:,frange)];  % standardized activity
  AR = [AR;ActHeatmap(:,frange)];    % original activity
end
Vel = Vel';  % time x sample
Act = Act';  % time x sample
AR = AR';    % time x sample
V_ref=V_ref';  % time x Vref ID

%--------------------------------------------------------------------------
% LABELING REFERENCE VELOCITY INDEX
%--------------------------------------------------------------------------
% Labeling velocity time-series data with reference velocity index
SampleNum=size(Vel,2); % Number of velocity time series data
RSS_Value=zeros(refVelNum, SampleNum);

% make RSS Matrix
for vel_idx=1:SampleNum
  for refVel_idx=1:refVelNum
    RS=(Vel(:,vel_idx)-V_ref(:,refVel_idx)).^2;
    RSMAX=max(RS);
    RSSUM=sum(RS);
    
    % maximum error is less than Tol per time & error sum is less than Tol
    if RSMAX < (Tol/refVelSize) && RSSUM <Tol
      RSS_Value(refVel_idx,vel_idx)=RSSUM;
    else
      % Large errors are made to RSS
      RSS_Value(refVel_idx,vel_idx)=Tol;
    end
  end
end

[minRSS, classidx] = min(RSS_Value); % Find the reference ID with the smallest RSS
classidx(minRSS == Tol) = 0; % If the minimum value is the tolerance parameter, select 0 as ID

%--------------------------------------------------------------------------
% PAIRING VELOCITY AND ACTIVITY
%--------------------------------------------------------------------------
% Extraction of velocity and activity time-series data
scMTA  = struct;  % single-cell MTA
param.dir = HeatmapDir;
param.refVelSize = refVelSize;
param.refVelNum = refVelNum;
param.usedRefVidx = usedRefVidx;

ref_id = 1;
Nsample = zeros(1,refVelNum1);

for ref_i = usedRefVidx
  if ref_i == refVelNum1  % Class 'Others'
    idx = find(classidx == 0);
  else
    idx = find(classidx == ref_i);
  end

  Nsample(ref_i) = length(idx);

  % Extract time series classified to 'class' label
  VI = Vel(:,idx)'; % time series idx x time
  AI = Act(:,idx)'; % time series idx x time

  scMTA(ref_i).V = VI;  % set of time series  (id x time)
  scMTA(ref_i).A = AI;  % set of time series  (id x time)
  scMTA(ref_i).Vavg = mean(VI);  % average of the cell
  scMTA(ref_i).Aavg = mean(AI);  % average of the cell

  ref_id = ref_id+1;
end

% SAVE scMTA results
outfile = sprintf('%s/scMTA.mat', HeatmapDir);
save(outfile, 'scMTA', 'param');

%--------------------------------------------------------------------------
% FIGURE DRAWING
%--------------------------------------------------------------------------

colors = [...
    242 255 157;...
    141 222 253;...
    224 228 255;...
    131 252 180;...
    255 140 120;...
    0 0 256;...
    ];
colors=colors./300;

Vmat = Vel;

r11=find(classidx==1|classidx==2);
Vmat(:,r11)=0;
r12=find(classidx==4|classidx==5);
Vmat(:,r12)=0.25;
r13=find(classidx==7|classidx==8);
Vmat(:,r13)=0.5;
r14=find(classidx==10|classidx==11);
Vmat(:,r14)=0.75;
r15=find(classidx==13|classidx==14|classidx==15);
Vmat(:,r15)=1;
r10=find(classidx==0|classidx==3|classidx==6|classidx==9|classidx==12);
Vmat(:,r10)=1.5;

% reconstruction of spatio-temporal data
VelHeatmap_re=[];
ActHeatmap_re=[];
VelHeatmap_reMTA=[];
for j =1:Ncycle
  srange=201*(j-1)+1:201*(j);
  VelHeatmap_re=[VelHeatmap_re,Vel(:,srange)'];
  ActHeatmap_re=[ActHeatmap_re,Act(:,srange)'];
  VelHeatmap_reMTA=[VelHeatmap_reMTA,Vmat(:,srange)'];
end

%--------------------------
% Draw histogram of sample number for each class
f1 = figure(1); clf;
f1.Position = [100 100 300 300]; 

bar(Nsample);
xticks([1:refVelNum1]);
xticklabels({[1:refVelNum],'Others'});
xlabel('Reference velocity index');
ylabel('Count');
title('Number of velocity samples classified');

figfile = sprintf('%s/SampleHistogram.png', HeatmapDir);
exportgraphics(gcf, figfile, 'Resolution', 600);

%--------------------------
% Draw velocity heatmap
f2 = figure(2); clf;
f2.Position = [200 200 400 400]; 

pcolor(VelHeatmap_re);
shading flat;
colormap( 'default' );
clim([-2 2]);
for xl=1:Ncycle
  xhline=line([7*(xl-1)+1 7*(xl-1)+1],[0 size(VelHeatmap_re,1)]);
  xhline.LineWidth=1;
  xhline.LineStyle='--';
end
for xl=1:10
  yhline=line([0 size(VelHeatmap_re,2)],[20*(xl-1)+1 20*(xl-1)+1]);
  xhline.LineWidth=1;
  xhline.LineStyle='--';
end
xlabel('Time (min)');
ylabel('Edge marker index');
title('Velocity heatmap');

figfile = sprintf('%s/VelHeatmap.png', HeatmapDir);
exportgraphics(gcf, figfile, 'Resolution', 600);

%--------------------------
% Draw activity heatmap
f3 = figure(3); clf;
f3.Position = [300 300 400 400]; 

pcolor(ActHeatmap_re);
shading flat;
colormap( 'default' );
clim([-2 2]);
for xl=1:Ncycle
  xhline=line([7*(xl-1)+1 7*(xl-1)+1],[0 size(ActHeatmap_re,1)]);
  xhline.LineWidth=1;
  xhline.LineStyle='--';
end
for xl=1:10
  yhline=line([0 size(ActHeatmap_re,2)],[20*(xl-1)+1 20*(xl-1)+1]);
  xhline.LineWidth=1;
  xhline.LineStyle='--';
end
xlabel('Time (min)');
ylabel('Edge marker index');
title('Activity heatmap');

figfile = sprintf('%s/ActHeatmap.png', HeatmapDir);
exportgraphics(gcf, figfile, 'Resolution', 600);

%--------------------------
% Draw label heatmap
f4 = figure(4); clf;
f4.Position = [400 400 400 400]; 

pcolor(VelHeatmap_reMTA);
shading flat;
colormap(colors);
for xl=1:Ncycle
  xhline=line([7*(xl-1)+1 7*(xl-1)+1],[0 size(VelHeatmap_re,1)]);
  xhline.LineWidth=1;
  xhline.LineStyle='--';
end
for xl=1:10
  yhline=line([0 size(VelHeatmap_re,2)],[20*(xl-1)+1 20*(xl-1)+1]);
  xhline.LineWidth=1;
  xhline.LineStyle='--';
end
xlabel('Time (min)');
ylabel('Edge marker index');
title('Label heatmap');

figfile = sprintf('%s/LabelHeatmap.png', HeatmapDir);
exportgraphics(gcf, figfile, 'Resolution', 600);

%--------------------------
% Draw time series
f5 = figure(5); clf;
f5.Position = [190 700 1700 300]; 

t = 1:refVelSize;
Ndraw = 30;
fprintf(1, 'Drawing randomly selected time series...   ');

for ref_i = 1:length(usedRefVidx)
  n = usedRefVidx(ref_i);

  name = sprintf('Ref V #%d', n);

  subplot(2,length(usedRefVidx), ref_i);
  yline(0); hold on;
  v = scMTA(n).V;
  
  if(size(v,1) >= Ndraw)
    rndidx = randperm(size(v,1),Ndraw); % pick up time series randomly
  else
    rndidx = 1:size(v,1);
  end
  
  for j = 1:length(rndidx)
    p = plot(t, v(rndidx(j),:)); p.Color(4) = 0.3;
  end
  plot(t, mean(v(rndidx,:)), 'k-', 'LineWidth', 2);
  title(name);
  xlabel('Time (min)'); ylabel('Velocity (um/min)');
  axis([t(1)-1 t(end)+1 -2 2]);

  subplot(2,length(usedRefVidx),ref_i+length(usedRefVidx));
  yline(0); hold on;
  a = scMTA(n).A;
  for j = 1:length(rndidx)
    p = plot(t, a(rndidx(j),:)); p.Color(4) = 0.3;
  end
  plot(t, mean(a(rndidx,:)), 'k-', 'LineWidth', 2);
  xlabel('Time (min)'); ylabel('Activity (a.u.)');
  axis([t(1)-1 t(end)+1 -3 3]);

end

figfile = sprintf('%s/scMTA_timeseries.png', HeatmapDir);
exportgraphics(gcf, figfile, 'Resolution', 600);

fprintf(1, 'Done.\n');
