%--------------------------------------------------------------------------
% Generate cell-common MTA
%
% Written by K. Kunida & Y. Sakumura
% Feb. 5th, 2023
% Ver. 1.0    
%
% Kunida et al. Cell Reports, 2023, https://doi.org/10.1016/j.celrep.2023.112071
%--------------------------------------------------------------------------

clear variables;
close all;

%---
% Specify the data dirctories
RootDir  = 'SampleData';
MTADirs  = {'rac-01', 'rac-03', 'rac-08'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following need not be edited
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract parameters
dir = sprintf('%s/%s', RootDir, MTADirs{1});
scMTAfile = sprintf('%s/scMTA.mat', dir);
load(scMTAfile,'scMTA', 'param');

refVelSize = param.refVelSize; % Reference velocity length
refVelNum  = param.refVelNum; % Reference velocity number
usedRefVidx = param.usedRefVidx;

%-------------------------------
% Load and format time series
Vts = cell(refVelNum, length(MTADirs));
Ats = cell(refVelNum, length(MTADirs));

for dir_i = 1:length(MTADirs)
  dir = sprintf('%s/%s', RootDir, MTADirs{dir_i});
  fprintf(1,'%s\n', dir);
  scMTAfile = sprintf('%s/scMTA.mat', dir);
  load(scMTAfile,'scMTA');

  for ref_i = 1:refVelNum

    Vts{ref_i, dir_i} = scMTA(ref_i).V;
    Ats{ref_i, dir_i} = scMTA(ref_i).Aavg;

  end
end

%-------------------------------
% Average over cells

ccMTA  = struct;  % cell-common MTA

Vavg = zeros(1,refVelSize); % ccMTA velocity
Aavg = zeros(1,refVelSize); % ccMTA activity

allVel = cell(1,1);
scMTAAct = cell(1,1);

minNvel = 10e5;

for ref_i = 1:refVelNum

  v = [];
  a = [];
  for dir_i = 1:length(MTADirs)
    v = [v; Vts{ref_i, dir_i}];
    a = [a; Ats{ref_i, dir_i}];
  end

  if(minNvel > size(v,1))
    minNvel = size(v,1);
  end

  allVel = v;
  scMTAAct = a;

  Vavg = mean(v,'omitnan');  % ccMTA velocity
  Aavg = mean(a,'omitnan');  % ccMTA activity

  ccMTA(ref_i).V = allVel;
  ccMTA(ref_i).A = scMTAAct;
  ccMTA(ref_i).Vavg = Vavg;
  ccMTA(ref_i).Aavg = Aavg;
end

%-------------------------------
% Save ccMTA
outfile = sprintf('%s/ccMTA.mat', RootDir);
save(outfile, 'ccMTA', 'param');


%-------------------------------
% Draw time series

f = figure(5); clf;
f.Position = [190 700 1700 300]; 

% Am = scMTA(d,:);
% V  = MT_Vel(d,:);
% A  = MT_Act(d,:);

t = 1:refVelSize;
Ndraw = 30;
fprintf(1, 'Drawing randomly selected time series...   ');

for ref_i = 1:length(usedRefVidx)
  n = usedRefVidx(ref_i);

  name = sprintf('Ref V #%d', n);

  subplot(2,length(usedRefVidx), ref_i);
  yline(0); hold on;
  v = ccMTA(n).V;
  
  if(size(v,1) >= Ndraw)
    rndidx = randperm(size(v,1),Ndraw); % pick up time series randomly
  else
    rndidx = 1:size(v,1);
  end
  
  for j = 1:length(rndidx)
    p = plot(t, v(rndidx(j),:)); p.Color(4) = 0.3;
  end
  plot(t, ccMTA(n).Vavg, 'k-', 'LineWidth', 2);
  title(name);
  xlabel('Time (min)'); ylabel('Velocity (um/min)');
  axis([t(1)-1 t(end)+1 -2 2]);

  subplot(2,length(usedRefVidx),ref_i+length(usedRefVidx));
  yline(0); hold on;
  a = ccMTA(n).A;
  
  p = plot(t, a); %p.Color(4) = 0.3;
  
  plot(t, ccMTA(n).Aavg, 'k-', 'LineWidth', 2);
  xlabel('Time (min)'); ylabel('Activity (a.u.)');
  axis([t(1)-1 t(end)+1 -3 3]);

end

figfile = sprintf('%s/ccMTA_timeseries.png', RootDir);
exportgraphics(gcf, figfile, 'Resolution', 600);

fprintf(1, 'Done.\n');


