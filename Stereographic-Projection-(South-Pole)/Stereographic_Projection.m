%% Stereographic Projection (South Pole)
% LPY 2023/09/05
close all;
clearvars;
% clc;
fclose all;
warning off all;
dbstop if error;
c = 299792458;

datapath = '..\DEM_50m\';
dataname = 'CE2_GRAS_DEM_50m_M001_77S158W_A';
dataname = 'CE2_GRAS_DEM_50m_M002_77S113W_A';
dataname = 'CE2_GRAS_DEM_50m_M003_77S068W_A';
dataname = 'CE2_GRAS_DEM_50m_M004_77S023W_A';
dataname = 'CE2_GRAS_DEM_50m_M005_77S023E_A';
dataname = 'CE2_GRAS_DEM_50m_M006_77S068E_A';
dataname = 'CE2_GRAS_DEM_50m_M007_77S113E_A';
dataname = 'CE2_GRAS_DEM_50m_M008_77S158E_A';
dataname = 'CE2_GRAS_DEM_50m_N001_87S000W_A';

savepath = ['..\MOON_50m_convert\',dataname,'.mat'];

fileID = fopen(['..\DEM_50m\',dataname,'.tfw']);
%% Load Pre Data
TFW = textscan(fileID,'%f');
TFW = TFW{1,1};

Eb = TFW(5);
Nb = TFW(6);
deltaE = TFW(1);
deltaN = TFW(4);

%
[DEM,R4] = readgeoraster([datapath,dataname,'.tif']);
% figure;imagesc(DEM);colorbar;caxis([-5e3,max(DEM,[],'all')*1.2]);
% axis equal

DEMsize = size(DEM);
numE = DEMsize(2);
numN = DEMsize(1);

DEM = reshape(DEM,[],1);

%% Para
R = 1737400;

%% Convert EN to lonlat -- Stereographic Projection (South Pole) 
E1 = Eb + (0:numE-1) * deltaE;
N1 = Nb + (0:numN-1) * deltaN;

[E,N] = meshgrid(E1,N1);
E = reshape(E,[],1);
N = reshape(N,[],1);
% figure; scatter3(E(1:300:end),N(1:300:end),DEM(1:300:end),[],DEM(1:300:end),"filled");
% zlim([-8e3 max(DEM,[],'all')*1.2]);caxis([-5e3,max(DEM,[],'all')*1.2]);

lon = rad2deg(atan2(E,N));       % 东西经
lat = 90 - 2*rad2deg(atan(vecnorm([E,N],2,2)./(2*R)));      % 南纬

Idx = find(DEM > -2e4);
lon = lon(Idx);
lat = -lat(Idx);
DEM = DEM(Idx);

scatter3(lon(1:300:end),lat(1:300:end),DEM(1:300:end),5,DEM(1:300:end),'filled');

save(string(savepath),"lon","lat","DEM",'-v7.3');