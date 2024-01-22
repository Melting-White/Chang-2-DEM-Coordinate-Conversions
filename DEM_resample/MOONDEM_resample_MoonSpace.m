clc;
clear;

%% DEM resample
latitude = [];      % phi
lontitude = [];     % lambda
DEM2 = [];

load ..\MOON_50m_convert\MOONDEM_50m_63S015E.mat
latitude = [latitude;phi];
lontitude = [lontitude;lambda];
DEM2 = [DEM2;DEM];

load ..\MOON_50m_convert\MOONDEM_50m_63S045E.mat
latitude = [latitude;phi];
lontitude = [lontitude;lambda];
DEM2 = [DEM2;DEM];

load ..\MOON_50m_convert\CE2_GRAS_DEM_50m_M005_77S023E_A.mat
latitude = [latitude;lat];
lontitude = [lontitude;lon];
DEM2 = [DEM2;DEM];

load ..\MOON_50m_convert\CE2_GRAS_DEM_50m_M005_77S023E_A.mat
latitude = [latitude;lat];
lontitude = [lontitude;lon];
DEM2 = [DEM2;DEM];

clear phi lambda DEM lat lon

% scatter3(lontitude(1:1000:end),latitude(1:1000:end),DEM2(1:1000:end),5,DEM2(1:1000:end));
%% Focusing area
resize_lon = [20,45];
resize_lat = [-80,-60];

Idx = latitude>resize_lat(1) & latitude<resize_lat(2) & lontitude>resize_lon(1) & lontitude<resize_lon(2);
resize_lontitude = lontitude(Idx);
resize_latitude = latitude(Idx);
resize_DEM = DEM2(Idx);

% scatter3(resize_lontitude(1:500:end),resize_latitude(1:500:end),resize_DEM(1:500:end),5,resize_DEM(1:500:end));
clear DEM2 latitude lontitude Idx

%% Resample
deltalon = 0.004;
deltalat = 0.0015;

[re_lontitude,re_latitude] = meshgrid(21:deltalon:42,-78:deltalat:-62);
re_DEM = griddata(resize_lontitude,resize_latitude,double(resize_DEM),re_lontitude,re_latitude,'natural');

DEM = flipud(re_DEM);
lontitude = 21:deltalon:42;
latitude = -78:deltalat:-62;

% save MOONDEM_50m_MoonSpace_21plus16deg_natural.mat DEM lontitude latitude 
% load MOONDEM_50m_MoonSpace_21plus16deg_natural.mat
DEM = flipud(DEM);
%%
viewAng = [-45,78];
ampAxis = [7e8, 17e8];
heiAxis = [-4e3, 6e3];
sceneAxis = [21,42,-78,-62,-0.7e4,0.7e4];

figure;s = surf(lontitude(1:10:end),latitude(1:10:end),DEM(1:10:end,1:10:end),'EdgeColor' ,'none'); 
xlabel('LON/deg');  ylabel('LAT/deg');  zlabel('H/m'); 
view(viewAng); 
axis(sceneAxis)
colorbar;
caxis(heiAxis);

figure;imagesc(DEM);