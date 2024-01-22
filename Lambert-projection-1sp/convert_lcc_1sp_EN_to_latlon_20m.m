%% convert_lcc_1sp_EN_to_latlon.m
% PL 01.06.2017
% Script to convert co-ordinates specified in easting, northing to latitude,
% longitude.
% This works with co-ordinates in LCC (Lambert conformal conic) projections
% with one standard parallel.
%
% Reverse procedure is in the function convert_lcc_1sp_latlon_to_EN.m
%%
%
% See: EPSG Guidance Note Number 7. European Petroleum Survey Group.
% POSC literature pertaining to Coordinate Conversions and Transformations including Formulas, p. 17-18.
%
%%
clc;

dataname1 = 'CE2_GRAS_DEM_20m_E009_35N010W_A';
dataname2 = 'CE2_GRAS_DEM_20m_E010_35N010E_A';
dataname3 = 'CE2_GRAS_DEM_20m_F010_21N009W_A';
dataname4 = 'CE2_GRAS_DEM_20m_F011_21N009E_A';

dataname = dataname2;
%% Moon Area
fileID = fopen(['..\DEM_50m\',dataname,'.tfw']);

tfw = textscan(fileID,'%f');
tfw = tfw{1,1};

[DEM,R4] = readgeoraster(['..\DEM_50m\',dataname,'.tif']);
% figure;imagesc(DEM4);colorbar;caxis([-5e3,max(DEM4,[],'all')*1.2]);
% axis equal
DEMsize = size(DEM);
numE = DEMsize(2);
numN = DEMsize(1);

tfw = [tfw;
numE;
numN;
];

x = tfw(5)+(0:tfw(7)-1)*tfw(1);
y = tfw(6)+(0:tfw(8)-1)*tfw(4);
E = reshape(repmat(x,tfw(8),1),[],1);
N = repmat(y',tfw(7),1);

% 转换 PARAMETER["Central_Meridian"], PARAMETER["Standard_Parallel_1"],
latitude_origin = str2double((dataname(23:24)));
longitude_origin = str2double((dataname(26:28)));
NS = (dataname(25));
EW = (dataname(29));

if NS == 'S'
    latitude_origin = latitude_origin*(-1);
end
if EW == 'W'
    longitude_origin = longitude_origin*(-1);
end

%% select projection to use
pmoon.r0=1737400; % value used by MERA for earth radius
pmoon.a=pmoon.r0; pmoon.b=pmoon.r0;  % ellipsoid semi-major and semi-minor axes. I think MERA uses a sphere.

pmoon.phi_0= latitude_origin   ;    % latitude of natural origin
pmoon.lambda_0= longitude_origin ;   % longitude of natural origin 
pmoon.k0= 1;            % scale factor at origin (1SP)
pmoon.FE=0 ;            % meters; this is the false easting of the MERA projection natural origin. it is probably zero (location of natural origin 15 E)
pmoon.FN=0;             % also probably zero.

pmoon.lambda_f= longitude_origin ;       % longitude of false origin. 

p=pmoon;

%% conversion of angles from degrees to rad
d2r=pi./180;
p.phi_0=p.phi_0.*d2r;
p.lambda_0=p.lambda_0.*d2r;
p.lambda_f=p.lambda_f.*d2r;

%% calculated  projection values - some of these are redundant in the 1sp case
% (see EPSG document section 14.1.2 and 1.4.1.1)
p.f=(p.a-p.b)./p.a;                 % flattening [VERIFIED]
p.e=sqrt(2*p.f-p.f^2);            % eccentricity [VERIFIED]
p.eprime=sqrt(p.e^2./(1-p.e^2));   % second eccentricity

%% iteration required for reverse conversion (EN to latlon):
% first set up the initial guess

%% n,m0,F,ro:
n=sin(p.phi_0); % [VERIFIED]
m0=cos(p.phi_0)./sqrt(1-(p.e^2).*(sin(p.phi_0)).^2); % [VERIFIED]
t0=(tan(pi./4 - p.phi_0./2))./( (1-p.e.*sin(p.phi_0))./(1+p.e.*sin(p.phi_0))).^(p.e./2);
F=m0./(n*(t0.^n)); % [VERIFIED]
r0=p.a.*F*(t0.^n); % [verified]

rprime=sign(n).*sqrt( (E-p.FE).^2 + (r0-(N-p.FN)).^2) ; % [VERIFIED]
tprime=(rprime./(p.a.*p.k0.*F)).^(1./n); % [VERIFIED]


thetaprime=atan( (E-p.FE)./(r0-(N-p.FN)) ); % [VERIFIED]
phi_trial=pi./2 - 2.*atan(tprime); % initial guess of phi


tol=0.0001; % convergence tolerance, degrees

clear E N rprime
%% iterate until phi  has  converged.
num_it=0; err=Inf; % intialisation of convergence variables
phi=phi_trial;
t=tprime;
while (abs(err)>tol)
    phinew=pi./2 - 2.*atan( t.*( (1-p.e.*sin(phi))./(1+p.e.*sin(phi)) ).^(p.e./2) );
    err=phinew-phi;
    phi=phinew; % update phi for next iteration
    disp(num_it);
    num_it=num_it+1;
end

clear err phinew tprime  phi_trial
%% lat , lon formulae as per 2SP case (EPSG document p. 17):
%% final version of t
% t=(tan(pi./4 - phi./2))./( (1-p.e.*sin(phi))./(1+p.e.*sin(phi))).^(p.e./2);
lambda=thetaprime./n  + p.lambda_f;

clear thetaprime
%% finish
disp('----Convert Finish----');
% disp(['phi = ',num2str(phi)]);
% disp(['lambda =',num2str(lambda)]);
% disp(['Degrees: ',num2str(phi./d2r),' ,',num2str(lambda./d2r)]);
phi = phi./d2r;         % 纬度
lambda = lambda./d2r;   % 经度

%%
DEM = reshape(DEM,[],1);

% remove invalid area
Idx = find(DEM > -2e4);
phi = phi(Idx);
lambda = lambda(Idx);
DEM = DEM(Idx);

scatter3(lambda(1:5000:end),phi(1:5000:end),DEM(1:5000:end),5,DEM(1:5000:end),'filled');

savepath = ['..\MOON_50m_convert\',dataname,'.mat'];

% save(string(savepath),"lon","lat","DEM",'-v7.3');

