% fm_bbl_setup
clear
global dz w K Erate

% Constants
vk = 0.41;    % ( ) von Karman's const.
g = 9.81;     % (m/s) grav. accel.
rhow = 1025;  % (kg/m3) water density
mu = 1e-3;    % (whatevertheseare) molecular viscosity
nu = mu/rhow; % (m^2/s) kinematic viscosity
lat = 45;     % (degrees) latitude
Omega = 7.292115e-5; % Earth rotation rate (rad/s) Kantha & Clayson p.842
f = 2*Omega .* sin( lat*pi/180 ); % (1/s) Corliolis parameter

% User input
h = 12;    % total water depth (m)
dzz = .5;  % vertical

% Switches to enable/disable various processes
do_settle = 1
do_floc = 1
do_sand = 1
do_wibergK=1
do_altG=1
do_floc=1

% Sand characteristics
nps = 1;        % number of sand classes
rhoss = 2650;   % sand particle density
d50s = 0.15e-3; % sand median size(m)
tauc = 0.16 % critical stress for erosion (Pa)
wss = -0.012 % sand settling velocity (m/s)
E0 = 5e-6; % erosion-rate parameter (kg m2 s-1)
Cb = 0.6;  % bed fraction (1-porosity) ()
kN = 2.5*d50s; % Soulsby (1997) p. 48
zoc = kN/30;

% Floc population parameters
npmud = 15;
Dp = 4e-6;
nf = 2.2;
rhos = 2650;
f_clim = 0.001;

% For winterwerp equilibrium diameter
ka=20;
%kb=ka/(1.6429*1e-3);
kb=ka/(1.22*1.55*1e-3);
% Floc model parameters
l_ADS=0;
l_ASH=1;
l_COLLFRAG=0;
f_dmax=0.001500;
f_nb_frag=2.;
alpha=0.35;
beta=0.15;
f_ater=0.;
f_ero_frac=0.0;
f_ero_nbfrag=2.;
f_ero_iv=1;
f_mneg_param=0.00;
f_collfragparam=0.00;
dfragmax=0.00003;
epsilon = 1e-10;

% Floc population distribution
Df = logspace(log10(4e-6),log10(1500e-6),npmud)';
% Df=Df(7:end);
% np=length(Df);
rhof = rhow+(rhos-rhow)*(Dp./Df).^(3-nf);
volf = (pi/6.0)*Df.^3;
mass = volf.*(rhof-rhow);
wsf =  -g*(rhof-rhow).*Df.^2 ./(18.*0.001);

% Initialize bbl
zf = (0:dzz:h)';         % faces
dz = diff(zf);
zc = (dzz/2:dzz:h-dzz/2)'; % centers
nzc = length(zc);

% Initial mud concentration
massconc = 0.10*ones(nzc,npmud)./npmud; % g/l = kg/m3

%% Time series of forcing
Tc = 4*3600; % tidal current period
omegac = 2*pi/(Tc);
Ttot = 1.75*Tc; % (s) overall time interval for solution
%Ttot = .2*Tc; % short run
CFLt = dzz/abs(wss)
dt = 10; % (s) for calculations
dtint = fix(CFLt/2) % (s) dt for settling solution

t = 0:dt:Ttot;
nt = length(t);

Tw = 10; % wave period (s)
omegaw = 2*pi/Tw;

% steady forcing, then off
if(0)
   uc = 0.15*ones(size(t));
   uw = 0.4*ones(size(t));
   uw(t>(2*3600))=0; % shut off waves after 2 hours.
   uc(t>(2*3600))=0;
end

% time-dependent forcing
if(1)
   Te = 2*3600; % wave event period
   uc = 0.2+0.2* sin( omegac*t );
   te = 0:dt:Te;
   uwa = 0.4 % uw amplitude
   
   uwe = uwa/2+(uwa/2)*cos(te*2*pi/Te - pi)
   tes = .75*Tc;
   ies = find(t>=tes,1,'first')
   uw = zeros(size(t));
   uw(ies:ies+length(uwe)-1)=uwe;
   figure(3); clf; plot(t/3600,uw,'-b',t/3600,uc,'-r')
end
%% Exectute main program
fm_bbl_main