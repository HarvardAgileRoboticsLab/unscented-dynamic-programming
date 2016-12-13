function p = yak_parameters()

% ----- Model Parameters for Yak RC Airplane ----- %

g = 9.81; %Gravitational acceleration (m/s^2)
rho = 1.2; %Air density at 20C (kg/m^3)
m = .075; %Mass of plane (kg)

Jx = 4.8944e-04; %roll axis inertia (kg*m^2)
Jy = 6.3778e-04; %pitch axis inertia (kg*m^2)
Jz = 7.9509e-04; %yaw axis inertia (kg*m^2)

J = diag([Jx Jy Jz]);
Jinv = diag([1/Jx 1/Jy 1/Jz]); %Assuming products of inertia are small

Jm = .007*(.0075)^2 + .002*(.14)^2/12; %motor + prop inertia (kg*m^2)

% All lifting surfaces are modeled as unsweapt tapered wings
b = 45/100; %wing span (m)
l_in = 6/100; %inboard wing length covered by propwash (m)
cr = 13.5/100; %root chord (m)
ct = 8/100; %tip chord (m)
cm = (ct + cr)/2; %mean wing chord (m)
S = b*cm; %planform area of wing (m^2)
S_in = 2*l_in*cr;
S_out = S-S_in;
%Ra = b^2/S; %wing aspect ratio (dimensionless)
Rt = ct/cr; %wing taper ratio (dimensionless)
r_ail = (b/6)*(1+2*Rt)/(1+Rt); %aileron moment arm (m)

ep_ail = 0.63; %flap effectiveness (Phillips P.41)
trim_ail = 106; %control input for zero deflection
g_ail = (15*pi/180)/100; %maps control input to deflection angle

b_elev = 16/100; %elevator span (m)
cr_elev = 6/100; %elevator root chord (m)
ct_elev = 4/100; %elevator tip chord (m)
cm_elev = (ct_elev + cr_elev)/2; %mean elevator chord (m)
S_elev = b_elev*cm_elev; %planform area of elevator (m^2)
Ra_elev = b_elev^2/S_elev; %wing aspect ratio (dimensionless)
r_elev = 22/100; %elevator moment arm (m)

ep_elev = 0.88; %flap effectiveness (Phillips P.41)
trim_elev = 106; %control input for zero deflection
g_elev = (20*pi/180)/100; %maps control input to deflection angle

b_rud = 10.5/100; %rudder span (m)
cr_rud = 7/100; %rudder root chord (m)
ct_rud = 3.5/100; %rudder tip chord (m)
cm_rud = (ct_rud + cr_rud)/2; %mean rudder chord (m)
S_rud = b_rud*cm_rud; %planform area of rudder (m^2)
Ra_rud = b_rud^2/S_rud; %wing aspect ratio (dimensionless)
r_rud = 24/100; %rudder moment arm (m)
z_rud = 2/100; %height of rudder center of pressure (m)

ep_rud = 0.76; %flap effectiveness (Phillips P.41)
trim_rud = 106; %control input for zero deflection
g_rud = (35*pi/180)/100; %maps from control input to deflection angle

trim_thr = 24; %control input for zero thrust (deadband)
g_thr = 0.006763; %maps control input to Newtons of thrust
g_mot = 3000*2*pi/60*7/255; %maps control input to motor rad/sec

% --- Pack everything into a struct --- %
p = struct();
p.g = g;
p.rho = rho;
p.m = m;
p.Jx = Jx;
p.Jy = Jy;
p.Jz = Jz;
p.J = J;
p.Jinv = Jinv;
p.Jm = Jm;
p.b = b;
p.l_in = l_in;
p.cr = cr;
p.ct = ct;
p.cm = cm;
p.S = S;
p.S_in = S_in;
p.S_out = S_out;
%P.Ra = Ra;
p.Rt = Rt;
p.r_ail = r_ail;
p.ep_ail = ep_ail;
p.trim_ail = trim_ail;
p.g_ail = g_ail;
p.b_elev = b_elev;
p.cr_elev = cr_elev;
p.cm_elev = cm_elev;
p.S_elev = S_elev;
p.Ra_elev = Ra_elev;
p.r_elev = r_elev;
p.ep_elev = ep_elev;
p.trim_elev = trim_elev;
p.g_elev = g_elev;
p.b_rud = b_rud;
p.cr_rud = cr_rud;
p.ct_rud = ct_rud;
p.cm_rud = cm_rud;
p.S_rud = S_rud;
p.Ra_rud = Ra_rud;
p.r_rud = r_rud;
p.z_rud = z_rud;
p.ep_rud = ep_rud;
p.trim_rud = trim_rud;
p.g_rud = g_rud;
p.trim_thr = trim_thr;
p.g_thr = g_thr;
p.g_mot = g_mot;

end
