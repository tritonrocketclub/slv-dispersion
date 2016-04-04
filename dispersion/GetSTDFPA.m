function [ gamma ] = GetSTDFPA( )
% GETSTDFPA
% 
% Objective: Calculate the standard deviation in the flight path angle
%   of the rocket, to be used in calculating the boost phase trajectory.
%
% input variables:
%
%
% output variables:
%   gamma - the standard deviation on the flight path angle, in rad
%
% functions called:
%   none
%

%Define Basic Variables

%Fin Variables
N = 4; %number of fins
c_R = .457; %fin root chord [m]
c_T = .0508; %fin tip chord [m]
fin_height = .152; %fin height [m]
S = N*((c_R+c_T)/2)*fin_height; %Aerodynamic reference area, or the total fin areas [m^2]


C_G = 2.27; %Center of mass [m] (From OpenRocket)
C_P = 2.41; %Center of pressure [m] (From OpenRocket)
T=7980; %thrust [N] (average thrust from OpenRocket O motor)
d = .1524; %rocket diameter [m]
A_xc = (d/2)^2*pi; %cross sectional area [m^2]
L_R = 3.2; %Rocket Length [m]
l_throat = L_R - C_G; %distance from nozzle throat to center of mass [m]
L = 4.8768; %Launcher Length [m]
m_0 = 56.2455; %wet rocket mass [kg]
m_b = 13.6; %dry rocket mass [kg]
Area_moment_inertia = (pi/4)*(d/2)^2; %Appromxation of the second moment of area of the rocket [m]
K_r = sqrt(Area_moment_inertia/A_xc);
delta_t = .0010908308; %thrust misalignment angle [radians] (from http://www.rsandt.com/media/Standardized%20Perturbation%20Values.pdf)
omega = .125*K_r; %Lateral offset of the center of mass [m] (from http://www.rsandt.com/media/Standardized%20Perturbation%20Values.pdf)
I_p = 1/12*m_b*L_R^2 + m_b*(.5*L_R - C_G)^2; %Pitch moment of inertia [kg*m^2] (assuming the rocket is a solid rod, using parallel axis theorem)
I_xx = .5*m_b*(d/2); %Roll moment of inertia (assuming the rocket is a solid cylinder)
U_burnout = 950; %Burnout velocity [m/s]
h_burnout = 13500; %Burnout altitude [m]
c_burnout = sqrt(1.4*287*216.7); %speed of sound at burnout
M_burnout = U_burnout/c_burnout; %Mach number at burnout
rho_burnout = .247205; %density of air at burnout altitude [kg/m^3] linearly interpolated from standard atmospheric conditions at 13500m
r_e = 6371000; %mean radius of the earth
g = 9.8*(r_e/(r_e + h_burnout)); %acceleration due to gravity [m/s^2]
a = 5*g; %acceleration (using assumption that g/2a = .1 from the document) [m/s^2]

%Define Variables to calculate pitch wave number (assume relationship
%between coefficient of normal force is linear with angle of attack)

%C_N_alpha = 8.59436; %Derivative of C_N with respect to alpha http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660004864.pdf
%C_M_alpha = C_N_alpha*(C_G - C_P)/L_R; %Derivative of C_M with respect to alpha from http://raceriv.com/arma2/AME02.pdf
C_N_alpha = 11.459; %http://naca.central.cranfield.ac.uk/reports/1952/naca-rm-a52c04.pdf
C_M_alpha = -28.648; %http://naca.central.cranfield.ac.uk/reports/1952/naca-rm-a52c04.pdf
lambda_p = sqrt(-rho_burnout*A_xc*d*C_M_alpha/(2*I_p)); %pitch wave number www.rsandt.com/media/The%20Pitch-Yaw%20Wave%20Number3.doc

%Define variables to calculate roll wave number

%http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19730012271.pdf
%delta_f = 0.00872665; %mean cant angle [rad]; rocket being built to 0 cant but assuming worst case scenario of .5 degrees
delta_f = .5*.00349066;
b = .152 + d/2; %fin semispan [m] (from rocket center line to tip of the fin)
A = (c_R*b - c_T*(d/2))/(b-d/2); %parameter A [m] www.rsandt.com/media/Fin%20Setting%20Theory[3].doc
C = (c_R - c_T)/(b-d/2); %slope of the chord
q = .5*rho_burnout*U_burnout^2; %dynamic pressure [Pa]
%p = .33069; %roll rate [rad/s] estimated by watching a rocket video https://www.youtube.com/watch?v=3U3GHqewb9Q
%Rolling_moment = q*C_N_alpha*(b-d/2)*delta_f*(A*(b+d/2)/2 - C*(b^2 + b*(d/2) + (d/2)^2)/3) + q*C_N_alpha*(b-d/2)*(p/U_burnout)*(-A*(b^2 b*(d/2) + (d/2)^2)/3 + C*(b^3 + b^2*(d/2) + b*(d/2)^2 +(d/2)^3)/4); %rolling moment
 

C_L_delta = ((C_N_alpha*(b-d/2))/(2*pi*(d/2)^3))*(A*(b+d/2)/2 - C*(b^2 + b*(d/2) + (d/2)^2)/3); %rolling stability derivative
C_L_p = ((C_N_alpha*(b-d/2)*N)/(2*pi*(d/2)^4))*(-A*(b^2 + b*(d/2) + (d/2)^2)/3 + C*(b^3 + b^2*(d/2) + b*(d/2)^2 + (d/2)^3)/4); %rolling damping derivative
lambda_R = -((C_L_delta*N*delta_f)/((d/2)*C_L_p))*(2*pi*rho_burnout*(h_burnout - L)*(d/2)^4*C_L_p)/(2*pi*rho_burnout*(h_burnout - L)*(d/2)^4*C_L_p - I_xx); %Roll wavenumber;

alpha_ss = T*(l_throat*delta_t + omega)/(I_p*U_burnout^2*(lambda_p^2 - lambda_R^2)); 

C_L_alpha = 1.35; %Lift coefficient slope http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930087519.pdf

L_A = q*S*C_L_alpha*alpha_ss;

gamma_T = (L_A/(m_b*g))*((h_burnout/L)^(g/(a) - 1)); %Standard deviation in dispersive flight angle
phi = lambda_R*(h_burnout - L);

%Weather Variables
sigma_G = 2.2352; %Standard deviation in gust velocity (need wind data at launch date) [m/s]
l_G = 300; %Longitudinal turbulence scale length (using assumed value from http://arc.aiaa.org/doi/pdf/10.2514/1.A32037) [m]

Northerly = gamma_T*sin(phi); %Northerly component of flight path angle
Easterly = gamma_T*cos(phi); %Easterly component of flight path angle

gamma_D = sqrt(Northerly^2 + Easterly^2); %Standard deviation in flight path angle due to thrust misalignment [rad]

gamma_G = FPAStdGust(rho_burnout, S, L, C_L_alpha, C_M_alpha, sigma_G, l_G, m_b, a, g, I_p, d, h_burnout); %Standard deviation in flight path angle due to gust [rad]

gamma = sqrt(gamma_D^2 + gamma_G^2); %standard deviation in flight path angle from the vertical down to the velocity vector [rad]

end