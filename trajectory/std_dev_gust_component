function [gammaG] = std_dev_gust_component(rho,S,L,Cla,Cma,sigmaG,lg,m,a,g,Ip,d,h)
% Objective:
% This function computes the standard deviation in one component of the
% flight path angle due to turbulence (wind gusts)
%
% Source: http://arc.aiaa.org/doi/pdf/10.2514/1.A32037
% input variables:
% rho - density of air [kg/m^3]
% S - aerodynamic reference area [m^2]
% L - launcher length [m]
% Cla - coefficient of lift derivative wrt angle of attack [1/rad]
% Cma - pitch moment coeff. derivative wrt angle of attack [1/rad]
% sigmaG - standard deviation in wind (gust) velocity
% lg - longitudinal turbulence scale length [m]
% m - rocket mass [kg]
% a - acceleration along the flight path [m/s^2]
% g - acceleration due to gravity [m/s^2]
% Ip - moment of inertia (of rocket) about the pitch or yaw axis[kg-m^2]
% d - aerodynamic reference length (body diameter) [m]
% h - altitude [m]
%
% output variables:
% gammaG - the standard deviation in one component of the flight path angle
% due to turbulence (wind gusts) 
% functions called:
% none
%
lambdaP = -0.5*rho*S*d*Cma/Ip;
lambdaG = 1/lg;
lambdaGstar = 1/(0.59*lg);


K=rho*S*L*Cla*sigmaG*lambdaP/(m*(2*a-g))*sqrt(a*lambdaGstar/(2*(lambdaP^2+lambdaGstar^2)));
gammaG = K*((h/L)-(h/L).^(g/2/a));

end
