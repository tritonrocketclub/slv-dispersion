%------------------------------------------------------
% TRITON ROCKET CLUB
% ROCKET ALTITUDE SIMULATION PROGRAM (MatRocket) v2.0
% (c)2016 TRITON ROCKET CLUB
%
% Authors: Trevor Irwin, Minh Nguyen, Toriana Dabkowski,
%   Robert Bacon, Patrick Lee, et al
%
% This script simulates 3DOF flight paths for a solid fuel rocket
% with a Monte Carlo approximation of dispersion probability.
% The script makes several key simplifying assumptions:
% --The rocket's boost flight path angle is constant, allowing
%   the rocket to be treated as a 3DOF particle during the boost phase
% --Due to the stabilizing nature of the rocket's fins, the rocket is
%   assumed to have zero angle of attack during the coast phase. This
%   allows the rocket to be treated as a 3DOF particle for purposes of drag
%   calculation. The only other force is gravity, acting on the CoM, so
%   rotation is neglected
% --During the recovery phase, the rocket travels always at terminal
%   velocity.
% --Wind adds constant velocity
% 
%
% Resources:
% http://www.rocketmime.com/rockets/rckt_eqn.html
% http://exploration.grc.nasa.gov/education/rocket/drageq.html
%------------------------------------------------------
%
% Clear the workspace of variables
%
clear all;
close all;
%
% Set default RNG for repeatability
%
rng default;
fprintf('\nTRITON ROCKET CLUB ALTITUDE SIMULATION (c)2016\n');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL VARIABLES
% 
% Establish global variables for use with all functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global rocketProp
global simuProp
rocketProp = struct();
simuProp = struct();

global logProp;
logProp = struct;

logProp.drag = zeros(1E6, 0);
logProp.drag_int = 1;
logProp.acc = zeros(1E6, 0);
logProp.acc_int = 1;
logProp.vel = zeros(1E6, 0);
logProp.vel_int = 1;
logProp.rho = zeros(1E6, 0);
logProp.rho_int = 1;
logProp.grav = zeros(1E6, 0);
logProp.grav_int = 1;
logProp.dm = zeros(1E6, 0);
logProp.dm_int = 1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE HANDLING
%
% Look for a local file, RocketProperties.txt
% If not found, ask for an input file of similar format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read the files
%
RProp = ReadFile('RocketProperties.txt', 'Rocket Properties', 0);
DragProp = ReadFile('RocketDrag.txt', 'Rocket Drag', 2);
AtmProp = ReadFile('StandardAtmosphere.txt', 'Standard Atmosphere', 6);
WindProp = ReadFile('WindProfile.txt', 'Winds Aloft', 3);
%
% Rocket Properties
%
% Set tburn to false. If burn time is specified in the
% text file, tburn will be overwritten with a value.
% If it remains false, it will be calculated later.
%
[m, n] = size(RProp{1});
tburn = false;
rocketProp.burnTime = 0;
simuProp.count = 0;
for i=1:m
    switch RProp{1}{i}
        
        % --Rocket Body Properties--
        case 'drymass'
            rocketProp.dryMass = RProp{2}(i); %[kg] mass of structure+propellant
        case 'wetmass'
            rocketProp.wetMass = RProp{2}(i); %[kg] mass of structure+propellant
        case 'radius'
            rocketProp.radius = RProp{2}(i); %[in]radius of body
        case 'numfins'
            rocketProp.numfins = RProp{2}(i); %number of fins on the rocket
        case 'rootchord'
            rocketProp.rootchord = RProp{2}(i); %[m] fin root chord
        case 'tipchord'
            rocketProp.tipchord = RProp{2}(i); %[m] fin tip chord
        case 'finheight'
            rocketProp.finheight = RProp{2}(i); %[m] fin height
        case 'C_G'
            rocketProp.C_G = RProp{2}(i); %[m]Center of mass (From OpenRocket)
        case 'C_P'
            rocketProp.C_P = RProp{2}(i); %[m] Center of pressure (From OpenRocket)
        case 'rocketlength'
            rocketProp.rocketlength = RProp{2}(i); %[m] Rocket length
        case 'launcherlength'
            rocketProp.launcherlength = RProp{2}(i); %[m] Launcher length
        case 'rho'
            simuProp.rho = RProp{2}(i); %[kg/m^3] air density
        case 'meanCantAngle'
            rocketProp.meanCantAngle = RProp{2}(i); %mean cant angle [rad]; rocket being built to 0 cant but assuming worst case scenario of .5 degrees
        case 'liftCoefficientSlope'
            rocketProp.liftCoefficientSlope = RProp{2}(i);
        case 'drag'
            simuProp.drag = RProp{2}(i); %drag coefficient
        case 'gravity'
            simuProp.gravity = RProp{2}(i); %[m/s^2] gravitational acceleration
        
        
        % --Motor Properties--
        case 'thrust'
            rocketProp.thrust = RProp{2}(i); %[N] thrust
        case 'Isp'
            rocketProp.isp = RProp{2}(i); %[s] specific impulse
        case 'Itot'
            rocketProp.iTotal = RProp{2}(i); %[N-s] total impulse
        case 'burn'
            rocketProp.burnTime = RProp{2}(i); %[s]An option to override the burn time
        case 'thrustma'
            rocketProp.thrustma = RProp{2}(i); %[rad] Thrust misalignmenet angle
        
        % --Chute Properties--
        case 'chutedrag'
            rocketProp.chuteCd = RProp{2}(i); %drag coefficient
        case 'chuteradius'
            rocketProp.chuteRadius = RProp{2}(i); %[m] chute radius
        case 'chutespillradius'
            rocketProp.chuteSpillRadius = RProp{2}(i); %[m] chute spill radius
        case 'chutealt'
            rocketProp.chuteAlt = RProp{2}(i); %[m] chute deployment altitude
        
        % --Simulation Properties--
        case 'res'
            simuProp.resolution = RProp{2}(i); %[N] thrust  
        case 'recoveryres'
            simuProp.recoveryRes = RProp{2}(i); % Resolution for recovery section
        case 'count'
            simuProp.count = RProp{2}(i); % Number of simulations to run
        case 'sigma_G'
            simuProp.sigma_G = RProp{2}(i); % [m/s]Standard deviation in gust velocity (need wind data at launch date)
        case 'l_G'
            simuProp.l_G = RProp{2}(i); % [m] Longitudinal turbulence scale length (using assumed value from http://arc.aiaa.org/doi/pdf/10.2514/1.A32037)
    end
end
%
% Convenience Properties
%
rocketProp.area = pi*rocketProp.radius^2;
if rocketProp.burnTime == 0
    rocketProp.burnTime = rocketProp.iTotal/rocketProp.thrust;
end
arg = (rocketProp.chuteSpillRadius*pi)/(2*rocketProp.chuteRadius);
rocketProp.chuteArea = pi*rocketProp.chuteRadius^2*cos(arg)^2;
rocketProp.chuteDrag = rocketProp.chuteArea*rocketProp.chuteCd;
if simuProp.count == 0  
    i = input('Input number of trajectories: ');
    simuProp.count = round(i);
end
%
% Rocket Drag
%
rocketProp.dragr    = DragProp{1};
rocketProp.drag     = DragProp{2};
%
% Atmospheric Properties
%
simuProp.tempr      = AtmProp{1};
simuProp.gravityr   = AtmProp{1};
simuProp.pressurer  = AtmProp{1};
simuProp.rhor       = AtmProp{1};
simuProp.viscosityr = AtmProp{1};
simuProp.temp       = AtmProp{2};
simuProp.gravity    = AtmProp{3};
simuProp.pressure   = AtmProp{4};
simuProp.rho        = AtmProp{5};
simuProp.viscosity  = AtmProp{6};
%
% Winds Aloft
%
simuProp.windsr = WindProp{1};
simuProp.windse = WindProp{2};
simuProp.windsn = WindProp{3};
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROCKET TRAJECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create zero initial conditions
%
zerosv = zeros(1, 7);
%
% Create trajectories array
%
traj = cell(1, simuProp.count);
%
% Calculate the characteristic trajectory with no FPA
%
%rocketProp.burnoutVelocity = 1000; %[m/s] hardcoded dummy variable
%rocketProp.burnoutAltitude = 13500; %[m] hardcoded dummy variable
%simuProp.burnoutAmbientTemp = 216.7; %[K] hardcoded dummy variable
%simuProp.burnoutAmbientDensity = .247205; %[kg/m^3] hardcoded dummy variable
fprintf('\nCalculating Characteristic Trajectory\n');
traj{1} = RocketTrajectory( zerosv, false );
rocketProp.burnoutVelocity = simuProp.burnOut(6); %[m/s] calculated from RocketTrajectory
rocketProp.burnoutAltitude = simuProp.burnOut(3); %[m] calculated from RocketTrajectory
simuProp.burnoutAmbientTemp = interp1(simuProp.tempr, simuProp.temp, rocketProp.burnoutAltitude); %[deg C] calculated using standard altitude values
simuProp.burnoutAmbientTemp = simuProp.burnoutAmbientTemp + 273.15; %convert from deg C to K
simuProp.burnoutAmbientDensity = interp1(simuProp.rhor, simuProp.rho, rocketProp.burnoutAltitude); %[kg/m^3] calculated based on altitude
%
% Draw trajectories actively
%
fig3 = figure(3);
itraj = traj{1};
time = itraj(:,7);
vel = itraj(:,6);
alt = itraj(:,3);
posy = itraj(:,2);
posx = itraj(:,1);
plot3(posx,posy,alt);
hold on;
xlabel( 'Easterly (m)' );
ylabel( 'Northerly (m)' );
zlabel( 'Altitude (m)' );
title( 'All Rocket Simulations' );
drawnow;
fig4 = figure(4);
plot(time,alt);
hold on;
drawnow;
%
% Calculate all other trajectories
%
for i = 2:simuProp.count
    fprintf('\nCalculating Trajectory #%d\n', i);
    traj{i} = RocketTrajectory( zerosv, true );
    itraj = traj{i};
    time = itraj(:,7);
    vel = itraj(:,6);
    alt = itraj(:,3);
    posy = itraj(:,2);
    posx = itraj(:,1);
    figure(3);
    plot3(posx,posy,alt);
    drawnow;
    figure(4);
    plot(time,alt);
    hold on;
    drawnow;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fig1 = figure(1);
grid on;
grid minor;
hold on;
itraj = traj{1};
time = itraj(:,7);
vel = itraj(:,6);
alt = itraj(:,3);
posy = itraj(:,2);
posx = itraj(:,1);
[hAx, hLine1, hLine2] = plotyy( time, vel, time, alt );
title( 'Characteristic Rocket Altitude and Velocity v. Time' );
xlabel( 'Time (s)');
ylabel(hAx(1), 'Velocity (m/s)');
ylabel(hAx(2), 'Altitude (m)');
set(hAx(1), 'Position',[0.14 0.18 0.72 0.72])

fig2 = figure(2);
for i = 1:size(traj,2)
    itraj = traj{i};
    time = itraj(:,7);
    vel = itraj(:,6);
    alt = itraj(:,3);
    posy = itraj(:,2);
    posx = itraj(:,1);
    plot3(posx,posy,alt);
    hold on;
end
xlabel( 'Easterly (m)' );
ylabel( 'Northerly (m)' );
zlabel( 'Altitude (m)' );
title( 'All Rocket Simulations' );
