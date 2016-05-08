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
        case 'rho'
            simuProp.rho = RProp{2}(i); %[kg/m^3] air density
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
        
        % --Chute Properties--
        case 'chutedrag'
            rocketProp.chuteDrag = RProp{2}(i); %drag coefficient
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
    end
end
%
% Convenience Properties
%
rocketProp.area = pi*rocketProp.radius^2;
if rocketProp.burnTime == 0
    rocketProp.burnTime = rocketProp.iTotal/rocketProp.thrust;
end
rocketProp.chuteArea = pi*(rocketProp.chuteRadius-rocketProp.chuteSpillRadius)^2;
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
fprintf('\nCalculating Characteristic Trajectory\n');
traj{1} = RocketTrajectory( zerosv, false );
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
