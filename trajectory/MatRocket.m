%------------------------------------------------------
% TRITON ROCKET CLUB
% ROCKET ALTITUDE SIMULATION PROGRAM (MatRocket) v0.1
% (c)2015 TRITON ROCKET CLUB
%
% Authors: Trevor Irwin, Yen Nguyen, Aidan Setran
%
% This script will calculate the maximum one-directional
% altitude and velocity of a rocket using three methods:
% (0) Tsiolkovsky's Rocket Equation, which calculates the rocket's
%       maximum delta-v.
% (1) The simplified "Rocketmime.com" method with constant mass
% (2) A modified, numeric integration of the "Rocketmime"
%       equations which allows for variable mass. This method
%       should be more accurate than (1) but less than (3).
% (3) Uses Euler's method to double integrate acceleration,
%       accounting for variation in Standard Atmosphere and drag
%       with mach number.
%
% Resources:
% http://www.rocketmime.com/rockets/rckt_eqn.html
% http://exploration.grc.nasa.gov/education/rocket/drageq.html
%
% TO DO: 
%   -Look into RockSim for secondary analysis
%   -Update drag formula to account for mach drag
%------------------------------------------------------

% Clear the workspace of variables
clear all;
close all;
%----------FILE HANDLING----------
% Look for a local file, RocketProperties.txt
% If not found, ask for an input file of similar format.
% TODO: Allow for input of multiple files for comparison.

fprintf('\nTRITON ROCKET CLUB ALTITUDE SIMULATION (c)2015\n');

fprintf('Loading rocket properties...');

% Open the default rocket properties file
% Define error message
fileID = fopen('RocketProperties.txt');
errmsg = 'No RocketProperties.txt found in local folder. Please specify a new file.';

% Check that the default file opened properly
% While no valid file is open, ask for a new one
% If the file we asked for is invalid, try appending .txt
while fileID < 0
    disp(errmsg);
    filename = input('Open file: ', 's');
    [fileID, errmsg] = fopen(filename);
    if fileID == -1
        filename = sprintf('%s.txt', filename);
        [fileID, errmsg] = fopen(filename);
    end
end

fprintf('done!\n');

% Scan the properties text file for contents, string = value
% Use comments styled like Matlab, %
RProp = textscan(fileID, '%s = %f', 'CommentStyle', '%');
fclose(fileID);

% Read the properties text file contents
% Set tburn to false. If burn time is specified in the
% text file, tburn will be overwritten with a value.
% If it remains false, it will be calculated later.
[m, n] = size(RProp{1});
tburn = false;
for i=1:m
    switch RProp{1}{i}
        
        % --Rocket Body Properties--
        case 'drymass'
            dryM = RProp{2}(i); %[kg] mass of structure+propellant
        case 'wetmass'
            wetM = RProp{2}(i); %[kg] mass of structure+propellant
        case 'radius'
            r = RProp{2}(i); %[in]radius of body
        case 'rho'
            rho = RProp{2}(i); %[kg/m^3] air density
        case 'drag'
            Cd = RProp{2}(i); %drag coefficient
        case 'gravity'
            g = RProp{2}(i); %[m/s^2] gravitational acceleration
        
        % --Motor Properties--
        case 'thrust'
            T = RProp{2}(i); %[N] thrust
        case 'Isp'
            Isp = RProp{2}(i); %[s] specific impulse
        case 'Itot'
            Itot = RProp{2}(i); %[N-s] total impulse
        case 'burn'
            tburn = RProp{2}(i); %[s]An option to override the burn time
            
        % --Simulation Properties--
        case 'res'
            res = RProp{2}(i); %[N] thrust  
   
    end
end

%----------VARIABLE PRE-CALCULATION----------
% Pre-calculate some commonly used terms for readibility
% These include the cross-sectional area, wind resistance
% factor used in the RocketMime equations, and the motor
% burn time.

% --Rocket Body Properties--
A = pi*r^2; %[in^2] area of rocket cross-section
k = 0.5*rho*Cd*A; %Wind Resistance factor

% --Motor Properties--
if ~tburn
    tburn = Itot/T; %[s] Burn Time
end

%----------ANONYMOUS FUNCTIONS----------
% Anonymous functions to simplify math and
% perform unit conversions.

% --Simplifying factors--
getz = @(T, M, g, k) sqrt((T-M*g)/k);
getx = @(k, z, M) 2*k*z/M;

% --Unit Conversions--
m2ft = @(m) 3.28084*m;
pa2psi = @(p) 0.000145038*p;


%% ----------Tsiolkovsky's Equation----------
% Assumes an ideal rocket. This section is included for
% the sake of interest and not meant as part of the simulation.

fprintf('[0] Calculating Tsiolkovsky rocket equation:\n');

% Compute delta-v
delv = Isp*g*log(wetM/dryM);

fprintf('done!\n');
fprintf('Delta-V: %2.2fm/s\n', delv);

%% ----------Rocketmime Calculations----------
% ----------ASSUMPTIONS----------
% +G-force constant, neglect change in gravity due to altitude
% +Mass constant, neglect change due to burn
% +Drag depends only on velocity, neglect air density changes
% +Thrust is constant throughout boost phase

fprintf('[1] Calculating Rocketmime equations:\n');
fprintf('Calculating...');

% Calculate all simplifying factors
z = getz(T, wetM, g, k);
x = getx(k, z, wetM);

% Calculate velocity at burnout
v = z*(1-exp(-x*tburn))/(1+exp(-x*tburn));

% Calculate distance travelled by end of burnout
yb = (-wetM/(2*k))*log((T-wetM*g-k*v^2)/(T-wetM*g));

% Calculate distance coasted
yc = (dryM/(2*k))*log((dryM*g+k*v^2)/(dryM*g));

% Combine boost and coast phases for total maximum altitude
rocketmime_max = yb+yc;

% Print to console
fprintf('done!\n');
fprintf('Rocketmime boost altitude: %2.2fm = %2.2fft\n', yb, m2ft(yb));
fprintf('Rocketmime coast distance: %2.2fm = %2.2fft\n', yc, m2ft(yc));
fprintf('Rocketmime maximum altitude: %2.2fm = %2.2fft\n', rocketmime_max, m2ft(rocketmime_max));

%% ----------Integrated Rocketmime Calculations----------
% ----------ASSUMPTIONS----------
% +G-force constant, neglect change in gravity due to altitude
% +Mass changes linearly with burn
% +Drag depends only on velocity, neglect air density changes
% +Thrust is constant throughout boost phase

fprintf('[2] Calculating integrated Rocketmime equations:\n');
fprintf('Calculating boost phase...');

% Set up variables for integration. Pre-allocating saves on runtime.
% Use a length defined by the resolution for each
time = linspace(0, tburn, res); %[s] Time between 0 and burnout, with specified resolution
mass = linspace(wetM, dryM, res); %[kg] Mass between wet and dry, with specified resolution
vel  = zeros(1, res); %[m/s] Velocity. 
alt  = zeros(1, res); %[m] Altitude.
t = tburn/res;

% Step through to the maximum resolution
% In each step, calculate a velocity differential and 
% add it to the previous velocity.
% Also in each step, use the previous velocity times the timestep
% as an altitude differential and add it to the altitude.

for i=1:res-1
    z = getz(T, mass(i), g, k);
    x = getx(k, z, mass(i));
    %TODO: Find the acceleration and integrate that instead?
    vel(i+1) = vel(i)+z*(1-exp(-x*t))/(1+exp(-x*t));
    alt(i+1) = alt(i)+vel(i)*t; 
end

% Define the height of the boost phase
rocketmime_int_boost = alt(res);
fprintf('done!\n');

fprintf('Calculating coast phase...');

% Calcualte the height of the coast phase
rocketmime_int_coast = (dryM/(2*k))*log((dryM*g+k*v^2)/(dryM*g)); %Altitude at end of coast
fprintf('done!\n');

% Define the total maximum altitude by combining boost and coast
rocketmime_int_max = rocketmime_int_boost+rocketmime_int_coast;

% Print to console
fprintf('Integrated Rocketmime boost altitude: %2.2fm = %2.2fft\n', rocketmime_int_boost, m2ft(rocketmime_int_boost));
fprintf('Integrated Rocketmime coast distance: %2.2fm = %2.2fft\n', rocketmime_int_coast, m2ft(rocketmime_int_coast));
fprintf('Integrated Rocketmime maximum altitude: %2.2fm = %2.2fft\n', rocketmime_int_max, m2ft(rocketmime_int_max));


%% ----------Euler's Method----------
% ----------ASSUMPTIONS----------
% +Changes in gravity, air density, temperature, cD will be modelled
% +Mass changes linearly with burn
% +Thrust is constant throughout boost phase
% +No other aerodynamic forces (other than drag) act on the body
% +Angle of attack is zero, motion is one-dimensional

fprintf('[3] Calculating Eulers method with resolution %d:\n', res)


%----------FILE HANDLING----------
% Look for local files, StandardAtmosphere.txt and RocketDrag.txt
% If not found, ask for input files of similar format.
fprintf('Loading atmospheric data...');

% Open the standard atmosphere file
% Define an error message
fileID = fopen('StandardAtmosphere.txt');
errmsg = 'No StandardAtmosphere.txt found in local folder. Please specify a new file.';

% Check that the default file opened properly
% While no valid file is open, ask for a new one
% If the file we asked for is invalid, try appending .txt
while fileID < 0
    disp(errmsg);
    filename = input('Open file: ', 's');
    [fileID, errmsg] = fopen(filename);
    if fileID == -1
        filename = sprintf('%s.txt', filename);
        [fileID, errmsg] = fopen(filename);
    end     
end

% Scan the properties text file for contents, where each
% column is a different atmospheric property.
% Use comments styled like Matlab, %
AtmProp = textscan(fileID, '%f\t%f\t%f\t%f\t%f\t%f', 'CommentStyle', '%');
fclose(fileID);

% Define vectors for temperature, gravity, pressure
% density, and viscosity.
ntemp = [AtmProp{1}, AtmProp{2}];
ng    = [AtmProp{1}, AtmProp{3}];
npres = [AtmProp{1}, AtmProp{4}];
nrho  = [AtmProp{1}, AtmProp{5}];
nvisc = [AtmProp{1}, AtmProp{6}];

fprintf('done!\n');

fprintf('Loading rocket drag data...');

% Open the rocket drag file
% Define an error message
fileID = fopen('RocketDrag.txt');
errmsg = 'No RocketDrag.txt found in local folder. Please specify a new file.';

% Check that the default file opened properly
% While no valid file is open, ask for a new one
% If the file we asked for is invalid, try appending .txt
while fileID < 0
    disp(errmsg);
    filename = input('Open file: ', 's');
    [fileID, errmsg] = fopen(filename);
    if fileID == -1
        filename = sprintf('%s.txt', filename);
        [fileID, errmsg] = fopen(filename);
    end
end

% Scan the properties text file for contents, where data is
% mach-number [tab] drag-coefficient
% Use comments styled like Matlab, %
DragProp = textscan(fileID, '%f\t%f', 'CommentStyle', '%');
fclose(fileID);

% Define a vector containing mach numbers and drag coefficients
ndrag = [DragProp{1}, DragProp{2}];

fprintf('done!\n');

% --Boost Phase--

% Ask if we want an altitude limit.
% If an altitude limit is defined, we will later iterate
% increasing amounts of dry mass to limit the altitude.
allow_alt = input('Enter altitude limit (or 0 for no limit): ');
if allow_alt == 0
    allow_alt = inf;
end

% While the current highest altitude is above the
% altitude limit, repeat all calculations until it is below
euler_max = allow_alt+1;
while euler_max >= allow_alt

    fprintf('Calculating boost phase...');
    
    % Record the current time. The difference with a later
    % current time will tell us the time taken to complete the task.
    tasktime=cputime;
    
    % Define a number of vectors. Pre-allocate them with
    % a number of indices corresponding to the resolution.
    ntime = linspace(0, tburn, res); %[s] Time between 0 and burnout, with specified resolution
    t     = tburn/res; %differential dt
    nmass = linspace(wetM, dryM, res); %[kg] Mass linearly interpolated between wet and dry
    nacc  = zeros(1, res); %[m/s^2] Acceleration. Zeroes used for memory pre-allocaiton (makes the script run faster ;) )
    nvel  = zeros(1, res); %[m/s] Velocity. 
    nalt  = zeros(1, res); %[m] Altitude.
    nfD  = zeros(1, res); %[m] Drag Force.
    nrhoc = zeros(1, res); %[kg/m^3] Calculated air density.
    nmach = zeros(1, res); % Mach number

    % For every index in the boost phase (1->the resolution)
    % Calculate an acceleration differential
    % To do so, determine mass, thrust, etc.
    % Interpolate any values which are not already defined
    % Double integrate with Euler's method to yield position
    for i=1:res-1
        imass = nmass(i); %Mass at time of integration
        fT = T; %[N] Force due to thrust
        aG = interp1(ng(:,1), ng(:,2), nalt(i)); %[m/s^2] Precise acceleration due to gravity
        fG = aG*imass;
        irho = interp1(nrho(:,1), nrho(:,2), nalt(i)); %[kg/m^3] Precise air density
        nrhoc(i) = irho;
        ivel = nvel(i); %Velocity of rocket at this point in time
        itemp = interp1(ntemp(:,1), ntemp(:,2), nalt(i)); %Temperature at this altitude
        itemp = itemp+273.15; %Convert *C to K
        ic = 20.05*sqrt(itemp);     %Speed of sound at this temperature
        imach = ivel/ic; %Mach number
        nmach(i+1) = imach;
        if imach <= ndrag(end, 1)
            iCd = interp1(ndrag(:,1), ndrag(:,2), imach); %Drag coefficient at mach number
        else
            iCd = ndrag(end, 2);
        end
        fD = iCd*irho*A/2*ivel^2; %Force due to drag
        fSum = fT - fG - fD; %Sum of forces
        aSum = fSum/imass;    %Acceleration, from newton's first law

        nfD(i+1)  = fD;
        nacc(i+1) = aSum; %Record acceleration
        nvel(i+1) = nvel(i) + aSum*t; %Numeric integration of velocity
        nalt(i+1) = nalt(i) + nvel(i+1)*t + 0.5*aSum*t^2; %Numeric integration of position
    end

    euler_boost = nalt(end);
    
    % Find the difference between the recorded time and current time
    tasktime = cputime-tasktime;
    % Output the time taken to complete the task
    fprintf('done! Took %2.2f seconds.\n', tasktime);

    % --Coast Phase--
    fprintf('Calculating coast phase...');
    
    % Record the current time. The difference with a later
    % current time will tell us the time taken to complete the task.
    tasktime = cputime;
    
    % While the rocket's altitude is above zero,
    % calculate an acceleration differential for the coast phase.
    % To do so, determine mass, thrust, etc. interpolate any values 
    % which are not already defined.
    % Double integrate with Euler's method to yield position
    % If the rocket has not fallen back to zero altitude yet,
    % pre-allocate memory by extending vectors by resolution.
    
    % Define the current index for vectors
    i = res;
    % Define the current "rollover". This is used to determine
    % how far to extend vectors when pre-allocating memory.
    iroll = 2;
    
    while nalt(end) >= 0
        %dryM
        if i+1 > length(ntime)
            ntime = [ntime, linspace(ntime(end)+t, ntime(res)*iroll, res)];
            nacc  = [nacc, zeros(1, res)];
            nvel  = [nvel, zeros(1, res)];
            nalt  = [nalt, zeros(1, res)];
            nfD   = [nfD, zeros(1, res)];
            nrhoc  = [nrhoc, zeros(1, res)];
            nmach = [nmach, zeros(1, res)];
            iroll = iroll+1;
        end

        aG = interp1(ng(:,1), ng(:,2), nalt(i)); %[m/s^2] Precise acceleration due to gravity
        fG = aG*dryM; %Force due to gravity, rocket now dry
        irho = interp1(nrho(:,1), nrho(:,2), nalt(i)); %[kg/m^3] Precise air density
        nrhoc(i) = irho;
        ivel = nvel(i); %Velocity of rocket at this point in time
        itemp = interp1(ntemp(:,1), ntemp(:,2), nalt(i)); %Temperature at this altitude
        itemp = itemp+273.15; %Convert *C to K
        ic = 20.05*sqrt(itemp);     %Speed of sound at this temperature
        imach = ivel/ic; %Mach number
        nmach(i+1) = imach;
        if abs(imach) < ndrag(end, 1);
            iCd = interp1(ndrag(:,1), ndrag(:,2), abs(imach)); %Drag coefficient at mach number
        elseif abs(imach) >= ndrag(end, 1)
            iCd = ndrag(end, 2);
        else
            iCd = ndrag(1, 2);
        end

        if imach < 0
            iCd = -1*iCd; %negate drag value if velocity is negative
        end

        fD = iCd*irho*ivel^2*A/2; %Force due to drag

        fSum = - fG - fD;
        aSum = fSum/dryM;

        nfD(i+1)  = fD;
        nacc(i+1) = aSum; %Record acceleration
        nvel(i+1) = nvel(i) + aSum*t; %Numeric integration of velocity
        nalt(i+1) = nalt(i) + nvel(i+1)*t + 0.5*aSum*t^2; %Numeric integration of position
        %TODO: Should acceleration term be included?

        i = i+1; %iterate looping variable
    end

    [euler_max, max_i] = max(nalt);
    euler_coast = euler_max-euler_boost;
    
    % Find the difference between the recorded time and current time
    tasktime = cputime-tasktime;
    % Output the time taken to complete the task
    fprintf('done! Took %2.2f seconds.\n', tasktime);

    % --Max Q--
    fprintf('Calculating Max Q...');
    tasktime = cputime;

    nq = 0.5.*nrhoc.*(nvel.^2); %[N/m^2] Dynamic Pressure
    [euler_maxq, max_nqi] = max(nq);
    euler_maxq_alt = nalt(max_nqi);
    [euler_maxv, max_nvi] = max(nvel);
    euler_maxv_alt = nalt(max_nvi);
    [euler_maxa, max_nai] = max(nacc);
    euler_maxa_alt = nalt(max_nai);
    [euler_maxm, max_nmi] = max(nmach);
    euler_maxm_alt = nalt(max_nmi);
    [euler_maxd, max_ndi] = max(nfD);
    euler_maxd_alt = nalt(max_ndi);

    tasktime = cputime-tasktime;
    fprintf('done! Took %2.2f seconds.\n', tasktime);

    fprintf('Euler boost altitude: %2.2fm = %2.2fft\n', euler_boost, m2ft(euler_boost));
    fprintf('Euler coast distance: %2.2fm = %2.2fft\n', euler_coast, m2ft(euler_coast));
    fprintf('Euler max altitude: %2.2fm = %2.2fft\n', euler_max, m2ft(euler_max));
    fprintf('Euler max Q: %2.2fN/m^2 = %2.2fPsi at %2.2fm = %2.2fft\n', euler_maxq, pa2psi(euler_maxq), euler_maxq_alt, m2ft(euler_maxq_alt));
    fprintf('Euler max velocity: %2.2fm/s at %2.2fm = %2.2fft\n', euler_maxv, euler_maxv_alt, m2ft(euler_maxv_alt));
    fprintf('Euler max acceleration: %2.2fm/s^2 at %2.2fm = %2.2fft\n', euler_maxa, euler_maxa_alt, m2ft(euler_maxa_alt));
    fprintf('Euler max mach: %2.2f at %2.2fm = %2.2fft\n', euler_maxm, euler_maxm_alt, m2ft(euler_maxm_alt));
    fprintf('Euler max drag: %2.2fN at %2.2fm = %2.2fft\n', euler_maxd, euler_maxd_alt, m2ft(euler_maxd_alt));
    
    dryM = dryM+5;
    wetM = wetM+5;

end

%TODO: Change 5 to a variable read from RocketProperties file!
fprintf('Optimal dry mass: %2.2fkg wet mass: %2.2fkg\n', dryM-5, wetM-5);

%% ----------Plots----------

fig1 = figure(1);
grid on;
grid minor;
hold on;
[hAx, hLine1, hLine2] = plotyy( ntime, nvel, ntime, nalt);
title( 'Rocket Altitude and Velocity v. Time' );
xlabel( 'Time (s)');
ylabel(hAx(1), 'Velocity (m/s)');
ylabel(hAx(2), 'Altitude (m)');
set(hAx(1), 'Position',[0.14 0.18 0.72 0.72])
%xmarkers = [ntime(max_nvi), ntime(max_nai)];
%ymarkers = [nvel(max_nvi), nalt(max_nai)];
line(ntime(max_nvi), nvel(max_nvi), 'linestyle', 'none', 'Marker', '*', 'Parent', hAx(1)); 
line(ntime(max_i), nalt(max_i), 'linestyle', 'none', 'Marker', '*', 'Parent', hAx(2)); 
vel_label=sprintf('Maximum Velocity = %2.2fm/s', euler_maxv);
alt_label=sprintf('Maximum Altitude = %2.2fm', euler_max);
text(ntime(max_nvi),nvel(max_nvi), vel_label, 'color', 'b', 'Parent', hAx(1));
text(ntime(max_i), nalt(max_i), alt_label, 'color', 'r', 'Parent', hAx(2));

fig2 = figure(2);
grid on;
grid minor;
hold on;
plot( nalt, nq )
title( 'Dynamic Pressure v. Altitude' )
xlabel( 'Altitude (m)' )
ylabel( 'Dynamic Pressure (kg*m/s^2)' )

% Open the SimResults file and be prepared to overwrite
fileID = fopen('SimResults.txt', 'w');
n = length(ntime);
step = round(n/100);
fprintf(fileID, 'Time (s)\tAlt (m)\tVel (m/s)\tQ (kg*m/s^2)\n');
for i = 1:step:n
    wtime = ntime(i);
    walt = nalt(i);
    wvel = nvel(i);
    wq = nq(i);
    fprintf(fileID, '%2.2f\t%2.2f\t%2.2f\t%2.2f\n', wtime, walt, wvel, wq);
end
fclose(fileID);
    

