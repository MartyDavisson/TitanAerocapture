%Titan Aerocapture Optimization Final Project: ME6103
%Martin Davisson, M.S. Aerospace Engineering, Georgia Institute of
%Technology
close all

%Step 0: Problem Constants and Parameters

%Titan's Gravitational Parameters(1):
G = 6.67408e-11; %Gravitational Constant [m^3/kg*s^2]
m_T = 1.3452e23; %Titan's Mass [kg]
mu_T = G * m_T; %Titan's Standard Gravitational Constant mu=GM [m^3/s^2]
r_T = 2574.7e3; %Titan's Equatorial Radius [m]
m_S = 5.683e26; %Saturn�s mass
r_ST = 1.22185e9; %Distance between Saturn and Titan [m](Mean Orbital Radius)
rp = 350e3; %Assume an altitude at periapsis of 350km (approx. necessary altitude for aerocapture) (2)
Vesc = sqrt(mu_T / r_T);
rSOI = r_ST * (m_T / m_S)^(2/5); %ref ch. 8 fund. astro

%Titan's Atmospheric Parameters(2):
H_T = 48.38e3; %Scale Height [km] based off Cassini Titan Flyby logarithmic approx. (full data available in ref 2
rho0_T = 1.006e-1; %Gravity at Titan's Surface [kg/m^3]
h0_T = 1000e3; %Atmospheric Interface Altitude
ks = 1.74153e-4; %Convective Heat Transfer coefficient for Titan's Atmosphere (assume same as earth for now)

%Vehicle Constants:
m = 420; %vehicle mass
LD = 0.25; %L/D Ratio
B = 90; %Ballistic Coefficient
rn = 0.91; %Effective Nose Radius
Isp = 312; %[s] vacuum Isp based off R-4D rocket engine
g0 = 9.80665; %[m/s^2] Standard Gravity
ve = Isp * g0; %[m/s]effective exhaust velocity

%Step 1: Define Optimization Problem

%Search Parameter Boundary Constraints:
Vmax = 12e3;
Vmin = 5.5e3;
gmax = -36;
gmin = -39.5;
smax = 90;
smin = 0;

%Exhaustive Search Domain:
Vatm_i = linspace(Vmin, Vmax, 25);
gamma0_i = linspace(gmin, gmax, 25);
sig_i = 15;

i = numel(Vatm_i);
j = numel(gamma0_i);
k = numel(sig_i);
num_iter = i * j * k;
%Solution Indexing:
valid = [];
invalid = [];
deviations = [];
n = 0; %Loop counter used for indexing
v = 0; %valid solution counter
iv = 0; %invalid solution counter
OptSol = NaN(3, 1); %Initialize optimum solutions
Opt_xyz = NaN(3, 1);
Opt_atmospheric = NaN(6, 1);
Opt_time = NaN;
Opt_qs = NaN;
Opt_hl = NaN;
% Target Values:
T_Vinf = 5000;

%Step 2: Begin Exhaustive Search
for Vatm = Vatm_i

    for gamma0 = gamma0_i

        for sig = sig_i

            % Step 2: Atmospheric Entry
            % Initial Conditions:
            gamma0 = deg2rad(gamma0); %EOM assume negative flight path angle [rad]
            sig = deg2rad(sig); %bank angle

            timerange = 0:10:6000;
            IC = [Vatm; gamma0; h0_T; 0; 0; 0];

            %Numerical Integration
            options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
            [t, x] = ode45(@(t, x) LiftingEntryTitan6DOF(t, x, B, H_T, r_T, mu_T, rho0_T, LD, m, sig), timerange, IC, options);

            %Identify Exit Conditions
            [val, exit] = min(abs(x(5:end, 3) - h0_T));
            ICexit = x(exit, :);

            %determine apogee radius to avoid integrating if invalid
            %solution
            r_vec = ICexit(1:3);
            v_vec = ICexit(4:6);
            angularmom = cross(r_vec, v_vec); %angular momentum of elliptic transfer orbit
            e_hat = (1 / mu_T) * ((norm(v_vec)^2 -mu_T / norm(r_vec)) .* r_vec - (dot(r_vec, v_vec)) .* v_vec);
            ecc = norm(e_hat);
            Energy = 0.5 * norm(v_vec)^2 - mu_T / norm(r_vec);
            semi_maj = -mu_T / Energy;
            r_apogee = semi_maj * (1 + ecc);

            t;
            x = x(1:exit, :);
            time = timerange(1:exit);
            v_atm_exit = x(exit, 1);
            dv_atm = Vatm - v_atm_exit;
            TOF_atm = timerange(exit); %Atmospheric Portion TOF
            atmospheric = x;
            %Heating Calculations:
            rho = rho0_T .* exp(-x(:, 3) ./ H_T); %recalculate the density from the height in numerical integration scheme
            qs = (ks .* sqrt(rho / rn) .* x(:, 1).^3) ./ 100^2; %Stagnation Point heating [w/cm^2]

            hl = [0]; %heatload

            for i = 2:size(qs, 1)
                hl = [hl; hl(i - 1) + trapz(time(i - 1:i), qs(i - 1:i))];
            end

            %Coordinate Transformation from Latitude/Longitude to Cartesian
            x1 = (r_T + x(:, 3)) .* cos(x(:, 6)) .* cos(x(:, 5));
            y1 = (r_T + x(:, 3)) .* cos(x(:, 6)) .* sin(x(:, 5));
            z1 = (r_T + x(:, 3)) .* sin(x(:, 6));
            atmospheric_xyz = [x1, y1, z1];

            r1 = [x1, y1, z1]; %radius vector array
            v1 = diff(r1); %velocity vector array

            %Step 3: Elliptic Exit Orbit
            %Identify Exit Conditions
            exit;

            if exit > 1 && r_apogee < rSOI
                r_i = r1(exit, :);
                v_i = x(exit, 1) * (v1(exit - 1, :) / norm(v1(exit - 1, :)));

                %Propagate Orbit
                t = 0:1:4 * 60000; %set time range for integration
                options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
                [t, r] = ode45(@orbdynTitan, t, [r_i v_i], options);

                r_mag = sqrt(r(:, 1).^2 + r(:, 2).^2 + r(:, 3).^2); %compute magnitude of radius [m]
                v_mag = sqrt(r(:, 4).^2 + r(:, 5).^2 + r(:, 6).^2); %compute magnitude of velocity [m/s]

                j = 1;
                hp = diff(r_mag);
                %Determine exit condition
                while ((hp(j) >= 0) && (j + 1 < size(r_mag, 1)))
                    j = j + 1;
                end

                exit2 = j;

                %Limit to range of interest
                r_e = r(1:exit2, :); %range of radii for plotting elliptic orbit
                TOF_e = t(exit); %TOF elliptic portion of Flight [s]

                rcs = r(exit2, [1, 2, 3]);

                vcs_hat = (r(exit2, [4, 5, 6]) - r(exit2, [4, 5, 6])) ./ (norm(r(exit2, [4, 5, 6]) - r(exit2, [4, 5, 6])));

                elliptic_xyz = [r_e(:, 1), r_e(:, 2), r_e(:, 3)];

                %Step 4: Circularization of Orbit
                r_cs = [r(exit2, 1), r(exit2, 2), r(exit2, 3)];
                v_cs = sqrt(mu_T / norm(r_cs));
                v_i = v_cs .* [r(exit2, 4), r(exit2, 5), r(exit2, 6)] ./ norm([r(exit2, 4), r(exit2, 5), r(exit2, 6)]); %

                t = 0:10:4 * 120000;
                options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
                [~, r] = ode45(@orbdynTitan, t, [r_cs v_i], options);
                r_c = r;
                circular_xyz = [r_c(:, 1), r_c(:, 2), r_c(:, 3)];

                TOF_c = 2 * pi * sqrt(r_T^3 / mu_T); %Circular Orbital Period [s]
                dv1 = v_cs - v_mag(exit2);

                %Hyperbolic entry
                %Initial Conditions
                t = 0:10:4 * 600000;
                r_i = [h0_T + r_T, 0, 0]; % program is integrated backwards from atmospheric entry to determin hyperbolic orbit
                v_i = [-sin(gamma0) * Vatm, -cos(gamma0) * Vatm, 0];

                options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
                [t_nom, r] = ode45(@orbdynTitan, t, [r_i v_i], options);

                r_mag = sqrt(r(:, 1).^2 + r(:, 2).^2 + r(:, 3).^2);

                [Vinf, enter1] = min(abs(r_mag - rSOI));
                %Vinf=2.982268906511366e+03; %What is this value doing? if anything?
                TOF_h = t(enter1);
                r_h = r(1:enter1, :);
                hyperbolic_xyz = [r_h(:, 1), r_h(:, 2), r_h(:, 3)];

                %Step 5: Trajectory Comparison Analysis
                %Direct Insertion
                E = Vinf^2/2; %Hyperbolic Orbit Energy
                Vp = sqrt((E + mu_T / norm(rcs)) * 2); %Velocity at Periapsis for Direct Insertion
                dv2 = Vp - v_cs;

                %Mass of Propellant
                mr1 = exp(dv1 / ve); %Mass ratio m0/mf from rocket equation
                mr2 = exp(dv2 / ve);

                mp1 = m * mr1 - m; %Mass of Propellant [kg]
                mp2 = m * mr2 - m;

                %Consolodating all coordinates of solution
                xyz = [hyperbolic_xyz; atmospheric_xyz; elliptic_xyz; circular_xyz];

                %Evaluating Constraints:
                if norm(r_cs) >= (r_T + h0_T) && norm(r_cs) <= rSOI && max(qs) <= 120 && max(hl) <= 10000 && min(sqrt(sum(circular_xyz.^2, 2))) >= (r_T + h0_T)%sigma>0 && gamma0>=Dmin && Vatm>0 && F>=Fmin && F<=Fmax && S<=Sy && W<1.5 && C>=4 && C<=25
                    v = v + 1;

                    %Deviations:
                    if T_Vinf - Vinf < 0%Goal 1 is to Maximize Arrival Velocity
                        d1p = 0;
                        d1m = abs(T_Vinf - Vinf) / T_Vinf;
                    elseif T_Vinf - Vinf > 0
                        d1p = abs(T_Vinf - Vinf) / T_Vinf;
                        d2m = 0;
                    else %Target Reached
                        d1p = 0;
                        d1m = 0;
                    end

                    deviations(1:2, v) = [d1p; d1m];

                    %Objective Function:
                    %Weighting Scheme 1: Equal weighting
                    w1 = 4/5; %Propellant Weighting
                    w2 = 1/5; %Arrival Velocity Target Weighting
                    f = w1 * mp1 / 100 + w2 * d1m;

                    valid(1:5, v) = [Vatm; gamma0; sig; mp1; f]; %index valid parameter values if meet constraints

                    if sum(isnan(OptSol)) > 0%If no optimum solution exists yet
                        OptSol = [Vatm; gamma0; sig; mp1; f];
                        Opt_xyz = xyz;
                        Opt_atmospheric = atmospheric;
                        Opt_time = time;
                        Opt_qs = qs;
                        Opt_hl = hl;
                        fBest = f;
                    elseif f < fBest%minimize fuel consumption
                        OptSol = [Vatm; gamma0; sig; mp1; f];
                        Opt_xyz = xyz;
                        Opt_atmospheric = atmospheric;
                        Opt_time = time;
                        Opt_qs = qs;
                        Opt_hl = hl;
                        fBest = f;
                    end

                else
                    iv = iv + 1;
                    invalid(:, iv) = [Vatm; gamma0; sig]; %invalid parameters
                end

            else %Case where integration fails due to impact with Titan's surface
                iv = iv + 1;
                invalid(:, iv) = [Vatm; gamma0; sig];
            end

            n = n + 1;
            percent_complete = round(100 * (n / num_iter));
            fprintf('percent complete: %d\n', percent_complete)
        end

    end

end

if sum(isnan(OptSol)) == 0
    %Plot of Optimum Solution Trajectory:
    figure(1)
    title('Optimal Trajectory: Fuel Consumption')
    npanels = 30;
    npanels2 = 30;
    r_atm = r_T + h0_T; % atmosphere radius (meters)
    axis(4.35e7 * [-1 1 -1 1 -1 1]);
    hold on;
    axis vis3d;
    [xx, yy, zz] = ellipsoid(0, 0, 0, r_T, r_T, r_T, npanels);
    globe = surf(xx, yy, -zz, 'FaceColor', 'g', 'EdgeColor', 'w');
    [xx, yy, zz] = ellipsoid(0, 0, 0, r_atm, r_atm, r_atm, npanels2);
    globe = surf(xx, yy, -zz, 'FaceColor', 'y', 'EdgeColor', 0.75 * [1 1 1]);
    alpha 0.35
    scatter3(Opt_xyz(:, 1), Opt_xyz(:, 2), Opt_xyz(:, 3), '.r')

    figure(2)
    plot(Opt_atmospheric(:, 1), Opt_atmospheric(:, 3), '-r');
    title('Atmospheric Entry')
    xlabel('Velocity [m/s]')
    ylabel('Altitude [m]')
    grid on
    hold on

    figure(3)
    plot(rad2deg(Opt_atmospheric(:, 2)), Opt_atmospheric(:, 3), '-r');
    title('Flight Path Angle vs Altitude')
    xlabel('Flight Path Angle [deg]')
    ylabel('Altitude [m]')
    grid on
    hold on

    figure(4)
    plot(Opt_time, rad2deg(Opt_atmospheric(:, 4)), '-r');
    title('Psi')
    xlabel('Time [s]')
    ylabel('Turning Angle [deg]')
    grid on
    hold on

    figure(5)
    plot(Opt_time, rad2deg(Opt_atmospheric(:, 5)), '-r');
    title('Longitude vs Time')
    xlabel('Time [s]')
    ylabel('Longitude [deg]')
    grid on
    hold on

    figure(6)
    plot(Opt_time, rad2deg(Opt_atmospheric(:, 6)), '-r');
    title('Latitude vs Time')
    xlabel('Time [s]')
    ylabel('Latitude [deg]')
    grid on
    hold on

    figure(7)
    plot(Opt_time, Opt_qs, '-r');
    title('Stagnation Point Heating vs Time')
    xlabel('Time [s]')
    ylabel('Heat Rate [W/cm^2]')
    grid on
    hold on

    figure(8)
    plot(Opt_time, Opt_hl, '-r');
    title('Integrated Heat Load vs Time')
    xlabel('Time [s]')
    ylabel('Heat Load [J/cm^2]')
    grid on
    hold on

    figure (9)
    hold on
    scatter(valid(1, :), rad2deg(valid(2, :)), 25, 'r', 'filled')
    scatter(OptSol(1), rad2deg(OptSol(2)), 45, 'b')
    xlabel('Atmospheric Entry Velocity, Vatm [m/s]')
    ylabel('Flight Path Angle, Gamma [rad]')
    title('Valid Range of Solutions: Vatm vs. Gamma')
    hold off

    figure (10)
    hold on
    scatter3(valid(1, :), rad2deg(valid(2, :)), valid(5, :))
    scatter3(OptSol(1), rad2deg(OptSol(2)), fBest, 45, 'b')
    xlabel('Atmospheric Entry Velocity, m/s')
    ylabel('Flight Path Angle, deg')
    title('Valid Range of Solutions: Vatm vs. Gamma')
    hold off

    %Interpolated surface plot
    x = valid(1, :);
    y = rad2deg(valid(2, :));
    z = valid(5, :);
    xv = linspace(min(x), max(x), 100);
    yv = linspace(min(y), max(y), 100);
    [X, Y] = meshgrid(xv, yv);
    Z = griddata(x, y, z, X, Y);
    figure(11)
    hold on
    title('Objective Function Surface Plot, Entry Velocity vs. Flight Path Angle')
    ylabel('Flight Path Angle, deg.')
    xlabel('Atmospheric Entry Velocity, m/s')
    zlabel('Mass of Propellant, kg')
    surf(X, Y, Z); %surface plot
    contour3(X, Y, Z, 'r')%contour lines
    scatter3(OptSol(1), rad2deg(OptSol(2)), fBest, 45, 'b')%optimal solution
    grid on
    set(gca, 'ZLim', [0 3])
    shading interp
else
end
