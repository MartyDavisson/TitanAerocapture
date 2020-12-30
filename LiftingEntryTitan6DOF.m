function rk = LiftingEntryTitan6DOF(t, x, B, H_T, r_T, mu_T, rho0_T, LD, m, sig)
    %NOTE: this function assumes a negative gamma value (must be in radians)
    %Additionally assumes T=0,eps=0,omega=0,
    format long

    rk = [-((rho0_T * exp(-x(3) / H_T)) * x(1)^2) / (2 * B) - (mu_T / (r_T + x(3))^2) * sin(x(2));

        (x(1) * cos(x(2))) / (r_T + x(3)) + ((rho0_T * exp(-x(3) / H_T)) * x(1) * LD) / (2 * B) - ((mu_T / (r_T + x(3))^2) * cos(x(2))) / x(1);

        x(1) * sin(x(2));

        ((rho0_T * exp(-x(3) / H_T)) * x(1) * LD) / (2 * B) * (sin(sig) / cos(x(2))) - x(1) * cos(x(2)) * cos(x(4)) * tan(x(6)) / (r_T + x(3));

        (x(1) * cos(x(2)) * cos(x(4))) / ((r_T + x(3)) * cos(x(6)));

        (x(1) * cos(x(2)) * sin(x(4))) / (r_T + x(3))];

    %y(1) = Velocity
    %y(2) = Gamma
    %y(3) = Height
    %y(4) = Psi
    %y(5) = Longitude
    %y(6) = Latitude
end
