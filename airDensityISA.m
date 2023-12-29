function rho = airDensityISA(altitude)
    % Constants for the International Standard Atmosphere (ISA)
    T0 = 288.15;    % Standard sea-level temperature in Kelvin
    L = 0.0065;     % Temperature lapse rate in K/m
    P0 = 101325;    % Standard sea-level pressure in Pascal
    R = 287.05;     % Specific gas constant for dry air in J/(kg·K)
    g0 = 9.80665;   % Standard acceleration of gravity in m/s^2

    % Calculate temperature at altitude
    temperature = T0 - L * altitude;

    % Calculate pressure at altitude
    pressure = P0 * (1 - L * altitude / T0) ^ (g0 / (R * L));

    % Calculate air density
    rho = pressure / (R * temperature);
end
