function [ttpp, ttanc] = tevent (mu, sma, ecc, argper, tanom)

% trajectory event times

% input

%  mu     = gravitational constant (km**3/sec**2)
%  sma    = semimajor axis (kilometers)
%  ecc    = orbital eccentricity
%  argper = argument of perigee (radians)
%  tanom  = true anomaly (radians)

% output

%  ttpp  = time until perigee passage (seconds)
%  ttanc = time until ascending node crossing (seconds)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi2 = 2 * pi;

% Keplerian period

period = pi2 * sma * sqrt(sma / mu);

a = sin(tanom) * sqrt(1 - ecc * ecc);
b = ecc + cos(tanom);

% eccentric anomaly

eanom2 = atan3(a, b);

c = -eanom2 + ecc * sin(eanom2);

% time until perigee passage
   
ttpp = c * period / pi2;
   
if (ttpp < 0)
   ttpp = period + ttpp;
end

% time until ascending node crossing

d = pi2 - argper;

a = sin(d) * sqrt(1 - ecc * ecc);
b = ecc + cos(d);

eanom1 = atan3(a, b);

c = eanom1 - eanom2 - ecc * (sin(eanom1) - sin(eanom2));

ttanc = c * period / pi2;

if (ttanc < 0)
   ttanc = period + ttanc;
end
