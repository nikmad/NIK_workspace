function [tnodal, tanomal] = tperiod (sma, ecc, inc, argper)

% orbital periods

% input

%  sma    = semimajor axis
%  ecc    = orbital eccentricity
%  inc    = orbital inclination (radians)
%  argper = argument of perigee (radians)

% output

%  tnodal  = nodal period (seconds)
%  tanomal = anomalistic period (seconds)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global req mu j2

pi2 = 2 * pi;

tkepler = pi2 * sma * sqrt(sma / mu);
    
e = sin(inc) * sin(inc);

ar = (sma / req) * (sma / req);
ep = 1 - ecc * ecc;
sw = sin(argper);
pp = 1 + ecc * cos(argper);

a = -0.75 * j2 * (4 - 5 * e) / (ar * sqrt(ep) * pp * pp);

b = -1.5 * j2 * pp * pp * pp / (ar * ep * ep * ep);

% nodal period
    
tnodal = tkepler * (1 + a + b);

% anomalistic period

tanomal = tkepler * (1 - 1.5 * j2 * (1 - 3 * e * sw * sw) ...
          / (ar * (1 - ecc) ^ 3));
