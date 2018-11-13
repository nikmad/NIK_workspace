function asrp = srp (t, sv)

% ECI acceleration vector due to solar radiation pressure

% input

%  t  = simulation time (seconds)
%  sv = current state vector (km and km/sec)

% output

%  asrp = eci acceleration vector (km/sec^2)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global req aunit csrp0 jdate0

asrp = zeros(3, 1);

% current julian date

jdate = jdate0 + t / 86400;

% compute solar ephemeris

[rasc, decl, rsun] = sun1(jdate);

% geocentric distance of the sun (kilometers)

rmsun = norm(rsun);

% geocentric unit vector of the sun

usun = rsun / norm(rsun);

% geocentric distance of the spacecraft

rscm = norm(sv(1:3));

% compute unit position vector of spacecraft

usat = sv(1:3) / rscm;

% determine shadow conditions

a = usat(2) * usun(3) - usat(3) * usun(2);
b = usat(3) * usun(1) - usat(1) * usun(3);
c = usat(1) * usun(2) - usat(2) * usun(1);

d = sqrt(a * a + b * b + c * c);

e = dot(usat, usun);

u = asin(0.00460983743 / (rmsun / aunit));

p = asin(0.0046951089 / (rmsun / aunit));

if (e > 0.0)
    q = -d;
else
    q = d;
end

% increase the Earth's radius by 90 km
% to account for the atmosphere

ratm = req + 90.0;

b = asin(ratm / rscm);

v = b - u;
w = b + p;

x = sin(v);
y = sin(w);

% determine shadow conditions

if (q <= y && q > x)
    
    % penumbra
    
    csrp = 0.0;
    
elseif (q <= x && q >= 0)
    
    % umbra
    
    csrp = 0.0;
    
else
    
    % sunlight
    
    csrp = csrp0;
    
end

% compute vector from the satellite to the sun

for i = 1:1:3
    
    rs2s(i) = sv(i) - rsun(i);
    
end

rs2sm = norm(rs2s);

rs2sm3 = rs2sm * rs2sm * rs2sm;

% compute acceleration vector

asrp = csrp * rs2s / rs2sm3;
