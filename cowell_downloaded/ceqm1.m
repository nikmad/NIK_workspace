function ydot = ceqm1 (t, y)

% first order form of Cowell's equations of orbital motion

% includes perturbations due to:

%   non-spherical earth gravity
%   aerodynamic drag (us 76 model)
%   solar radiation pressure
%   point mass sun and moon gravity

% input

%  t = current simulation time (seconds)
%  y = current eci state vector (km and km/sec)

% output

%  ydot = eci acceleration vector (km/sec^2)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jdate0 smu mmu

global isun imoon idrag isrp

asun = zeros(3, 1);

amoon = zeros(3, 1);

adrag = zeros(3, 1);

asrp = zeros(3, 1);

ydot = zeros(1, 6);

% compute geopotential perturbation

agrav = gravity(t, y);

if (isun == 1)
    
    % solar perturbations

    jdate = jdate0 + t / 86400;

    [rasc, decl, rsun] = sun1(jdate);

    % compute heliocentric position vector of the spacecraft

    for i = 1:3
        
        rs2sc(i) = y(i) - rsun(i);
        
    end

    % f(q) formulation

    for i = 1:1:3
        
        vtmp(i) = y(i) - 2.0d0 * rsun(i);
        
    end

    dot1 = dot(y(1:3), vtmp);

    dot2 = dot(rsun, rsun);

    qsun = dot1 / dot2;

    fsun = qsun * ((3.0d0 + 3.0d0 * qsun + qsun * qsun) ...
        / (1.0d0 + (1.0d0 + qsun)^1.5d0));

    d3sun = norm(rs2sc)^3;

    % point-mass gravity of the sun

    for i = 1:1:3
        
        asun(i) = -smu * (y(i) + fsun * rsun(i)) / d3sun;
        
    end
    
end

if (imoon == 1)

    % lunar perturbations

    jdate = jdate0 + t / 86400;

    [rasc, decl, rmoon] = moon(jdate);

    % compute selenocentric position vector of the spacecraft

    for i = 1:1:3
        
        rm2sc(i) = y(i) - rmoon(i);
        
    end

    % f(q) formulation

    for i = 1:1:3
        
        vtmp(i) = y(i) - 2.0d0 * rmoon(i);
        
    end

    dot1 = dot(y(1:3), vtmp);

    dot2 = dot(rmoon, rmoon);

    qmoon = dot1 / dot2;

    fmoon = qmoon * ((3.0d0 + 3.0d0 * qmoon + qmoon * qmoon) ...
        / (1.0d0 + (1.0d0 + qmoon)^1.5d0));

    d3moon = norm(rm2sc)^3;

    % point-mass gravity of the moon

    for i = 1:1:3
        
        amoon(i) = -mmu * (y(i) + fmoon * rmoon(i)) / d3moon;
        
    end
end

if (idrag == 1)
    
    % compute atmospheric drag perturbation
    
    adrag = drag1(y);
    
end

if (isrp == 1)
    
    % compute solar radiation pressure perturbation
    
    asrp = srp(t, y);
    
end

% total acceleration vector

ydot(1) = y(4);
ydot(2) = y(5);
ydot(3) = y(6);

for i = 1:1:3
    
    ydot(i + 3) = agrav(i) + asun(i) + amoon(i) + adrag(i) + asrp(i);
    
end

