% cowell1.m       May 6, 2008

% orbital motion of Earth satellites
% using Cowell's method and equations
% of motion in rectangular coordinates

% includes perturbations due to:

%   non-spherical earth gravity
%   aerodynamic drag (us 76 model)
%   solar radiation pressure
%   sun and moon

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global dtr rtd atr
global flat req mu smu mmu aunit omega

global j2 lgrav mgrav
global jdate0 gst0 isun imoon idrag isrp
global rkcoef ccoef scoef ad76 bcoeff csrp0

% astrodynamic and utility constants

om_constants;

dsun = 696000;             % radius of the sun
dmoon = 1738;              % radius of the moon
ps = 0.00456;              % solar constant

atr = dtr / 3600;

% initialize rkf78 integrator

rkcoef = 1;

% define number of differential equations

neq = 6;

% read gravity model coefficients

gmfile = 'egm96.dat';

[ccoef, scoef] = readegm(gmfile);

% extract j2 coefficient

j2 = -ccoef(3, 1);

if (gmfile == 'egm96.dat')
    % egm96 constants

    mu = 398600.4415;
    req = 6378.1363;
    omega = 7.292115e-5;
end

% begin simulation

clc; home;

fprintf('\n                program cowell1\n');

fprintf('\n  < Earth orbital motion - Cowell''s method >\n');

% request initial calendar date and time

fprintf('\n\ninitial calendar date and time\n');

[month, day, year] = getdate;

[uthr, utmin, utsec] = gettime;

while (1)
    fprintf('\nplease input the simulation period (days)\n');

    ndays = input('? ');

    if (ndays > 0)
        break;
    end
end

tsim = 86400 * ndays;

while (1)
    fprintf('\nplease input the algorithm error tolerance');
    fprintf('\n(a value between 1.0e-8 and 1.0e-10 is recommended)\n');

    tetol = input('? ');

    if (tetol > 0)
        break;
    end
end

% request degree and order of gravity model

fprintf('\n\ngravity model inputs \n');

while(1)
    fprintf('\nplease input the degree of the gravity model (zonals)\n');
    fprintf('(0 <= zonals <= 18)\n');

    lgrav = input('? ');

    if (lgrav >= 0 && lgrav <= 18)
        break;
    end
end

while(1)
    fprintf('\nplease input the order of the gravity model (tesserals)\n');
    fprintf('(0 <= tesserals <= 18)\n');

    mgrav = input('? ');

    if (mgrav >= 0 && mgrav <= 18)
        break;
    end
end

fprintf('\n\norbital perturbations \n');

% include solar perturbation?

while(1)
    fprintf('\nwould you like to include solar perturbations (y = yes, n = no)\n');

    yn = lower(input('? ', 's'));

    if (yn == 'y' || yn == 'n')
        break;
    end
end

if (yn == 'y')
    isun = 1;
else
    isun = 0;
end

% include lunar perturbation?

while(1)
    fprintf('\nwould you like to include lunar perturbations (y = yes, n = no)\n');

    yn = lower(input('? ', 's'));

    if (yn == 'y' || yn == 'n')
        break;
    end
end

if (yn == 'y')
    imoon = 1;
else
    imoon = 0;
end

% include drag perturbation?

while(1)
    fprintf('\nwould you like to include drag perturbations (y = yes, n = no)\n');

    yn = lower(input('? ', 's'));

    if (yn == 'y' || yn == 'n')
        break;
    end
end

if (yn == 'y')
    idrag = 1;
else
    idrag = 0;
end

% include solar radiation pressure perturbation?

while(1)
    fprintf('\nwould you like to include srp perturbations (y = yes, n = no)\n');

    yn = lower(input('? ', 's'));

    if (yn == 'y' || yn == 'n')
        break;
    end
end

if (yn == 'y')
    isrp = 1;
else
    isrp = 0;
end

% create and display graphics?

while (1)

    fprintf('\n\nwould you like to create and display graphics (y = yes, n = no)\n');

    yn = lower(input('? ', 's'));

    if (yn == 'y' || yn == 'n')
        break;
    end
end

if (yn == 'y')
    igraph = 1;
else
    igraph = 0;
end

if (igraph == 1)

    % request graphics step size

    while (1)
        fprintf('\n\nplease input the graphics step size (minutes)\n');

        dtstep = input('? ');

        if (dtstep > 0)
            break;
        end
    end

    dtstep = 60 * dtstep;

end

% request initial orbital elements

while (1)
    fprintf('\n\n orbital elements menu\n');

    fprintf('\n   <1> user input\n');

    fprintf('\n   <2> data file\n\n');

    slct = input('? ');

    if (slct == 1 || slct == 2)
        break;
    end
end

if (slct == 1)
    % interactive request

    fprintf('\n\ninitial orbital elements \n');

    oev1 = getoe([1;1;1;1;1;1]);
else
    % data file
    fprintf('\nplease input the orbital elements data file name');
    fprintf('\n(be sure to include the file name extension)\n');

    filename = input('? ', 's');

    [fid, oev1] = readoe1(filename);
end

if (idrag == 1)
    % atmospheric drag perturbations

    fprintf('\naerodynamic drag inputs \n');

    % read atmosphere data file

    [fid, ad76] = read76;

    while(1)
        fprintf('\nplease input the drag coefficient (non-dimensional)\n');

        cd = input('? ');

        if (cd > 0)
            break;
        end
    end

    while(1)
        fprintf('\nplease input the cross-sectional area (square meters)\n');

        areadrag = input('? ');

        if (areadrag > 0)
            break;
        end
    end

    while(1)
        fprintf('\nplease input the spacecraft mass (kilograms)\n');

        scmass = input('? ');

        if (scmass > 0)
            break;
        end
    end

    bcoeff = 1.0e-6 * areadrag * cd /scmass;
end

if (isrp == 1)
    % solar radiation pressure perturbations

    fprintf('\nsolar radiation pressure inputs \n');

    while(1)
        fprintf('\nplease input the reflectivity constant (non-dimensional)\n');

        reflect = input('? ');

        if (reflect > 0)
            break;
        end
    end

    while(1)
        fprintf('\nplease input the cross-sectional area (square meters)\n');

        areasrp = input('? ');

        if (areasrp > 0)
            break;
        end
    end

    while(1)
        fprintf('\nplease input the spacecraft mass (kilograms)\n');

        scmass = input('? ');

        if (scmass > 0)
            break;
        end
    end

    csrp0 = 0.000001 * reflect * ps * aunit^2 * areasrp / scmass;
end

% determine initial eci state vector

[ri, vi] = orb2eci(mu, oev1);

% load initial position and velocity vectors

for i = 1:1:3
    yi(i) = ri(i);

    yi(i + 3) = vi(i);
end

% compute initial julian date

jdate0 = julian(month, day, year) ...
    + uthr / 24 + utmin / 1440 + utsec / 86400;

% compute initial greenwich sidereal time

gst0 = gast1(jdate0);

fprintf('\n\n  working ...\n');

if (igraph == 1)
    ti = -dtstep;
end

npts = 0;

if (igraph == 1)
    % create initial graphics data

    npts = npts + 1;

    t = 0;

    oev = eci2orb2 (mu, gst0, omega, t, ri, vi);

    xdata(npts) = t / 86400;

    for i = 1:1:11
        switch i
            case {1, 2}
                ydata(i, npts) = oev(i);
            case {3, 4, 5, 6}
                ydata(i, npts) = rtd * oev(i);
            case 7
                ydata(i, npts) = rtd * oev(21);
            case 8
                ydata(i, npts) = oev(22);
            case 9
                ydata(i, npts) = oev(16);
            case 10
                ydata(i, npts) = rtd * oev(15);
            case 11
                ydata(i, npts) = rtd * oev(14);
        end
    end
end

while(1)
    % step size guess

    h = 30;

    if (igraph == 1)
        ti = ti + dtstep;
        tf = ti + dtstep;
    else
        ti = 0;
        tf = tsim;
    end

    % integrate from ti to tf

    yfinal = rkf78('ceqm1', neq, ti, tf, h, tetol, yi);

    if (igraph == 1)

        % create graphics data

        npts = npts + 1;

        % compute current state vector

        for i = 1:1:3
            rf(i) = yfinal(i);
            vf(i) = yfinal(i + 3);
        end

        % compute current orbital elements

        oev2 = eci2orb2(mu, gst0, omega, tf, rf, vf);

        % altitude check

        alt = oev2(1) * (1.0 - oev2(2)) - req;

        if (alt <= 90)
            break;
        end

        xdata(npts) = tf / 86400;

        for i = 1:1:11
            switch i
                case {1, 2}
                    ydata(i, npts) = oev2(i);
                case {3, 4, 5, 6}
                    ydata(i, npts) = rtd * oev2(i);
                case 7
                    ydata(i, npts) = rtd * oev2(21);
                case 8
                    ydata(i, npts) = oev2(22);
                case 9
                    ydata(i, npts) = oev2(16);
                case 10
                    ydata(i, npts) = rtd * oev2(15);
                case 11
                    ydata(i, npts) = rtd * oev2(14);
            end
        end
    end

    yi = yfinal;

    % check for end of simulation

    if (tf >= tsim)
        break;
    end
end

% compute final state vector and orbital elements

for i = 1:1:3
    rf(i) = yfinal(i);
    
    vf(i) = yfinal(i + 3);
end

oev2 = eci2orb1(mu, rf, vf);

% print results

[cdstr0, utstr0] = jd2str(jdate0);

jdatef = jdate0 + tf / 86400;

[cdstrf, utstrf] = jd2str(jdatef);

fprintf('\n                     program cowell1\n');

fprintf('\n    < Earth orbital motion - Cowell''s method >\n');

fprintf('\ninitial calendar date       ');

disp(cdstr0);

fprintf('initial universal time      ');

disp(utstr0);

fprintf('\nfinal calendar date         ');

disp(cdstrf);

fprintf('final universal time        ');

disp(utstrf);

oeprint1(mu, oev2);

fprintf('\ndegree of gravity model    %2i', lgrav);

fprintf('\norder of gravity model     %2i \n', mgrav);

if (isun == 1)
    fprintf('\nsimulation includes solar perturbations');
end

if (imoon == 1)
    fprintf('\nsimulation includes lunar perturbations');
end

if (idrag == 1)
    fprintf('\nsimulation includes drag perturbations');
end

if (isrp == 1)
    fprintf('\nsimulation includes srp perturbations');
end

if (igraph == 1)

    % request item to plot

    while (1)

        while (1)
            fprintf('\n\nplease select the item to plot\n');

            fprintf('\n  <1> semimajor axis');

            fprintf('\n  <2> eccentricity');

            fprintf('\n  <3> orbital inclination');

            fprintf('\n  <4> argument of perigee');

            fprintf('\n  <5> right ascension of the ascending node');

            fprintf('\n  <6> true anomaly');

            fprintf('\n  <7> geodetic perigee altitude');

            fprintf('\n  <8> geodetic apogee altitude');

            fprintf('\n  <9> geodetic altitude');

            fprintf('\n  <10> east longitude');

            fprintf('\n  <11> geodetic latitude\n\n');

            oetype = input('? ');

            if (oetype >= 1 && oetype <= 11)
                break;
            end
        end

        % create and label plot

        plot(xdata, ydata(oetype, :));

        switch oetype
            case 1
                ylabel('Semimajor Axis (kilometers)', 'FontSize', 12);
            case 2
                ylabel('Eccentricity', 'FontSize', 12);
            case 3
                ylabel('Inclination (degrees)', 'FontSize', 12);
            case 4
                ylabel('Argument of Perigee (degrees)', 'FontSize', 12);
            case 5
                ylabel('RAAN (degrees)', 'FontSize', 12);
            case 6
                ylabel('True Anomaly (degrees)', 'FontSize', 12');
            case 7
                ylabel('Perigee Altitude (kilometers)', 'FontSize', 12');
            case 8
                ylabel('Apogee Altitude (kilometers)', 'FontSize', 12);
            case 9
                ylabel('Altitude (kilometers)', 'FontSize', 12);
            case 10
                ylabel('East Longitude (degrees)', 'FontSize', 12);
            case 11
                ylabel('Geodetic Latitude (degrees)', 'FontSize', 12);
        end

        title('Long Term Orbit Evolution - Cowell''s Method', 'FontSize', 16);

        xlabel('Simulation Time (days)', 'FontSize', 12);

        grid;

        zoom on;

        % create eps file

        print -depsc -tiff -r300 cowell1.eps

        % create another plot?

        while (1)

            fprintf('\n\nwould you like to create another plot (y = yes, n = no)\n');

            yn = lower(input('? ', 's'));

            if (yn == 'y' || yn == 'n')
                break;
            end
        end

        if (yn == 'n')
            break;
        end
    end
end
