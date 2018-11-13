function [ccoef, scoef] = readegm(fname)

% read gravity model data file

% input

%  fname = name of gravity data file

% output

%  ccoef, scoef = gravity model coefficients

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open gravity model data file

fid = fopen(fname, 'r');

% check for file open error

if (fid == -1)
    clc; home;
    fprintf('\n\n  error: cannot find file!!');
    keycheck;
    return;
end

% read the data file

a = fscanf(fid, '%i %i %E %E', [inf]);

status = fclose(fid);

% initialize coefficients

ccoef = zeros(18, 18);

scoef = zeros(18, 18);

% extract gravity model coefficients

kk = -3;

for i = 1:1:187

    kk = kk + 4;

    ii = a(kk);
    jj = a(kk + 1);

    ccoef(ii + 1, jj + 1) = a(kk + 2);

    scoef(ii + 1, jj + 1) = a(kk + 3);
end


