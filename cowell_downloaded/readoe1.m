function [fid, oev] = readoe1(filename)

% read orbital elements data file

% required by cowell.m

% input

%  filename = name of orbital element data file

% output

%  fid = file id

%  oev(1) = semimajor axis
%  oev(2) = orbital eccentricity
%  oev(3) = orbital inclination
%  oev(4) = argument of perigee
%  oev(5) = right ascension of the ascending node
%  oev(6) = true anomaly

% NOTE: all angular elements are returned in radians

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtr = pi / 180;

% open data file

fid = fopen(filename, 'r');

% check for file open error

if (fid == -1)
    
    clc; home;
    
    fprintf('\n\n  error: cannot find this file!!');
    
    pause;
    
    return;
    
end

% read 23 lines of data file

for i = 1:1:23
    
    cline = fgetl(fid);
    
    switch i
        
        case 3
            
            oev(1) = str2double(cline);
            
        case 7
            
            oev(2) = str2double(cline);
            
        case 11
            
            oev(3) = dtr * str2double(cline);
            
        case 15
            
            oev(4) = dtr * str2double(cline);
            
        case 19
            
            oev(5) = dtr * str2double(cline);
            
        case 23
            
            oev(6) = dtr * str2double(cline);
            
    end
    
end

status = fclose(fid);

