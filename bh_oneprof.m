% Madeline Sprague (2/1/2024) - mjsprague@uri.edu
%
% This function is similar to hycom_bh except it runs Bellhop using a
% single non-range-dependent sound speed profile with no specified real
% bathymetry. However, this function does not use prompts; all variables
% are specified in the argument. 
%
% USAGE: bh_oneprof(c, z, r, sz, freq, runtype, fr_min, fr_max, nrays, filename, title)
% 
% INPUTS: 
%   c: double array (Nx1)
%       Vertical sound speed profile in meters/second
%   z: double array (Nx1)
%       Corresponding depth at each sound speed value in meters, positive
%       downwards
%   r: double scalar 
%       Desired maximum distance of propagation to be modeled in kilometers
%   sz: double scalar 
%       Depth of source in meters
%   freq: double scalar
%       Frequency of source in Hz
%   runtype: char 
%       Enter 'R' for ray tracing, 
%             'A' for arrivals, 
%             'I' for incoherent transmission loss, 
%             'C' for coherent transmission loss 
%   fr_min, fr_max: double scalar 
%       The minimum and maximum elevation angles for the ray fan range 
%   nrays: double integer 
%       The number of rays to be used by Bellhop (set to 0 to let Bellhop
%       choose). 
%   filename: char
%       The desired filename for the environment file and the Bellhop
%       output files. 
%   title: char 
%       Desired title of outputted plots

function bh_oneprof(c, z, r, sz, freq, runtype, fr_min, fr_max, nrays, filename, title)

% create folder for output files 

    filepath = cd; 
    mkdir(filename)
    addpath([filepath '\' filename])
    cd([filepath '\' filename]) 
    
% prep variables 

    runtype = sprintf('''%s''', runtype);
    title   = sprintf('''%s''', title);
    maxdep = max(z); 
    
% define receiver parameters using specified max distance

    n_rec_z = 100; 
    z_rec_min = ceil(min(z)); 
    z_rec_max = floor(max(z)); 
    n_rec_x = 1000; 
    x_rec_min = 0; 
    x_rec_max = floor(r); 
    
% write .env file 

    fid_env = fopen([filename '.env'], 'w'); 
    fprintf(fid_env, [title '\r\n']);                 % title 
    fprintf(fid_env, '%3.1f\r\n', freq);              % ray frequency
    fprintf(fid_env, '%d\r\n', 1);                    % NMEDIA
    fprintf(fid_env, '''SVM''\r\n');                  % interpolation/parameters 
    fprintf(fid_env, '%d %1.1f %4.1f\r\n', ...
            [0 0.0 maxdep]);                          % default flat bottom at 5000 m
    fprintf(fid_env, '%3.2f %5.4f /\r\n', ...
            [z'; c']);                                % c profile
    fprintf(fid_env, '''A_'' 0.0 \r\n');              % acousto-elastic bottom condition 
    fprintf(fid_env, [char(string(maxdep)) ...
            ' 1600.0 0.0 1.0 /\r\n']);                % bottom depth, sound speed of bottom
    fprintf(fid_env, [char(string(1)) '\r\n']);       % number of sound sources 
    fprintf(fid_env, '%3.1f /\r\n', sz);              % depth of sound source (approximate)
    fprintf(fid_env, [char(string(n_rec_z)) '\r\n']); % number of receiver depths (z)
    fprintf(fid_env, [char(string(z_rec_min)) ' ' ...
                      char(string(z_rec_max)) ...
                      ' /\r\n']);                     % z range of receivers 
    fprintf(fid_env, [char(string(n_rec_x)) '\r\n']); % number of receiver ranges (x)
    fprintf(fid_env, [char(string(x_rec_min)) ' ' ...
                      char(string(x_rec_max)) ...
                      ' /\r\n']);                     % x range of receivers 
    fprintf(fid_env, [runtype '\r\n']);               % RunType ('R' for ray trace)
    fprintf(fid_env, [char(string(nrays)) '\r\n']);   % Number of rays to plot 
    fprintf(fid_env, [char(string(fr_min)) ' ' ...
                      char(string(fr_max)) ...
                      ' /\r\n']);                     % Angular limits of beam fan   
    fprintf(fid_env, ['0.0 ' char(string(maxdep)) ...
                      ' ' char(string(r))]);          % ray step size, z limit, x limit 
    fclose(fid_env); 

% now run bellhop 

    bellhop(filename); 

% plot results 

    % plot SSP in background
    
        [r_grid, z_grid] = meshgrid([0 r*1000], z); 
        c_grid = [c c]; 

    if contains(runtype, 'R')
        
        figure; hold on 
        pcolor(r_grid, z_grid, c_grid); shading('flat'); 
        plotray([filename '.ray'])
        cbar = colorbar; cbar.Label.String = 'Sound Speed (m/s)'; 
        a = gca; 
            xticklabels(string(a.XTick/1000)); 
        xlabel('Range (km)')
        ylabel('Depth (m)')
        labelformat([14 16])
        
    elseif strcmp(runtype, '''C''') | strcmp(runtype, '''I''')
        
        figure; hold on 
        plotshd([filename '.shd'])
        cbar = colorbar; cbar.Label.String = 'Transmission Loss (dB/mkHz)'; 
        a = gca; 
            xticklabels(string(a.XTick/1000)); 
        xlabel('Range (km)')
        ylabel('Depth (m)')
        labelformat([14 16])

    elseif contains(runtype, 'A')

            plotarr([filename '.arr'], 1, 1, 1)     
            
    end
    
    cd(filepath)

end

