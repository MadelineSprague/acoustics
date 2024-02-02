function [traj, raychars] = getray( rayfil )

% This is an amended version of plotray that stores the individual ray
% trajectories allowing for analysis and plotting
%
% Output arguments: traj is a cell array, each cell containing the range
% (traj{i}(1,:)) and the depth (traj{i}(2,:)) variables pertaining to the
% ith range trajectory. raychars is a struct containing the characteristics
% (launch angle, number of top bounces, number of bottom bounces)

global units jkpsflag

if ( strcmp( rayfil, 'RAYFIL' ) == 0 && ~contains( rayfil, '.ray' ) )
   rayfil = [ rayfil '.ray' ]; % append extension
end

% plots a BELLHOP ray file

fid = fopen( rayfil, 'r' );   % open the file
if ( fid == -1 )
   disp( rayfil );
   error( 'No ray file exists; you must run BELLHOP first (with ray ouput selected)' );
end

% read header stuff

TITLE       = fgetl(  fid );
FREQ        = fscanf( fid, '%f', 1 );
Nsxyz       = fscanf( fid, '%f', 3 );
NBeamAngles = fscanf( fid, '%i', 2 );

DEPTHT      = fscanf( fid, '%f', 1 );
DEPTHB      = fscanf( fid, '%f', 1 );

Type        = fgetl( fid );
Type        = fgetl( fid );

Nsx    = Nsxyz( 1 );
Nsy    = Nsxyz( 2 );
Nsz    = Nsxyz( 3 );

Nalpha = NBeamAngles( 1 );
Nbeta  = NBeamAngles( 2 );

% Extract letters between the quotes
nchars = strfind( TITLE, '''' );   % find quotes
TITLE  = [ TITLE( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 7 - ( nchars( 2 ) - nchars( 1 ) ) ) ];
TITLE  = deblank( TITLE );  % remove whitespace

nchars = strfind( Type, '''' );   % find quotes
Type   = Type( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
%Type  = deblank( Type );  % remove whitespace

% read rays

traj = cell([Nalpha Nsz]); 
raychars = struct; 

% this could be changed to a forever loop
for isz = 1 : Nsz
   for ibeam = 1 : Nalpha
      raychars(ibeam, isz).alpha0    = fscanf( fid, '%f', 1 );
      nsteps    = fscanf( fid, '%i', 1 );

      raychars(ibeam, isz).NumTopBnc = fscanf( fid, '%i', 1 );
      raychars(ibeam, isz).NumBotBnc = fscanf( fid, '%i', 1 );

      

      if isempty( nsteps ); break; end
      switch Type
         case 'rz'
            ray = fscanf( fid, '%f', [2 nsteps] );
            
            r = ray( 1, : );
            z = ray( 2, : );
         case 'xyz'
            ray = fscanf( fid, '%f', [3 nsteps] );
            
            xs = ray( 1, 1 );
            ys = ray( 2, 1 );
            r = sqrt( ( ray( 1, : ) - xs ).^2 + ( ray( 2, : ) - ys ).^2 );
            z = ray( 3, : );
      end
      
      if ( strcmp( units, 'km' ) )
         r = r / 1000;   % convert to km
      end
      
      traj{ibeam, isz} = [r; z]; 
      
   end	% next beam
end % next source depth

fclose( fid );

end