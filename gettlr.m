% Amended version of plottlr that gets the plotted data without plotting
% it. 
% 
% usage:
% [r, tl] = gettlr( filename, rdt )
% where
%   filename is the shadefile (with extension)
%   rdt is the receiver depth in m
%   if rdt is a vector then one plot is generated for each element


function [r, tl] = gettlr( filename, rdt )



global units

disp( 'PlotTLr uses the first bearing and source depth in the shade file; check OK' )
itheta = 1;
isd    = 1;

% read

[ PlotTitle, ~, freq, ~, ~, Pos, pressure ] = read_shd( filename );
rkm = Pos.r.r / 1000.0;         % convert to km

pressure = pressure( itheta, isd, :, : );

tlt = abs( pressure );	            % this is really the negative of TL
tlt( tlt == 0 ) = max( max( tlt ) ) / 1e10;      % replaces zero by a small number
tlt = -20.0 * log10( tlt );          % so there's no error when we take the log

% interpolate the TL field at the receiver depth
% note: interp1 won't interpolate a vector with only 1 element

if ( length( Pos.r.z ) == 1 )
  tlslice = squeeze( tlt( 1, : ) );
  rdt     = Pos.r.z( 1 );
else
  TLtemp  = squeeze( tlt );   % need to avoid dimensional problems when TLT has 1, 2, or 3 dimensions
  tlslice = interp1( Pos.r.z, TLtemp, rdt );
end

r = rkm; 
tl = tlslice'; 

end