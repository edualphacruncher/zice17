function [ coefficients ] = tsgEstimateAnisotropicCoefficients( lGrid, sType, iMinNew, iOut )
%
% [ new_points ] = tsgRefineIsotropic( lGrid, sType, iMinNew )
%
% computes anisotropic weights corresponding to sType and refines the grid
% so that at least iMinNew points are added to the grid while respecting
% the computed anisotropic weights
%
% NOTE: this can be called only for global and sequence grids
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% sType: (string giving the tensor selection strategy)
%        'level'       'curved'         'tensor'         'iptensor'
%        'iptotal'     'ipcurved'       'qptotal'        'qpcurved'
%
% iMinNew: (integer strictly positive)
%          minimum number of new points to include in the new grid
%
% OUTPUT:
%
% coefficients: vector
%               contains the anisotropic coefficients
%

[ sFiles, sTasGrid ] = tsgGetPaths();
[ sFileG, sFileX, sFileV, sFileO, sFileW, sFileC ] = tsgMakeFilenames( lGrid.sName );

sCommand = [sTasGrid,' -refineaniso'];

sCommand = [ sCommand, ' -gridfile ', sFileG];

sCommand = [ sCommand, ' -type ', sType];
sCommand = [ sCommand, ' -ming ', num2str(iMinNew)];

if ( exists( 'iOut' ) )
    sCommand = [ sCommand, ' -refout ', num2str(iOut)];
end

% read the points for the grid
sCommand = [ sCommand, ' -of ',sFileO];
lClean.sFileO = 1;

[status, cmdout] = system(sCommand);

if ( max( size( findstr( 'ERROR', cmdout ) ) ) ~= 0 )
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if ( ~isempty(cmdout) )
        fprintf(1,['WARNING: Command had non-empty output:\n']);
        disp(cmdout);
    end
    [ coefficients ] = tsgReadMatrix( sFileO );
end

if ( exist( 'lClean' ) )
    tsgCleanTempFiles( lGrid, lClean );
end

end