function tsgCleanTempFiles( lGrid, lFiles )
%
% tsgCleanTempFiles( lGrid, bFlags )
%
% Cleans the temporary files associated with lGrid, i.e., all files other
% than the _FileG. If lFlags is [], then all files are cleaned. Otherwise
% only the files with entries in lFlags are cleaned, this is done to limit
% unecessary calls to the 'rm' command.
%
% You may call this function, but it exists mostly to be called by other
% functions. In fact, other functions should clean their temp files so you
% should never have to worry about this.
%
% INPUT:
%       lGrid:  a grid list created by a tsgMakeXXX(...) command (i.e., an
%               existing grid)
%
%       lFiles: object with fields sFileX, sFileV, sFileO, sFileW, sFileC
%               corresponding to the files to be deleted. If the field
%               exist, then the file is deleted. If lFiles is omitted, then
%               all files are deleted (except the permanent grid file FileG
%

% generate filenames
[ sFiles, sTasGrid ] = tsgGetPaths();
[ sFileG, sFileX, sFileV, sFileO, sFileW, sFileC ] = tsgMakeFilenames( lGrid.sName );

if ( isfield( lFiles, 'all' ) )
    sCommand = [ 'rm -fr ',sFileX,' ',sFileV,' ',sFileO,' ',sFileW,' ',sFileC];
else
    sCommand = [ 'rm -fr ' ];
    if ( isfield( lFiles, 'sFileX' ) )
        sCommand = [ sCommand, sFileX, ' ' ];
    end
    if ( isfield( lFiles, 'sFileV' ) )
        sCommand = [ sCommand, sFileV, ' ' ];
    end
    if ( isfield( lFiles, 'sFileO' ) )
        sCommand = [ sCommand, sFileO, ' ' ];
    end
    if ( isfield( lFiles, 'sFileW' ) )
        sCommand = [ sCommand, sFileW, ' ' ];
    end
    if ( isfield( lFiles, 'sFileC' ) )
        sCommand = [ sCommand, sFileC, ' ' ];
    end
end

[status, cmdout] = system(sCommand);
if ( ~isempty(cmdout) )
    disp(cmdout);
end

end