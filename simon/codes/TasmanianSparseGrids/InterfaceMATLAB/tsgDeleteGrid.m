function tsgDeleteGrid( lGrid )
%
% tsgDeleteGrid( lGrid )
%
% deletes all of the background files used by the grid
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeGrid(...)
%
% NOTE: lGrid gets deleted and it can no longer be used
%

[ sFiles, sTasGrid ] = tsgGetPaths();
[ sFileG, sFileX, sFileV, sFileO, sFileW, sFileC ] = tsgMakeFilenames( lGrid.sName );

sCommand = [ 'rm -fr ',sFileG,' ',sFileX,' ',sFileV,' ',sFileO,' ',sFileW,' ',sFileC];

[status, cmdout] = system(sCommand);
if ( ~isempty(cmdout) )
    fprintf(1,['WARNING: Command had non-empty output:\n']);
    disp(cmdout);
end

end