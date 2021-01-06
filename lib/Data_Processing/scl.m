function scl()
    % Set [path to] Current Location
    %
    % DESCRIPTION:
    %   Changes the directory to the path of the current file in editor.
    %
    % SETUP:
    %   Place scl.m in a directory and run the following line once.
    %
    %       addpath('<path to the directory of scl.m>'); savepath;
    %
    % Eg:
    %   addpath('/Users/totallySomeoneElse/Documents/MATLAB/'); savepath;
    try
        editor_service = com.mathworks.mlservices.MLEditorServices;
        editor_app = editor_service.getEditorApplication;
        active_editor = editor_app.getActiveEditor;
        storage_location = active_editor.getStorageLocation;
        file = char(storage_location.getFile);
        path = fileparts(file);
        cd(path);
    catch e
        error('scl failed.\n\n%s',e.message);
    end
end