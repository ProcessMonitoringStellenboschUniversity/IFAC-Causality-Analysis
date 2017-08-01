%Remove project folders from path
rmpath(projectRootDir);
rmpath(genpath(fullfile(projectRootDir,'data')));
rmpath(genpath(fullfile(projectRootDir,'documents')));
rmpath(genpath(fullfile(projectRootDir,'libraries')));
rmpath(genpath(fullfile(projectRootDir,'models')));
rmpath(genpath(fullfile(projectRootDir,'work')));

%Reset file generation folder to its default
Simulink.fileGenControl('reset');