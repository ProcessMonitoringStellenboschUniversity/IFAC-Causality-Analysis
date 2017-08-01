%Get root project directory
projectRootDir = fileparts(mfilename('fullpath'));

%Add project directories to path
addpath(projectRootDir);
addpath(genpath(fullfile(projectRootDir,'data')),'-end');
addpath(genpath(fullfile(projectRootDir,'documents')),'-end');
addpath(genpath(fullfile(projectRootDir,'libraries')),'-end');
addpath(genpath(fullfile(projectRootDir,'models')),'-end');
addpath(genpath(fullfile(projectRootDir,'work')),'-end');

% direct generatedd files to "work" folder
Simulink.fileGenControl('set','CacheFolder',fullfile(projectRootDir,'work'),'CodeGenFolder',fullfile(projectRootDir,'work'));
