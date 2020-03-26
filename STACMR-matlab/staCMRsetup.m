function staCMRsetup(filename)
% function staCMRsetup(filename)
% If filename is specified, add it to the java classpath
% Otherwise, add the newest version of fxMR to the java classpath from the
% higher director. Filename must be in the format "fxMR-X.X.X.jar".
bookinfo = 'Dunn, J. C. & Kalish, M. L. (2018). State-Trace Analysis. Springer.';
fprintf('STACMR program library Version 26.03.2020\n');
fprintf (['Utility programs for use with the book:\n', bookinfo, '\n\n'])

if nargin==0
    currdir = pwd;
    pathstr = fileparts(which('/staCMRsetup.m')); eval(['cd ''', pathstr,'''']); % go to STACMR folder
    cd ('java');  % change directory
    files=dir('fxMR-*.jar');
    versions=zeros(size(files, 1),3);
    for i=1:size(files)
       fi=files(i).name;
       vit=regexp(fi, 'fxMR-(\d*)\.(\d*)\.(\d*)\.jar','tokens');
       if ~isempty(vit)
           vis=vit{1};
           versions(i,:)=[str2double(vis{1}),str2double(vis{2}),str2double(vis{3})];
       end
    end
    versions=sortrows(versions); 
    v=versions(size(versions,1),:);
    filename=strcat('fxMR-',num2str(v(1)),'.',num2str(v(2)),'.',num2str(v(3)),'.jar');
end
disp(['STACMR linked to java library ', filename])
%which(filename);
javaclasspath (which(filename));
eval(['cd ''', currdir,'''']); % return to working directory

