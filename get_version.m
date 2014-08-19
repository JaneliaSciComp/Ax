function VERSION=get_version()

if isdeployed
  importdata('version.txt');
  VERSION=ans{1};
else
  tmp=fullfile(fileparts(mfilename('fullpath')),'.git');
  if ispc
    [s,VERSION]=system(['"c:\\Program Files (x86)\Git\bin\git" --git-dir ' tmp ' log -1 --pretty=format:"%ci %H"']);
  else
    [s,VERSION]=system(['TERM=xterm git --git-dir ' tmp ' log -1 --pretty=format:"#%ci %H#"']); 
    tmp=strfind(VERSION,'#');
    VERSION=VERSION((tmp(1)+1):(tmp(2)-1));
  end
  if s
    warning('can''t find git.  to save version info, git-bash must be installed.');
    VERSION='don''t know';
  end
end
