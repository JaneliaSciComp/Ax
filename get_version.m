function VERSION=get_version(repo)

if isdeployed
  importdata('version.txt');
  VERSION=ans{1};
else
  if ispc
    [s,VERSION]=system('"c:\\Program Files (x86)\Git\bin\git" log -1 --pretty=format:"%ci %H"');
  else
    [s,VERSION]=system(['cd $(dirname $(which ' repo ')) && git log -1 --pretty=format:"#%ci %H#"']);  % TERM=xterm
    tmp=strfind(VERSION,'#');
    VERSION=VERSION((tmp(1)+1):(tmp(2)-1));
  end
  if s
    warning('can''t find git.  to save version info, git-bash must be installed.');
    VERSION='don''t know';
  end
end