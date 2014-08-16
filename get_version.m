function VERSION_AX=get_version()

if ispc
  [s,VERSION_AX]=system('"c:\\Program Files (x86)\Git\bin\git" log -1 --pretty=format:"%ci %H"');
else
  [s,VERSION_AX]=system('TERM=xterm git log -1 --pretty=format:"#%ci %H#"');
  tmp=strfind(VERSION_AX,'#');
  VERSION_AX=VERSION_AX((tmp(1)+1):(tmp(2)-1));
end
if s
    warning('cant''t find git.  to save version info, git-bash must be installed.');
end