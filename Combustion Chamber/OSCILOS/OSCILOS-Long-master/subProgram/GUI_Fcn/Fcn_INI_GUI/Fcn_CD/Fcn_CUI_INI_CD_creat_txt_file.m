% --- creat one txt file including the combustor shape
function Fcn_CUI_INI_CD_creat_txt_file 
% This program is used to create a example file including a cambridge
% combustor information
global CI
% ---------------------------------
l_sample    = [200 300 350]./1000;
ro_sample   = [50 17.5 35]./1000;
ri_sample   = [4 12.5 0]./1000; 
r_sample    = (ro_sample.^2-ri_sample.^2).^0.5;
x_sample    = [0,cumsum(l_sample)];
r_sample    = [r_sample,r_sample(end)];
SectionIndex= [0 0 11 0];
TubeIndex   = [0 0 0 0];
data        = cat(1,x_sample,r_sample,SectionIndex,TubeIndex);
% ---------------------------------
currentFolder   = pwd;
currentFolder   = fullfile(currentFolder,CI.SD.name_program,'CD_example.txt');
fid             = fopen(currentFolder,'wt');
data_title      = {'x[m]','r[m]','InterfaceIndex','ModuleIndex'};
fprintf(fid,'%s\b',data_title{1:end});
fprintf(fid,'\n');
fprintf(fid,'%6.5f\b%6.5f\b%6.0f\b%6.0f\b\n',data);
fclose(fid);
% -----------------------------end-----------------------------------------