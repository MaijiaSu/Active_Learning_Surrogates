function G = ansys_project(u,fun_par)
% 设置工作路径，在该路径下三个文件：Project_APDL.mac、DataIn.txt、DataOut.txt
Ns = size(u,1);
OriginFolder = pwd;
t_start = tic;
disp(['Running ANSYS Program','......'])

Ansys_Path = 'F:\ANSYS\ANSYS Inc\v192\ansys\bin\winx64\ANSYS192.exe';
Ansys_Worksapce = 'F:\matlab2016b\bin\dsmj\ToolBox\LSF_Library\ANSYS_PROJECT';
Jobname = 'MAT2ANS_project';
Project_APDL = 'command_input';
OutputFile = 'process';

cd(Ansys_Worksapce)

for iii = 1:Ns
    
    X = u(iii,:);
    
    % 
    fid = fopen([Ansys_Worksapce,'\','DataIn.txt'],'w');   
    fprintf(fid,'%f\n',X);
    fclose(fid);
       
    % 记录当前文件
    
    
    cmd =  ['SET KMP_STACKSIZE=1000M & ',...
            '"' Ansys_Path '"' ...
            ' -b' ...
            ' -j ' Jobname ...
            ' -i "' Project_APDL  '.mac"' ...
            ' -o "' OutputFile '.txt"' ...
            ];
        
    [o, txt] = system(cmd, '-echo');

    disp(['Ending(',num2str(iii),'/',num2str(Ns),')',':',...
        'elapsed time with ',num2str(toc(t_start),'%.1f'),'s'])
    
     % 输出
     data = textread('DataOut.txt');
     G(iii,1) = 0.4-max(abs(data(:,2)));
               
end

cd(OriginFolder);

end