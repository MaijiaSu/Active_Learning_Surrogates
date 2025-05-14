function G = MAT2ANS_project(u,fun_par)
% ���ù���·�����ڸ�·���������ļ���Project_APDL.mac��DataIn.txt��DataOut.txt
% ����Project_APDL.mac�ǵ���ģ�ͣ��������ANSYSģ���ļ�.db
if strcmp(fun_par.StadNrom,'Yes')    
    [Ns,Ndim] = size(u);  
    muX = fun_par.muX;
    sigmaX = fun_par.sigmaX; 
    X = zeros(Ns,Ndim);
    for n = 1:Ndim
        X(:,n) = u(:,n)*sigmaX(n) + muX(n);
    end
else  % ����ת�� 
    X = u;
end

Ns = size(u,1);
OriginFolder = pwd;
t_start = tic;
disp(['Running ANSYS Program','......'])
 
% fun_par.Ansys_solver_Path = 'C:\Software\Ansys 19.2\ANSYS Inc\v192\ansys\bin\winx64\ANSYS192.exe';
% fun_par.Ansys_Worksapce = 'E:\MAT2ANS_project\project1_isoAP100_spectrum';
Ansys_Path = fun_par.Ansys_solver_Path;
Ansys_Worksapce = fun_par.Ansys_Worksapce;
Project_APDL = fun_par.Project_APDL;
Jobname = 'MAT2ANS_project';
OutputFile = 'process';

cd(Ansys_Worksapce)

for iii = 1:Ns
    
    x = X(iii,:);
    
    % 
    fid = fopen([Ansys_Worksapce,'\','DataIn.txt'],'w');   
    fprintf(fid,'%f\n',x);
    fclose(fid);
       
    % ��¼��ǰ�ļ�
    
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
    
     % ���
     dataout = textread('DataOut.txt');
     G(iii,1) = fun_par.threshold-dataout;
               
end

cd(OriginFolder);

end

