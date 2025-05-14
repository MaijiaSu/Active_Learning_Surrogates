function G = MResponceMDF_OS_smj(CandidatePool,fun_par)

Ns = size(CandidatePool,1);
for jj = 1:Ns
    %% uncertainty variables
    u = CandidatePool(jj,:);
    %%
    fun_par = [];
    t = 0.01:0.01:10;
    Mat.t = t;
    x0=1/6*(u(1)*sin(Mat.t*2*pi/0.125)+u(2)*cos(Mat.t*2*pi/0.25)+u(3)*sin(Mat.t*2*pi/0.5)+u(4)*sin(Mat.t*2*pi/1));
    fid1 = fopen(['accel.txt'],'w');
    fprintf(fid1,'%6.10f\n ',x0');
    fclose(fid1);
    !OpenSees.exe "M".tcl
    HistVar.u(1,:) = load(['output/disp_1.out']);
    HistVar.u(2,:) = load(['output/disp_2.out']);
    HistVar.u(3,:) = load(['output/disp_3.out']);
    G(jj,1) = max(max(abs(HistVar.u)));
    if mod(jj,100)==0
        disp(['Running:MResponceMDF_OS_smj--',num2str(jj),'/',num2str(Ns)])
    end
end

end