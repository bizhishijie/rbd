save_address=['D:\Research\rbd_contact\'];
wx_loop=[1.02 1.04 1.06 1.08 1.1:0.05:1.5];
load criticalvalue.mat
wy=5;N=100000;N_loop=100000;
for in=1:length(N_loop)
    Nt=N_loop(in);
    %     for ix=1:length(wx_loop)
    for ix=10% 10:13
        wx=wx_loop(ix)
        eps_multpile=eps_critical_value_loop(ix);
        cal_s_all=zeros(1,6);
        for irep=1:500
            irep
            file_name=[num2str(wx) '_' num2str(wy) '_' num2str(N) '_' num2str(irep)];
            [Rc,tole_degree_loop]=Sequential_deposition3D_(wx,wy,N);
            save([save_address file_name '.mat'],'Rc','tole_degree_loop','wx','wy')
            
            tmp_all=touching_n_(Rc(:,1:Nt),Rc(:,1:Nt),eps_multpile);
            cal_s_all(irep,1:6)=tmp_all./sum(tmp_all);
            
            file_name_save=[num2str(wx) '_' num2str(wy) '_' num2str(N) '_' num2str(Nt) 'contractn'];
            save([save_address file_name_save '.mat'],'cal_s_all')
        end
    end
end
