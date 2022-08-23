%cacluate the short wx

save_address='./save_addr/1.16/';
wx=1.16;wy=5;N=20000;
xs=[];dxs=[];
de42=[];dz42=[];dy42=[];de30=[];dz30=[];dy30=[];de56=[];dz56=[];dy56=[];
pll30=[];pll42=[];pll56=[];pln30=[];pln42=[];pln56=[];
pl30m=[];pl30a=[];pl42m=[];pl42a=[];pl56m=[];pl56a=[];
tl30=[];tl42=[];tl56=[];tn30=[];tn42=[];tn56=[  ];

for irep=394:1000
    try
        irep
        [Rc,tole_degree_loop]=Sequential_deposition3D_(wx,wy,N);
        wall_y=wy/2;wall_x=wx/2;
        file_name=[num2str(wx) '_' num2str(wy) '_' num2str(N) '_' num2str(irep)];
        save([save_address file_name '.mat'],'Rc','tole_degree_loop','wx','wy')
        Rc_old=Rc;
        Period_random_deposition%%%
        file_name=[num2str(wx) '_' num2str(wy) '_period_staff_' num2str(irep)];
        save([save_address file_name '.mat'],'x_save','dx_save','de42t','dz42t','dy42t','de30t','dz30t'...
            ,'dy30t','de56t','dz56t',['d' .../
            'y56t'],'pll30t','pll42t','pll56t','pln30t','pln42t','pln56t','pl30mt','pl30at'...
            ,'pl42mt','pl42at','pl56mt','pl56at','tl30t','tl42t','tl56t','tn30t','tn42t','tn56t')

        xs=[xs x_save];dxs=[dxs dx_save];
        de42=[de42 de42t];dz42=[dz42 dz42t];dy42=[dy42 dy42t];
        de30=[de30 de30t];dz30=[dz30 dz30t];dy30=[dy30 dy30t];
        de56=[de56 de56t];dz56=[dz56 dz56t];dy56=[dy56 dy56t];
        pll30=[pll30 pll30t];pll42=[pll42 pll42t];pll56=[pll56 pll56t];
        pln30=[pln30 pln30t];pln42=[pln42 pln42t]; pln56=[pln56 pln56t];
        pl30m=[pl30m pl30mt];pl30a=[pl30a pl30at];pl42m=[pl42m pl42mt];
        pl42a=[pl42a pl42at];pl56m=[pl56m pl56mt];pl56a=[pl56a pl56at];
        tl30=[tl30 tl30t];tl42=[tl42 tl42t];tl56=[tl56 tl56t];
        tn30=[tn30 tn30t];tn42=[tn42 tn42t];tn56=[tn56 tn56t];
    catch
        irep
        [Rc,tole_degree_loop]=Sequential_deposition3D_(wx,wy,N);
        wall_y=wy/2;wall_x=wx/2;
        file_name=[num2str(wx) '_' num2str(wy) '_' num2str(N) '_' num2str(irep)];
        save([save_address file_name '.mat'],'Rc','tole_degree_loop','wx','wy')
        Rc_old=Rc;
        Period_random_deposition%%%
        file_name=[num2str(wx) '_' num2str(wy) '_period_staff_' num2str(irep)];
        save([save_address file_name '.mat'],'x_save','dx_save','de42t','dz42t','dy42t','de30t','dz30t'...
            ,'dy30t','de56t','dz56t',['d' .../
            'y56t'],'pll30t','pll42t','pll56t','pln30t','pln42t','pln56t','pl30mt','pl30at'...
            ,'pl42mt','pl42at','pl56mt','pl56at','tl30t','tl42t','tl56t','tn30t','tn42t','tn56t')

        xs=[xs x_save];dxs=[dxs dx_save];
        de42=[de42 de42t];dz42=[dz42 dz42t];dy42=[dy42 dy42t];
        de30=[de30 de30t];dz30=[dz30 dz30t];dy30=[dy30 dy30t];
        de56=[de56 de56t];dz56=[dz56 dz56t];dy56=[dy56 dy56t];
        pll30=[pll30 pll30t];pll42=[pll42 pll42t];pll56=[pll56 pll56t];
        pln30=[pln30 pln30t];pln42=[pln42 pln42t]; pln56=[pln56 pln56t];
        pl30m=[pl30m pl30mt];pl30a=[pl30a pl30at];pl42m=[pl42m pl42mt];
        pl42a=[pl42a pl42at];pl56m=[pl56m pl56mt];pl56a=[pl56a pl56at];
        tl30=[tl30 tl30t];tl42=[tl42 tl42t];tl56=[tl56 tl56t];
        tn30=[tn30 tn30t];tn42=[tn42 tn42t];tn56=[tn56 tn56t];
    end
end
file_name=[num2str(wx) '_' num2str(wy) '_all'];
save([save_address file_name '.mat'],'xs','dxs','de42','dz42','dy42','de30','dz30'...
    ,'dy30','de56','dz56','dy56','pll30','pll42','pll56','pln30','pln42','pln56','pl30m','pl30a'...
    ,'pl42m','pl42a','pl56m','pl56a','tl30','tl42','tl56','tn30','tn42','tn56')
