clear
load_address=['D:\Research\rd_zdf\'];
wx_loop=[1.1 1.15 1.2 1.25 1.3 1.4 1.5];
share_loop=[1 2 5 10 15 20 25];
xi_loop=zeros(size(wx_loop));

figure(1);clf;hold on
ave_window=30;

color_loop=[linspace(1,0,length(wx_loop));...
    linspace(0.2,0.8,length(wx_loop));...
    linspace(0,1,length(wx_loop))]';


for ix=1:length(wx_loop)
    
%     figure(1);clf;hold on
%     for is=1:length(share_loop)
is=6;

        file_name=[num2str(wx_loop(ix)) '_5_100000_share' num2str(share_loop(is)) '_zdf_staff'];
        load([load_address file_name '.mat'])
        z=mean(z_all,1);
        g=mean(g_all,1);
        gm=zeros(size(g));
        for ii=1:length(g)
            gm(ii)=mean(g(max(1,ii-ave_window):min(length(g),ii+ave_window)));
        end
        
        [p_height,p_idx]=findpeaks(gm,'MinPeakDistance',800);
        p_height=[p_height inf];
        idx_cut=find(p_height(2:end)>p_height(1:end-1),1);% this value needs some modification or chosen manually
        p_idx=p_idx(1:idx_cut);
        P=polyfit(z(p_idx),log(gm(p_idx)-1),1);
        xi_loop(ix)=-1./P(1);
        
        plot(z,abs(gm-1),'.-','Color',color_loop(ix,:))
        plot(z(p_idx),abs(gm(p_idx)-1),'ko')
        plot(z,exp(polyval(P,z)),'-','Color',max(0,color_loop(ix,:)-0.2),'LineWidth',1.5)
%     end

end
set(gca,'yscale','log')
xlabel('z')
ylabel('|zdf-1|')

%%
figure(2);clf
plot(wx_loop-1,1./xi_loop,'o-')



