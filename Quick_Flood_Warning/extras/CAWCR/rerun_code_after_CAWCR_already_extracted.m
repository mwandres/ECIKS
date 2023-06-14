clear all
close all
%%
%output_point_fl = 'C:\Users\moritzw\OneDrive - SPC\Documents\Projects\ECIKS\Quick_Flood_Warning\extras\output_points.csv';
output_point_fl = '\\GSD-MORITZW-PC\D_Mead_Rd\ECIKS\Quick_Flood_Warning\extras\output_points.csv';
pts = readtable(output_point_fl);
pts_no = pts.point_no;

cawcr_folder = '\\spc\shared\GEM\FJ_NAB\O_M\Oceans\Original_Sources\PACCSAP_WAVE_HINDCAST_DATA\GRIDDED\GLOB\';
ORAS5_folder = 'C:\Users\moritzw\OneDrive - SPC\Documents\Projects\ECIKS\Quick_Flood_Warning\extras\ORAS5\';
ARI_1X = [];
ARI_2X = [];
ARI_5X = [];
ARI_10X= [];

for i = 1:length(pts_no)
    target_lon = pts.lon(i);
    target_lat = pts.lat(i);
    out_fl_name = [pts.atoll{i} '_' num2str(pts.point_no(i)) '_waves.csv'];
    disp(out_fl_name)
    tab = readtable(out_fl_name);
    time_all = tab.Var1;
    hs_all = tab.Var2;
    t0m1_all = tab.Var3;
    dp_all = tab.Var4;
    
    %writematrix([time_all, hs_all, t0m1_all, dp_all],out_fl_name)
    oras5_data = readtable([ORAS5_folder pts.atoll{i} '_sossheig.txt']);
    oras5_data.dates(1) = datetime(time_all(1),'ConvertFrom','datenum');
    oras5_data.dates(end) = datetime(time_all(end),'ConvertFrom','datenum');
    oras5_hourly = interp1(datenum(oras5_data.dates),oras5_data.sossheig,time_all);
    [tide_hourly,conList]=tmd_tide_pred('C:\Users\moritzw\OneDrive - SPC\Documents\MATLAB\TMD_Matlab_Toolbox_v2.5\TMD_Matlab_Toolbox_v2.5\TMD\DATA\Model_atlas',time_all,target_lat,target_lon,'z');

    %% Merrifield
    gamma = 1.0;
    b1 = 0.33;
    b0 = -0.1;
    theta_N = pts.shore_normal_angle(i);
    Hs = hs_all;
    Tp = t0m1_all;
    angles = cosd(dp_all - theta_N);
    angles(angles<0) = 0;
    Hb = (Hs.^2.* Tp * (4*pi)^(-1).* angles.* sqrt(gamma.*9.81)).^(2/5);
    runup = b1 * Hb + (b0);
    twl = runup + tide_hourly + oras5_hourly;
    %% EVA
    
    ixx = find(isnan(twl));
    twl(ixx) = min(twl);
    lambda = length(twl)/40;
    [zor,zor_ix] = sort(twl);%real
    Promedioz = (((1:length(zor)))/(length(zor)+1))';
    Tamedioz = 1./(1-Promedioz);
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 1 scrsz(3) scrsz(4).*.9],'Visible','On')
    semilogx(Tamedioz/lambda,zor,'o','MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',2);
    xlim([0.05 50])
    xlabel('ARI (years)')
    ylabel('Estimated TWL nearshore')
    fileprint=[pts.atoll{i} '_' num2str(pts.point_no(i)) '_twl_EVA.png'];
    print('-dpng','-r200',fileprint)
    close
    
	ARI_0p5_ix = find(abs(Tamedioz/lambda-0.5) == min(abs(Tamedioz/lambda-0.5)));
    ARI_0p5 = zor(ARI_0p5_ix);
    ARI_1_ix = find(abs(Tamedioz/lambda-1) == min(abs(Tamedioz/lambda-1)));
    ARI_1 = zor(ARI_1_ix);
    ARI_2_ix = find(abs(Tamedioz/lambda-2) == min(abs(Tamedioz/lambda-2)));
    ARI_2 = zor(ARI_2_ix);
    ARI_5_ix = find(abs(Tamedioz/lambda-5) == min(abs(Tamedioz/lambda-5)));
    ARI_5 = zor(ARI_5_ix);
    ARI_10_ix = find(abs(Tamedioz/lambda-10) == min(abs(Tamedioz/lambda-10)));
    ARI_10 = zor(ARI_10_ix);
    
    fileprint_twl_csv=[pts.atoll{i} '_' num2str(pts.point_no(i)) '_twl_EVA.csv'];
    writematrix([1,2,5,10;ARI_1,ARI_2,ARI_5,ARI_10],fileprint_twl_csv)
    ARI_1X = [ARI_1X; ARI_1]; 
    ARI_2X = [ARI_2X; ARI_2];
    ARI_5X = [ARI_5X; ARI_5];
    ARI_10X= [ARI_10X; ARI_10];
    
    pts.threshold_1(i) = ARI_0p5;
    pts.threshold_2(i) = ARI_5;

end
fileprint_twl_csv=['All_Islands' '_twl_EVA.csv'];
writematrix([1,2,5,10;ARI_1X,ARI_2X,ARI_5X,ARI_10X],fileprint_twl_csv)


writetable(pts,'\\GSD-MORITZW-PC\D_Mead_Rd\ECIKS\Quick_Flood_Warning\extras\output_points_new.csv')


