clear all
close all
%%
output_point_fl = 'C:\Users\moritzw\OneDrive - SPC\Documents\Projects\ECIKS\Quick_Flood_Warning\extras\output_points.csv';
pts = readtable(output_point_fl);
pts_no = pts.point_no;
%%
cnt = 0;
figure()
hold on
for i = 1:8
for j=1:4
    cnt = cnt+1;
    plot(pts.lon(cnt),pts.lat(cnt),'kx')
    label = [pts.atoll{cnt} ' ' num2str(pts.point_no(cnt))];
    text(pts.lon(cnt),pts.lat(cnt),label)
    
    
end
end