vertices=34835;
coord=importdata('bunny_watermarking_similaritytransform_number1.off');
all_coord=coord.data(2:34836,:);


mod_x_coord=[all_coord(:,1)]; %ｘ座標を格納
mod_y_coord=[all_coord(:,2)]; %ｙ座標を格納
mod_z_coord=[all_coord(:,3)]; %ｚ座標を格納

mod_x_mean=sum(mod_x_coord)/vertices; %x座標の平均値
mod_y_mean=sum(mod_y_coord)/vertices; %y座標の平均値
mod_z_mean=sum(mod_z_coord)/vertices; %z座標の平均値

for n=1:vertices
     p(n,1)=sqrt(mod_x_coord(n,1)^2+mod_y_coord(n,1)^2+mod_z_coord(n,1)^2);
end

order=[1:vertices];
order=order';
order_p=[p order];
sort_p=sortrows(order_p);
cut_p=sort_p(1742:33093,1:2);

i=1;
for n=1:64
    B(1:489,n)=cut_p(i:i+488,1);
    i=i+489;
end

for n=1:64
    nor_B(:,n)=(B(:,n)-min(B(:,n)))/(max(B(:,n))-min(B(:,n)));
end

mean_nor_B=sum(nor_B)/489;

for n=1:64
    if mean_nor_B(1,n)<0.5
        w(1,n)=0;
    else
        w(1,n)=1;
    end
end