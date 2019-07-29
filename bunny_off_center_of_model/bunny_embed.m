vertices=34835;
alpha=0.3
%埋め込み透かしを定義しておく
embedding_watermarking=[1 0 1 0 1 1 1 1 0 1 1 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 1 1 0 1 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 1 1];


coord=importdata('bunny.off'," ",2);
all_coord=coord.data(1:34835,:);

%(a)calculate and sort eigenvalues.
x_coord=[all_coord(:,1)]; %ｘ座標を格納
y_coord=[all_coord(:,2)]; %ｙ座標を格納
z_coord=[all_coord(:,3)]; %ｚ座標を格納

x_mean=sum(x_coord)/vertices; %x座標の平均値
y_mean=sum(y_coord)/vertices; %y座標の平均値
z_mean=sum(z_coord)/vertices; %z座標の平均値

for n=1:vertices
    p(n,1)=sqrt(x_coord(n,1)^2+y_coord(n,1)^2+z_coord(n,1)^2);
    psi(n,1)=acos(z_coord(n,1)/sqrt(x_coord(n,1)^2+y_coord(n,1)^2+z_coord(n,1)^2));
end

for n=1:vertices
    if y_coord(n,1)>=0
        sgn_y=1;
     theta(n,1)=sgn_y*acos(x_coord(n,1)/sqrt(x_coord(n,1)^2+y_coord(n,1)^2));
    else
        sgn_y=-1;
        theta(n,1)=sgn_y*acos(x_coord(n,1)/sqrt(x_coord(n,1)^2+y_coord(n,1)^2));
    end
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


copy_nor_B=nor_B; %copy_Bは不変の値として、Bの値を変えて行くことにする
mean_copy_nor_B=sum(copy_nor_B)/489;
for i=1:64
    if embedding_watermarking(1,i)==0
            if mean_nor_B(1,i)<0.5-alpha
                mod_B(:,i)=nor_B(:,i);
            else
                n=1;
                while mean_nor_B(1,i) >= 0.5-alpha
                    for j=1:489
                        if copy_nor_B(j,i)>=0 && copy_nor_B(j,i)<0.3
                            mod_B(j,i)=copy_nor_B(j,i).^(1+n*0.002);
                        elseif copy_nor_B(j,i)>=0.3 && copy_nor_B(j,i)<0.9
                            mod_B(j,i)=copy_nor_B(j,i).^(1+n*0.001);
                        else
                            mod_B(j,i)=copy_nor_B(j,i).^(1+n*0.0015);
                        end
                    end
                    mean_nor_B(1,i)=sum(mod_B(:,i))/489;
                    n=n+1;
                end 
            end
    elseif embedding_watermarking(1,i)==1
        if mean_nor_B(1,i)>0.5+alpha
            mod_B(:,i)=nor_B(:,i);
        else
            n=1;
            while mean_nor_B(1,i) <= 0.5+alpha
                for j=1:489
                    if copy_nor_B(j,i)>=0 && copy_nor_B(j,i)<0.6
                            mod_B(j,i)=copy_nor_B(j,i).^(1-n*0.001);
                        elseif copy_nor_B(j,i)>=0.6 && copy_nor_B(j,i)<0.8
                            mod_B(j,i)=copy_nor_B(j,i).^(1-n*0.0015);
                        else
                            copy_nor_B(j,i)>=0.8 && copy_nor_B(j,i)<1
                            mod_B(j,i)=copy_nor_B(j,i).^(1-n*0.002);
                        end
                    end
                    mean_nor_B(1,i)=sum(mod_B(:,i))/489;
                    n=n+1;
            end
        end
    end
end

for n=1:64
    inv_pmj(:,n)=mod_B(:,n)*(max(B(:,n))-min(B(:,n)))+min(B(:,n));
end

i=1;
for n=1:64
    vertical_inv_pmj(i:i+488,1)=inv_pmj(1:489,n);
    i=i+489;
end

final_p=sort_p;
final_p(1742:33037)=vertical_inv_pmj;
final_p_sort=sortrows(final_p,2);


for n=1:34835
    mod_x(n,1)=final_p_sort(n,1)*cos(theta(n,1))*sin(psi(n,1));
    mod_y(n,1)=final_p_sort(n,1)*sin(theta(n,1))*sin(psi(n,1));
    mod_z(n,1)=final_p_sort(n,1)*cos(psi(n,1));
end

mod_coord=[mod_x mod_y mod_z];


total=0;
 for n=1:489
 total=total+sum(abs(mod_B(n,1)-nor_B(n,1)));
 end
