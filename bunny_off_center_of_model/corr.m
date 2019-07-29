ext_mean=sum(w)/64;
ori_mean=sum(embedding_watermarking)/64;

sigma_bunsi=0;
sigma_bunbo_hidari=0;
sigma_bunbo_migi=0;
for i=1:64
    sigma_bunsi=sigma_bunsi+(w(1,i)-ext_mean)*(embedding_watermarking(1,i)-ori_mean);
    sigma_bunbo_hidari=sigma_bunbo_hidari+(w(1,i)-ext_mean).^2;
    sigma_bunbo_migi=sigma_bunbo_migi+(embedding_watermarking(1,i)-ori_mean).^2;
end

correlation=sigma_bunsi/sqrt(sigma_bunbo_hidari*sigma_bunbo_migi);
