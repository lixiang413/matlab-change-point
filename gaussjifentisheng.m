function res = gaussjifentisheng(sigma, XX,m)
    %m = 10;
    mu = 0;
    bound = 10;
    dz = 0.01;
    z = -bound:dz:bound;
    % 计算 normcdf 和 normpdf 的值，减少循环中的重复计算
    A = normcdf(z, mu, sigma);
    B = normpdf(z, mu, sigma);
    % 为了减少内存分配，提前预分配 Int_t_all
    N_z = length(XX);
   % Int_t_all = zeros(N_z, 1);

%     for i_z = 1:N_z
%         % 计算 normcdf 的值
%         Int_n = normcdf(XX(i_z) + z, mu, sigma) - A;
%         
%         % 计算 Int_n2 的值
%         Int_n2 = Int_n.^(m - 2);
% 
%         % 计算 pdf 的值
%         pdf = normpdf(z + XX(i_z), mu, sigma);
% 
%         % 计算 all 的值
%         all = Int_n2 .* B .* pdf;
% 
%         % 避免使用 trapz 函数，可以通过矩阵乘积和 sum 函数来计算积分
%         Int_t_all(i_z) = log(m * (m - 1) * dz * sum(all));
%     end
PP= length(z);
XX_matrix = repmat(XX, 1, PP);
z_matrix = repmat(z, length(XX), 1);
Int_n = normcdf(XX_matrix + z_matrix, mu, sigma) - A;
Int_n2 = Int_n.^(m-2);
pdf = normpdf(z_matrix + XX_matrix, mu, sigma);
all = dz*(Int_n2 .* B .* pdf);
LOGfx=log(sum(all,2));
Int_t_all=N_z*log(m*(m-1))+sum(LOGfx);
    % 对 Int_t_all 求和得到最终结果
    res = Int_t_all;
end