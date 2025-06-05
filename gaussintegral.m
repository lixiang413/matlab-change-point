function res = gaussintegral(sigma, XX,m)
    %m = 10;
    mu = 0;
    bound = 5;
    dz = 0.01;
    z = -bound:dz:bound;
    % Calculate the values of normcdf and normpdf to reduce repetitive calculations in the loop
    A = normcdf(z, mu, sigma);
    B = normpdf(z, mu, sigma);
    N_z = length(XX);

PP= length(z);
XX_matrix = repmat(XX, 1, PP);
z_matrix = repmat(z, length(XX), 1);
Int_n = normcdf(XX_matrix + z_matrix, mu, sigma) - A;
Int_n2 = Int_n.^(m-2);
pdf = normpdf(z_matrix + XX_matrix, mu, sigma);
all = dz*(Int_n2 .* B .* pdf);
LOGfx=log(sum(all,2));
Int_t_all=N_z*log(m*(m-1))+sum(LOGfx);
    % Sum Int_t_all to obtain the final result
    res = Int_t_all;
end