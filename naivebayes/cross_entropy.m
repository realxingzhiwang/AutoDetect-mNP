function [H12,H21] = cross_entropy(data1,data2)
%Cross entropy between two distributions
%   [H12,H21] = cross_entropy(data1,data2) computes cross entropies between
%   data1 and data2. H12 = -sum(P11(i)*log(P12(i)), H21 =
%   -sum(P22(i)*log(P21(i)).

mu1 = mean(data1);
mu2 = mean(data2);
sigma1 = diag(sqrt(var(data1)));
sigma2 = diag(sqrt(var(data2)));

P11 = mvnpdf(data1, mu1, sigma1);
P12 = mvnpdf(data1, mu2, sigma2);
P22 = mvnpdf(data2, mu2, sigma2);
P21 = mvnpdf(data2, mu1, sigma1);

H12 = -P11(P11>1e-3)'*log(P12(P11>1e-3));
H21 = -P22(P22>1e-3)'*log(P21(P22>1e-3));

%H12 = -log(P11(P11>1e-3))'*P12(P11>1e-3);
%H21 = -log(P22(P22>1e-3))'*P21(P22>1e-3);




end

