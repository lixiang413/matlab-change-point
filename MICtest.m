%%m=10, sample size n=50, cycle number N=1000,ERm10=3.07753; ERm20 = 3.73492, ERm50 = 4.49815; ERm30 = 4.08552;
%% Critical value n=100; Alpha = 8.626; n=200; alpha= 8.088;
% % ERm16 = 3.53198

n1=50;
n2=50;
n=(n1+n2)/2;
N=1000;
%The number of j is NJ
NJ=1;
j=1;
outDATA=zeros(N,4*NJ);
outDATAMIC=zeros(N,4*NJ);
m=10;
ERm=3.07753;

%%The number of cycles N=10 and the sample size n
sigma1=0.01;
sigma2=0.01;

for p=1:N
H1=sigma1*randn(n1,m);
DATA11=max(H1,[],2);
DATA12=min(H1,[],2);
DATA1=DATA11-DATA12;
H2=sigma2*randn(n2,m);
DATA21=max(H2,[],2);
DATA22=min(H2,[],2);
DATA2=DATA21-DATA22;
DATA=[DATA1;DATA2];

%%Null hypothesis SICH0, parameter estimation sigmahat
sigmahat=(1/ERm)*mean(DATA);
SICH0=-2*gaussintegral(sigmahat,DATA,m)+log(2*n);

MICH0=-2*gaussintegral(sigmahat,DATA,m)+log(2*n);


%%Alternative hypothesis SICH1, parameter estimation sigmahat
SICH1=zeros(2*n,1);
MICH1=zeros(2*n,1);
%%parfor i=2:(2*n-1)
parfor i=2:(2*n-1)
    A1=DATA(1:i);
    A2=DATA(i+1:2*n);
sigmahat1=(1/ERm)*mean(A1);
sigmahat2=(1/ERm)*mean(A2);
SICH1(i)=-2*(gaussintegral(sigmahat1,A1,m)+gaussjifentisheng(sigmahat2,A2,m))+2*log(2*n);

MICH1(i)=-2*(gaussintegral(sigmahat1,A1,m)+gaussjifentisheng(sigmahat2,A2,m))+(2+((2*i)/(2*n)-1)^2)*log(2*n);

end 

%%SICH1k=min(SICH1(2:2*n-1));
[value,location]=min(SICH1(2:2*n-1));

outDATA(p,4*j-2)=value;
outDATA(p,4*j-1)=location+1;
outDATA(p,4*j-3)=SICH0;
logical_array = SICH0-value +log(2*n) > 10.3663646423061;
%logical_array = outDATA(p,4*j-2) +8.626 < outDATA(p,4*j-3);
result_array = zeros(1);
result_array(logical_array) = 1;
outDATA(p,4*j)=result_array;

[value,location]=min(MICH1(2:2*n-1));

outDATAMIC(p,4*j-2)=value;
outDATAMIC(p,4*j-1)=location+1;
outDATAMIC(p,4*j-3)=MICH0;
logical_array = MICH0 - value +log(2*n) > 8.47645432305867;
%logical_array = outDATAMIC(p,4*j-2) +2.5 < outDATAMIC(p,4*j-3);
result_array = zeros(1);
result_array(logical_array) = 1;
outDATAMIC(p,4*j)=result_array;

p

end