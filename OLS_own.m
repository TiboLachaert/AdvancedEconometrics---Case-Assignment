function out = OLS_own(y,x)

%% preamble
n  = size(y,1);
k  = size(x,2);
df = n-k;

%% OLS estimation
xy     = x'*y;
xxi    = (x'*x)^(-1);
coefs  = xxi*xy;

yhat   = x*coefs;
res    = y-yhat;
sigma2 = (res'*res)/df;

stdvs  = sqrt(sigma2)*sqrt(diag(xxi));
tstats = coefs./stdvs;
pvals  = 2*(1-tcdf(abs(tstats),df));

%% Save output
out.estimates = [coefs,stdvs,tstats,pvals];

end
