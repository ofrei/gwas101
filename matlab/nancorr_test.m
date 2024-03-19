X = randn(30,4);
Y = randn(30,4);
 
Y(:,4) = Y(:,4)+X(:,2);
Y(1,2)= inf;
 
[rho,pval] = corr(X,Y);
 
[coef, t, n] = nancorr(X,Y);
    
a=py.nancorr.nancorry(py.numpy.array(X),py.numpy.array(Y));
 
coefpy=double(a{1});
tpy=double(a{2});
npy=double(a{3});