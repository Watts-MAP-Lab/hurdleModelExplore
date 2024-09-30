# true coefficient
bbeta = π;
# random x data
x = randn(100,1)*5 .+ 3;
# OLS data generating process
y = bbeta.*x .+ randn(100,1)*10;
# OLS estimation
bbeta_hat = inv(x'x)x'y;
println("β-hat is $(round(bbeta_hat[1],digits=3)) and the true β is $(round(bbeta,digits=3)).")