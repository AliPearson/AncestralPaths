#!/usr/bin/env Rscript

exp_fit<-function(res_chrom, path1, path2){
	res_comb<-res_chrom[,path1,,]+res_chrom[,path2,,]
	std<-apply((res_comb[,path1,]+res_comb[,path2,])/apply(res_comb, 3, rowSums), 1, function(x) var(x, na.rm=TRUE))
	res<-rowSums(res_chrom, dims=3)
	res_comb<-res[,path1,]+res[,path2,]
	y<-(res_comb[,path1]+res_comb[,path2])/apply(res_comb,1,sum)
	#print(y)
	std[which(std==0)]=1e-10
	#print(std)
	
	x=length(y)
	#print(x)
	data.df<-data.frame(x=((1:x)/100), y=y, std=std)
	theta.0 <- min(data.df$y, na.rm=TRUE) * 0.5
	model.0 <- lm(log(y - theta.0) ~ x, data=data.df, na.action=na.exclude)  
	alpha.0 <- exp(coef(model.0)[1])
	beta.0 <- coef(model.0)[2]
	start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
	nlc <- nls.control(maxiter = 5000, warnOnly=T)
	model <- nls(y ~ (alpha)*exp(beta * x) + theta , data = data.df, start = start, control=nlc, algorithm="port", weights=(1/std))
	estimate_beta<-as.numeric(coef(model)[2])
	estimate_theta<-as.numeric(coef(model)[3])
	std_error_beta<-summary(model)$parameters[2,2]
	std_error_theta<-summary(model)$parameters[3,2]
	rm(model)
	rm(model.0)
	return(c(estimate_beta, std_error_beta, estimate_theta, std_error_theta))
}
