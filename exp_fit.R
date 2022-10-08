#!/usr/bin/env Rscript

exp_fit<-function(res_chrom, path){
	
	std<-apply(res_chrom[,path,path,]/apply(res_chrom[,path,,], 3, rowSums), 1, function(x) var(x, na.rm=TRUE))
	res<-rowSums(res_chrom, dims=3)
	y<-(res[,path,path]/apply(res[,path,], 1, sum))	
	#print(y)
	std[which(std==0)]=1e-10
	#print(std)
	
	x=length(y)	
	data.df<-data.frame(x=((1:x)/100), y=y, std=std)
	theta.0 <- min(data.df$y, na.rm=TRUE) * 0.5
	model.0 <- lm(log(y - theta.0) ~ x, data=data.df, na.action=na.exclude)  
	alpha.0 <- exp(coef(model.0)[1])
	beta.0 <- coef(model.0)[2]
	start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
		  #nlc <- nls.control(maxiter = 2000, minFactor=0.000244140, warnOnly=T)
	nlc <- nls.control(maxiter = 5000, warnOnly=T)
	model <- nls(y ~ (alpha)*exp(beta * x) + theta , data = data.df, start = start, algorithm="port", control=nlc, weights=(1/std))
	estimate_beta<-summary(model)$parameters[2,1]
	estimate_theta<-summary(model)$parameters[3,1]
	std_error_beta<-summary(model)$parameters[2,2]
	std_error_theta<-summary(model)$parameters[3,2]
	rm(model)
	rm(model.0)
	return(c(estimate_beta, std_error_beta, estimate_theta, std_error_theta))
}
