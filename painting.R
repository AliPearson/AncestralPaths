paint<-function(ints, sample, main=""){
	chunks<-t(ints[sample,,])
	start_path<-1
	path_counter<-0
	rows_remove<-c()
	while(start_path+path_counter <= nrow(chunks)){
		if(chunks[start_path,3]!=chunks[(start_path+path_counter),3]){
			if(path_counter==1){
				start_path<-start_path+path_counter
				path_counter=0
			}else{
				rows_remove<-c(rows_remove, start_path:(start_path+(path_counter-2)))
				start_path<-start_path+path_counter
				path_counter=0
			}
		}else{
			path_counter<-path_counter+1
		}
	}
	rows_remove<-c(rows_remove, start_path:(nrow(chunks)-1))
	chunks<-chunks[-rows_remove,]

	palette(c("red", "purple", "black", "orange", "cyan", "magenta", "grey"))
	#par(mar=c(2,2,2,2))
	if(is.null(dim(chunks))){
		plot(0, type="n",xlim=c(0,chunks[1]), ylim=c(0,5), bty='L', yaxt="n", yaxs = "i", xaxs="i", main=main,xlab="Position (bp)")
		rect(0,0,chunks[1],5,col=chunks[3], border=NA)
	}else{
		plot(0, type="n",xlim=c(0,chunks[nrow(chunks),4]), ylim=c(0,5), bty='L', yaxt="n", yaxs = "i", xaxs="i", main=main,xlab="Position (bp)")
		for (i in 1:nrow(chunks)){
			if (i==1){
				rect(0,0,chunks[i,4],5,col=chunks[i,3], border=NA)
			}else{
				rect(chunks[(i-1),4],0,chunks[i,4],5, col=chunks[i,3],border=NA)
			}
		}

	}
	legend("topright",5, fill=1:6, legend=1:6, cex=0.75, inset=c(1,0), xpd=TRUE)
	return(chunks)
}
