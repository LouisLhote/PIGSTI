# Rscript sexing_residual_method_ovis.R input.idx output.pdf output.tsv

args = commandArgs(trailingOnly=TRUE)

#header should be "chr length sample"

#take the idx file
d<-read.table(args[1],header=T)
pdf(args[2])

critt<-qt(0.05,df=24)
probs<-0
ratios<-0
calls<-character(length(colnames(d)) - 2)

for( i in 3:length(colnames(d)) )
{
	# Start with plot
	plot(d[c(1:27),i]~d$length[c(1:27)],xlab="Chromsome length (Mb)",ylab="Number of reads", main=colnames(d)[i],pch=16,col="slategrey",xaxt="n")
	points(d$length[27],d[27,i],col="red",pch=16)
	text(d$length[27],d[27,i],"X",col="red",pos=1)
	text(d$length[1],d[1,i],"1",col="slategrey",pos=1, cex=0.5)
	text(d$length[2],d[2,i],"2",col="slategrey",pos=1, cex=0.5)
	text(d$length[3],d[3,i],"3",col="slategrey",pos=1, cex=0.5)
	text(d$length[4],d[4,i],"4",col="slategrey",pos=1, cex=0.5)
	text(d$length[5],d[5,i],"5",col="slategrey",pos=1, cex=0.5)
	text(d$length[6],d[6,i],"6",col="slategrey",pos=1, cex=0.5)
	text(d$length[7],d[7,i],"7",col="slategrey",pos=1, cex=0.5)
	text(d$length[8],d[8,i],"8",col="slategrey",pos=1, cex=0.5)
	text(d$length[9],d[9,i],"9",col="slategrey",pos=1, cex=0.5)
	text(d$length[10],d[10,i],"10",col="slategrey",pos=1, cex=0.5)
	text(d$length[11],d[11,i],"11",col="slategrey",pos=1, cex=0.5)
	text(d$length[12],d[12,i],"12",col="slategrey",pos=1, cex=0.5)
	text(d$length[13],d[13,i],"13",col="slategrey",pos=1, cex=0.5)
	text(d$length[14],d[14,i],"14",col="slategrey",pos=1, cex=0.5)
	text(d$length[15],d[15,i],"15",col="slategrey",pos=1, cex=0.5)
	text(d$length[16],d[16,i],"16",col="slategrey",pos=1, cex=0.5)
	text(d$length[17],d[17,i],"17",col="slategrey",pos=1, cex=0.5)
	text(d$length[18],d[18,i],"18",col="slategrey",pos=1, cex=0.5)
	text(d$length[19],d[19,i],"19",col="slategrey",pos=1, cex=0.5)
	text(d$length[20],d[20,i],"20",col="slategrey",pos=1, cex=0.5)
	text(d$length[21],d[21,i],"21",col="slategrey",pos=1, cex=0.5)
	text(d$length[22],d[22,i],"22",col="slategrey",pos=1, cex=0.5)
	text(d$length[23],d[23,i],"23",col="slategrey",pos=1, cex=0.5)
	text(d$length[24],d[24,i],"24",col="slategrey",pos=1, cex=0.5)
	text(d$length[25],d[25,i],"25",col="slategrey",pos=1, cex=0.5)
	text(d$length[26],d[26,i],"26",col="slategrey",pos=1, cex=0.5)
	points(d$length[27],(2*d[27,i]),col="lightgrey",pch=16)
	text(d$length[27],(2*d[27,i]),"2X",col="lightgrey",pos=3)
	axis(side=1,at=axTicks(1),labels=axTicks(1)/1000000)

	# Fit linear regression
	model<-lm(d[c(1:26),i]~d$length[c(1:26)])
	abline(model)

	# 95% CI lines for regression line
	prd <- predict( model, data.frame(d$length[c(1:26)]), interval="confidence", level=0.95)
	confidence <- data.frame( cbind( d$length[c(1:26)], prd[,2], prd[,3] ) )
	confordered <- confidence[order(confidence[,1]),]
	lines(confordered$X1,confordered$X2,lty=3)
	lines(confordered$X1,confordered$X3,lty=3)

	# Externally Studentized estimate of variance of residuals
	studres<-rstudent(model)
	modelwithX<-lm(d[c(1:27),i]~d$length[c(1:27)])
	studreswithX<-rstudent(modelwithX)
	femaleprob<-pt(studreswithX[27],df=24)
	modelwithtwiceX<-lm(c(d[c(1:26),i],(d[27,i]*2))~d$length[c(1:27)])
	studreswithtwiceX<-rstudent(modelwithtwiceX)
	maleprob<-pt(studreswithtwiceX[27],df=24,lower.tail=F)
	ratio<-femaleprob/maleprob
	probs[i-2]<-femaleprob

#	R2<-format(summary(model)$adj.r.squared, digits=4)


#	lm_coef <- round(coef(model), 10) 
	#extract coefficients
#	mtext(bquote(y == .(lm_coef[2])*x + .(lm_coef[1])), adj=1, padj=0) # display equation

#	m_coef <- round(coef(model), 10) 
	#extract coefficients
#	eq = mtext(bquote(y == .(lm_coef[2])*x + .(lm_coef[1]))

 # 	eq <- substitute(italic(d[c(1:22),i]) == a + b %.% italic(d$length[c(1:22)])*", 
  #       list(a = format(coef(model)[1], digits = 2), 
   #           b = format(coef(model)[2], digits = 2))



	lm_eqn = function(df){
   		 model<-lm(d[c(1:26),i]~d$length[c(1:26)]);
   		 eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
        	 list(a = format(coef(m)[1], digits = 2), 
           	  b = format(coef(m)[2], digits = 2), 
          	 r2 = format(summary(m)$r.squared, digits = 3)))
    	as.character(as.expression(eq));    
	}

	# Add descriptive legend
	if( ratio >= 1 )
	{
		lrt<--2*log(maleprob/femaleprob)		# <- stab at likelihood ratio test
		lrtp<-1-(pchisq(lrt,df=1))			# <- stab at computing probability based on LRT
		legend( "topleft",c(
					paste("Female (",signif(ratio,digits=3),"X more likely)",sep=""),
					paste("p = ", signif(lrtp,digits=3))
				   ), cex = 0.6
		      )
		ratios[i-2]<-ratio
		calls[i-2]<-"Female"
	}
	else
	{
		lrt<--2*log(femaleprob/maleprob)		# <- stab at likelihood ratio test
		lrtp<-1-(pchisq(lrt,df=1))			# <- stab at computing probability based on LRT
		legend( "topleft",c(
					paste("Male (",signif((1/ratio),digits=3),"X more likely)",sep=""),
					paste("p = ", signif(lrtp,digits=3))
				   ), cex = 0.6
		      )
		ratios[i-2]<-1/ratio
		calls[i-2]<-"Male"
	}
}

dev.off()

if (length(args) >= 3) {
	samples <- colnames(d)[3:ncol(d)]
	outdf <- data.frame(
		sample=samples,
		call=calls,
		female_prob=probs,
		likelihood_ratio=ratios,
		stringsAsFactors=FALSE
	)
	write.table(outdf, args[3], sep="\t", row.names=FALSE, quote=FALSE)
}
