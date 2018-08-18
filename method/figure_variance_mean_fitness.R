## In file data-cleaning/Yfitness.R
Y<-readRDS('data/Y.rds')

# pdf("figs/mean_to_sd_ratio.pdf",useDingbats = FALSE,width = 5,height = 5)
yforvariance=Y
toplot<-data.frame(y=sqrt(Vy(yforvariance$oFitness,yforvariance$row)),
       x=My(yforvariance$oFitness,yforvariance$row)
       )
p1<-ggdotscolor(toplot$y,
						x=toplot$x) %>%
	addggregression(se=F) + ylim(c(0,2))+
	ylab("Standard Deviation of fitness")+
	xlab("Mean fitness")

yforvariance=filter(Y, oFitness!=0)
toplot<-data.frame(y=sqrt(Vy(yforvariance$oFitness,yforvariance$row)),
       x=My(yforvariance$oFitness,yforvariance$row)
       )
p2<-ggdotscolor(toplot$y,
						x=toplot$x) %>%
  addggregression(se=F) + ylim(c(0,2))+
	ylab("Standard Deviation of fitness")+
	xlab("Mean fitness")
# dev.off()


save_plot(filename = "figs/mean_to_sd_ratio.pdf",useDingbats = FALSE,base_width = 10,base_height = 5,
          plot=plot_grid(p1,p2,ncol = 2,labels=c('a','b')))
