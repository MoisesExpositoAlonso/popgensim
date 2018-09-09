# Kondrasove
pk<-readRDS('figsfiles/kondrasov.rda')

# GWA
pgwa<-readRDS('figsfiles/gwa.rda')

pgwa
dev.off()
save_plot(plot=plot_grid(myplot,pgwa[[1]],pgwa[[2]],pgwa[[3]],ncol=2), filename = 'figs/kondrasov_gwa_fitness_oriented.pdf',base_height = 8,base_width = 8)
