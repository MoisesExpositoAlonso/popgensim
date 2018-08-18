# argumentsgws<-function(fitness='rFitness',extrafitness='Fitness',exp='mlp',mod='lm', thr='FDR', predi='geoenv',window=5e5, nTraining=1e3,nTesting=10000, clim='pres'  ){
#
# #####**********************************************************************#####
#
# require('argparser')
#
# ####  set up arguments
# cat("\n")
# p <- arg_parser("Get arguments for command line jobs:")
# p <- add_argument(p, "--fitness", help="Fitness component", default=fitness)
# p <- add_argument(p, "--extrafitness", help="Fitness component on other environmenst", default=extrafitness)
# p <- add_argument(p, "--mod", help="The gwa model", default=mod)
# p <- add_argument(p, "--exp", help="The code of experiment", default=exp)
# p <- add_argument(p, "--threshold", help="Code of threshold. Possible options: BAY, FDR, BONF, ALL", default=thr)
# p <- add_argument(p, "--predi", help="Code of predictors Possible options: env, geo, geoenv", default=predi)
# p <- add_argument(p, "--rfname", help="Which trained random Forest model is wanted", default='spager')
# p <- add_argument(p, "--window", help="Select best SNP within window", default=window)
# p <- add_argument(p, "--nTraining", help="Number of observations to train models", default=nTraining)
# p <- add_argument(p, "--nTesting", help="Number of observations to test models", default=nTesting)
# p <- add_argument(p, "--clim", help="Climate layer: present or climate change", default=clim)
#
# argv<-parse_args(p)
#
# # Check correct
# stopifnot(argv$rfname %in% c('spager','spa','ger','sixenv') )
# stopifnot(argv$fitness %in% c('rFitness','rSeeds','rSurvival_fruit') )
# stopifnot(argv$mod %in% c('lm','bslmm') )
# stopifnot(argv$thr %in% c('ALL','BAY','FDR','BONF','RANDOM'))
# stopifnot(argv$predi %in% c('env','geo','geoenv') )
# stopifnot(argv$clim %in% c('pres',
#                            'pres2000',
#                            'MP2650',
#                            'MP2670',
#                            'MP4550',
#                            'MP4570',
#                            'MP6050',
#                            'MP6070',
#                            'MP8550',
#                            'MP8570'
#                            ) )
#
#
# # Save and pass to environment
# fvar=argv$fitness
# fvarextra=argv$extrafitness
# mod=argv$mod
# exp=argv$exp
# threshold=argv$thr
# predi=argv$predi
# window=argv$window
# nTraining=argv$nTraining
# nTesting=argv$nTesting
# clim=argv$clim
# if(mod =='bslmm'){
#   threshold<-'BAY'
# }
# myrfname=argv$rfname
#
# assign("myrfname", myrfname, envir=globalenv())
# assign("fvar", fvar, envir=globalenv())
# assign("fvarextra", fvarextra, envir=globalenv())
# assign("mod", mod, envir=globalenv())
# assign("exp", exp, envir=globalenv())
# assign("threshold", threshold, envir=globalenv())
# assign("predi", predi, envir=globalenv())
# assign("window", window, envir=globalenv())
# assign("nTraining", nTraining, envir=globalenv())
# assign("nTesting", nTesting, envir=globalenv())
# assign("clim", clim, envir=globalenv())
#
# # print(p)
#
#
# return(p)
# #####**********************************************************************#####
# }
