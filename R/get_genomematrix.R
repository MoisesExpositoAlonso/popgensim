get_genomematrix<-function (size = 1000, positions = NULL, type = c("random", "window",
    "top"), start = 1)
{
    genomes <- readRDS("databig/genomes.rda")
    G <- attachgenomes(genomes)
    Map <- genomes$map
    if (!is.null(positions)) {
        message("subsetting the provided positions")
        stopifnot(is.numeric(positions))
        X = inputena.mat(G[, positions])
    }
    else if (type == "random") {
        message("subsetting ", size, " random positions ")
        X = inputena.mat(G[, sample(1:ncol(G), size = size)])
    }
    else if (type == "window") {
        message("subsetting a window of ", size, " base pairs starting at ",
            start)
        X = inputena.mat(G[, start:(start + size)])
    }
    else if (type == "top") {
        X <- get_topSNPs()
    }
    else {
        stop("None of the possibilities were given")
    }
    return(X)
}
