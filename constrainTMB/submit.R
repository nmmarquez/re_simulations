rm(list=ls())
library(dplyr)

Nsets <- 20
Sizes <- c(10^3, 10^4, 10^5, 10^6)

configDF <- data.frame(
    N=rep(Sizes, each=Nsets * 2),
    const=rep(0:1, Nsets * length(Sizes)),
    seed=rep(1:(Nsets * length(Sizes)), each=2)) %>%
    mutate(M=.1*N)

for(i in 1:nrow(configDF)){
    N <- configDF$N[i]
    const <- configDF$const[i]
    seed <- configDF$seed[i]
    call_ <- paste0(
        "qsub -pe multi_slot 20 ",
        "-e /share/temp/sgeoutput/$USER/errors/ ",
        "-o /share/temp/sgeoutput/$USER/output/ ",
        "-N constrain_", N, "_", const, "_", seed,
        " /share/singularity-images/lbd/shells/singR.sh ",
        "-m 1 -o 1 -e s ",
        "$HOME/Documents/TMBtests/constrain/constrain.R ",
        N, " ", const, " ", seed)
    print(call_)
    print("")
    system(call_)
}

