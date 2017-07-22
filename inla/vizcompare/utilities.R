library(ggplot2)
library(data.table)

load_data <- function(sigma, range, rho, N){
    save_folder <- "/share/scratch/users/nmarquez/sim2dresults/"
    save_file_data <- paste0(save_folder, "sigma_", sigma, "_range_", range, 
                             "_rho_", rho, "_N_", N, "_data.Rda")
    load(save_file_data)
}


plot_latent <- function(sigma, range, rho, N){
    load_data(sigma, range, rho, N)
    ggplot(DataList$latent[time %in% 1:3 & !is.na(obs)], aes(x, y, z=obs)) + 
        geom_tile(aes(fill=obs)) + theme(plot.title=element_text(hjust=0.5)) +
        facet_grid(model~time) +
        scale_fill_gradientn(colors=heat.colors(8)) + labs(title="Latent Field")
}


plot_var <- function(sigma, range, rho, N){
    load_data(sigma, range, rho, N)
    ggplot(DataList$variance[time %in% 1:3 & !is.na(obs)], aes(x, y, z= obs)) + 
        geom_tile(aes(fill=obs)) + theme(plot.title=element_text(hjust=0.5)) + 
        facet_grid(model~time) + labs(title="Variance Estimates")
        scale_fill_gradientn(colors=rev(terrain.colors(8)))
}

plot_params <- function(sigma, range, rho, N){
    load_data(sigma, range, rho, N)
    ggplot(data=DataList$fixed, aes(x=value, fill=method, group=method)) + 
        geom_density(alpha=.6) +  
        geom_vline(aes(xintercept=value), data=DataList$params) + 
        facet_wrap(~par, scales="free")
}