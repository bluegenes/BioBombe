
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(ggplot2)

set.seed(123)

`%>%` <- dplyr::`%>%`

source(file.path('scripts/viz_utils.R'))

dataset_name <- snakemake@params[["dataset_name"]][[1]]
dataset_name

##################################################
                #    VAE     #
##################################################

# Load Data
tybalt_file <- file(snakemake@input[["tybalt"]])
tybalt_data <- readParamSweep(tybalt_file, algorithm = 'Tybalt')

tybalt_data$select_df

# Set output directories and files
tybalt_param_z_png <- file.path(snakemake@output[["vae_loss_png"]])
tybalt_param_z_pdf <- file.path(gsub(pattern = ".png", replacement = ".pdf",  tybalt_param_z_png))
output_fig <- dirname(tybalt_param_z_png) 


# plot final loss
p <- plotFinalLoss(tybalt_data$select_df, algorithm = 'Tybalt', dataset = dataset_name)
p


ggsave(tybalt_param_z_png, plot = p, height = 3, width = 5.5)
ggsave(tybalt_param_z_pdf, plot = p, height = 3, width = 5.5)

tybalt_one_model <- tybalt_data$full_df %>% dplyr::filter(learning_rate == "0.0005", batch_size == "10", epochs == 100)
tybalt_one_model

tybalt_one_model$num_components <-
  dplyr::recode_factor(tybalt_one_model$num_components, 
                       `5` = "Latent Dim: 5",
                       `25` = "Latent Dim: 25", 
                       `50` = "Latent Dim: 50",
                       `75` = "Latent Dim: 75",
                       `100` = "Latent Dim: 100",
                       `125` = "Latent Dim: 125")


tybalt_one_model_png <- file.path(snakemake@output[["vae_training_png"]])
tybalt_one_model_pdf <- gsub(pattern = ".png", replacement = ".pdf",  tybalt_one_model_png) 

p <- plotOneModel(tybalt_data$one_model_df, algorithm = 'Tybalt', dataset = dataset_name)
p

ggsave(tybalt_one_model_png, plot = p, height = 2.5, width = 5)
ggsave(tybalt_one_model_pdf, plot = p, height = 2.5, width = 5)

tybalt_best_params <- tybalt_data$best_params
tybalt_best_params

best_param_file <- file.path(snakemake@output[["vae_best_params"]])
readr::write_tsv(tybalt_best_params, best_param_file)

tybalt_good_training_df <- tybalt_data$melt_df %>%
  dplyr::filter(batch_size == 50, epochs == 100, kappa == "0.0") %>%
  dplyr::filter(
    (learning_rate == "0.002" & num_components == 5) |
      (learning_rate == "0.0015" & num_components == 25) |
      (learning_rate == "0.0015" & num_components == 50) |
      (learning_rate == "0.0015" & num_components == 75) |
      (learning_rate == "0.001" & num_components == 100) |
      (learning_rate == "0.0005" & num_components == 125)
  )

# Reorder the latent space dimensionality for plotting
num_com <- tybalt_good_training_df$num_components
num_com <- factor(num_com, levels = sort(as.numeric(paste(unique(num_com)))))
tybalt_good_training_df$num_components <- num_com

tybalt_best_model_png <- file.path(snakemake@output[["vae_best_model_png"]]) #file.path(output_fig, "z_parameter_best_model_tybalt.png")
tybalt_best_model_pdf <- gsub(pattern = ".png", replacement = ".pdf",  tybalt_best_model_png) 
#tybalt_best_model_pdf <- file.path(output_fig, "z_paraeter_best_model_tybalt.pdf")

p <- plotBestModel(tybalt_good_training_df, algorithm = 'Tybalt', dataset = dataset_name, output_fig_dir=output_fig)
p

ggsave(tybalt_best_model_png, plot = p, height = 2.5, width = 4)
ggsave(tybalt_best_model_pdf, plot = p, height = 2.5, width = 4)

print("end VAE section")
##################################################
                #    DAE     #
##################################################

# Load Data
adage_file <- file(snakemake@input[["adage"]])
adage_data <- readParamSweep(adage_file, algorithm = "ADAGE")

# set output dirs and files
adage_param_z_png <- file.path(snakemake@output[["dae_loss_png"]]) 
adage_param_z_pdf <- gsub(pattern = ".png", replacement = ".pdf",  adage_param_z_png)
output_fig <- dirname(adage_param_z_png)

# Specify that the model has tied weights
p <- plotFinalLoss(adage_data$select_df, algorithm = 'ADAGE', dataset = dataset_name)
p

ggsave(adage_param_z_png, plot = p, height = 2.5, width = 5.5)
ggsave(adage_param_z_pdf, plot = p, height = 2.5, width = 5.5)

# Several hyperparameter combinations did not converge
# This was particularly a result of the low learning rates - filter and replot

### this actually didn't work for my data, needed to modify filters
#adage_data$select_df <- adage_data$select_df %>%
#    dplyr::filter(end_loss < 0.01, learning_rate != 'Learn: 1e-05')

adage_filtered_png <- gsub(pattern = ".png", replacement = "_filtered.png",  adage_param_z_png)


filter_by_learning_rate_end_loss <- function(dataDF, dataset_name, pngfile, LR_filter = "Learn: 1e-05", Loss_filter=0.01){
    dataDF$select_df <- dataDF$select_df %>%
    dplyr::filter(end_loss < Loss_filter, learning_rate != LR_filter)
    p <- plotFinalLoss_ADAGE(dataDF$full_df, algorithm = 'ADAGE', dataset = dataset_name, output_fig_dir = dirname(pngfile))
    ggsave(adage_filtered_png, plot = p, height = 2.5, width = 5.5)
    adage_filtered_pdf <-  gsub(pattern = ".png", replacement = ".pdf",  adage_filtered_png) 
    ggsave(adage_filtered_pdf, plot = p, height = 2.5, width = 5.5)
    return(adage_data$select_df)
}

orig_select_df <- adage_data$select_df
orig_select_df

adage_data$select_df <- tryCatch(
    {
    filter_by_learning_rate_end_loss(adage_data, dataset_name, adage_filtered_png)
    },
    error = function(e){
         # restore original adage_data$select_df
        adage_data$select_df <- orig_select_df
    }
)

adage_best_params <- adage_data$best_params
adage_best_params

best_param_file <- file.path(snakemake@output[["dae_best_params"]])
readr::write_tsv(adage_best_params, best_param_file)

adage_tied_good_training_df <- adage_data$melt_df %>%
  dplyr::filter(sparsity == "0.0",
                epochs == 100,
                batch_size == 50,
                noise == "0.0") %>%
  dplyr::filter(
    (num_components == 5 & learning_rate == "0.0015") |
      (num_components == 25 & learning_rate == "0.0015") |
      (num_components == 50 & learning_rate == "0.0005") |
      (num_components == 75 & learning_rate == "0.0005") |
      (num_components == 100 & learning_rate == "0.0005") |
      (num_components == 125 & learning_rate == "0.0005"))

num_com <- adage_tied_good_training_df$num_components
num_com <- factor(num_com, levels = sort(as.numeric(paste(unique(num_com)))))
adage_tied_good_training_df$num_components <- num_com

adage_best_model_png <- file.path(snakemake@output[["dae_best_model_png"]])
adage_best_model_pdf <- gsub(pattern = ".png", replacement = ".pdf",  adage_best_model_png)

p <- plotBestModel(adage_tied_good_training_df, algorithm = 'ADAGE', dataset = dataset_name, output_fig_dir=output_fig)
p

ggsave(adage_best_model_png, plot = p, height = 2.5, width = 4)
ggsave(adage_best_model_pdf, plot = p, height = 2.5, width = 4)
