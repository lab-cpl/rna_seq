#################################################################
##                          Libraries                          ##
#################################################################

pacman::p_load(
    tidyverse,
    tidymodels,
    ggplot2,
    bannerCommenter,
    edgeR,
    doFuture,
    finetune,
    tictoc,
    see,
    ggpubr,
    beepr,
    themis,
    annotables,
    furrr
)

##################################################################
##                       Helper functions                       ##
##################################################################

# https://stackoverflow.com/questions/47044068/get-the-path-of-current-script
# get path of source file
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
    if (length(this_file)==0)
    {
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}

##################################################################
##                     Data sets for models                     ##
##################################################################

# sets path on source file location
script_path <- getCurrentFileLocation()
setwd(script_path)

# prepare data
# raw counts per gene
# ID 6 is left out
# both directions are summed
import <- list.files("~/data/rna-seq/", pattern = "*.tab", full.names = TRUE) %>% 
    map_dfr(., function(x){
        read_delim(x, col_names = c(
            "ensembl_gene_id",
            "val",
            "val1",
            "val2"
        ), delim = "\t") %>% 
            mutate(
                ID = str_extract(x, "[0-9]"),
                ensembl_gene_id = str_extract(ensembl_gene_id, "ENSMUSG[0-9]*")
                ) %>% 
            tail(-4) %>% 
            select(-c("val1", "val2")) %>% 
            group_by(ensembl_gene_id, ID) %>% 
            summarise(
                val_sum = sum(val)
            )
    }) %>% filter(ID != 6) # this sample does not cluster



# library(annotables)
#inner_join(grcm38, by = c("ensembl_gene_id" = "ensgene"))

raw_data <- import %>% 
    group_by(ensembl_gene_id) %>% 
    mutate(s = sum(val_sum) == 0) %>% 
    filter(s == FALSE) %>% 
    select(-s) %>% 
    ungroup() %>% 
    pivot_wider(
        names_from = ID,
        values_from = val_sum
    )

#################################################################
##                        Normalization                        ##
#################################################################


# normalize data by library length
# this is the most typical normalization
# sum all gene variants, so we get count per gene ignoring isoforms
raw_data_sum <- raw_data %>%
    group_by(
        ensembl_gene_id
    ) %>%
    summarise(
        across(where(is.numeric), list(sum))
    )
# write back col names
colnames(raw_data_sum) <- colnames(raw_data)

# set data as matrix
raw_data_sum_mat <- raw_data_sum %>%
    select(-ensembl_gene_id) %>%
    as.matrix()
# set data into DGE matrix format
rownames(raw_data_sum_mat) <- raw_data_sum$ensembl_gene_id


# just to check number of 0 counts
raw_data_sum %>% 
    filter(across(is.numeric, ~ .x < 0))

# create DGE object
raw_data_sum_dge <- DGEList(raw_data_sum)
norm_data_dge <- calcNormFactors(raw_data_sum_dge)
# TMM + log2 transform DGE object
cpm_log2_data <- cpm(
    norm_data_dge,
    log = TRUE
    ) %>%
    as_tibble() %>%
    mutate(
        ensembl_gene_id = raw_data_sum$ensembl_gene_id,
        .before = `1`
    )


# create data for models

# sample info
# this specifies treatment for each ID
sample_info <- read_csv("info_sample.csv") %>%
    rename(
        ID = SampleName
    ) %>%
    mutate(
        ID = as.character(ID)
    )

# full model data
full_model_data <- cpm_log2_data %>%
    pivot_longer(
        cols = where(is.numeric),
        names_to = "ID"
    ) %>%
    left_join(
        ., sample_info, by = c('ID')
    ) %>%
    pivot_wider(
        names_from = ensembl_gene_id,
        values_from = value
    ) %>% 
    mutate(
        ID = as.factor(ID),
        group = as.factor(group),
    )

# plot cpm distribution across all genes
# Chen Y, Lun AT, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Res. 2016 Jun 20;5:1438. doi: 10.12688/f1000research.8987.2. PMID: 27508061; PMCID: PMC4934518.

# get cpm gene distribution
cpm_dist <-
    full_model_data %>% 
    pivot_longer(
        cols = starts_with("ENS")
    ) %>% 
    inner_join(grcm38, by = c("name"="ensgene")) %>% 
    filter(
        !is.na(entrez),
        !grepl("predicted", description),
           grepl("protein", biotype)
        )

# rule of thumb cpm is above 0.5 or 10 / L
# where L is the minimum library size
lib_sizes <- 
    raw_data_sum_dge$samples$lib.size %>% 
    {10 / (min(.)/1000000)}
lib_sizes

cpm_dist_plot <-
    cpm_dist %>% 
    ggplot(aes(
        value, group = ID, color = ID
    )) +
    geom_density() +
    geom_vline(xintercept = 0.5) +
    theme(legend.position = "none")
cpm_dist_plot
    

# this is the actual filter
# select only protein coding genes
# cpm > 0.5
# genes in entrez data base
# after previous filters, genes that are present in more than 3 animals
highly_expressed_genes <- full_model_data %>% 
    pivot_longer(
        cols = starts_with("ENS")
    ) %>% 
    inner_join(grcm38, by = c("name" = "ensgene")) %>% 
    filter(!is.na(entrez), !grepl("predicted", description),
           grepl("protein", biotype)) %>% 
    filter(value > 0.5) %>% 
    group_by(name) %>% 
    summarise(n = n()) %>% 
    filter(n > 3) %>%
    pull(name)

# DO filter
full_model_data <- full_model_data %>% 
    select(
        ID, group, all_of(highly_expressed_genes)
    )

# variance ranked data
# set genes in a descending order y standard deviation
# the idea is to think of gene in two axis: deviation and expression
ord_var_rank_data <- full_model_data %>% 
    pivot_longer(
        cols = contains("ENSMUSG"),
        names_to = "gene_name",
        values_to = "cpm"
    ) %>% 
    group_by(gene_name) %>% 
    summarise(
        sd_val = sd(cpm)
    ) %>% 
    ungroup() %>% 
    arrange(desc(sd_val)) %>% 
    pull(gene_name)

# TODO
# plot var pre var post

# get the first 100 partitions, from a 1000 division, 10% most variable
partitions <- floor(seq(1, length(ord_var_rank_data), length.out = 1000)) %>% 
    head(n=100)

# generate a list with partitions
full_model_data_parts <- partitions %>% 
    map(., function(part){
        sel_vec <- c("ID", "group", ord_var_rank_data[1:part])
        out <- full_model_data %>% 
            select(all_of(sel_vec))
        return(out)
    })




##################################################################
##                        Model workflow                        ##
##################################################################

# TODO
# plot var pre var post


# optimization of parameters is done through grid-search
# objective function ~ cross validation
# hard-coded are the results of optimization
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1228-x#MOESM17
tree_grid <- full_model_data_recipe %>% 
    map(., function(x){
        # optimize hyper-parameters
        # create grid
        # get total number of predictor genes
        n_genes <- summary(x)$variable %>% 
            str_subset("ENSMUSG") %>% 
            { print(length(.)) }
        sqrt_n_genes <- round(sqrt(n_genes))
        tree_grid <- expand.grid(
            mtry = c(2, 10, 20),#527
            min_n = c(2, 4, 6), #15
            trees = c(5, 50, 100) #66
        ) %>% 
            unique
    })


# this is the tune process
tuned_models <- 1:5 %>%
    future_map(., function(x){
        tictoc::tic()
        full_model_data_optim <- full_model_data_workflow[[x]] %>%
            tune_grid(
                full_model_data_fold[[x]],
                grid = tree_grid[[x]],
                metrics = metric_set(accuracy, roc_auc, kap)
                )
        tictoc::toc()
        fn <- paste("tuned_model_uncertainty", x, sep = "_")
        saveRDS(full_model_data_optim, file = fn)
        return(full_model_data_optim)
    })

# load tuned models into a tibble
tuned_model_list <- list.files(pattern = "tuned_model_uncertainty_[0-9]+") %>% 
    str_sort(., numeric = TRUE) %>% 
    map_dfr(., function(x){
        data_perc <- str_extract(x, pattern = "[0-9]+")
        tuned_model <- readRDS(x) %>% 
            unnest(.metrics) %>% 
            mutate(data_perc = as.numeric(data_perc))
        return(tuned_model)
    })

##################################################################
##                  Optimization visualization                  ##
##################################################################

# data for modeling contributions of hyperparameters
params_results <-
tuned_model_list %>% 
    group_by(data_perc, .metric, mtry, trees, min_n) %>% 
    drop_na() %>% 
    summarise(
        mean_estimate = mean(.estimate) * 100,
    ) %>% 
    pivot_wider(
        names_from = .metric,
        values_from = mean_estimate
    ) %>% 
    arrange()
params_results

#################################################################
##                   Optimized model fitting                   ##
#################################################################

# # load tuned models objects
# tuned_models_obj <- list.files(pattern = "tuned_model_uncertainty_no_sample_6_10") %>% 
#     str_sort(., numeric = TRUE) %>% 
#     map(., function(x){
#         tuned_model <- readRDS(x)
#         return(tuned_model)
#     })

# after optimization this were the optimal parameters
rf_params <- tibble(
    mtry = 20,
    trees = 10,
    min_n = 5,
    .config = "Preprocessor1_Model1"
)

# partition 3 (30% sd) was the optimal partition
# re-fit model with optimal parameters
full_model_data_optim <- full_model_data_workflow[[3]] %>%
    tune_grid(
        full_model_data_fold[[3]],
        grid = data.frame(
            mtry = 20,
            trees = 10,
            min_n = 5
            ),
        metrics = metric_set(accuracy, roc_auc, kap)
        )

# show accuracy of re-fit model 
# rf are highly stochastic so this could change considerably every time
# acc is ~0.7
full_model_data_optim %>% 
    show_best(metric = "accuracy")

# here the model is finalized
# using optimal param and optimal partition
# obtain a more robust estimate of performance
# change vector from 1:10 to compute all relevant partitions
fitted_models <- 1:10 %>%
    map(., function(x){
        tic()
        folds <- vfold_cv(full_model_data_parts[[x]], v = 6, repeats = 10)
        fitted_model <- finalize_workflow(
            full_model_data_workflow[[x]],
            rf_params
        ) %>%
            fit_resamples(folds, metrics = metric_set(accuracy, roc_auc, kap),
                       control = control_resamples(save_pred = TRUE))
        toc()
        print(paste("Model", x, "DONE!", sep = " "))
        return(fitted_model)
    })

# accuracy over partitions
1:10 %>% 
    map_dfr(
        ., function(x){
            collect_metrics(fitted_models[[x]]) %>% 
                mutate(n = x)
        }
    ) %>% 
    filter(.metric == "accuracy") %>% 
    ggplot(aes(
        n, mean
    )) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = 1:10)




#saveRDS(fitted_models, "fitted_models_full_data")
fitted_models <- readRDS("fitted_models_full_data")

    
# test <- 1:10 %>% 
#  map(., function(x){
#      tic()
#      best_params <- select_best(tuned_models_obj[[1]], metric = "kap")
#      folds <- vfold_cv(full_model_data_parts[[x]], v = 6)
#      fitted_model <- finalize_workflow(
#          full_model_data_workflow[[x]],
#          best_params
#      ) %>% 
#          fit_resamples(folds, metrics = metric_set(accuracy, roc_auc, kap),
#                        control = control_resamples(save_pred = TRUE))
#      toc()
#      print(paste("data model", x, "DONE!", sep = " "))
#      return(fitted_model)
#  })

# tic()
# best_params <- rf_params
# out <- finalize_workflow(
#         full_model_data_workflow[[3]],
#         best_params
#     ) %>% 
#         fit(full_model_data_parts[[100]]) %>% 
#         extract_fit_engine()
# toc()

# compute variable importance
# for optimal partition with optimal parameters
# due to stochastic rf process this should be computed multiple times
variable_importance <- finalize_workflow(
    full_model_data_workflow[[3]],
    rf_params
) %>% 
    fit(full_model_data_parts[[3]]) %>% 
    extract_fit_engine()
variable_importance

var_imp_est <-
    1:1000 %>% 
    future_map_dfr(
        ., function(x){
            var_imp_obj <-
                finalize_workflow(
                    full_model_data_workflow[[3]],
                    rf_params
                ) %>% 
                fit(full_model_data_parts[[3]]) %>% 
                extract_fit_engine()
            s <-
                var_imp_obj %>% 
                # key function for estimated gini importance
                randomForest::importance() %>% 
                as.table() %>% 
                as.data.frame() %>% 
                as_tibble() %>% 
                filter(Var2 == "MeanDecreaseGini") %>% 
                arrange(desc(Freq)) %>% 
                mutate(rank = row_number()) %>% 
                rename(ensembl_gene_id = Var1) %>% 
                inner_join(grcm38, by = c("ensembl_gene_id" = "ensgene")) %>% 
                filter(
                    !is.na(entrez),
                    !grepl("predicted", description)
                    ) %>% 
                mutate(iteration = x)
            return(s)
        }
    )
var_imp_est

# mean variable importance
mean_var_imp <-
    var_imp_est %>% 
        group_by(description) %>% 
        summarise(
            mdg = mean(Freq),
            mrank = mean(rank)
        ) %>% 
    arrange(desc(mdg))
mean_var_imp

# variable importance plot
gini_decrease_plot <-
    var_imp_est %>% 
        ggplot(aes(
            forcats::fct_reorder(factor(description,), Freq, mean, .desc = TRUE), Freq
        )) +
        stat_summary(fun.data = "mean_cl_boot", geom = "point") +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar") +
        theme_pubr() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
gini_decrease_plot

# Rank plot
rank_plot <-
    var_imp_est %>% 
        ggplot(aes(
            forcats::fct_reorder(factor(description,), rank, mean, .desc = FALSE), rank
        )) +
        stat_summary(fun.data = "mean_cl_boot", geom = "point") +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar") +
        theme_pubr() +
    scale_y_reverse() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
rank_plot

var_imp_list <-
    1:100 %>% 
    future_map(
        ., function(x){
            var_imp_obj <-
                finalize_workflow(
                    full_model_data_workflow[[3]],
                    rf_params
                ) %>% 
                fit(full_model_data_parts[[3]]) %>% 
                extract_fit_engine()
        }
    )
var_imp_list

pdp_plots_data <-
    var_imp_list %>% 
    future_map_dfr(
        ., function(x){
            # get gene names
            gn <-
                var_imp_est %>% 
                filter(Freq > 0) %>% 
                pull(ensembl_gene_id) %>% 
                unique()
            out <-
                gn %>% 
                map_dfr(., function(g){
                    p <- pdp::partial(
                    x,
                    pred.var = c(g),
                    grid.resolution = 100,
                    which.class = "treatment",
                    prob = TRUE,
                    plot = FALSE, 
                    train = as.data.frame(full_model_data_parts[[3]])
                ) %>% 
                    as_tibble() %>% 
                    mutate(ensembl_gene_id = colnames(.)[1]) %>% 
                    left_join(grcm38, by = c("ensembl_gene_id" = "ensgene"))
                    colnames(p)[1] <- "x"
                    return(p)
                })
            return(out)
        }
    )
pdp_plots_data

pdp_plots_data$symbol <- factor(
    pdp_plots_data$symbol,
    levels = c(
            "Mup6", "Pmch", "Hcrt",
            "Pvalb", "Cplx3", "Chrna2"
    )
)
pdp_gam <-
    pdp_plots_data %>% 
    drop_na() %>% 
    filter(
        symbol %in% c(
            "Mup6", "Pmch", "Hcrt",
            "Pvalb", "Cplx3", "Chrna2"
            )
    ) %>% 
        ggplot(aes(
            x, yhat
        )) +
        geom_smooth(method = "gam") +
        facet_wrap(~symbol, scales = "free") +
    theme_pubr() +
    theme(
        text = element_text(size = 30)
    ) +
    xlab("CPM") +
    ylab("P(class = Uncertainty)")
pdp_gam
    

# get the slopes
pdp_gini_slopes <-
    pdp_plots_data %>% 
    group_by(symbol) %>% 
    group_split() %>% 
    map_dfr(
        ., function(gene){
            mdl <- lm(yhat ~ x, data = gene)
            out <- broom::tidy(mdl, conf.int = TRUE)
            slope <- out %>%
                pull(estimate) %>% 
                {.[2]}
            conf.low <- out %>% 
                pull(conf.low) %>% 
                {.[2]}
            conf.high <- out %>% 
                pull(conf.high) %>% 
                {.[2]}
            r <- tibble(
                slope = slope,
                conf.low = conf.low,
                conf.high = conf.high,
                gene = gene$symbol %>% unique()
            )
            return(r)
        }
    ) %>% 
    mutate(gene = as.factor(gene), sign = sign(slope)) %>% 
    arrange(desc(abs(slope)))
pdp_gini_slopes

pdp_gini_slopes$gene <- factor(pdp_gini_slopes$gene,
                               levels = pdp_gini_slopes$gene[order(desc(pdp_gini_slopes$slope))])

pdp_gini_slopes_plot <-
    pdp_gini_slopes %>% 
    drop_na() %>% 
        ggplot(aes(
            gene, slope, ymin = conf.low, ymax = conf.high, fill = as.factor(sign)
        )) +
        geom_hline(yintercept = 0, color = "gray", linewidth = 1) +
        geom_col(shape = 15, size = 2) +
        geom_errorbar(width = 0.5, linewidth = 1) +
        theme_pubr() +
        theme(
                text = element_text(size = 30),
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none"
        ) +
        ylab("P(class = Uncertainty) slope") +
        xlab("")
pdp_gini_slopes_plot
    
