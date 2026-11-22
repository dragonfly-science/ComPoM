# Helper functions shared by the simulation vignettes and scenario scripts.
# Not part of the ComPoM package. Source with source("vignette_helpers.R") from the vignettes
# directory.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(tidybayes); library(cowplot); library(patchwork)})

# One theme for all vignette figures: the package theme (cowplot with faded major gridlines)
ggplot2::theme_set(ComPoM::theme_compom())

# Logistic selectivity: 50% selected at L50, 95% selected at L50 + (L95 - L50)
Selectivity <-function (bins,
                        L50,
                        L95) {

  Selex=1/(1+exp(-log(19)*( (L50-bins)/(L95-L50) )))

  return(Selex)

}

# Relative error of modelled and scaled compositions against the known true composition.
# Used to compare estimation methods in the simulation vignettes.
scaled_comp_comp <- function(scaled_comp=NULL,
                             grps = NULL,
                             scaled=T,
                             scales = "free_y",
                             comp_are = NULL,
                             comp_are2 = NULL,
                             cvar='tcatch'){

  if(scaled == T) {

    preda <- scaled_comp %>%
      group_by(across(all_of(grps) ), bin, .draw) %>%
      summarise(tot = sum(tot_by_bin)) %>%
      group_by(across(all_of(grps) ), .draw) %>%
      mutate(prop = tot/sum(tot)) %>%
      group_by(across(all_of(grps) ), bin) %>%
      mean_qi(prop) %>%
      filter(!is.na(prop)) %>%
      mutate(bin = as.numeric(as.character(bin)))

    if(!is.null(comp_are)){
      comp <- comp_are %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(tot = sum(!!sym(cvar))) %>%
        group_by(across(all_of(grps) )) %>%
        mutate(prop = tot/sum(tot)) %>%
        filter(!is.na(prop),prop>0,prop<1) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }

    if(!is.null(comp_are2)){
      comp2 <- comp_are2 %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(tot = sum(!!sym(cvar))) %>%
        group_by(across(all_of(grps) )) %>%
        mutate(prop = tot/sum(tot)) %>%
        filter(!is.na(prop),prop>0,prop<1) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }

  } else {

    preda <- scaled_comp %>%
      group_by(across(all_of(grps) ), bin, .draw) %>%
      summarise(tot = sum(tot_by_bin)) %>%
      group_by(across(all_of(grps) ), bin) %>%
      mean_qi(prop = tot)  %>%
      filter(!is.na(prop)) %>%
      mutate(bin = as.numeric(as.character(bin)))

    if(!is.null(comp_are)){

      comp <- comp_are %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(prop = sum(!!sym(cvar))) %>%
        filter(!is.na(prop),sum(prop)>100) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }

    if(!is.null(comp_are2)){

      comp2 <- comp_are2 %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(prop = sum(!!sym(cvar))) %>%
        filter(!is.na(prop),sum(prop)>100) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }
  }


  preda %>% rename(model=prop) %>%
    inner_join(comp %>% rename(real=prop) %>% select(-tot)) %>%
    inner_join(comp2 %>% rename(scaled=prop)) %>%
    pivot_longer(c(model, scaled), values_to = 'value', names_to = 'method') %>%
    mutate(err = abs((value-real)/real),
           method = factor(method, c('model', 'scaled'), c('Model', 'Raw samples'))) %>%

  ggplot( ) +
    geom_boxplot(aes(x=method,y=err, col=method), outliers = F, width = 0.4) +
    scale_color_manual(values = c('skyblue','orange')) +
    xlab(NULL) +
    ylab('Relative error') +
    theme_compom() +
    theme(legend.position = 'none')


}
