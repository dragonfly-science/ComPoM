#' Post predictive
#' @param mod An object from fit_model
#' @param grp group to do PPC over
#' @param xlab label for composition plot
#'
#' @import cowplot
#' @export
#'

post_pred_group <- function(mod, grp=NULL, xlab = 'Length (cm)'){

  pds <- mod$preds
  if(!is.null(grp)) grps <- sym(grp)
  pps <- pds  %>%
    group_by(!!grps,bin, .draw) %>%
    summarise(tot = sum(tot_by_bin), .pred = sum(.prediction)) %>%
    group_by(!!grps,bin, tot) %>%
    median_qi(.pred)

  p <- ggplot(pps) +
    geom_pointrange(aes(x=bin,y=.pred,ymin=.lower,ymax=.upper)) +
    geom_point(aes(x=bin,y=tot), col='skyblue') +
    xlab(xlab) +
    ylab('Prediction') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))

  if(!is.null(grp)) print(p+facet_wrap(vars(!!grps),scales = 'free_y')) else print(p)
}

#' Effects
#' @param comp_data Compositional dataframe
#' @param mod An object from fit_model
#' @param grp group to do PPC over
#' @param form formula for effect
#'
#' @import cowplot
#' @importFrom posterior draws_df
#' @export
#'
Fx_plot <- function(comp_data, mod, grp='Year', form='~(1|bin:Year)', grid=NULL){

  if(!is.null(grp)) grps <- as.symbol(grp)

  int <- as_draws_df(mod$mod, "b_Intercept") %>% select(-.chain,-.iteration)

  preda <-
    comp_data %>% ungroup() %>% select(bin,!!grp) %>% distinct() %>% mutate(n=100) %>%
    complete(nesting(!!!syms(grp)),bin,fill=list(n=100)) %>%
    add_linpred_draws(mod$mod, re_formula = form, allow_new_levels=T) %>%
    inner_join(int) %>%
    group_by(bin,!!!syms(grp)) %>%
    median_qi(pred = exp((.linpred-b_Intercept)-mean(.linpred-b_Intercept)))


  p <- ggplot(preda)  + geom_hline(yintercept = 1,linetype=2,alpha=0.5) +
    geom_pointrange(aes(x=bin,y=pred,ymin=.lower,ymax=.upper)) +
    xlab('Length bin (cm)') +
    ylab('Prediction') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))
  #browser()
  if(is.null(grid)) print(p+facet_wrap(facets = vars(!!grps),scales = 'free_y'))
  if(!is.null(grid)) print(p+facet_grid(rows = vars(!!as.symbol(grid[1])),cols = vars(!!as.symbol(grid[2])),scales = 'free_y'))
}

#' scale compositions
#' @param scale_df data frame used for scaling
#' @param predvar Column in the scale_df used for scaling compositions
#' @param grps groups to do scaling
#' @param mod An object from fit_model
#'
#' @import cowplot
#' @export
#'
scale_comps <- function(scale_df, predvar='catch', mod, grps){

  scale_df %>%
    group_by(across(all_of(grps) )) %>%
    summarize(n=sum(!!sym(predvar),na.rm=T)) %>%
    mutate(bin = factor(min(unique(mod$mod$data$bin)), levels=unique(mod$mod$data$bin))) %>%
    ungroup() %>% # need to ungroup before augmenting
    complete(nesting(!!!syms(grps)),bin,fill = list(n=0))  %>%
    group_by(across(all_of(grps) )) %>%
    mutate(n = sum(n)) %>%
    filter(n>0) %>%
    add_predicted_draws(mod$mod, allow_new_levels=T, value = 'tot_by_bin')

}

#' Plot scaled compositions
#' @param scaled_comp data frame used for scaling
#' @param grps groups to do scaling
#' @param scaled proportions or catch-at
#' @param scales free y axis?
#'
#' @import cowplot
#' @export
#'
scaled_comp_plot <- function(scaled_comp=NULL,
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
      median_qi(prop) %>%
      filter(!is.na(prop),prop>0,prop<1) %>%
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
      median_qi(prop = tot)  %>%
      filter(!is.na(prop),sum(prop)>100) %>%
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


  p <- ggplot( preda ) +
    geom_pointrange(aes(x=bin,y=prop,ymin=.lower,ymax=.upper)) +
    {if(!is.null(comp_are)) geom_point(aes(x=bin,y=prop), col='skyblue',data=comp)} +
    {if(!is.null(comp_are2)) geom_point(aes(x=bin,y=prop), col='orange',data=comp2)} +
    xlab('Length bin (cm)') +
    ylab('Prediction') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))

  if(length(grps)==1) print(p+facet_wrap(facets = vars(!!as.symbol(grps)),scales = scales,drop = T,ncol = 4))
  if(length(grps)>=2) print(p+facet_grid(rows = vars(!!as.symbol(grps[1])),cols = vars(!!!syms(grps[2:length(grps)])),scales = scales))

}