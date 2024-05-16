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

  if(!is.null(grp)) p <- p+facet_wrap(vars(!!grps),scales = 'free_y')

  return(p)
}

#' Effects
#' @param comp_data Compositional dataframe
#' @param mod An object from fit_model
#' @param grp group to do PPC over
#' @param form formula for effect
#' @param grid optional plotting with facet grid along two variables (need to be in grp)
#' @param cvar continous variable if used in model
#'
#' @import cowplot
#' @import dplyr
#' @importFrom posterior draws_df
#' @export
#'
Fx_plot <- function(comp_data, mod, grp='Year', form='~(1|bin:Year)', grid=NULL, cvar=NULL){

  if(!is.null(cvar)) grps <- as.symbol(grp[(!grp==cvar)])

  int <- as_draws_df(mod$mod, "b_Intercept") %>% dplyr::select(-.chain,-.iteration)

  preda <-
    comp_data %>% ungroup() %>% dplyr::select(bin,!!grp) %>% distinct() %>% mutate(n=100) %>%
    complete(nesting(!!!syms(grp)),bin,fill=list(n=100)) %>%
    add_linpred_draws(mod$mod, re_formula = form, allow_new_levels=T) %>%
    inner_join(int) %>%
    group_by(bin,!!!syms(grp)) %>%
    median_qi(pred = exp((.linpred-b_Intercept)-mean(.linpred-b_Intercept)))


  p <- ggplot(preda)  + geom_hline(yintercept = 1,linetype=2,alpha=0.5) +
    geom_pointrange(aes(x=bin,y=pred,ymin=.lower,ymax=.upper)) +
    xlab('Length bin (cm)') +
    ylab('Multiplier') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))
  #browser()
  if(is.null(grid))  p <- p+facet_wrap(facets = vars(!!grps),scales = 'free_y')
  if(!is.null(grid)) p <- p+facet_grid(rows = vars(!!as.symbol(grid[1])),cols = vars(!!as.symbol(grid[2])),scales = 'free_y')

  return(p)
}

#' scale compositions
#' @param scale_df data frame used for scaling
#' @param predvar Column in the scale_df used for scaling compositions
#' @param grps groups to do scaling
#' @param fit An object from fit_model
#' @param form an alternative random effect structure to the fitted model
#'
#' @import cowplot
#' @export
#'
scale_comps <- function(scale_df, predvar='catch', fit = NULL, grps, iters=NULL, form=NULL){

  if(!is.null(form)){
    form_parts <- stringr::str_remove_all(stringr::str_split(form, '\\+')[[1]],pattern = ' ')
    form <- paste('~ (1|bin) +',paste0('(1|bin:',form_parts,')', collapse = ' + '))
  } else if(!is.null(fit)){
    form = paste("~",paste(paste('(1|',fit$model$ranef$group,')'), collapse='+'))
  } else {stop("Must provide either a formula or a model fit from fit_model")}

  #browser()
  scale_df %>%
    group_by(across(all_of(grps) )) %>%
    summarize(n=sum(!!sym(predvar),na.rm=T)) %>%
    mutate(bin = factor(min(unique(fit$mod$data$bin)), levels=unique(fit$mod$data$bin))) %>%
    ungroup() %>% # need to ungroup before augmenting
    complete(nesting(!!!syms(grps)),bin,fill = list(n=0))  %>%
    group_by(across(all_of(grps) )) %>%
    mutate(n = sum(n),
           bin = as.numeric(as.character(bin))) %>%
    filter(n>0) %>%
    add_predicted_draws(fit$mod, allow_new_levels=T, ndraws = iters, value = 'tot_by_bin', re_formula = form)

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

  if(length(grps)==1) p <- p+facet_wrap(facets = vars(!!as.symbol(grps)),scales = scales,drop = T,ncol = 4)
  if(length(grps)>=2) p <- p+facet_grid(rows = vars(!!as.symbol(grps[1])),cols = vars(!!!syms(grps[2:length(grps)])),scales = scales)

  p
}


#' Plot scaled compositions as ridgeline plot
#' @param scaled_comp data frame used for scaling
#' @param grps groups to do scaling
#' @param scales free y axis?
#' @param plot_y y axis for ridgeline
#' @param plot_facet faceting ridgelines (e.g., monthly or quarterly ridges within years)
#'
#' @import cowplot
#' @export
#'
scaled_ridge_plot <- function(scaled_comp = NULL,
                              grps = NULL,
                              scales = "free_y",
                              plot_y = NULL,
                              plot_facet = NULL) {
  sdat <- scaled_comp %>%
    ungroup() %>%
    group_by(across(all_of(grps)), .draw, bin) %>%
    summarise(n=sum(tot_by_bin)) %>%
    mutate(p=n/sum(n)) %>%
    group_by(across(all_of(grps)), bin) %>%
    median_qi(p)

  sd <- sdat %>%
    mutate(bin = as.numeric(as.character(bin))) %>%
    filter(bin < 80, bin > 40,
           !is.na(quarter)
           # ,!quarter %in% c(3)
    ) %>%
    mutate(cat = cut(bin,breaks = c(0,55,67,100)))

  p <- ggplot(sd, aes(x=bin, y=!!!syms(plot_y))) +
    geom_density_ridges_gradient(aes(x=bin, y=!!sym(plot_y), height=p, fill=bin),
                                 scale = 2, rel_min_height = 0.01, stat='identity') +
    geom_density_ridges_gradient(aes(x=bin, y=!!sym(plot_y), height=p),
                                 scale = 2, stat='identity',fill=NA) +
    scale_fill_viridis_c(name = "", guide='none') +
    facet_grid(as.formula(paste(plot_facet, "~ .")),as.table = F) +
    xlab('Length (cm)') +
    ylab(plot_y) +
    theme_bw()

  return(p)
}
