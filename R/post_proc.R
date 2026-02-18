#' Post predictive
#' @param mod An object from fit_model
#' @param grp group to do PPC over
#' @param xlab label for composition plot
#'
#' @import cowplot
#' @import ggplot2
#' @export
#'
post_pred_group <- function(mod, grp=NULL, xlab = 'Length (cm)'){

  pds <- mod$preds
  if(!is.null(grp)) grps <- sym(grp)

  #browser()
  pps <- pds  %>%
    group_by(!!grps,bin, .draw) %>%
    summarise(tot = sum(tot_by_bin), .pred = sum(.prediction)) %>%
    group_by(!!grps, .draw) %>%
    mutate( .pred = (.pred/sum(.pred))*sum(tot)) %>%
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
#' @param grp group to do PPC over, needs to include cvar if present
#' @param form formula for effect
#' @param grid optional plotting with facet grid along two variables (need to be in grp)
#' @param cvar continuous variable if used in model
#'
#' @import cowplot
#' @import dplyr
#' @importFrom posterior draws_df
#' @export
#'
Fx_plot <- function(mod, grp='Year', form='~(1|bin:Year)', grid=NULL, cvar=NULL){

  if(!is.null(cvar)) grps <- syms(grp[(!grp %in% cvar)]) else grps <- syms(grp)

  #int <- as_draws_df(mod$mod, "b_Intercept") %>% dplyr::select(-.chain,-.iteration)
  helpers = any(grp %in% mod$data$bin)
  helper_vars = grp[grp %in% mod$data$bin]

  # Columns to exclude from factor conversion (coordinates for spatial models)
  exclude_cols <- paste0('^', grp, '$')

  if(any(grepl('TMB', attr(mod, "call")))) {
   #if(!is.null(coords)) exclude_col <- paste(exclude_col, paste0('^', coords, '$', collapse='|'), sep='|')
   exclude_cols <- paste(exclude_cols, paste0('^', mod$time, '$'), sep='|')
  }

  mean_or_mode <- function(x) {
    if(is.numeric(x)) mean(x) else factor(names(sort(-table(x)))[1])
  }

  #browser()
  preda_data <-
    mod$model$data %>%
    ungroup() %>%
    dplyr::select(-n, -tot_by_bin) %>%
    group_by(bin) %>%
    mutate(across(!matches(exclude_cols), ~mean_or_mode(.x))) %>%
    distinct() %>%
    mutate(n=100) %>% ungroup() %>%
    {if(!is.null(cvar)) mutate(., !!cvar:=mean(mod$data[[cvar]])) else .} %>%
    complete(nesting(!!!syms(grp),bin),fill=list(n=100))

  if(class(mod$model)=='brmsfit'){

    preda <- preda_data %>%
      add_linpred_draws(mod$model, allow_new_levels=T) %>%
      ungroup()

  } else if(class(mod$model)=='sdmTMB'){

    # For spatial models, add coordinate columns using mean coordinates
    if(mod$spatial_varying && !is.null(mod$coords)){
      for(coord in mod$coords){
        if(!coord %in% names(preda_data) && coord %in% names(mod$data)){
          preda_data[[coord]] <- mean(mod$data[[coord]], na.rm=TRUE)
        }
      }
    }

    preda <- predict(mod$model,
                     newdata = preda_data,
                     offset = preda_data$n,
                     type = c("link"),
                     nsim = mod$nsim) %>%
      as.data.frame() %>%
      cbind(preda_data) %>%
      mutate(.row = 1:n()) %>%
      pivot_longer(cols=matches("V[0-9]+"), names_to = "iter", values_to = ".linpred") %>%
      mutate(.draw = as.numeric(gsub("V","",iter))) %>%
      ungroup()

  } else {stop("Model must be either brmsfit or sdmTMB")}

  # this is needed for when there are helper variables to run the splines over categories. We only want the combos of helpers
  # if(helpers){
  #
  #   preds <- c()
  #   for (v in helper_vars)  preds <- bind_rows(preds, preda[preda$bin==v, ] %>% filter(!!sym(v) == 1))
  #
  #   preds <- bind_rows(preds, preda[!preda$bin %in% helper_vars, ] %>% filter(across(helper_vars, ~. != 1)))
  #   preda <- preds
  # }
  #browser()
  preda <- preda %>%
    group_by(!!!syms(grps),.draw) %>%
    mutate(lp = exp(.linpred)/sum(exp(.linpred))) %>%
    group_by(bin) %>%
    mutate(ml=gmean(lp)) %>%
    group_by(bin,!!!syms(grps)) %>%
    median_qi(pred = lp/ml)


  p <- ggplot(preda)  + geom_hline(yintercept = 1,linetype=2,alpha=0.5) +
    geom_pointrange(aes(x=bin,y=pred,ymin=.lower,ymax=.upper)) +
    xlab('Length bin (cm)') +
    ylab('Multiplier') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))
  #browser()
  if(is.null(grid))  p <- p+facet_wrap(facets = vars(!!grps[[1]]),scales = 'free_y')
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
scale_comps <- function(scale_df, 
  predvar='catch', 
  fit = NULL, 
  grps, 
  iters=NULL, 
  form=NULL, 
  pgrid = NULL # To do
){

  if(!is.null(form)){
    form_parts <- stringr::str_remove_all(stringr::str_split(form, '\\+')[[1]],pattern = ' ')
    form <- paste('~ (1|bin) +',paste0('(1|bin:',form_parts,')', collapse = ' + '))
  } else if(!is.null(fit)){
    form = paste("~",paste(paste0('(1|bin:',grps,')'), collapse=' + '))
  } else {stop("Must provide either a formula or a model fit from fit_model")}

  # Get bin levels from the model data
  if(fit$spatial_varying) grps <- c(grps, fit$coords)
  bin_levels <- unique(fit$data$bin)
  bin_is_factor <- is.factor(fit$data$bin)

  scale_dfs <- scale_df %>%
    group_by(across(all_of(grps) )) %>%
    summarize(n=sum(!!sym(predvar),na.rm=T), .groups = "drop") %>%
    mutate(bin = if(bin_is_factor) factor(bin_levels[1], levels=bin_levels) else bin_levels[1]) %>%
    ungroup() %>% # need to ungroup before augmenting
    complete(nesting(!!!syms(grps)),bin,fill = list(n=0))  %>%
    group_by(across(all_of(grps) )) %>%
    mutate(n = sum(n),
           bin = if(bin_is_factor) bin else as.numeric(as.character(bin))) %>%
    filter(n>0)


  if(class(fit$model)=='brmsfit'){

    add_predicted_draws(scale_dfs, fit$mod, allow_new_levels=T, ndraws = iters, value = 'tot_by_bin', re_formula = form)

  } else if(class(fit$model)=='sdmTMB'){

   ff <- gsub("\\(1\\|","",form)
    ff <- gsub("\\)","",ff)
    form_nrfx <- as.formula(paste("tot_by_bin",ff))
    lt <- labels(terms(form_nrfx))

    for(col in lt){
      scale_dfs[,col] <- paste(scale_dfs$bin, scale_dfs[[gsub('bin:','',col)]], sep = '_')
    }

    # Convert to factors but preserve coordinates, time, and numeric bin
    exclude_pattern_tmb <- '.*_by_.*|^n$'
    if(!bin_is_factor) exclude_pattern_tmb <- paste(exclude_pattern_tmb, '^bin$', sep='|')
    if(!is.null(fit$coords)) exclude_pattern_tmb <- paste(exclude_pattern_tmb, paste0('^', fit$coords, '$', collapse='|'), sep='|')
    if(!is.null(fit$time)) exclude_pattern_tmb <- paste(exclude_pattern_tmb, paste0('^', fit$time, '$'), sep='|')
    scale_dfs <- scale_dfs %>%
      mutate(across(-matches(exclude_pattern_tmb), as.factor))
  
    if(fit$spatial_varying){
      svc_setup <- sdmTMB::make_category_svc(
        data = scale_dfs,
        category_column = "bin",
        time_column = fit$time,
        share_spatial_sd = TRUE,
        share_spatiotemporal_sd = TRUE
      )
    }

    predict(fit$model, type = "response",
             newdata = if(fit$spatial_varying) svc_setup$data_expanded else scale_dfs,
             offset = log(scale_dfs$n),
             nsim = fit$nsim) %>%
      as.data.frame() %>%
      cbind(scale_dfs) %>%
      mutate(.row = 1:n()) %>%
      pivot_longer(cols=matches("V[0-9]+"), names_to = "iter", values_to = 'tot_by_bin') %>%
      mutate(.draw = as.numeric(gsub("V","",iter))) %>%
      ungroup()
  }
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
      mean_qi(prop) %>%
      filter(!is.na(prop)) %>%
      mutate(bin = as.numeric(as.character(bin)))

    if(!is.null(comp_are)){
      comp <- comp_are %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(tot = sum(!!sym(cvar)), .groups = "drop") %>%
        group_by(across(all_of(grps) )) %>%
        mutate(prop = tot/sum(tot)) %>%
        filter(!is.na(prop), prop > 0) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }

    if(!is.null(comp_are2)){
      comp2 <- comp_are2 %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(tot = sum(!!sym(cvar)), .groups = "drop") %>%
        group_by(across(all_of(grps) )) %>%
        mutate(prop = tot/sum(tot)) %>%
        filter(!is.na(prop), prop > 0) %>%
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
