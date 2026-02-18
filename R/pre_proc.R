#' Preparing data
#'
#' This function takes a dataset of counts (number of samples) per bin
#' (length/age/etc) across a range of strata (vars_for_grouping) and expands
#' it for use in the poisson-multinomial model
#' @param comp A compositional dataset with columns for counts (can be names differently) per bin (can be named differently), and all `vars_for_grouping`
#' @param vars_for_grouping A character vector of variables used for grouping the compositional dataframe
#' @param bin_lab A string with the column label for the bins if not "bin"
#' @param count_lab A string with the column label for the counts by bins if not "counts"
#' @export
#' @import dplyr
#' @import rlang

data_prep <- function(
    comp = NULL,
    vars_for_grouping = NULL,
    bin_lab = "bin",
    count_var = 'count'
){

  counts <- sym(count_var)
  bin <- sym(bin_lab)
  ## bin counts by area and year (0 counts are included)
  comp %>% group_by(across(all_of(c(vars_for_grouping, bin_lab)))) %>%
    summarize(tot_by_bin=sum(!!counts,na.rm=T)) %>%
    ungroup() %>% # need to ungroup before augmenting
    complete(nesting(!!!syms(vars_for_grouping)),!!bin,fill = list(tot_by_bin=0))  %>%
    group_by(across(all_of(vars_for_grouping) )) %>%
    mutate(n = sum(tot_by_bin)) %>%
    rename(bin = bin_lab) %>%
    filter(n>0)
}

#' Parse a formula
#'
#' Appends a formula by attaching a bin interaction term
#' @param form A RHS formula for the poisson-multinomial model, without explicit mention of bins
#' @param backend the backend used for model fitting
#' @param cvars Continuous variable to be fitted as 2D spline with bin
#' @param knots knots (1D or 2D) for cvar spline; as character for now
#' @param spl 1D spline over cvar per bin, or 2D spline over bins and cvar; as character for now
#' @param coords Character vector of coordinate column names to preserve (TMB backend only)
#' @param time Time column name to preserve (TMB backend only)
#' @import stringr
#' @import recipes
#' @export
#'
parse_form <- function(data,
                       backend = 'brms',
                       form = 'gear + area + area:yy',
                       cvars=NULL,
                       spl="2D",
                       knots='c(5,10)',
                       coords=NULL,
                       time=NULL
                       ){

  form_parts <- stringr::str_remove_all(stringr::str_split(form, '\\+')[[1]],pattern = ' ')
  newform <- paste('0 + bin +',paste0('(1|bin:',form_parts,')', collapse = ' + '))

  if(is.numeric(data$bin)) data$bin <- factor(data$bin, levels = sort(unique(data$bin)))
  if(!is.null(cvars) & spl=='2D') {
    data$cbin <- as.numeric(as.character(data$bin))
    for(i in 1:length(cvars)) newform <- paste(newform, paste("t2(",cvars[i], ", cbin, k = ",knots,")"), sep = ' + ')
  }
  if(!is.null(cvars) & spl=='by') {
    mm <- model.matrix(data=data, ~bin)
    binvars <- unique(data$bin)
    colnames(mm) <- binvars
    lb <- length(binvars)
    data <- cbind(data, mm[,2:lb])

    for(i in 1:length(cvars)) {

      for (b in binvars[2:lb]) newform <- paste(newform,
                                                paste("t2(",cvars[i], sprintf(", by=%s, k = ",b),knots,")"), sep = ' + ')

    }
  }

  print(newform)

  form <- switch(backend,
                 brms = as.formula(paste("tot_by_bin~offset(log(n)) + ",newform)),
                 TMB  = as.formula(paste("tot_by_bin~",newform))
  )

  # make sure all non-cvars are factors (preserve coordinates and time for TMB)
  exclude_pattern <- paste('.*_by_.*|^n$', cvars, sep='|')
  if(!is.null(coords)) exclude_pattern <- paste(exclude_pattern, paste0('^', coords, '$', collapse='|'), sep='|')
  exclude_pattern <- paste(exclude_pattern, paste0('^', time, '$'), sep='|')
  data <- data %>% ungroup() %>%
    mutate(across(-matches(exclude_pattern), as.factor))


  if(backend == 'TMB') {

    # to make interaction columns

    ff <- gsub("\\(1\\s\\|\\s","",form[3])
    ff <- gsub("\\)","",ff)
    form_nrfx <- as.formula(paste("tot_by_bin~",ff))
    lt <- labels(terms(form_nrfx))[-1]

    for(col in lt){
      data[,col] <- paste(data$bin, data[[gsub('bin:','',col)]], sep = '_')
    }

    # Convert to factors but preserve coordinates and time
    exclude_pattern_tmb <- '.*_by_.*|^n$'
    if(!is.null(coords)) exclude_pattern_tmb <- paste(exclude_pattern_tmb, paste0('^', coords, '$', collapse='|'), sep='|')
    if(!is.null(time)) exclude_pattern_tmb <- paste(exclude_pattern_tmb, paste0('^', time, '$'), sep='|')
    data <- data %>%
      mutate(across(-matches(exclude_pattern_tmb), as.factor))
  }

  list(form=form, data=data)

}

#' Fit the poisson multinomial model; see brms doc for brm() other parameters
#' @param form A RHS formula for the poisson-multinomial model, without explicit mention of bins
#' @param backend the backend used for model fitting
#' @param cvars Continuous variable to be fitted as 2D spline with bin
#' @param knots knots (1D or 2D) for cvar spline; as character for now
#' @param spl 1D spline over cvar per bin, or 2D spline over bins and cvar; as character for now
#' @param add_preds add predictions? (Can be very slow - only use with brms)
#' @param mesh An sdmTMB mesh object created with make_mesh() (TMB backend only). When provided, a spatially varying coefficient model is always used.
#' @param time Column name for time (year) variable (TMB backend only, required when mesh is provided)
#' @param share_spatial_sd Share spatial SD across bins (TMB backend only, default TRUE)
#' @param share_spatiotemporal_sd Share spatiotemporal SD across bins (TMB backend only, default TRUE)
#' @param coords Character vector of length 2 specifying coordinate column names (TMB backend only, default c("x", "y"))
#' @export
#' @import brms
#' @import sdmTMB
#' @import tidybayes
#'
fit_model <- function(form,
                      cvars=NULL,
                      spl='2D',
                      knots='c(5,10)',
                      data = NULL,
                      dist = 'poisson()',
                      chains=4,
                      threads = 16,
                      cores=32,
                      backend = 'brms',
                      brms_backend = 'rstan',
                      iter = 1000,
                      warmup = 500,
                      mtd=10,
                      ad=0.8,
                      thin = 2,
                      refresh=1000,
                      add_preds = FALSE,
                      mesh = NULL,
                      time = NULL,
                      share_spatial_sd = TRUE,
                      share_spatiotemporal_sd = TRUE,
                      coords = c("x", "y")){

  # Determine if spatial model based on mesh

  spatial_varying <- !is.null(mesh)

  # get formula and data in the right format for the model backend
  df <- parse_form(data = data,
                   backend = backend,
                   form = form,
                   cvars = cvars,
                   spl = spl,
                   knots = knots,
                   coords = if(backend == 'TMB' && spatial_varying) coords else NULL,
                   time = if(backend == 'TMB') time else NULL)

  form <- df$form
  data <- df$data

  fixed <- paste0("bin",unique(data$bin)[1])
  pp <-  set_prior("student_t(10, 0, 1)", class = "b") #+
   # set_prior("constant(0)", class = "b", coef = fixed)

  if(backend == 'brms'){
    mod <- brm(bf(form),
               family = dist,
               prior = pp,
               data = data,
               backend = brms_backend,
               chains=chains,
               threads = threads,
               cores=cores,
               iter = iter,
               warmup = warmup,
               thin = thin,
               refresh=refresh,
               control = list(max_treedepth=mtd, adapt_delta=ad)
    )


  cat('Converged: ', all(brms::rhat(mod)<1.05))

  } else if(backend == 'TMB'){

    if(spatial_varying){
      # Spatial model: always uses spatial varying coefficients
      if(is.null(time)) stop("'time' must be provided for spatial models")

      svc_setup <- sdmTMB::make_category_svc(
        data = data,
        category_column = "bin",
        time_column = time,
        share_spatial_sd = share_spatial_sd,
        share_spatiotemporal_sd = share_spatiotemporal_sd
      )

      control_sdmTMB <- sdmTMBcontrol(
        map = svc_setup$svc_map,
        parallel = parallel::detectCores(),
        multiphase = FALSE,
        newton_loops = 0,
        profile = "b_j",
        getsd = TRUE)

      mod <- sdmTMB(
        data = svc_setup$data_expanded,
        offset = log(data$n),
        formula = form,
        family = poisson(),
        mesh = mesh,
        spatial = "off",
        spatiotemporal = "off",
        time = time,
        spatial_varying = svc_setup$svc_formula,
        control = control_sdmTMB
      )
    } else {
      # Non-spatial model
      mod <- sdmTMB(
        data = data,
        offset = log(data$n),
        formula = form,
        family = poisson(),
        spatial = "off"
      )
    }
  }

  if(add_preds) pp <- switch(backend,
                             brms = add_predicted_draws(data, mod),
                             TMB = simulate(mod,
                                            nsim = (iter-warmup)/thin,
                                            type = "mle-mvn",
                                            mle_mvn_samples = "multiple"
                                            ) %>%
                               as.data.frame() %>%
                               cbind(data) %>%
                               mutate(.row = 1:n()) %>%
                               pivot_longer(cols=matches("V[0-9]+"), names_to = "iter", values_to = ".prediction") %>%
                               mutate(.draw = as.numeric(gsub("V","",iter)))
  )   else pp <- NULL

  out <- list(model = mod,
       data=data,
       nsim = (iter-warmup)/thin,
       preds = pp,
       mesh = mesh,
       spatial_varying = spatial_varying,
       time = time,
       coords = if(spatial_varying) coords else NULL)

  class(out) <- 'ComPoM'
  attr(out, 'call') <- match.call()
  out
}



