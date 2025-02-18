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
#' @param cvars Continuous variable to be fitted as 2D spline with bin
#' @param knots knots (1D or 2D) for cvar spline; as character for now
#' @import stringr
#'
parse_form <- function(data, form = 'gear + area + area:yy', cvars=NULL, spl="2D", knots='c(5,10)'){

  form_parts <- stringr::str_remove_all(stringr::str_split(form, '\\+')[[1]],pattern = ' ')
  newform <- paste('0 + bin +',paste0('(1|bin:',form_parts,')', collapse = ' + '))

  if(!is.null(cvars) & spl=='2D') for(i in 1:length(cvars)) newform <- paste(newform, paste("t2(",cvars[i], ", bin, k = ",knots,")"), sep = ' + ')
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

  form <- as.formula(paste("tot_by_bin~offset(log(n)) + ",newform))
  list(form=form, data=data)

}

#' Fit the poisson multinomial model; see brms doc for brm() other parameters
#' @param form A RHS formula for the poisson-multinomial model, without explicit mention of bins
#' @param cvars Continuous variable to be fitted as 2D spline with bin
#' @param knots knots (1D or 2D) for cvar spline; as character for now
#' @param add_preds add predictions? (Can be very slow)
#' @export
#' @import brms
#' @import tidybayes
#'
fit_model <- function(form,
                      cvars=NULL,
                      spl='2D',
                      knots='c(5,10)',
                      data = NULL,
                      dist = 'poisson',
                      chains=4,
                      threads = 16,
                      cores=32,
                      backend = 'cmdstanr',
                      iter = 1000,
                      warmup = 500,
                      mtd=10,
                      ad=0.8,
                      thin = 2,
                      refresh=10,
                      add_preds = TRUE){

  df <- parse_form(data, form, cvars, spl, knots)

  form <- df$form
  data <- df$data

  fixed <- paste0("bin",unique(data$bin)[1])
  pp <-  set_prior("student_t(10, 0, 1)", class = "b") +
    set_prior("constant(0)", class = "b", coef = fixed)

  mod <- brm(bf(form),
             family = dist,
             prior = pp,
             data = data,
             backend = backend,
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

  if(add_preds) pp <- add_predicted_draws(data, mod) else pp <- NULL

  list(model = mod,
       data=data,
       preds = pp)

}



