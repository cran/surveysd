#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data with rotating panel
#'   design. Survey information, like ID, sample weights, strata and population
#'   totals per strata, should be specified to ensure meaningfull survey
#'   bootstraping.
#'
#' @param dat either data.frame or data.table containing the survey data with
#'   rotating panel design.
#' @param REP integer indicating the number of bootstrap replicates.
#' @param hid character specifying the name of the column in `dat` containing
#'   the household id. If `NULL` (the default), the household structure is not
#'   regarded.
#' @param weights character specifying the name of the column in `dat`
#'   containing the sample weights.
#' @param period character specifying the name of the column in `dat` containing
#'   the sample periods. If `NULL` (the default), it is assumed that all
#'   observations belong to the same period.
#' @param strata character vector specifying the name(s) of the column in `dat`
#'   by which the population was stratified. If `strata` is a vector
#'   stratification will be assumed as the combination of column names contained
#'   in `strata`. Setting in addition `cluster` not NULL stratification will be
#'   assumed on multiple stages, where each additional entry in `strata`
#'   specifies the stratification variable for the next lower stage. see Details
#'   for more information.
#' @param cluster character vector specifying cluster in the data. If not
#'   already specified in `cluster` household ID is taken es the lowest level
#'   cluster.
#' @param totals character specifying the name of the column in `dat` containing
#'   the the totals per strata and/or cluster. Is ONLY optional if `cluster` is
#'   `NULL` or equal `hid` and `strata` contains one columnname! Then the
#'   households per strata will be calcualted using the `weights` argument. If
#'   clusters and strata for multiple stages are specified `totals` needs to be
#'   a vector of `length(strata)` specifying the column on `dat` that contain
#'   the total number of PSUs at each stage. `totals` is interpreted from left
#'   the right, meaning that the first argument corresponds to the number of
#'   PSUs at the first and the last argument to the number of PSUs at the last
#'   stage.
#' @param single.PSU either "merge" or "mean" defining how single PSUs need to
#'   be dealt with. For `single.PSU="merge"` single PSUs at each stage are
#'   merged with the strata or cluster with the next least number of PSUs. If
#'   multiple of those exist one will be select via random draw. For
#'   `single.PSU="mean"` single PSUs will get the mean over all bootstrap
#'   replicates at the stage which did not contain single PSUs.
#' @param boot.names character indicating the leading string of the column names
#'   for each bootstrap replica. If NULL defaults to "w".
#' @param split logical, if TRUE split households are considered using `pid`,
#'   for more information see Details.
#' @param pid column in `dat` specifying the personal identifier. This
#'   identifier needs to be unique for each person throught the whole data set.
#' @param new.method logical, if TRUE bootstrap replicates will never be
#'   negative even if in some strata the whole population is in the sample.
#'   WARNING: This is still experimental and resulting standard errors might be
#'   underestimated! Use this if for some strata the whole population is in the
#'   sample!
#'
#' @return the survey data with the number of REP bootstrap replicates added as
#'   columns.
#'
#' @details `draw.bootstrap` takes `dat` and draws `REP` bootstrap replicates
#' from it.
#' `dat` must be household data where household members correspond to multiple
#' rows with the same household
#' identifier. For most practical applications, the following columns should be
#' available in the dataset
#' and passed via the corresponding parameters:
#'
#' * Column indicating the sample period (parameter `period`).
#' * Column indicating the household ID (parameter `hid`)
#' * Column containing the household sample weights (parameter `weights`);
#' * Columns by which population was stratified during the sampling process
#'   (parameter: `strata`).
#'
#' For single stage sampling design a column the argument `totals` is optional,
#' meaning that a column of the number of PSUs at the first stage does not need
#' to be supplied.
#' For this case the number of PSUs is calculated and added to `dat` using
#' `strata` and `weights`. By setting `cluster` to NULL single stage sampling
#' design is always assumed and
#' if `strata` contains of multiple column names the combination of all those
#' column names will be used for stratification.
#'
#' In the case of multi stage sampling design the argument `totals` needs to be
#' specified and needs to have the same number of arguments as `strata`.
#'
#' If `cluster` is `NULL` or does not contain `hid` at the last stage, `hid`
#' will automatically be used as the final cluster. If, besides `hid`,
#' clustering in additional stages is specified the number of column names in
#' `strata` and `cluster` (including `hid`) must be the same. If for any stage
#' there was no clustering or stratification one can set "1" or "I" for this
#' stage.
#'
#' For example `strata=c("REGION","I"),cluster=c("MUNICIPALITY","HID")` would
#' speficy a 2 stage sampling design where at the first stage the municipalities
#' where drawn stratified by regions
#' and at the 2nd stage housholds are drawn in each municipality without
#' stratification.
#'
#' Bootstrap replicates are drawn for each survey period (`period`) using the
#' function [rescaled.bootstrap].
#' Afterwards the bootstrap replicates for each household are carried forward
#' from the first period the household enters the survey to all the censecutive
#' periods it stays in the survey.
#'
#' This ensures that the bootstrap replicates follow the same logic as the
#' sampled households, making the bootstrap replicates more comparable to the
#' actual sample units.
#'
#' If `split` ist set to `TRUE` and `pid` is specified, the bootstrap replicates
#' are carried forward using the personal identifiers instead of the houshold
#' identifier.
#' This takes into account the issue of a houshold splitting up.
#' Any person in this new split household will get the same bootstrap replicate
#' as the person that has come from an other household in the survey.
#' People who enter already existing households will also get the same bootstrap
#' replicate as the other households members had in the previous periods.
#'
#' @return Returns a data.table containing the original data as well as the
#'   number of `REP` columns containing the bootstrap replicates for each
#'   repetition.\cr
#'   The columns of the bootstrap replicates are by default labeled "w*Number*"
#'   where *Number* goes from 1 to `REP`. If the column names of the bootstrap
#'   replicates should start with a different character or string the parameter
#'   `boot.names` can be used.
#'
#' @seealso [`data.table`][data.table::data.table] for more information on
#'   data.table objects.
#'
#' @author Johannes Gussenbauer, Alexander Kowarik, Statistics Austria
#'
#' @examples
#' \dontrun{
#' eusilc <- demo.eusilc(prettyNames = TRUE)
#'
#' ## draw sample without stratification or clustering
#' dat_boot <- draw.bootstrap(eusilc, REP = 10, weights = "pWeight",
#'                            period = "year")
#'
#' ## use stratification w.r.t. region and clustering w.r.t. households
#' dat_boot <- draw.bootstrap(
#'   eusilc, REP = 10, hid = "hid", weights = "pWeight",
#'   strata = "region", period = "year")
#'
#' ## use multi-level clustering
#' dat_boot <- draw.bootstrap(
#'   eusilc, REP = 10, hid = "hid", weights = "pWeight",
#'   strata = c("region", "age"), period = "year")
#'
#'
#' # create spit households
#' eusilc[, pidsplit := pid]
#' year <- eusilc[, unique(year)]
#' year <- year[-1]
#' leaf_out <- c()
#' for(y in year) {
#'   split.person <- eusilc[
#'     year == (y-1) & !duplicated(hid) & !(hid %in% leaf_out),
#'     sample(pid, 20)
#'   ]
#'   overwrite.person <- eusilc[
#'     (year == (y)) & !duplicated(hid) & !(hid %in% leaf_out),
#'     .(pid = sample(pid, 20))
#'   ]
#'   overwrite.person[, c("pidsplit", "year_curr") := .(split.person, y)]
#'
#'   eusilc[overwrite.person, pidsplit := i.pidsplit,
#'          on = .(pid, year >= year_curr)]
#'   leaf_out <- c(leaf_out,
#'                 eusilc[pid %in% c(overwrite.person$pid,
#'                                   overwrite.person$pidsplit),
#'                 unique(hid)])
#' }
#'
#' dat_boot <- draw.bootstrap(
#'   eusilc, REP = 10, hid = "hid", weights = "pWeight",
#'   strata = c("region", "age"), period = "year", split = TRUE,
#'   pid = "pidsplit")
#' # split households were considered e.g. household and
#' # split household were both selected or not selected
#' dat_boot[, data.table::uniqueN(w1), by = pidsplit][V1 > 1]
#' }
#'
#' @export draw.bootstrap
#'


draw.bootstrap <- function(
  dat, REP = 1000, hid = NULL, weights, period = NULL, strata = NULL,
  cluster = NULL, totals = NULL, single.PSU = c("merge", "mean"), boot.names =
    NULL, split = FALSE, pid = NULL, new.method = FALSE) {

  occurence_first_period <- NULL

  ##########################################################
  # INPUT CHECKING
  if (is.data.frame(dat)) {
    dat <- as.data.table(dat)
  }else if (!is.data.table(dat)) {
    stop("dat must be a data.frame or data.table")
  }
  dat <- copy(dat)

  c.names <- colnames(dat)

  removeCols <- c()
  # check REP
  if (length(REP) != 1) {
    stop("REP must have length 1")
  }
  if (!is.numeric(REP)) {
    stop("REP must contain one numeric value")
  }
  if (REP %% 1 != 0) {
    stop("REP cannot have a decimal part")
  }

  # check hid
  hidNULL <- is.null(hid)
  if (hidNULL) {
    hid <- generateRandomName(20, colnames(dat))
    dat[, c(hid) := 1:.N]
    removeCols <- c(removeCols, hid)
  }

  if (length(hid) != 1) {
    stop("hid must have length 1")
  }
  if (!hid %in% c.names) {
    stop(paste0(hid, " is not a column in dat"))
  }

  # check weights
  if (length(weights) != 1) {
    stop("weights must have length 1")
  }
  if (!weights %in% c.names) {
    stop(paste0(weights, " is not a column in dat"))
  }
  if (!is.numeric(dt.eval("dat[,", weights, "]"))) {
    stop(paste0(weights, " must be a numeric column"))
  }

  # check period
  periodNULL <- is.null(period)
  if (periodNULL) {
    period <- generateRandomName(20, colnames(dat))
    dat[, c(period) := 1]
    removeCols <- c(removeCols, period)
  }

  if (length(period) != 1) {
    stop("period must have length 1")
  }
  if (!period %in% c.names) {
    stop(paste0(period, " is not a column in dat"))
  }
  if (!is.numeric(dat[[period]]) & !is.integer(dat[[period]])) {
    stop(paste0(period, " is not an integer or numeric column"))
  }

  # check design
  if (is.null(strata)) {
    strata <- "I"
  }

  if (is.null(cluster)) {
    cluster <- hid
  } else {
    if (length(cluster) == 1) {
      if (cluster %in% c("1", "I")) {
        cluster <- hid
      }
    }
    if (!hid %in% cluster) {
      cluster <- c(cluster, hid)
    }
  }

  if (!all(strata[!strata %in% c("1", "I")] %in% c.names)) {
    stop("Not all elements in strata are column names in dat")
  }
  if (any(!cluster[!cluster %in% c("1", "I")] %in% c.names)) {
    stop("Not all names in cluster are column names in dat")
  }
  if (any(!totals %in% c.names)) {
    stop("Not all names in totals are column names in dat")
  }

  # check for missing values
  spec.variables <- c(hid, weights, period, strata, cluster, totals, pid)
  spec.variables <- spec.variables[!spec.variables %in% c("1", "I")]
  dat.na <- dat[, mget(spec.variables)]
  dat.na <- sapply(dat.na, function(z) {
    any(is.na(z))
  })
  if (any(dat.na)) {
    stop("Missing values found in column(s): ",
         paste(names(dat.na[dat.na == TRUE]), collapse = ", "))
  }

  if (length(cluster) > 1) {
    if (length(cluster) != length(strata)) {
      stop("strata and cluster need to have the same number of stages!\n ",
           "Please use either '1' or 'I' if there was no clustering or ",
           "stratification in one of the stages.")
    }
  }else{

    if (length(strata) > 1) {
      if (any(c("1", "I") %in% strata)) {
        stop("When defining multiple strata variables for single stage",
             " sampling design\n none of them can be '1' or 'I'.")
      }
      strata_var_help <- generateRandomName(20, colnames(dat))
      dat[, c(strata_var_help) := do.call(paste, c(.SD, sep = "-")), ]
      strata <- strata_var_help
      removeCols <- c(removeCols, strata)
    }
  }

  # check single.PSUs
  single.PSU <- single.PSU[1]
  if (is.null(single.PSU) || !single.PSU %in% c("merge", "mean")) {
    warning("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap",
            " replicates for single PSUs cases will be missing!")
    single.PSU <- FALSE
  }

  # check boot.names
  if (!is.null(boot.names)) {
    if (!grepl("^[[:alpha:]]", boot.names)) {
      stop("boot.names must start with an alphabetic character")
    }
  }

  # check split and pid
  if (is.null(split)) {
    stop("split needs to be logical")
  }
  if (!is.logical(split)) {
    stop("split needs to be logical")
  }
  if (split) {
    if (!is.character(pid)) {
      stop("when split is TRUE pid needs to be a string")
    } else {
      if (length(pid) > 1) {
        stop("pid can only have length 1")
      } else {
        if (!pid %in% c.names) {
          stop(pid, "is not a column of dat")
        }
      }
    }
    # check if pid is unique in each household and period
    unique.pid <- dt.eval("dat[,uniqueN(", pid,
                          ")==.N,by=c(period,hid)][V1==FALSE]")
    if (nrow(unique.pid) > 0) {
      stop("pid is not unique in each household for each period")
    }
  }

  # check totals
  # if clusters are specified the finite population correction factors must
  #   be user specified
  # check input for totals
  # if no totals are specified then leave them NULL
  if (is.null(totals)) {
    if (length(cluster) == 1) {
      # if no clusters are specified calculate number of households in each
      #   strata
      totals <- generateRandomName(20, existingNames = colnames(dat))
      fpc.strata <- strata[!strata %in% c("I", "1")] # nolint
      dt.eval("dat[,", totals, ":=sum(", weights, "[!duplicated(",
              hid, ")]),by=c(fpc.strata,period)]")
      removeCols <- c(removeCols, totals)
    } else {

      stop("For multistage sampling the number of PSUs at each level needs to ",
           "be specified!")
    }

  } else {

    if (length(totals) != length(strata)) {
      stop("totals must be specified for each stage")
    }
    if (any(!totals %in% c.names)) {
      stop("Not all elements in totals are column names in dat")
    }
    if (!any(unlist(dat[, lapply(.SD, is.numeric), .SDcols = c(totals)]))) {
      stop("Not all elements in totals are numeric columns in dat")
    }
  }
  
  # check for each stage that PSUs are not in mutiple strata
  for(i in seq_along(strata)){
    if(!strata[i]%in%c("1","I") & !cluster[i]%in%c("1","I")){
      countMultiple <- dt.eval("dat[,uniqueN(",strata[i],"),by=c(cluster[i],period)][V1>1]")
      if(nrow(countMultiple)>0){
        stop("Some sampling units in ",cluster[i]," occur in multiple strata of ",strata[i])
      }
    }
  }
  
  ##########################################################

  # define sample design
  strata <- paste(strata, collapse = ">")
  cluster <- paste(cluster, collapse = ">")
  totals <- paste(totals, collapse = ">")

  if (is.null(boot.names)) {
    w.names <- paste0("w", 1:REP)
  } else{
    w.names <- paste0(boot.names, 1:REP)
  }

  # calculate bootstrap replicates
  dat[, c(w.names) := rescaled.bootstrap(
    dat = copy(.SD), REP = REP, strata = strata, cluster = cluster,
    fpc = totals, single.PSU = single.PSU, return.value = "replicates",
    check.input = FALSE, new.method = new.method), by = c(period)]

  # respect split households
  if (split) {
    dat <- generate.HHID(dat, period = period, pid = pid, hid = hid)
  }

  # keep bootstrap replicates of first period for each household
  # if households drops out of survey and re-enters survey at a later stage
  # -> treat as different instances
  help_survey_grouping <- function(x){
    x_holes <- which(diff(x)>1)
    x_groups <- diff(c(0,x_holes,length(x)))
    x_groups <- rep(1:length(x_groups),times=x_groups)
    return(x_groups)
  }
  dt.eval("dat[,hid_survey_segment_help:=help_survey_grouping(", period, "),by=c(hid)]")
  dt.eval("dat[,occurence_first_period :=min(", period, "),by=c(hid,'hid_survey_segment_help')]")
  select.first.occurence <- paste0(c(hid,"hid_survey_segment_help", w.names), collapse = ",")
  dat.first.occurence <- unique(
    dt.eval("dat[", period, "==occurence_first_period,.(",
            select.first.occurence, ")]"
    ), by = c(hid,"hid_survey_segment_help"))
  dat[, c(w.names) := NULL]
  dat <- merge(dat, dat.first.occurence, by = c(hid, "hid_survey_segment_help"), all.x = TRUE)
  dat[, c("occurence_first_period","hid_survey_segment_help") := NULL]


  # remove columns
  if (split) {
    dt.eval("dat[,", hid, ":=", paste0(hid, "_orig"), "]")
    dat[, c(paste0(hid, "_orig")) := NULL]
  }
  if (length(removeCols) > 0) {
    dat[, c(removeCols) := NULL]
  }


  if (periodNULL) {
    period <- NULL
  }
  if (hidNULL) {
    hid <- NULL
  }

  setattr(dat, "weights", weights)
  setattr(dat, "period", period)
  setattr(dat, "b.rep", w.names)
  setattr(dat, "hid", hid)

  return(dat)
}
