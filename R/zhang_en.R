#' Zhang DNA methylation age calculation
#'
#' Calculate the \insertCite{zhang2019;textual}{methylAge} DNA methylation ages with the
#' Elastic Net (EN) and Best Linear Unbiased Prediction (BLUP) clocks.
#' Please cite the referenced article if using this function.
#'
#' @inheritParams generic_clock
#' @param en_out Column name for the Elastic Net DNA methylation age estimate.
#' @param blup_out Column name for the Best Linear Unbiased Prediction (BLUP) DNA methylation age estimate.
#' @param clock The clock/s to calculate from Zhang et al. The Elastic Net will be calculated with "en" and the BLUP will be calculated with "blup". Both may be included as a vector to calculate both.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
zhang_clock <- function(x, id_out = "ID",
                         en_out = "zhang_en_mage", blup_out = "zhang_blup_mage",
                         clock = c("en", "blup"), dim_warning = TRUE) {
  # Calculate Zhang Elastic Net Methylation Age
  check_methylation_data(x, dim_warning = dim_warning)

  data <- t(x) # Temp

  message("1.1 Replacing missing values with mean value")
  if(anyNA(data)) {
    data<-apply(data,2,function(x) addna(x))   ###############  replace the NA with mean value for each probe
  }

  message("1.2 Standardizing")
  dataNona.norm<- apply(data,1,scale)        ############### standardize the DNA methylation within each individual, remove the mean and divided by the SD of each individual     Probe * IND
  rownames(dataNona.norm)<-colnames(data)
  rm(data)


  ############# 3. get the coefficients of each probe from Elastic Net/BLUP method, !!!!WE HAVE TWO PREDICTORS!!!#############

  message("2. Loading predictors")
  encoef <- zhang_en_coef
  blupcoef <- zhang_blup_coef

  en_int<-encoef[1,2]
  blup_int<-blupcoef[1,2]

  encoef<-encoef[-1,]
  blupcoef<-blupcoef[-1,]

  rownames(encoef)<-encoef$marker
  rownames(blupcoef)<-blupcoef$marker

  ############# 4. get common probes between predictors and data ##############
  message("3. Checking misssing probes")

  encomm<- intersect(rownames(encoef),rownames(dataNona.norm))
  blupcomm<- intersect(rownames(blupcoef),rownames(dataNona.norm))

  endiff<- nrow(encoef) - length(encomm)
  blupdiff<- nrow(blupcoef) - length(blupcomm)

  message(paste0(endiff," probe(s) in Elastic Net predictor is(are) not in the data"))
  message(paste0(blupdiff," probe(s) in BLUP predictor is(are) not in the data"))
  message("BLUP can perform better if the number of missing probes is too large!")

  ############# 5. extract the common probes and do age prediction ###############
  message("4. Predicting")

  encoef<-encoef[encomm,]
  blupcoef<-blupcoef[blupcomm,]
  encoef$coefficient%*%dataNona.norm[encomm,]+en_int->enpred
  blupcoef$coefficient%*%dataNona.norm[blupcomm,]+blup_int->blupred


  ## # Load coefficients
  ## coefs_clock <- setNames(zhang_en_coef$Coefficient, zhang_en_coef$Marker)
  ## # Calculate Methylation Age
  ## x_mat <- as.matrix( x[names(coefs_clock), ] )
  ## m_age <- coefs_clock %*% x_mat
  ## # Cleanup
  ## tibble(!!id_out := colnames(m_age), !!age_out := m_age[1,])

  list(en = enpred, blup = blupred)
  if("en" %in% clock) {
    age_pred <- tibble::tibble(!!id_out := colnames(enpred), !!en_out := enpred[1,])
  } else {
    age_pred <- tibble::tibble(!!id_out := colnames(blupred))
  }
  if("blup" %in% clock) {
    age_pred <- tibble::add_column(age_pred, !!blup_out := blupred[1,])
  }
  age_pred
}


addna<-function(methy){
  methy[is.na(methy)]<-mean(methy,na.rm=T)
  return(methy)
}
