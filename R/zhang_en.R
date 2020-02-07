#' Zhang DNA methylation age calculation
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
zhang_clocks <- function(x, id_col = "ID",
                         en_col = "zhang_en_mAge", blup_col = "zhang_blup_mAge",
                         clocks = c("en", "blup"), dim_warning = TRUE) {
  # Calculate Zhang Elastic Net Methylation Age
  check_methylation_data(x, dim_warning = dim_warning)

  data <- t(x) # Temp

  message("1.1 Replacing missing values with mean value")
  if(anyNA(data)) {
    dataNona<-apply(data,2,function(x) addna(x))   ###############  replace the NA with mean value for each probe
  } else {
    dataNona <- data
  }

  message("1.2 Standardizing")
  dataNona.norm<- apply(dataNona,1,scale)        ############### standardize the DNA methylation within each individual, remove the mean and divided by the SD of each individual     Probe * IND
  rownames(dataNona.norm)<-colnames(dataNona)


  ############# 3. get the coefficients of each probe from Elastic Net/BLUP method, !!!!WE HAVE TWO PREDICTORS!!!#############

  message("2. Loading predictors")
  encoef <- zhang_en_coef
  blupcoef <- zhang_blup_coef

  en_int<-encoef[1,2]
  blup_int<-blupcoef[1,2]

  encoef<-encoef[-1,]
  blupcoef<-blupcoef[-1,]

  rownames(encoef)<-encoef$probe
  rownames(blupcoef)<-blupcoef$probe

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
  encoef$coef%*%dataNona.norm[encomm,]+en_int->enpred
  blupcoef$coef%*%dataNona.norm[blupcomm,]+blup_int->blupred


  ## # Load coefficients
  ## coefs_clock <- setNames(zhang_en_coef$Coefficient, zhang_en_coef$Marker)
  ## # Calculate Methylation Age
  ## x_mat <- as.matrix( x[names(coefs_clock), ] )
  ## m_age <- coefs_clock %*% x_mat
  ## # Cleanup
  ## tibble(!!id_col := colnames(m_age), !!age_col := m_age[1,])

  list(en = enpred, blup = blupred)
  if("en" %in% clocks) {
    age_pred <- tibble::tibble(!!id_col := colnames(enpred), !!en_col := enpred[1,])
  } else {
    age_pred <- tibble::tibble(!!id_col := colnames(blupred))
  }
  if("blup" %in% clocks) {
    age_pred <- tibble::add_column(age_pred, !!blup_col := blupred[1,])
  }
  age_pred
}


addna<-function(methy){
  methy[is.na(methy)]<-mean(methy,na.rm=T)
  return(methy)
}
