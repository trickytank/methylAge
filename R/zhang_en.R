#' Zhang DNA methylation age calculation
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
zhang_en <- function(x, id_col = "ID", age_col = "zhang_en_mAge", dim_warning = TRUE) {
  # Calculate Zhang Elastic Net Methylation Age
  check_methylation_data(x, dim_warning = dim_warning)

  print("1. Data loading and QC")

  print("1.1 Reading the data")
  readRDS(infile)-> data        ########## IND * Probe, each row represents one individual, it should be "RAW BETA" DNA methylation value

  if(nrow(data) > ncol(data)){
    print("I guess you are using Probe in the row, data will be transformed!!!")
    data<-t(data)
  }

  print("1.2 Replacing missing values with mean value")
  dataNona<-apply(data,2,function(x) addna(x))   ###############  replace the NA with mean value for each probe


  print("1.3 Standardizing")
  dataNona.norm<- apply(dataNona,1,scale)        ############### standardize the DNA methylation within each individual, remove the mean and divided by the SD of each individual     Probe * IND
  rownames(dataNona.norm)<-colnames(dataNona)


  ############# 3. get the coefficients of each probe from Elastic Net/BLUP method, !!!!WE HAVE TWO PREDICTORS!!!#############

  print("2. Loading predictors")
  read.table("en.coef",stringsAsFactor=F,header=T)->encoef
  read.table("blup.coef",stringsAsFactor=F,header=T)->blupcoef

  en_int<-encoef[1,2]
  blup_int<-blupcoef[1,2]

  encoef<-encoef[-1,]
  blupcoef<-blupcoef[-1,]

  rownames(encoef)<-encoef$probe
  rownames(blupcoef)<-blupcoef$probe

  ############# 4. get common probes between predictors and data ##############
  print("3. Checking misssing probes")

  encomm<- intersect(rownames(encoef),rownames(dataNona.norm))
  blupcomm<- intersect(rownames(blupcoef),rownames(dataNona.norm))

  endiff<- nrow(encoef) - length(encomm)
  blupdiff<- nrow(blupcoef) - length(blupcomm)

  print(paste0(endiff," probe(s) in Elastic Net predictor is(are) not in the data"))
  print(paste0(blupdiff," probe(s) in BLUP predictor is(are) not in the data"))
  print("BLUP can perform better if the number of missing probes is too large!")

  ############# 5. extract the common probes and do age prediction ###############
  print("4. Predicting")

  encoef<-encoef[encomm,]
  blupcoef<-blupcoef[blupcomm,]
  encoef$coef%*%dataNona.norm[encomm,]+en_int->enpred
  blupcoef$coef%*%dataNona.norm[blupcomm,]+blup_int->blupred


  ############# 6. Save the predicted result ###########
  read.table(agefile,header=T,stringsAsFactor=F)->age.raw
  enpred<-enpred[,age.raw$ID]
  blupred<-blupred[,age.raw$ID]

  age.raw$enpred<-as.double(enpred)
  age.raw$blupred<-as.double(blupred)


  # Load coefficients
  coefs_clock <- setNames(zhang_en_coef$Coefficient, zhang_en_coef$Marker)
  # Calculate Methylation Age
  x_mat <- as.matrix( x[names(coefs_clock), ] )
  m_age <- coefs_clock %*% x_mat
  # Cleanup
  tibble(!!id_col := colnames(m_age), !!age_col := m_age[1,])
}


addna<-function(methy){
  methy[is.na(methy)]<-mean(methy,na.rm=T)
  return(methy)
}
