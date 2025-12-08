#' @title Sort pedigree information
#' @description \code{orderped}  Sort pedigree information in ancestral order
#'
#' @importFrom dplyr semi_join anti_join pull %>%
#' @importFrom checkmate assert_names
#' @importFrom utils capture.output
#' @importFrom rlang .data
#' @param ped A data frame of pedigree information, comprised of three columns: IID, Seed, Pollen. IID column corresponds the individual name, and Seed and Pollen columns correspond to the mother and father name respectively. The idividuals with unknown parents are allowed, and unknown should be assigned a clear name (e.g. NA, "0", "unknown).
#' @param missingCode A character assigned for unknown parent. Default value is NA.
#' @return A data frame of pedigree information, sorted in ancestral order.
#' @export

orderPed <- function( ped, missingCode = NA ){

  assert_names(names(ped), must.include = c("IID", "Seed", "Pollen"))
  ped <- ped[ ! duplicated( ped ), ]

  ancestorCode = "0"

  if( is.na(missingCode) ){
    ped <- replace( ped, is.na(ped), ancestorCode )
  } else {
    ped <- replace( ped, ped == missingCode, ancestorCode )
  }

  parents <- ancestorCode
  exist_cross <- unique( data.frame( Var1 = ped$Seed, Var2 = ped$Pollen ) )
  cross_history <- data.frame( Var1 = character(0), Var2 = character(0) )

  n <- 1
  repeat{
    parents0grid <- expand.grid( parents, parents )
    parents0grid <- semi_join( parents0grid, exist_cross, by = c( "Var1", "Var2" ) )
    if( n > 1 ){
      parents0grid <- anti_join( parents0grid, cross_history, by = c( "Var1","Var2" ) )
    }
    offspring <- character(0)
    if( nrow( parents0grid ) == 0 ) break
    offspring <- ped %>%
      semi_join(parents0grid, by = c("Seed" = "Var1", "Pollen" = "Var2")) %>%
      pull(.data$IID) %>% unique()
    offspring <- unique( offspring )

    if( length( offspring ) == 0 ) break
    parents <- c( parents, offspring )
    cross_history <- rbind( cross_history, parents0grid )

    n <- n + 1

  }
  parents <- parents[ -1 ]
  ped <- ped[ ( ped$IID %in% parents ) , ]
  ped <- ped[ match( parents, ped$IID ), ]
  rownames( ped ) <- 1:nrow( ped )
  ped <- replace( ped, ped == ancestorCode, missingCode )
  return( ped )
}

#' @title Extract pedigree information
#' @description \code{selectPed} Extract pedigree information of specific individuals.
#'
#' @importFrom dplyr semi_join anti_join pull %>%
#' @importFrom checkmate assert_names
#' @importFrom utils capture.output
#' @param ped A data frame of pedigree information, comprised of three columns: IID, Seed, Pollen. IID column corresponds the individual name, and Seed and Pollen columns correspond to the mother and father name respectively. The idividuals with unknown parents are allowed, and unknown should be assigned a clear name (e.g. NA, "0", "unknown).
#' @param nameList A vector of the target individuals whose pedigree information are to be extracted.
#' @param missingCode A character assigned for unknown parent. Default value is NA.
#' @return A data frame of pedigree information for specific individuals
#' @export

selectPed <- function( ped, nameList, missingCode = NA ){

  assert_names(names(ped), must.include = c("IID", "Seed", "Pollen"))
  ped <- ped[ ! duplicated( ped ), ]

  ## Inner function to search ancestors of specific individual from pedigree file
  tracePed <- function( name ){

    if( is.na( missingCode ) ){
      is.missing <- is.na( name )
    } else {
      is.missing <- name == missingCode
    }

    if( is.missing ){
      return( NULL )
    } else {
      Ind <- ped[ ped$IID == name, ]
      if( nrow( Ind ) == 0 ){
        Ind <- data.frame( IID = name, Seed = missingCode, Pollen = missingCode, stringsAsFactors = FALSE )
      } else if( nrow( Ind ) > 1 ) {
        errorInd <- Ind[1,1]
        stop(sprintf(
          "Several crosses have matched to the IID '%s'.\nDetails:\n%s. \nPlease remove one of them.",
          errorInd,
          paste(capture.output(print(Ind)), collapse = "\n")
        ))
      }
      Male <- Ind$Pollen
      Female <- Ind$Seed

      return( rbind( Ind, tracePed( Male ), tracePed( Female ) ) )
    }
  }

  pedList <- lapply( nameList, tracePed )
  newPed <- do.call( rbind, pedList )
  newPed <- orderPed( newPed, missingCode = missingCode )
  rownames( newPed ) <- NULL
  return( newPed )

}


#' @title Calculate selection intensity
#' @description \code{calcSI}  Calculate selection intensity for all individuals
#'
#' @importFrom checkmate assert_data_frame assert_set_equal assert_character assert_numeric assert_matrix assert_subset
#' @importFrom furrr future_map
#' @importFrom future plan availableCores multisession
#' @importFrom stats setNames
#' @param Marker A data frame with 5 columns; MarkerID, Chr, Map, MarkerEff.a, and MarkerEff.d. MarkerID specifies the ID of each marker, and Chr and Map are corresponding to the chromosome number and map position, respectively. MarkerEff.a and MarkerEff.d are the effect of each markers calculated by genomic prediction model. When not considering dominance effect, please fill MarkerEff.d column with 0.
#' @param Pedigree A data frame of pedigree information, comprised of three columns: IID, Seed, Pollen. IID column corresponds the individual name, and Seed and Pollen columns correspond to the mother and father name respectively. The idividuals with unknown parents are allowed, and unknown should be assigned a clear name (e.g. NA, "0", "unknown).
#' @param genoPhased A matrix of phased genotypic data. Row and columns are corresponding to individuals and markers, respectively. Rownames must be a subset of IID in Pedigree, and Colnames must be a subset of MarkerID in Marker. Each cell must contain two numbers separated by a vertical bar (e.g. "0|1"). In this package, the left number are regarded as mathernal allele, and the right as paternal allele.
#' @param nCore Integer indicating the number of cores used for parallel calculation. If not specified, the maximum number of cores minus two is used.
#' @return A list of length three. First element of the list contains the result attributed to the additive effect, and second and third are the result attributed to the dominance and total effects, respectively. When considering only additive effect, length of the list become one. The result in each element are summarized in data frame and has following eight columns;
#' \describe{
#' \item{gEffect_s}{Genetic effect of seed parent}
#' \item{gEffect_p}{Genetic effect of pollen parent}
#' \item{gEffect_mid}{Genetic effect of mid-parent}
#' \item{gEffect_o}{Genetic effect of target individual}
#' \item{Vg_s}{Variance of expected progeny distribution, attributed to the seed parent}
#' \item{Vg_p}{Variance of expected progeny distribution, attributed to the pollen parent}
#' \item{Vg}{Variance of expected progeny distribution}
#' \item{Select_Intens}{Estimated selection intensity of target individual}
#' }
#' @export

calcSI <- function( Marker, Pedigree, genoPhased, nCore = NULL ){

  assert_data_frame(Marker)
  assert_set_equal( colnames( Marker ), c( "MarkerID", "Chr", "Map", "MarkerEff.a", "MarkerEff.d" ) )
  assert_character( Marker$MarkerID )
  assert_numeric( Marker$Chr )
  assert_numeric( Marker$Map )
  assert_numeric( Marker$MarkerEff.a )
  assert_numeric( Marker$MarkerEff.d )

  assert_data_frame(Pedigree)
  assert_set_equal( colnames( Pedigree ), c( "IID", "Seed", "Pollen" ) )
  assert_character( Pedigree$IID )
  assert_character( Pedigree$Seed )
  assert_character( Pedigree$Pollen )

  assert_matrix( genoPhased )
  assert_character( colnames( genoPhased ) )
  assert_character( rownames( genoPhased ) )

  assert_subset( colnames(  genoPhased ), Marker$MarkerID )

  if( all( Marker$MarkerEff.d == 0 ) ){
    Model <- "Add"
  } else {
    Model <- "Dom"
  }

  genoPhased <- genoPhased[ , match( Marker$MarkerID, colnames( genoPhased ) ) ]

  crossCache <- list()

  Indiv <- list()
  nIndiv <- length( Pedigree$IID )

  Effect <- if(Model == "Dom") c("Add","Dom","Tot") else c("Tot")
  Summary <- setNames( vector("list",length(Effect)), Effect )
  for( i in 1:length(Effect) ){
    Summary[[i]] <- list()
  }

  if( is.null( nCore ) ){
    nCore <- max( 1L, availableCores() - 2L )
  }
  plan(multisession, workers = nCore)

  for( i in 1:nIndiv ){

    cross <- Pedigree[ i, ]

    if( sum( unlist( cross )  %in% rownames( genoPhased ) ) == 3 ){

      key <- paste( cross$Seed, cross$Pollen, sep = "_")
      vgList <- crossCache[[key]]

      offspring <- genoPhased[ cross$IID, ]
      seed <- genoPhased[ cross$Seed, ]
      pollen <- genoPhased[ cross$Pollen, ]

      oDf <- data.frame( m = as.numeric( substr( offspring, 1, 1 ) ),
                         f = as.numeric( substr( offspring, 3, 3 ) ) )
      sDf <- data.frame( m = as.numeric( substr( seed, 1, 1 ) ),
                         f = as.numeric( substr( seed, 3, 3 ) ) )
      pDf <- data.frame( m = as.numeric( substr( pollen, 1, 1 ) ),
                         f = as.numeric( substr( pollen, 3, 3 ) ) )

      oDf$zygote <- oDf$m + oDf$f
      sDf$zygote <- sDf$m + sDf$f
      pDf$zygote <- pDf$m + pDf$f

      rownames( oDf ) <- rownames( sDf ) <- rownames( pDf ) <-  colnames( genoPhased )

      if( is.null(vgList) ){
        vgListNew <- list()

        nMarker <- dim( genoPhased )[2]
        map <- Marker$Map
        Chr <- Marker$Chr

        chrList <- split( seq_len(nMarker), Marker$Chr )
        combList <- lapply( chrList, function( ind ) {
          l <- length(ind)
          comb <- expand.grid( 1:l, 1:l )
          comb <- comb[ ( comb[,1] <= comb[,2] ), ]
          return( comb )
        } )

        varCovList <- future_map( as.numeric(names(chrList)), function(chr){

          calcCovMat_OneChr = function( idx, comb, map, SNPs) {
            i <- idx[comb[, 1]]
            j <- idx[comb[, 2]]
            d <- abs(map[i] - map[j])
            r <- 0.5 * (1 - exp(-2 * d))
            z <- SNPs[i, 1] * SNPs[j, 1] + SNPs[i, 2] * SNPs[j, 2] -
              SNPs[i, 1] * SNPs[j, 2] - SNPs[i, 2] * SNPs[j, 1]
            return( 0.25 * (1 - 2 * r) * z )
          }

          Sym_mat = function( upper ) {
            l <- ( - 1 + sqrt( 1 + 8 * length( upper ) ) ) / 2
            Mat <- matrix( 0, ncol = l, nrow = l )
            Mat[ upper.tri( Mat, diag = TRUE ) ] <- upper
            Mat <- Mat + t( Mat )
            diag(Mat) <- 1/2 * diag(Mat)
            return( Mat )
          }

          idx <- chrList[[chr]]
          comb <- combList[[chr]]

          varCov_s <- Sym_mat( calcCovMat_OneChr(idx, comb, map, sDf) )
          varCov_p <- Sym_mat( calcCovMat_OneChr(idx, comb, map, pDf) )

          zyg_s <- sDf$zygote[idx]
          zyg_p <- pDf$zygote[idx]

          varCov1 <- varCov_s + varCov_p
          varCov2 <- varCov_s * varCov_p + 1/4 * ( (zyg_s %o% zyg_s) * varCov_s + (zyg_p %o% zyg_p) * varCov_p )
          varCov3 <- 1/2 * (varCov_s * zyg_s[col(varCov_s)] + varCov_p * zyg_p[col(varCov_p)])

          markerName <- list( colnames( genoPhased )[idx], colnames( genoPhased )[idx] )
          dimnames( varCov_s ) <- dimnames( varCov_p ) <- dimnames( varCov1 ) <-  dimnames( varCov2 ) <- dimnames( varCov3 ) <- markerName

          return( list( varCov_s = varCov_s, varCov_p = varCov_p, varCov1 = varCov1, varCov2 = varCov2, varCov3 = varCov3 ) )
        } )

        vgListNew <- Reduce( `+`, future_map( as.numeric(names(chrList)), function(chr){
          idx <- chrList[[chr]]

          mEffectAdd <- Marker$MarkerEff.a[ idx, ]
          mEffectDom <- Marker$MarkerEff.a[ idx, ]
          mEffectSum <- mEffectAdd + mEffectDom

          gEffect_s <- sDf[ idx, ]$zygote %*% mEffectAdd
          gEffect_p <- pDf[ idx, ]$zygote %*% mEffectAdd
          gEffect_mid <- mean( c( gEffect_s, gEffect_p ) )
          gEffect_o <- oDf[ idx, ]$zygote %*% mEffectAdd

          gEffect_s.d <- ( sDf[ idx, ]$zygote %% 2 ) %*% mEffectDom
          gEffect_p.d <- ( pDf[ idx, ]$zygote %% 2 ) %*% mEffectDom
          gEffect_mid.d <- ( ( pDf[ idx, ]$zygote + sDf[ idx, ]$zygote - pDf[ idx, ]$zygote * sDf[ idx, ]$zygote )  %*% mEffectDom ) / 2
          gEffect_o.d <- ( oDf[ idx, ]$zygote %% 2 ) %*% mEffectDom

          gEffect_s.t <- gEffect_s + gEffect_s.d
          gEffect_p.t <- gEffect_p + gEffect_p.d
          gEffect_mid.t <- gEffect_mid + gEffect_mid.d
          gEffect_o.t <- gEffect_o + gEffect_o.d

          varCovSelect1 <- varCovList[[chr]]$varCov1[ idx, idx ]
          varCovSelect2 <- varCovList[[chr]]$varCov2[ idx, idx ]
          varCovSelect3 <- varCovList[[chr]]$varCov3[ idx, idx ]
          varCovSelect_s <- varCovList[[chr]]$varCov_s[ idx, idx ]
          varCovSelect_p <- varCovList[[chr]]$varCov_p[ idx, idx ]

          Vga <- t( mEffectAdd ) %*% varCovSelect1 %*% mEffectAdd
          Vgd <- t( mEffectDom ) %*% ( varCovSelect1 + 4 * varCovSelect2 - 4 * varCovSelect3 ) %*% mEffectDom
          Vgad <- 2 * t( mEffectAdd ) %*% varCovSelect1 %*% mEffectDom - 4 * t( mEffectAdd ) %*% varCovSelect3 %*% mEffectDom
          Vg <- Vga + Vgd + Vgad

          Vg_s <- t( mEffectAdd ) %*% varCovSelect_s %*% mEffectAdd
          Vg_p <- t( mEffectAdd ) %*% varCovSelect_p %*% mEffectAdd

          Vg_s.t <- t( mEffectSum) %*% varCovSelect_s %*% mEffectSum
          Vg_p.t <- t( mEffectSum ) %*% varCovSelect_p %*% mEffectSum

          return( c( gEffect_s, gEffect_p, gEffect_mid, gEffect_o, gEffect_s.d, gEffect_p.d, gEffect_mid.d, gEffect_o.d, gEffect_s.t, gEffect_p.t, gEffect_mid.t, gEffect_o.t, Vga, Vgd, Vgad, Vg, Vg_s, Vg_s.t, Vg_p, Vg_p.t ) )
        } )
        )
        crossCache[[key]] <- vgListNew

      } else {

        mEffectAdd <- Marker$MarkerEff.a
        mEffectDom <- Marker$MarkerEff.d

        gEffect_o <- oDf$zygote %*% mEffectAdd
        gEffect_o.d <- ( oDf$zygote %% 2 ) %*% mEffectDom
        gEffect_o.t <- gEffect_o + gEffect_o.d

        vgList[4] <- c( gEffect_o )
        vgList[8] <- c( gEffect_o.d )
        vgList[12] <- c( gEffect_o.t )
      }

      selectInt <- ( vgList[4] - vgList[3] ) / sqrt( vgList[13] )

      AddVgList <- c( vgList[1], vgList[2], vgList[3], vgList[4], vgList[17], vgList[19], vgList[13], selectInt )

      if( Model == "Add" ){
        Summary$Tot[[i]] <- AddVgList
      } else {
        selectInt.d <- ( vgList[8] - vgList[7] ) / sqrt( vgList[14] )
        selectInt.t <- ( vgList[12] - vgList[11] ) / sqrt( vgList[16] )

        DomVgList <- c( vgList[5], vgList[6], vgList[7], vgList[8], NA, NA, vgList[14], selectInt.d )
        TotVgList <- c( vgList[9], vgList[10], vgList[11], vgList[12], vgList[18], vgList[20], vgList[16], selectInt.t )

        Summary$Add[[i]] <- AddVgList
        Summary$Dom[[i]] <- DomVgList
        Summary$Tot[[i]] <- TotVgList
      }

    } else {
      if( Model == "Add" ){
        Summary$Tot[[i]] <- rep( NA, 8 )
      } else {
        Summary$Add[[i]] <- rep( NA, 8 )
        Summary$Dom[[i]] <- rep( NA, 8 )
        Summary$Tot[[i]] <- rep( NA, 8 )
      }
    }
    cat( paste( i, "/", nIndiv, " Completed\n", sep = "" ) )
  }

  for( i in 1:length(Effect) ){
    Summary[[i]] <- do.call( rbind, Summary[[i]] )
    rownames( Summary[[i]] ) <- Pedigree$IID
    colnames( Summary[[i]] ) <- c( "gEffect_s", "gEffect_p", "gEffect_mid", "gEffect_o", "Vg_s", "Vg_p", "Vg", "Select_Intens" )
  }

  return( Summary )

}


#' @title Calculate marker effects
#' @description \code{calcMarkerEff}  Calculate marker effects from either "Ridge", "Lasso", "BRR", "BayesB", or "BayesC". Useder can also choose the effects of marker to be included (additive only or additive + dominance ).
#'
#' @importFrom checkmate assert_subset
#' @importFrom BGLR BGLR
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats complete.cases cor
#' @param Phenotype A matrix of phenotypic data. Rownames must contain individual ID but rownames is not necessarily required. The missing value should be coded as NA.
#' @param Genotype A matrix of genotypic data. Both phased and unphased genotypic data are allowed. Rownames are individual ID, and colnames are marker ID (number is not allowed). In case of phased data, each cell should contain two numbers seperated by a vertical var (e.g. 1|0), and in case of unphased genotype, each cell should be either 0, 1, or 2.
#' @param phased A logical value indicating the genotypic data give is phased or unphased (TRUE corresponds to phased). The devault value is TRUE.
#' @param Model A character string specifying the model to apply. It should be either "Ridge", "Lasso", "BRR" (Bayesian Ridge Redgression), "BayesB", or "BayesC". The former two model is implemented by glmnet and the latter three models are implemented by BGLR.
#' @param Effect A character string specifying marker effect applied. It should be either "A" or "AD". "A" corresponds to the model with only additive effect, and "AD" corresponds to the model with both additive and dominance effect.
#' @param ... Additional parameters for cv.glmnet glmnet, and BGLR can be specified here.
#' @return A data frame of pedigree information, sorted in ancestral order.
#' \describe{
#' \item{model}{Model applied. Either "Ridge", "Lasso", "BRR", "BayesB", or "BayesC" }
#' \item{coefDetermin}{Coefficient of determination}
#' \item{rFit}{Pearson's correlation coefficient between phenotype and estimated genotypic effect}
#' \item{mseFit}{Mean square error}
#' \item{yEst}{Estimated genetic effect}
#' \item{mEffect}{A list of length two. One contains the additive marker effect, and another contains the dominance marker effect}
#' }
#' @export

calcMarkerEff <- function( Phenotype, Genotype, phased = TRUE, Model = "Ridge", Effect = "A",... ){

  assert_matrix( Phenotype )
  assert_matrix( Genotype )
  assert_character( rownames(Phenotype) )
  assert_character( rownames(Genotype) )
  assert_character( colnames(Genotype) )
  assert_subset( Model, c("Ridge", "Lasso", "BRR", "BayesB", "BayesC") )
  assert_subset( Effect, c("A", "AD"))

  dots <- list( ... )

  MarkerEstim <- vector( "list", 6 )
  names( MarkerEstim ) <- c( "model", "coefDetermin", "rFit", "mseFit", "yEst", "mEffect" )

  Phenotype <- Phenotype[complete.cases(Phenotype), , drop = FALSE]

  if( phased ){
    genome <- matrix( sapply( Genotype, function( x ) {
      if( is.na(x) ){
        return( x )
      } else {
        fst <- as.numeric( substr( x, 1, 1 ) )
        snd <- as.numeric( substr( x, 3, 3 ) )
        return( fst + snd )
      }
    }), ncol = ncol( Genotype ) )
    colnames( genome ) <- colnames( Genotype )
    rownames( genome ) <- rownames( Genotype )
  } else {
    genome <- Genotype
  }

  commonID <- intersect( rownames( genome ), rownames( Phenotype ) )
  if( length(commonID) == 0 ){
    stop( "There are no individuals with both phenotype and genotype\n Please comfirm that rownames of Phenotype and Genotype are correctly specified." )
  }
  genoSelectAdd <- genome[ match( commonID, rownames( genome ) ) ,]
  phenoSelect <- Phenotype[ match( commonID, rownames( Phenotype ) ), ]

  if( Effect == "AD" ){
    genoSelectDom <- genoSelectAdd %% 2
    genoSelectAll <- cbind( genoSelectAdd, genoSelectDom )
  } else {
    genoSelectAll <- genoSelectAdd
  }

  if( Model == "Ridge" ){
    param <- 0
  } else if ( Model == "Lasso" ){
    param <- 1
  } else if ( Model %in% c( "BRR", "BayesB", "BayesC" ) ){
    param <- Model
  }

  if ( Model %in% c("Ridge", "Lasso") ){
    args1 <- dots[ names(dots) %in% names(formals(cv.glmnet)) ]
    lasso.cv <- do.call( cv.glmnet,
                         c( list(x = genoSelectAll, y = phenoSelect, alpha = param), args1 ) )
    lambda <- lasso.cv$lambda.min
    args2 <- dots[ names(dots) %in% names(formals(glmnet)) ]
    result <- do.call( glmnet,
                       c( list(x = genoSelectAll, y = phenoSelect, lambda = lambda, alpha = param), args2 )  )
    b <- as.matrix( result$beta )
    mu <- result$a0
  } else if( Model %in% c( "BRR", "BayesB", "BayesC" ) ){
    args3 <- dots[ names(dots) %in% names(formals(BGLR)) ]
    if( is.null( args3[["saveAt"]] ) ){
      tmpdir <- tempdir()
      args3[["saveAt"]] <- paste0( tmpdir, "/" )
    }
    Eta <- list( list( X = genoSelectAll, model = param ) )
    result <- do.call( BGLR,
                       c( list( y = phenoSelect, ETA = Eta ), args3 ) )
    b <- as.matrix( result$ETA[[1]]$b )
    mu <- result$ETA[[1]]$S0
  }

  yEst <- as.vector( genoSelectAll %*% b ) + mu
  names( yEst ) <- rownames( genoSelectAdd )
  MarkerEstim$model <- Model
  MarkerEstim$coefDetermin <- 1 - sum(( phenoSelect - yEst)^2) / sum((phenoSelect - mean(phenoSelect))^2)
  MarkerEstim$rFit <- cor( yEst, phenoSelect )
  MarkerEstim$mseFit <- mean( ( yEst - phenoSelect )^2 )
  MarkerEstim$yEst <- yEst
  MarkerEstim$mEffect <- vector( "list", 2 )
  if( Effect == "AD" ){
    MarkerEstim$mEffect$Add <- b[ 1:(nrow(b)/2), , drop = FALSE ]
    MarkerEstim$mEffect$Dom <- b[ -( 1:(nrow(b)/2) ), , drop = FALSE ]
  } else {
    MarkerEstim$mEffect$Add <- b
    b[] <- 0
    MarkerEstim$mEffect$Dom <- b
  }
  rownames( MarkerEstim$mEffect$Add ) <- rownames( MarkerEstim$mEffect$Dom ) <- colnames( Genotype )

  return( MarkerEstim )

}
