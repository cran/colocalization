summary.colocal <- function(object, ...){
  UseMethod("summary")
}

summary.colocal <- function(object = NULL, ...){
  # results is the results returned from nsinc.d() or nsinc.z(),
  # which is a list of length 11 containing information of the input data,
  # calculation of the index,
  # and other useful information for user's reference.
  # object <- results

  if(is.null(object)){
    stop(paste("There are no results returned from nsinc.d or nsinc.z!"))
  }

  method <- object$method
  K <- object[["post.data.summary"]]["membership.levels"]
  dim <- object[["post.data.summary"]]["dim"]
  membership <- as.data.frame(object[["post.data.summary"]]$membership)
  r.summary <- object[["r.summary"]]

  study.region <- object[["study.region"]]

  r <- object[["r"]]
  r <- data.frame(t(r))
  colnames(r) <- 1:ncol(r)
  rownames(r) <- "r:"

  edge.effect <- object[["edge.effect"]]
  edge.effect.message <- ifelse(edge.effect,""," 'NO'")
  strata <- object[["strata"]]

  strata.choice <- ifelse(strata$strata, "single","double")
  base.member <- strata$base.member


  index.all <- t(object[["index.all"]])
  index.all[c(1,2),] <-index.all[c(2,1),]
  rownames(index.all) <- c(rownames(index.all)[2], rownames(index.all)[1])
  index <- object$index

  message0 <- ifelse(dim==2,paste("[",study.region$xmin,",",
                                  study.region$xmax,"] X [",study.region$ymin,",",
                                  study.region$ymax,"]", sep=""),
                     paste("[",study.region$xmin,",",
                           study.region$xmax,"] X [",study.region$ymin,",",
                           study.region$ymax,"] X [",study.region$zmin,",",
                           study.region$zmax,"]",sep=""))

  message1 <- ifelse(strata.choice == "single",
                     paste("Method '",method,"' for '", strata.choice, "'-direction colocalization
                           to base channel - '", base.member, "' with",edge.effect.message,
                           " edge effect correction is used for a '",dim,
                           "D' dataset over the study region '",message0, "' with '",K,"' channels of signals: \n\n",
                           sep=""),
                     paste("Method '",method,"' for '",  strata.choice, "'-direction colocalization with",
                           edge.effect.message," edge effect correction is used for a '",dim,
                           "D' dataset over the study region '",message0, "' with '",K,"' channels of signals: \n\n",
                           sep=""))


  message2 <- paste("Setting of proximity sizes:\n r.model = '",
                    r.summary$r.model, "', r.min = ",r.summary$r.min, ", r.max = ",r.summary$r.max,
                    ", r.count= ", r.summary$r.count,", r.adjust=",r.summary$r.adjust,
                    "\n\n",
                    sep="")
  message3 <- paste("The r series is:\n")

  message4 <- paste("The separate index at each proximity size: \n",paste="")

  message5 <- paste("The average index across the whole proximity range for the whole image = ",index,"!\n\n",sep="" )


  cat(message1)
  print(membership,row.names = FALSE)
  cat("\n")
  cat(message2)
  #cat(message3)
  #print(r)
  cat(message4)
  print(index.all)
  cat("\n\n")
  cat(message5)
}
