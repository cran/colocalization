
nsinc.z <- function(data, membership, dim=2, r.min=NULL,r.max=NULL, r.count=NULL,
                             r.adjust = NULL, box=NULL, edge.effect=TRUE,
                             strata=FALSE, base.member = NULL,r.model="full",...){

  # data:
  # data must be a dataframe and contain at least 3 or 4 columns including "membership" and
  # coordinates "x", "y", and "z" matching with "dim";

  # membership:
  # specification of the columne containing membership of signals;

  # dim:
  # if dim = 2, x,y columns are dealth with;
  # if dim = 3, x, y z columns are dealth with;

  # r.min, r.max, r.adjust must be numerical;
  # r.min <= r.max if specified;
  # r.min + r.adjust <= r.max - r.adjust if specified;

  # r.count must be integers;
  # r.count >= 1;
  # if r.max - r.adjust  = r.max - r.adjust, then r.count = 1;

  # r.adjust must be a nonnegative numeric;
  # r.adjust is used to form the interval for the r series: [r.min + r.adjust, r.max - r.adjust]
  # to avoid zero standard deviation of average proportion densities;

  # box:
  # must be a data.frame showing xmin, xmax, ymin, ymax, zmin, zmax if specified
  # xmin < xmax, ymin < ymax, zmin < zmax;

  # strata:
  # = FALSE if colocalization is for bi-direction;
  # = TRUE if colocalization is for single directions and must have a base.member either specified
  #           or designated by default;

  # r.model:
  # can be any of "full","r.med","other"

  ##########################################################
  # Step 1. preprocess the data
  # rename column names

  data <- as.data.frame(data)
  if (membership %in% colnames(data)){
    colnames(data)[colnames(data) == membership] <- "membership"
  } else{
    stop(paste("There is no column names in the data called ",membership, "!",sep=""))
  }

  if (dim != 2 & dim !=3){
    stop("dim must be either 2 or 3!")
  }

  if (!(any(c("x","xc","X","Xc") %in% colnames(data)))) {
    stop("Data must contain a 'x' column!")
  }
  colnames(data)[colnames(data) == "x" | colnames(data) == "X"
                 | colnames(data) == "xc" | colnames(data) == "Xc"] <- "x"

  if (!(any(c("y","yc","Y","Yc") %in% colnames(data)))) {
    stop("Data must contain a 'y' column!")
  }

  colnames(data)[colnames(data) == "y" | colnames(data) == "Y"
                 | colnames(data) == "yc" | colnames(data) == "Yc"] <- "y"
  if (dim==3){
    if (!(any(c("z","zc","Z","Zc") %in% colnames(data)))) {
      stop("Data must contain a 'z' column!")
    }
    colnames(data)[colnames(data) == "z" | colnames(data) == "Z"
                   | colnames(data) == "zc"| colnames(data) == "Zc"] <- "z"
  }

  cat(paste("The method is 'nsinc.z'!\n"))
  cat(paste("The input data has dimension of ", dim, "!\n",sep=""))
  ##########################################################
  # Step 2: functions to calculate v.r for each signal for edge effect corrections later

  # subfunction for the area of each quarter in 2D contained inside the box
  calculate.v.r.2d.sub <- function(d,edge.effect){
    # d is a vector of length 0,1 or 2, with values between 0 and 1
    if (edge.effect) {
      quad.area <- ifelse(length(d) == 1, 0.5*( asin(1-d) + (1-d)*sqrt(1-(1-d)^2) ),
                          ifelse(length(d) == 2 &  sum((1-d)^2) > 1,  0.5*sum(asin(1-d) + (1-d)*sqrt(1-(1-d)^2)) - (pi/4),
                                 ifelse(length(d) == 2 &  sum((1-d)^2) <= 1,
                                        prod(1-d), pi/4 ) ) )
      return(quad.area) # with edge effect correction
    } else {
      return(pi/4) # with no edge effect correction
    }
  }

  # function for the area of each circle in 2D contained inside the box
  calculate.v.r.2d <-function(d.vector,edge.effect){
    # d.vector is a vector of length 4 with values either <=0 or 0< <=1

    dx <- d.vector[c(1,2)] # left, right distance out of the box
    dy <- d.vector[c(3,4)] # top, bottom distance out of the box

    quad <- expand.grid(dx,dy) # get all 4 quarters: 2,1,3,4 quadrants respectively

    colnames(quad) <- c("x","y")

    quad.list <- split(quad, seq(nrow(quad))) # get the x,y dimensions of each quarter

    quad.positive.list <- lapply(quad.list, function(x) {x[x> 0]})
    # get which dimension of each quarter is outside of the box and how far

    quad.area <- lapply(quad.positive.list,
                        function(x,edge.effect2 = edge.effect){
                          calculate.v.r.2d.sub(x,edge.effect = edge.effect2)
                        })
    circle.area <- sum(unlist(quad.area))
    return(circle.area) # this is the enclosed area with normalization by r^2
  }

  # subfunction for the volume of each octant in 3D contained inside the box
  calculate.v.r.3d.sub <- function(d,edge.effect){
    # d is a vector of length 0,1,2 or 3, with values between 0 and 1

    # length(d) == 3 &  length(which(colSums(combn((1-d)^2,2)) <= 1))==1
    intersec.one <- function(d){
      select <- which(colSums(combn((1-d)^2,2)) <= 1)
      combine <- combn(d,2) [,select]

      integ.combine <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-combine[2], sqrt(1-(x^2)))$value}),
        1-combine[1], sqrt(1-(1-combine[2])^2)) $value

      return(2-sum((d^2)*(3-d))+integ.combine)
    }
    # length(d) == 3 &  length(which(colSums(combn((1-d)^2,2)) <= 1))==2
    intersec.two <- function(d){
      select <- which(colSums(combn((1-d)^2,2)) <= 1)
      combine <- combn(d,2) [,select]

      integ.combine1 <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-combine[2,1], sqrt(1-(x^2)))$value}),
        1-combine[1,1], sqrt(1-(1-combine[2,1])^2)) $value
      integ.combine2 <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-combine[2,2], sqrt(1-(x^2)))$value}),
        1-combine[1,2], sqrt(1-(1-combine[2,2])^2)) $value

      return(2-sum((d^2)*(3-d)) + integ.combine1+integ.combine2 )
    }

    #length(d) == 3 & length(which(colSums(combn((1-d)^2,2)) > 1))==0  &  sum((1-d)^2) > 1
    intersec.three <- function(d){

      integ.combine1 <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-d[2], sqrt(1-(x^2)))$value}),
        1-d[1], sqrt(1-(1-d[2])^2)) $value

      integ.combine2 <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-d[3], sqrt(1-(x^2)))$value}),
        1-d[1], sqrt(1-(1-d[3])^2)) $value

      integ.combine3 <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-d[3], sqrt(1-(x^2)))$value}),
        1-d[2], sqrt(1-(1-d[3])^2)) $value

      return(2-sum((d^2)*(3-d)) + integ.combine1+integ.combine2 + integ.combine3)
    }

    #length(d) == 3 &  sum((1-d)^2) <= 1
    intersec.three.union <- function(d){

      integ.combine1 <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-d[2], sqrt(1-(x^2)))$value}),
        1-d[1], sqrt(1-(1-d[2])^2)) $value

      integ.combine2 <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-d[3], sqrt(1-(x^2)))$value}),
        1-d[1], sqrt(1-(1-d[3])^2)) $value

      integ.combine3 <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-d[3], sqrt(1-(x^2)))$value}),
        1-d[2], sqrt(1-(1-d[3])^2)) $value

      integ.union <-  (12/pi) * integrate(Vectorize(function(x) {
        integrate(function(y){sqrt(1-(x^2)-(y^2))-(1-d[3])}, 1-d[2], sqrt(1-(x^2)-(1-d[3])^2))$value}),
        1-d[1], sqrt(1-(1-d[2])^2-(1-d[3])^2)) $value

      return(2-sum((d^2)*(3-d)) + integ.combine1+integ.combine2 + integ.combine3 - integ.union)
    }

    if (edge.effect) {
      octant.volume <- ifelse(length(d) ==0, 2,
                              ifelse(length(d) == 1, 2 - (d^2)*(3-d),
                                     ifelse(length(d) == 2 &  sum((1-d)^2) > 1,  2 - sum((d^2)*(3-d)),
                                            ifelse(length(d) == 2 &  sum((1-d)^2) <= 1,
                                                   2 - sum((d^2)*(3-d)) + (12/pi) * integrate(Vectorize(function(x) {
                                                     integrate(function(y){sqrt(1-(x^2)-(y^2))}, 1-d[2], sqrt(1-(x^2)))$value}),
                                                     1-d[1], sqrt(1-(1-d[2])^2))$value,
                                                   ifelse(length(d) == 3 & all(colSums(combn((1-d)^2,2)) > 1), 2 -  sum((d^2)*(3-d)),
                                                          ifelse(length(d) == 3 &  length(which(colSums(combn((1-d)^2,2)) <= 1))==1,
                                                                 intersec.one(d),
                                                                 ifelse(length(d) == 3 &  length(which(colSums(combn((1-d)^2,2)) <= 1))==2,
                                                                        intersec.two(d),
                                                                        ifelse(length(d) == 3 & length(which(colSums(combn((1-d)^2,2)) > 1))==0  &  sum((1-d)^2) > 1,
                                                                               intersec.three(d),
                                                                               intersec.three.union(d)
                                                                        ) ) ) ) ) ) ) )
      return(octant.volume) # with edge effect correction
    } else {
      return(2) # with no edge effect correction
    }
  }

  # function for the volume of each ball in 3D contained inside the box
  calculate.v.r.3d <-function(d.vector, edge.effect){
    # d.vector is a vector of length 6 with values either <=0 or 0< <=1

    dx <- d.vector[c(1,2)]
    dy <- d.vector[c(3,4)]
    dz <- d.vector[c(5,6)]

    octant <- expand.grid(dx,dy,dz) # get all 8 octants

    colnames(octant) <-c("x","y","z")

    octant.list <- split(octant, seq(nrow(octant)))  # split into 8 octants

    octant.positive.list <- lapply(octant.list, function(x) {x[x > 0]})

    octant.volume <- lapply(octant.positive.list,
                            function(x,edge.effect2 = edge.effect){
                              calculate.v.r.3d.sub(x,edge.effect = edge.effect2)
                            })

    ball.volume <- sum(unlist(octant.volume))

    return(ball.volume)
  }

  # function to calculate v.r for each signal
  calculate.v.r <- function(data,box.post,r,dim,A,edge.effect){
    if (dim == 2){

      # find the 4 vertices of each signal
      data.vertices <- data.frame(r-data$x, data$x + r, r-data$y, data$y + r )
      d <- as.matrix(data.vertices) + as.matrix(box.post)
      d.positive <- (1/(r))* (d > 0)*d
      # normalized by r
      # negative d's are replaced with 0
      d.positive.list <- split(d.positive,seq(nrow(d.positive))) # split rows of d.positive to be a list
      v.r <- as.vector(unlist(lapply(d.positive.list,
                                     function(x, edge.effect1=edge.effect) {
                                       calculate.v.r.2d(x,edge.effect=edge.effect1)
                                     })))*(r^2)/A # add the scale r back
    }
    else if (dim == 3){

      # find the 6 vertices of each signal
      data.vertices <- data.frame(r-data$x, data$x + r, r-data$y, data$y + r , r-data$z, data$z + r )
      d <- as.matrix(data.vertices) + as.matrix(box.post)
      d.positive <- (1/(r))* (d > 0)*d # normalized by r
      # negative d's are replaced with 0
      d.positive.list <- split(d.positive,seq(nrow(d.positive)))
      v.r <- as.vector(unlist(lapply(d.positive.list,
                                     function(x, edge.effect1=edge.effect) {
                                       calculate.v.r.3d(x,edge.effect=edge.effect1)
                                     })))*(r^3*pi/12)/A # add the scale r back
    }
    return(v.r)
  }

  ##########################################################
  # Step 3. calculate distance.matrix, box.post,
  # and the area of the whole region: A
  # for edge effect corrections later

  if (dim == 2){

    if (is.null(box)){

      data.post <- data[,c("membership","x","y")]
      distance.matrix <- as.matrix(dist(data.post[,c("x","y")]))

      distance.matrix.x <- as.matrix(dist(data.post$x))
      distance.matrix.y <- as.matrix(dist(data.post$y))
      diag(distance.matrix.x) <- NA
      diag(distance.matrix.y) <- NA

      distance.x.min <- median(apply(distance.matrix.x,2,function(x) min(x, na.rm=TRUE)))
      distance.y.min <- median(apply(distance.matrix.y,2,function(x) min(x, na.rm=TRUE)))
      #Minimum and maximum of every column:  apply(a,2,min)    apply(a,2,max)

      box <- data.frame(xmin=min(data.post$x)-distance.x.min,
                        xmax=max(data.post$x)+distance.x.min,
                        ymin=min(data.post$y)-distance.y.min,
                        ymax=max(data.post$y)+distance.y.min)
      study.region <- box
      study.region$buffer.x <- distance.x.min
      study.region$buffer.y <- distance.y.min

      cat(paste("The input study region is 'NULL'!\n"))
      cat(paste("The observed study region is [", # without buffer edges added
                min(data.post$x),",",max(data.post$x),"] X [",min(data.post$y),",",
                max(data.post$y),"]!\n", sep=""))
      #cat(paste("The buffer edges are added to x and y dimensions with width ",
      #          distance.x.min," and ",distance.y.min,", respectively!\n",sep=""))
    } else {
      if (any(!(c("xmin","xmax","ymin","ymax") %in% colnames(box)))) {
        stop("'box' must be a dataframe containing columns 'xmin','xmax','ymin' and 'ymax'!")
      }
      if (box$xmax <= box$xmin | box$ymax <= box$ymin){
        stop("'xmax' or 'ymax' must be larger than 'xmin' or 'ymin' in 'box'!")
      }


      data.post <- data[,c("membership","x","y")]
      data.post <- data.post[data.post$x >= box$xmin & data.post$x <= box$xmax &
                             data.post$y >= box$ymin & data.post$y <= box$ymax,]

      cat(paste("The observed study region is [",
                min(data$x),",",max(data$x),"] X [",min(data$y),",",
                max(data$y),"]!\n", sep=""))
      cat(paste("The input study region is [",box$xmin,",",box$xmax, "] X [",box$ymin,",",
                box$ymax,"]!\n",sep=""))

      distance.matrix <- as.matrix(dist(data.post[,c("x","y")]))
      study.region <- box
      study.region$buffer.x <- 0
      study.region$buffer.y <- 0
    }

    box.post <- box
    box.post$xmax<- -box$xmax
    box.post$ymax <- -box$ymax
    box.post <- matrix(rep(as.numeric(box.post),each=nrow(data.post)),ncol=4,byrow=FALSE)
    A <- (box$xmax-box$xmin)*(box$ymax-box$ymin)

  } else if (dim==3){

    if (is.null(box)){
      data.post <- data[,c("membership","x","y","z")]
      distance.matrix <- as.matrix(dist(data.post[,c("x","y","z")]))

      distance.matrix.x <- as.matrix(dist(data.post$x))
      distance.matrix.y <- as.matrix(dist(data.post$y))
      distance.matrix.z <- as.matrix(dist(data.post$z))
      diag(distance.matrix.x) <- NA
      diag(distance.matrix.y) <- NA
      diag(distance.matrix.z) <- NA

      distance.x.min <- median(apply(distance.matrix.x,2,function(x) min(x, na.rm=TRUE)))
      distance.y.min <- median(apply(distance.matrix.y,2,function(x) min(x, na.rm=TRUE)))
      distance.z.min <- median(apply(distance.matrix.z,2,function(x) min(x, na.rm=TRUE)))

      #Minimum and maximum of every column:  apply(a,2,min)    apply(a,2,max)

      box <- data.frame(xmin=min(data.post$x)-distance.x.min,
                        xmax=max(data.post$x)+distance.x.min,
                        ymin=min(data.post$y)-distance.y.min,
                        ymax=max(data.post$y)+distance.y.min,
                        zmin=min(data.post$z)-distance.z.min,
                        zmax=max(data.post$z)+distance.z.min)
      study.region <- box
      study.region$buffer.x <- distance.x.min
      study.region$buffer.y <- distance.y.min
      study.region$buffer.z <- distance.z.min

      cat(paste("The input study region is 'NULL'!\n"))
      cat(paste("The observed study region is [", # without buffer edges added
                min(data.post$x),",",max(data.post$x),"] X [",min(data.post$y),",",
                max(data.post$y),"] X [",min(data.post$z),",",max(data.post$z),"]!\n", sep=""))
      # cat(paste("The buffer edges are added to x, y and z dimensions with width ",
      #          distance.x.min,",", distance.y.min," and ", distance.z.min,", respectively!\n",sep=""))

    }else{
      if (any(!(c("xmin","xmax","ymin","ymax","zmin","zmax") %in% colnames(box)))) {
        stop("'box' must be a dataframe containing columns'xmin','xmax','ymin', 'ymax','zmin',
             and 'zmax'!")
      }
      if (box$xmax<=box$xmin | box$ymax<=box$ymin | box$zmax<=box$zmin){
        stop("'xmax','ymax' or 'zmax' must be larger than 'xmin', 'ymin' or 'zmin' in 'box'!")
      }
      data.post <- data[,c("membership","x","y","z")]
      data.post <- data.post[data.post$x >= box$xmin & data.post$x <= box$xmax &
                             data.post$y >= box$ymin & data.post$y <= box$ymax &
                             data.post$z >= box$zmin & data.post$z <= box$zmax,]

      cat(paste("The observed study region is [",
                min(data$x),",",max(data$x),"] X [",min(data$y),",",
                max(data$y),"] X [",min(data$z),",",max(data$z),"]!\n", sep=""))
      cat(paste("The input study region is [",box$xmin,",",box$xmax, "] X [",box$ymin,",",
                box$ymax,"] X [",box$zmin,",",box$zmax,"]!\n",sep=""))

      distance.matrix <- as.matrix(dist(data.post[,c("x","y","z")]))
      study.region <- box
      study.region$buffer.x <- 0
      study.region$buffer.y <- 0
      study.region$buffer.z <- 0
      }

    box.post <- box
    box.post$xmax<- -box$xmax
    box.post$ymax <- -box$ymax
    box.post$zmax <- -box$zmax
    box.post <- matrix(rep(as.numeric(box.post),each=nrow(data.post)),ncol=6,byrow=FALSE)
    A <- (box$xmax-box$xmin)*(box$ymax-box$ymin)*(box$zmax-box$zmin)
}

  ##########################################################
  # Step 4. find K, n, and membership names

  # summarize membership information of the input data
  data$membership <- as.factor(data$membership)
  droplevels(data$membership)
  data$membership <- factor(data$membership)

  membership.names.pre <- unique(data$membership)
  K.pre <- length(membership.names.pre)
  if (K.pre < 2){
    stop("There must be at least two memberships of signals in the input data!")
  }
  n.pre <- nrow(data)

  # summarize membership information of the data after removing signals out of the specified study region
  data.post$membership <- as.factor(data.post$membership)
  droplevels(data.post$membership)
  data.post$membership <- factor(data.post$membership)

  membership.names <- unique(data.post$membership)
  K <- length(membership.names)
  if (K < 2){
    stop("There must be at least two memberships of signals enclosed in the study region!")
  }
  n <- nrow(data.post)

  ##########################################################
  # Step 5. create membership matrix of dimension: n by K
  # M, M.norm
  # M.norm is M normalized by number of signals in each membership: n.k

  # membership matrix for the input data
  M.pre <- matrix(0,nrow=n.pre,ncol=K.pre)
  for (k in 1:K.pre){
    M.pre[,k] <- data$membership == membership.names.pre[k]
  }
  colnames(M.pre) <- membership.names.pre

  # membership matrix for the data after the removal of signals outside of study region
  M <- matrix(0,nrow=n,ncol=K)
  for (k in 1:K){
    M[,k] <- data.post$membership == membership.names[k]
  }
  colnames(M) <- membership.names
  M.norm <- M %*% (diag(1/colSums(M)))
  colnames(M.norm) <- membership.names

  # find the number of signals in each channel: n.1,n.2,...n.k
  # and their reciprocals: n.k.recipr, n.k.recipr.squar

  n.k.pre <- colSums(M.pre)

  n.k <- colSums(M)
  n.k.recipr <- vector()
  for (k in 1:K){
    n.k.recipr[k] <-  mean(M.norm[M.norm[,k]>0,k])
  }
  n.k.recipr.squar <- n.k.recipr * n.k.recipr


  ##########################################################
  # Step 6. set up the r sequence

  r.min.post.data <-  min(distance.matrix[lower.tri(distance.matrix)])
  r.max.post.data <- max(distance.matrix[lower.tri(distance.matrix)])
  r.med.post.data <- median(distance.matrix[lower.tri(distance.matrix)])
  r.max.post.data.half <- 0.5*r.max.post.data

  cat(paste("\n The default full r range for colocalization index of type z is [",r.min.post.data,",",
            r.max.post.data.half,"]!\n\n",sep=""))

  if (max(r.min.post.data,r.med.post.data) >= r.max.post.data.half){
    stop(paste("The smallest or median of interpoint distances is no less than half of the largest interpoint distance! ",
               "Hence, the z-type colocalization method is not appropriate for calculating the colocalization degree! ",
               "It is suggested to use the d-type colocalization index function: 'colocalization.d.r'!",sep= ""))
  } else {
    if (r.model == "full"){
      cat("The r.model = 'full'!\n")
      r.min <- r.min.post.data
      cat(paste("r.min = ",r.min," is used!\n",sep=""))

      r.max <- r.max.post.data.half
      cat(paste("r.max = ",r.max," is used!\n",sep=""))

      r.count.default <- 30

      if (is.null(r.count)){
        r.count <- r.count.default
        cat(paste("r.count = ",r.count," is used!\n",sep=""))
      } else if (round(r.count,digits=0) <= 0){
        r.count <- r.count.default
        warning(paste("r.count must be a positive integer and r.count = ", r.count," is used instead!",sep=""))
        cat(paste("r.count must be a positive integer and r.count = ",r.count," is used for r.count!\n",sep=""))
      } else {
        r.count <- round(r.count,digits=0)
        cat(paste("r.count = ",r.count," is used!\n",sep=""))
      }

      r.adjust.default <- (r.max-r.min)/(r.count + 1)
      r.adjust.middle <- (r.max-r.min)/2

      if (is.null(r.adjust)){
        r.adjust <- r.adjust.default
        cat(paste("r.adjust = ",r.adjust," is used!\n",sep=""))
      } else if (r.adjust > r.adjust.middle | r.adjust <= 0){
        warning("The r.adjust must be a positive number no larger than half of the difference between r.max and r.min!
                We suggest to use the default value by assigning r.adjust = NULL!")
        cat(paste("The r.adjust must be a positive number no larger than half of the difference between r.max and r.min! \n",
                  "The input r.adjust = ",r.adjust," is invalid!\n ",
                  "r.adjust = ",r.adjust.default, " is used instead!\n",sep=""))
        r.adjust <- r.adjust.default
      } else {
        cat(paste("r.adjust = ",r.adjust," is used!\n",sep=""))
      }

    } else if (r.model == "r.med"){
      cat("The r.model = 'r.med'!\n")

      r.min = r.med.post.data
      cat(paste("r.min = ",r.min," is used!\n",sep=""))

      r.max = r.med.post.data
      cat(paste("r.max = ",r.max," is used!\n",sep=""))

      r.count <- 1
      r.adjust <- 0
      cat(paste("Because of r.min = r.max, r.count = 1 and r.adjust = 0 are used!\n"))

    } else if (r.model == "other"){
      cat("The r.model = 'other'!\n")

      if (is.null(r.min)){
        stop("If choose the 'other' for r.model, then r.min must be specified by the user!")
      } else if (is.null(r.max)){
        stop("If choose the 'other' for r.model, then r.max must be specified by the user!")
      } else if(r.min < r.min.post.data | r.min > r.max.post.data.half) {
        stop(paste("r.min must be between the smallest and half of the largest interpoint distances: [",
                   r.min.post.data,",",r.max.post.data.half,"]!",sep=""))
      } else if (r.max < r.min.post.data | r.max > r.max.post.data.half) {
        stop(paste("r.max must be between the smallest and half of the largest interpoint distances: [",
                   r.min.post.data,",",r.max.post.data.half,"]!",sep=""))
      } else if (r.min > r.max){
        stop("The r.min must be smaller than r.max!")
      } else if (r.min == r.max){
        r.count <- 1
        r.adjust <- 0
        cat(paste("Because of r.min = r.max, r.count = 1 and r.adjust = 0 are used!\n"))
      } else { # r.min* <= r.min < r.max <= 0.5*r.max*
        cat(paste("r.min = ",r.min," is used!\n",sep=""))
        cat(paste("r.max = ",r.max," is used!\n",sep=""))

        r.count.default <- 30

        if (is.null(r.count)){
          r.count <- r.count.default
          cat(paste("r.count = ",r.count," is used!\n",sep=""))
        } else if (round(r.count,digits=0) <= 0){
          r.count <- r.count.default
          warning(paste("r.count must be a positive integer and r.count = ", r.count," is used instead!",sep=""))
          cat(paste("r.count must be a positive integer and r.count = ",r.count," is used for r.count!\n",sep=""))
        } else {
          r.count <- round(r.count,digits=0)
          cat(paste("r.count = ",r.count," is used!\n",sep=""))
        }

        r.adjust.default <- (r.max-r.min)/(r.count + 1)
        r.adjust.middle <- (r.max-r.min)/2

        if (is.null(r.adjust)){
          if (r.min == r.min.post.data | r.max == r.max.post.data) {
            r.adjust <- r.adjust.default
            cat(paste("Since [r.min,r.max] = [",r.min,",",r.max,
                      "] is on the boundary of the default full r range: [",
                      r.min.post.data,",",r.max.post.data,"], r.adjust = ",
                      r.adjust," is used!\n",sep=""))
          } else {
            r.adjust <- 0
            cat(paste("Since [r.min,r.max] = [",r.min,",",r.max,"] is within the default full r range: [",
                      r.min.post.data,",",r.max.post.data,"], r.adjust = ",r.adjust,
                      " is used!\n",sep=""))
          }
        } else if (r.adjust > r.adjust.middle | r.adjust < 0){
          stop (paste("The r.adjust must be a nonnegative number smaller than half of the difference between r.max and r.min!\n",
                      "It is suggested to try r.adjust = 0 or r.adjusts = NULL!"))
        } else {
          if ((r.min == r.min.post.data | r.max == r.max.post.data) & r.adjust==0) {
            r.adjust <- r.adjust.default
            cat(paste("Since [r.min,r.max] = [",r.min,",",r.max,
                      "] is on the boundary of the default full r range: [",
                      r.min.post.data,",",r.max.post.data,"], r.adjust = ",r.adjust,
                      " is used!\n",sep=""))
          }else{
            cat(paste("r.adjust = ",r.adjust," is used!\n",sep=""))
          }
        }
      }
    } else {
      stop("r.model must be one of 'full' (default),'r.med','other'!")
    }
  }

  r.seq <- seq(r.min + r.adjust, r.max - r.adjust,length.out = r.count)


  ##########################################################
  # Step 7. finally calculate the index at each r

  # set weights for coefficient coefficients to get average index at each r
  # weight <- (matrix(rep(n.k,each=length(n.k)),ncol=length(n.k),byrow=FALSE) +
  #             matrix(rep(n.k,each=length(n.k)),ncol=length(n.k),byrow=TRUE) ) / (n*(K-1))

  if (strata==TRUE){
    if (is.null(base.member)){
      base.member <- membership.names[1]
      cat(paste("Single-direction colocalization is considered and '",base.member,"' is chosen as the base membership!\n", sep=""))
    } else if (base.member %in% membership.names) {
      cat(paste("Single-direction colocalization is considered and '",base.member,"' is chosen as the base membership!\n"))
    } else{
      stop(paste0("The specified base membership '",base.member, "' is not found in the provided data!"))
    }
  } else{
    cat("Bi-direction colocalization is considered!\n")
    base.member <- NULL
  }

  index <- data.frame()
  P.orig <- data.frame()
  P.z <- data.frame()
  for (r in r.seq){
    # cat(paste("Calculating index for r = ",r,"!", sep=""))
    # calculate v.r = ratio of area of the local circle over the area of the whole box
    v.r <- calculate.v.r(data=data.post,box.post=box.post,r=r,dim=dim,A=A,edge.effect=edge.effect)
    v.r <- ifelse(v.r<=1,v.r,1)
    # make duplicates of v.r
    v.r.post <- matrix(rep(v.r,each=K), ncol=K, byrow=TRUE)
    mu <-  matrix(1,ncol=K,nrow=n)  * v.r.post + (M.norm %*% (diag(K)) ) * (1-v.r.post)
    tau.squared <- (matrix(rep(n.k.recipr,each=n), ncol=K, byrow=FALSE) - (M.norm^2)) * v.r.post * (1-v.r.post)
    # find the number of points enclosed in the r neighborhoold of every data point
    D <- (distance.matrix <= r) * 1
    # normalize the number of points to get the proportion in the r neighborhood of each data point
    # in order to find the index matrix, whose rows are points and columns are P.1,...P.K
    P <- D %*% M.norm
    # normalize the proportion by the proportion of area of r neighborhood over the whole region
    #P.v.r <- P / v.r.post

    P.normal <- (P-mu) / sqrt(tau.squared)

    # find the Pearson correlation for each pair of memberships
    if (strata==TRUE){
      P.normal.base <- P.normal[data$membership==base.member,]
      if ( 0 %in% apply(P.normal.base, 2, sd)){
        stop("The standard deviation of P.v.r.base is zero at r = ", r,"! Choose appropriate proximity size r!")
      }else{
        correlation <- cor(P.normal.base,...)
      }
    } else{
      if ( 0 %in% apply(P.normal, 2, sd)){
        stop("The standard deviation of P.v.r.base is zero at r = ", r,"! Choose appropriate proximity size r!")
      }else{
        correlation <- cor(P.normal,...)
      }
    }

    # find the mean of all pairs and store the corresponding r
    #index.r <- data.frame(index.z = sum(weight[lower.tri(weight)]*correlation[lower.tri(correlation)]),r=r)
    index.r <- data.frame(index.z = mean(correlation[lower.tri(correlation)]),r=r)
    cat(paste("index.z for r = ",r, " is ", index.r$index.z," ;\n",sep=""))
    index <- rbind(index,index.r)
    P.orig<- rbind(P.orig,data.frame(P,r=r,type="P"))
    if (strata==TRUE){
      P.z <- rbind(P.z, data.frame(P.normal.base,r=r,type="P.z"))
    } else{
      P.z <- rbind(P.z, data.frame(P.normal,r=r,type="P.z"))
    }
  }

  index.ave <- mean(index$index.z)
  cat(paste("The average coolocalization index of type z over all choosen r's is = ",index.ave,"!\n",sep=""))

  results <- list()
  results[["method"]] <- "nsinc.z"
  results[["input.data.summary"]] <- list("membership.levels" = K.pre,
                                          "membership" = data.frame("membership" = membership.names.pre,
                                                                    "signal.count" = n.k.pre,row.names = NULL),
                                          "dim" = dim)
  results[["post.data.summary"]] <- list("membership.levels" = K,
                                         "membership" = data.frame("membership" = membership.names,
                                                                   "signal.count" = n.k,row.names = NULL),
                                         "dim" = dim)
  results[["r.summary"]] <- data.frame("r.min" = r.min,"r.max" = r.max,
                                       "r.count"= r.count,"r.adjust"=r.adjust,
                                       "r.model" = r.model,
                                       "r.min.post.data" =  min(distance.matrix[distance.matrix!=0]),
                                       "r.max.post.data.half" = 0.5* max(distance.matrix[distance.matrix!=0]),
                                       "r.med.post.data" =  median(distance.matrix[distance.matrix!=0]))
  results[["strata"]] <- list("strata" = strata, "base.member" = base.member)
  results[["edge.effect"]] <- data.frame("edge.effect" = edge.effect)
  results[["index.all"]] <- index
  results[["index"]] <- index.ave # store the average index
  results[["post.data"]] <- data.post # return the data
  results[["study.region"]] <- study.region # store the study region
  results[["P.all"]] <- rbind(P.orig,P.z) # return all P values involved at different r
  results[["r"]] <- r.seq # return all r's
  attr(results, "class") <- "colocal"
  return(results)
}
