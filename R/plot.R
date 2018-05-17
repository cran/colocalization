#####################
## P-plot function ##
#####################
plot.colocal <- function(x, ...){
  UseMethod("plot")
}

plot.colocal <- function(x = NULL, ...){

  # get the membership for all observations in the data
  data.color.rgb <- x[["post.data"]]$membership

  if (length(unique(data.color.rgb)) > 2){
    stop(paste("The P plots can only be shown for two members!"))
  }

  # read r in from returned results
  allP <- x[["P.all"]]
  ra <- unique(allP$r)

  # get the names of the membership from the data
  membership.names <- x[["post.data.summary"]]$membership$membership
  name1<-as.character(as.vector(membership.names)[1])
  name2<-as.character(as.vector(membership.names)[2])

  # assign colors to the membership
  # if the memberships are red or green colors, then we use the same color as the showing color for each observation in the P plot
  # otherwise, we use red and green for the first and second detected memberships

  data.color.rgb <- ifelse(data.color.rgb=="red","red",
                           ifelse(data.color.rgb=="green","green",
                                  ifelse(data.color.rgb==name1,"red","green")))

  color1 <- ifelse(name1=="red","red",
                   ifelse(name1=="green","green","red"))

  color2 <- ifelse(name2=="red","red",
                   ifelse(name2=="green","green","green"))

  color.rgb <- unique(data.color.rgb)

  pplot <- list()

  for (i in 1:length(ra))
    local({
      i <- i

      # read P.all at r from returned result
      P.all <- allP[allP$r==ra[i],]

      if (x[["method"]]=="nsinc.d"){
        levels(P.all$type) <- c(expression(paste(italic(P),sep="")), expression(paste(italic(P^{d}),sep="")))

        P.all$data.color.rgb <- data.color.rgb

        data.refline <- data.frame(vline=c(NA,1),hline=c(NA,1),type=levels(P.all$type))

        p.x.label <- bquote({italic(P)^{italic(d)}}[.(paste(name1)),"=",.(paste(color1))][paste(",r=",.(paste(round(ra[i], digits = 4))),sep="")])
        p.y.label <- bquote({italic(P)^{italic(d)}}[.(paste(name2)),"=",.(paste(color2))][paste(",r=",.(paste(round(ra[i], digits = 4))),sep="")])
      } else if (x[["method"]]=="nsinc.z"){
        levels(P.all$type) <- c(expression(paste(italic(P),sep="")), expression(paste(italic(P^{z}),sep="")))

        P.all$data.color.rgb <- data.color.rgb

        p.x.label <- bquote({italic(P)^{italic(z)}}[.(paste(name1)),"=",.(paste(color1))][paste(",r=",.(paste(round(ra[i], digits = 4))),sep="")])
        p.y.label <- bquote({italic(P)^{italic(z)}}[.(paste(name2)),"=",.(paste(color2))][paste(",r=",.(paste(round(ra[i], digits = 4))),sep="")])
      }

      # plot P and P^d or P^z together
      P.all.plot<- ggplot(P.all, aes(x=P.all[,1],y=P.all[,2]))+
        geom_point(size=1,aes(color=data.color.rgb),show.legend = FALSE,alpha = 0.5)+
        geom_smooth(aes(group=P.all$type),method="lm")+
        scale_color_manual(name='',values = color.rgb)+
        facet_wrap(~P.all$type,scales="free",labeller=label_parsed) +
        labs(title="",x=p.x.label,y=p.y.label)

      pplot[[i]] <<- P.all.plot
    })

  return(pplot)
}
