require(stringr)
require(dplyr)
require(data.table)
require(ggplot2)

make_plot_table <- function(x,y,z) {
  data <- data.table(melt(x))
  setDT(data)[, val.mean := mean(value), by = variable]
  setDT(data)[, val.length := length(value), by = variable]
  setDT(data)[, val.sd := sd(value), by = variable]
  data$se <- data$val.sd/sqrt(data$val.length)
  data$se.min <- data$val.mean - data$se
  data$se.max <- data$val.mean + data$se
  data$color <- as.character(y)
  data$grp <- z
  data <- data.frame(data)
  return(data)
}

PatternFlow <- function(plist, pfg, nFactor, nIter, Repeats, ColsNames, PatternOccurence, ProteinOccurence, OutputDir) {

  #plist - Pattern List; pfg - Protein List
  
  ftab <- data.frame(first=numeric(), second=numeric(), xpat=numeric(), ypat=numeric(), distpat=numeric())
  dd <- expand.grid(1:Repeats, 1:nFactor, 1:Repeats, 1:nFactor)
  colnames(dd)<- c("first", "xpat", "second", "ypat")
  distpat <- c()
  for ( j in 1:nrow(dd))
  {
    m <- dist(rbind(plist[[dd$first[j]]][dd$xpat[j],],plist[[dd$second[j]]][dd$ypat[j],]))
    distpat <- c(distpat, m)
  }
  ftab <- data.frame(cbind(dd, distpat))
  colnames(ftab)<- c("first","xpat","second","ypat","distpat")
  ftabmod <- ftab[ftab$distpat!=0,]

  fsel <- data.frame(ftabmod %>%
                       group_by(first,second,xpat) %>%
                       dplyr::slice(which.min(distpat)))

  png(paste0(OutputDir,"/disthist.png",sep=""))
  print(plot(density(fsel$distpat)))
  dev.off()

  threshold <- sqrt(sum(rep(0.01^2,6)))

  fsel_sign <- ftabmod[ftabmod$distpat<threshold,]
  fsel_round <- fsel_sign

  pattab <- list()
  patnames <- list()
  vnames <- c()
  track <- list()
  for (i in 1:Repeats)
  {
    for (j in 1:nFactor)
    {
      list_name <- paste(i, "_", j,sep="")
      fsub <- fsel_round[fsel_round$first==i & fsel_round$xpat==j,]
      if (nrow(fsub)>0)
      {
        vtab <- data.frame()
        vnames <- pfg[[i]][[j]]
        fsub_new <- fsub
        track_ids <- c()
        for (k in 1:nrow(fsub))
        {
          vtab <- rbind(vtab,plist[[fsub$second[k]]][fsub$ypat[k],])
          #here I merge all patterns bound to the one determined and so define all connected patterns
          fsub_new <- rbind(fsub_new, fsel_round[fsel_round$first==fsub$second[k] & fsel_round$xpat==fsub$ypat[k],])
          vnames <- c(vnames, pfg[[fsub$second[k]]][[fsub$ypat[k]]])
          track_name <- paste(fsub$second[k], "_", fsub$ypat[k], sep="")
          track_ids <- c(track_ids, track_name)
        }
        vtab <- rbind(plist[[i]][j],vtab)
        colnames(vtab) <- ColsNames
        track_ids <- c(list_name, track_ids)
        track[[list_name]] <- track_ids
        if (nrow(vtab)/Repeats>PatternOccurence)
        {
          pattab[[list_name]] <- vtab
          nstat <- table(vnames)/nrow(vtab)
          patnames[[list_name]] <- nstat[order(nstat,decreasing = T)]
        }
        fsel_round <- fsetdiff(data.table(fsel_round),data.table(fsub_new), all=T)
      }
    }
  }

  patnames_sel <-lapply(patnames,function(x) x[x>=ProteinOccurence])
  cond <- lapply(patnames_sel, function(x) length(x) > 0)
  patnames_sign <- patnames_sel[unlist(cond)]

  pattab_sign <- pattab[names(patnames_sign)]

  cols <- str_order(rainbow(length(pattab_sign)))
  colpat <- cbind(names(pattab_sign), cols)
  print("Patterns colors:")
  print(colpat)
  plot <- ggplot()+ xlab('')+ylab('relative strength') +  theme_bw() +
    theme(legend.position="none", axis.text.y = element_text(size=15),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size=15),
          axis.title.y = element_text(size=20, vjust=5, margin = ggplot2::margin(l=20, unit = "pt")))

  pat <- list()

  patdistribution <- list()
  print(head(pattab_sign))
  if (length(pattab_sign)>0) {
  for (i in 1:length(pattab_sign))
  {
    pat[[i]] <- make_plot_table(pattab_sign[[i]],colpat[,1][i],i)

    patdistribution[[colpat[,1][i]]] <- unique(pat[[i]][,c("variable", "val.mean")])

    plot <- plot + geom_ribbon(data=pat[[i]], aes(x=variable, ymin = se.min, ymax = se.max, group=grp, fill=color), alpha = .25) +
      geom_line(data=pat[[i]], aes(x=variable, y=val.mean, group=grp), colour="grey20", alpha=.5)

  }
  plot <- plot

  tiff(paste0(OutputDir,"/Plot_", PatternOccurence, ".tiff", sep=""), width = 300, height=400, units="px")
  par(mar=c(6,6,2,2))
  print(plot)
  dev.off()
  }
  else {print("No patterns found")}
  
 save(patnames_sign,file=paste0(OutputDir,"/Patterns_", PatternOccurence, ".Rdata",sep=""))
 save(patdistribution,file=paste0(OutputDir,"/Patterns_", PatternOccurence, "_Annotation.Rdata",sep=""))
}
