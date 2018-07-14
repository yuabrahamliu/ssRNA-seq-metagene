###########################
#Config
getwd()
configs <- readLines('config')
idx <- grep('drawfig', configs)
drawfig <- configs[idx]
pointer <- unlist(strsplit(drawfig, split = ' = '))[2]

if(toupper(pointer) != 'YES'){
  q(save = 'no')
}else{
  idx <- grep('fig_tag', configs)
  fig_tag <- configs[idx]
  if(length(fig_tag) == 0){
    idx <- grep('result_tag', configs)
    fig_tag <- configs[idx]
  }
  idx <- grep('project', configs)
  project <- configs[idx]
  idx <- grep('sampletags', configs)
  sampletags <- configs[idx]
  idx <- grep('halfwindow', configs)
  halfwindows <- configs[idx]
  idx <- grep('qcutoff', configs)
  qcutoff <- configs[idx]
  idx <- grep('allcut', configs)
  allcut <- configs[idx]
}

fig_tag <- unlist(strsplit(fig_tag, split = ' = '))[2]
project <- unlist(strsplit(project, split = ' = '))[2]
sampletags <- unlist(strsplit(sampletags, split = ' = '))[2]
halfwindows <- unlist(strsplit(halfwindows, split = ' = '))[2]
qcutoff <- as.numeric(unlist(strsplit(qcutoff, split = ' = '))[2])
allcut <- as.logical(toupper(unlist(strsplit(allcut, split = ' = '))[2]))


result_dir <- paste0(getwd(), '/', project, '/result')
setwd(result_dir)

tssfile <- paste0('sense', fig_tag, 
                  '_tss.htseq.txt')
ttsfile <- paste0('sense', fig_tag, 
                  '_tts.htseq.txt')

tssfile_anti <- paste0('antisense', fig_tag, 
                       '_tss.htseq.txt')

tsspdf_ori <- paste0(fig_tag, '_meta_tss.htseq.pdf')
ttspdf_ori <- paste0(fig_tag, '_meta_tts.htseq.pdf')


sampletags <- unlist(strsplit(sampletags, split = ','))
halfwindows <- as.numeric(unlist(strsplit(halfwindows, split = ',')))


##############################

generateplot_ttss <- function(filenames = tssfile, 
                              position = 'around_tss', 
                              halfwidth = 500, 
                              qcutoff = 0.99, allcut = FALSE){
  kind <- unlist(strsplit(filenames, split = '_', fixed = TRUE))
  kind <- kind[length(kind)-1]
  
  datafile <- read.table(filenames, sep = '\t', header = TRUE, 
                         stringsAsFactors = FALSE, check.names = FALSE, 
                         quote = '')
  
  genelen <- nrow(datafile)

  rowstart <- (nrow(datafile)/2 - halfwidth) + 1
  rowend <- nrow(datafile)/2 + halfwidth
  datafile <- datafile[rowstart:rowend,]
  
  for(j in seq(1, ncol(datafile), 2)){
    subfile <- datafile[c(j, j+1)]
    subfile$group <- NA
    subfile$subgroup <- NA
    subfile$genetype <- NA
    for(k in 1:nrow(subfile)){
      
      groupcontent <- unlist(strsplit(subfile[k, 2], split = '.', 
                                      fixed = TRUE))[1]
      grouplab <- paste0(groupcontent, '_', kind)
      
      subfile[k, 3] <- groupcontent
      subfile[k, 4] <- grouplab
      subfile[k, 5] <- kind
    }
    subfile$xcord <- 1:nrow(subfile)
    names(subfile)[1] <- 'rpm'
    
    if(allcut == FALSE){
      if(nrow(subfile) > 2000){
        midpoint <- round(nrow(subfile)/2)
        point1 <- midpoint - 1000
        point2 <- midpoint + 1000
        subfile1 <- subfile[1:point1,]
        subfile2 <- subfile[(point1+1):point2,]
        subfile3 <- subfile[(point2+1):nrow(subfile),]
        
        names(subfile)[1] <- 'rpm'
        cutoff <- as.vector(quantile(subfile[,1], qcutoff))
        low <- subset(subfile, rpm < cutoff)
        lowmax <- max(low[,1])
        
        
        names(subfile1)[1] <- 'rpm'
        subfile1$rpm[subfile1$rpm >= cutoff] <- lowmax
        
        names(subfile3)[1] <- 'rpm'
        subfile3$rpm[subfile3$rpm >= cutoff] <- lowmax
        
        names(subfile2)[1] <- 'rpm'
        
        subfile <- do.call(rbind, list(subfile1, subfile2, subfile3))
      }
      
    }else{
      
      names(subfile)[1] <- 'rpm'
      cutoff <- as.vector(quantile(subfile[,1], qcutoff))
      low <- subset(subfile, rpm < cutoff)
      lowmax <- max(low[,1])
      
      subfile$rpm[subfile$rpm >= cutoff] <- lowmax
      
    }
    
    
    
    if(j == 1){
      totalfile <- subfile
    }else{
      totalfile <- rbind(totalfile, subfile)
    }
  }  
  totalfile$range <- position
  
  samplenames <- rev(names(table(totalfile$group)))
  
  samplelist <- list()
  for(l in 1:length(samplenames)){
    samplelist[[l]] <- subset(totalfile, group == samplenames[l])
  }
  
  totalfile <- do.call(rbind, samplelist)
  
  totalfile$GeneType <- factor(totalfile$genetype, levels = kind, 
                               ordered = TRUE)
  
  return(totalfile)
  
  
}

processtab <- function(frame = tss_total, 
                       samplenames = sampletags){
  frame$condition <- NA
  groupnames <- names(table(frame$group))
  for(i in 1:length(samplenames)){
    frame$condition[frame$group == groupnames[i]] <- samplenames[i]
  }
  
  frame <- frame[-c(2, 3)]
  
  frame$rpm[frame$strand == 'antisense'] <- 
    -frame$rpm[frame$strand == 'antisense']
  frame$group <- paste0(frame$genetype, '\n', frame$strand)
  
  
  
  return(frame)
  
}

for(i in 1:length(halfwindows)){
  halfwindow <- halfwindows[i]
  
  tss <- generateplot_ttss(halfwidth = halfwindow, qcutoff = qcutoff, 
                           allcut = allcut)
  tts <- generateplot_ttss(ttsfile, 'around_tts', 
                           halfwidth = halfwindow, 
                           qcutoff = qcutoff, allcut = allcut)
  
  tss$strand <- 'sense'
  tts$strand <- 'sense'
  
  tss_anti <- generateplot_ttss(tssfile_anti, 'around_tss', 
                                halfwidth = halfwindow, 
                                qcutoff = qcutoff, allcut = allcut)
  
  tss_anti$strand <- 'antisense'
  
  
  tss_total <- rbind(tss, tss_anti)
  tts_total <- tts
  
  
  tss_total$Strand <- factor(tss_total$strand, levels = c('sense', 
                                                          'antisense'), 
                             ordered = TRUE)
  tts_total$Strand <- factor(tts_total$strand, levels = 'sense', 
                             ordered = TRUE)
  #########
  temptss <- tss_total
  temptts <- tts_total
  
  tss_total <- processtab()
  tts_total <- processtab(tts_total)
  
  kind <- unlist(strsplit(tssfile, split = '_', fixed = TRUE))
  kind <- kind[length(kind)-1]
  
  tss_total$group <- factor(tss_total$group, 
                            levels = paste0(kind, c('\nsense', '\nantisense')), 
                            ordered = TRUE)
  
  kind <- unlist(strsplit(ttsfile, split = '_', fixed = TRUE))
  kind <- kind[length(kind)-1]
  
  tts_total$group <- factor(tts_total$group, 
                            levels = paste0(kind, '\nsense'), 
                            ordered = TRUE)
  
  ################################
  
  library(ggplot2)
  library(scales)
  library(grid)
  
  pdf(sub('htseq', halfwindow, tsspdf_ori), height = 6, width = 7)
  grid.newpage()
  p <- ggplot(tss_total, mapping = aes(x = xcord, y = rpm, color=condition))
  
  p <- p + geom_line() + 
    scale_x_continuous(breaks = c(0, halfwindow, 2*halfwindow), 
                       labels = c(paste0('-', halfwindow), 'TSS', 
                                  paste0('+', halfwindow))) + 
    xlab('') + ylab('RPM') + 
    geom_vline(xintercept = halfwindow, linetype = 2, color = 'red') + 
    facet_grid(group~., scales = 'free_y') + 
    ggtitle(paste0('Range around TSS ', halfwindow, 'bp')) + 
    scale_color_discrete(name = 'condition') + 
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    theme(axis.text.x = element_text(size = 10)) + 
    theme(axis.text.x=element_text(hjust=1, angle = 90)) + 
    theme(panel.spacing = unit(0, 'lines'))
  
  library(grid)
  
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-r', g$layout$name))
  fills <- hue_pal()(2)
  k <- 1
  for(i in stripr){
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  grid.draw(g)
  
  dev.off()
  
  pdf(sub('htseq', halfwindow, ttspdf_ori), height = 4, width = 7)
  grid.newpage()
  p <- ggplot(tts_total, mapping = aes(x = xcord, y = rpm, color=condition))
  
  p <- p + geom_line() + 
    scale_x_continuous(breaks = c(0, halfwindow, 2*halfwindow), 
                       labels = c(paste0('-', halfwindow), 'TTS', 
                                  paste0('+', halfwindow))) + 
    xlab('') + ylab('RPM') + 
    geom_vline(xintercept = halfwindow, linetype = 2, color = 'red') + 
    facet_grid(group~., scales = 'free_y') + 
    ggtitle(paste0('Range around TTS ', halfwindow, 'bp')) + 
    scale_color_discrete(name = 'condition') + 
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    theme(axis.text.x = element_text(size = 10)) + 
    theme(axis.text.x=element_text(hjust=1, angle = 90))
  
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-r', g$layout$name))
  fills <- hue_pal()(1)
  k <- 1
  for(i in stripr){
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  grid.draw(g)
  dev.off() 
  
  ######################################
  write.table(tss_total, paste0('tss_total.', halfwindow, '.txt'), sep = '\t', 
              quote = FALSE, row.names = FALSE)
  write.table(tts_total, paste0('tts_total.', halfwindow, '.txt'), sep = '\t', 
              quote = FALSE, row.names = FALSE)
  
}
