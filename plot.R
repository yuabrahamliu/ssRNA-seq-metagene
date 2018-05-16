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
  idx <- grep('project', configs)
  project <- configs[idx]
  idx <- grep('sampletags', configs)
  sampletags <- configs[idx]
  idx <- grep('halfwindow', configs)
  halfwindows <- configs[idx]
}

fig_tag <- unlist(strsplit(fig_tag, split = ' = '))[2]
project <- unlist(strsplit(project, split = ' = '))[2]
sampletags <- unlist(strsplit(sampletags, split = ' = '))[2]
halfwindows <- unlist(strsplit(halfwindows, split = ' = '))[2]


result_dir <- paste0(getwd(), '/', project, '/result')
setwd(result_dir)

tssfiles <- paste0('sense', fig_tag, c('_cds', '_lnc', '_nonpaf', '_paf'), 
                   '_tss.htseq.txt')
ttsfiles <- paste0('sense', fig_tag, c('_cds', '_lnc', '_nonpaf', '_paf'), 
                   '_tts.htseq.txt')

tssfiles <- paste0('antisense', fig_tag, c('_cds', '_lnc', '_nonpaf', '_paf'), 
                   '_tss.htseq.txt')
ttsfiles <- paste0('antisense', fig_tag, c('_cds', '_lnc', '_nonpaf', '_paf'), 
                   '_tts.htseq.txt')

tsspdf_ori <- paste0(fig_tag, '_meta_tss.htseq.pdf')
ttspdf_ori <- c(fig_tag, '_meta_tts.htseq.pdf')


sampletags <- unlist(strsplit(sampletags, split = ','))
halfwindows <- as.numeric(unlist(strsplit(halfwindows, split = ',')))


##############################

generateplot_ttss <- function(filenames = tssfiles, 
                              position = 'around_tss', 
                              halfwidth = 500, 
                              qcutoff = 0.99, allcut = FALSE){
  file1 <- read.table(filenames[1], sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE, 
                      quote = '')
  file2 <- read.table(filenames[2], sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE, 
                      quote = '')
  file3 <- read.table(filenames[3], sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE, 
                      quote = '')
  file4 <- read.table(filenames[4], sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, check.names = FALSE, 
                      quote = '')
  
  genelen <- nrow(file1)
  
  filelist <- list(file1, file2, file3, file4)
  
  for(i in 1:length(filelist)){
    
    datafile <- filelist[[i]]
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
        if(i == 1){
          kind <- 'P_coding'
          grouplab <- paste0(groupcontent, '_', kind)
        }else if(i == 2){
          kind <- 'lnc'
          grouplab <- paste0(groupcontent, '_', kind)
        }else if(i == 3){
          kind <- 'nonPaf'
          grouplab <- paste0(groupcontent, '_', kind)
        }else if(i == 4){
          kind <- 'Paf'
          grouplab <- paste0(groupcontent, '_', kind)
        }
        
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
        reshapefile <- subfile
      }else{
        reshapefile <- rbind(reshapefile, subfile)
      }
    }
    if(i == 1){
      totalfile <- reshapefile
    }else{
      totalfile <- rbind(totalfile, reshapefile)
    }
  }
  totalfile$range <- position
  
  samplenames <- rev(names(table(totalfile$group)))
  
  samplelist <- list()
  for(l in 1:length(samplenames)){
    samplelist[[l]] <- subset(totalfile, group == samplenames[l])
  }
  
  totalfile <- do.call(rbind, samplelist)
  
  totalfile$GeneType <- factor(totalfile$genetype, levels = c('lnc', 'Paf', 
                                                              'P_coding', 
                                                              'nonPaf'), 
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
  
  tss <- generateplot_ttss(halfwidth = halfwindow)
  tts <- generateplot_ttss(ttsfiles, 'around_tts', 
                           halfwidth = halfwindow, 
                           qcutoff = 0.95, allcut = TRUE)
  
  tss$strand <- 'sense'
  tts$strand <- 'sense'
  tts <- subset(tts, genetype == 'Paf')
  
  
  tss_anti <- generateplot_ttss(tssfiles_anti, 'around_tss', 
                                halfwidth = halfwindow, 
                                qcutoff = 0.95)
  
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
  
  tss_total$group <- factor(tss_total$group, 
                            levels = c('lnc\nsense', 'lnc\nantisense', 
                                       'Paf\nsense', 'Paf\nantisense', 
                                       'P_coding\nsense', 
                                       'P_coding\nantisense', 
                                       'nonPaf\nsense', 'nonPaf\nantisense'), 
                            ordered = TRUE)
  tts_total$group <- factor(tts_total$group, 
                            levels = c('Paf\nsense'), 
                            ordered = TRUE)
  
  ################################
  
  library(ggplot2)
  library(scales)
  library(grid)
  
  pdf(sub('htseq', halfwindow, tsspdf_ori[1]), height = 7, width = 8)
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
    theme(axis.text.x = element_text(size = 10)) + 
    theme(axis.text.x=element_text(hjust=1, angle = 90)) + 
    theme(panel.spacing = unit(c(0, 0.375, 0, 0.375, 0, 0.375, 0), 'lines'))
  
  library(grid)
  
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-r', g$layout$name))
  fills <- hue_pal()(4)
  fills <- rep(fills, each = 2)
  k <- 1
  for(i in stripr){
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  grid.draw(g)
  
  dev.off()
  
  pdf(sub('htseq', halfwindow, ttspdf_ori[1]), height = 4, width = 7)
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
    theme(axis.text.x = element_text(size = 10)) + 
    theme(axis.text.x=element_text(hjust=1, angle = 90))
  
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-r', g$layout$name))
  fills <- hue_pal()(4)[2]
  fills <- rep(fills, each = 2)
  k <- 1
  for(i in stripr){
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  grid.draw(g)
  dev.off() 
  
  ######################################
  temptss <- tss_total
  temptts <- tts_total
  
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
  
  tss_total <- processtab()
  tts_total <- processtab(tts_total)
  
  tss_total$group <- factor(tss_total$group, 
                            levels = c('lnc\nsense', 'lnc\nantisense', 
                                       'Paf\nsense', 'Paf\nantisense', 
                                       'P_coding\nsense', 
                                       'P_coding\nantisense', 
                                       'nonPaf\nsense', 'nonPaf\nantisense'), 
                            ordered = TRUE)
  tts_total$group <- factor(tts_total$group, 
                            levels = c('Paf\nsense'), 
                            ordered = TRUE)
  
  
  ##############
  write.table(tss_total, paste0('tss_total.', halfwindow, '.txt'), sep = '\t', 
              quote = FALSE, row.names = FALSE)
  write.table(tts_total, paste0('tts_total.', halfwindow, '.txt'), sep = '\t', 
              quote = FALSE, row.names = FALSE)
  
}
