library(data.table)
library(ggplot2)

raw_dir <- "results/raw_chunks/power"
fig_dir <- "results/figures/size_power"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

compute_curve <- function(A) {
  alpha <- seq(0,1,length.out=100)
  
  data.table(
    alpha = rep(alpha,3),
    power = c(
      sapply(alpha, function(a) mean(A$p_proposed <= a)),
      sapply(alpha, function(a) mean(A$p_schott <= a)),
      sapply(alpha, function(a) mean(A$p_cheng <= a))
    ),
    method = rep(c("Proposed","Schott","Cheng"), each=100)
  )
}

for(dtag in c("d20","d40","d60")) {
  
  for(s in sprintf("S%02d",1:40)) {
    
    files <- list.files(raw_dir,
                        pattern=paste0("^",s,"_",dtag),
                        full.names=TRUE)
    
    A <- rbindlist(lapply(files, function(f) readRDS(f)$results))
    
    curve <- compute_curve(A)
    
    p <- ggplot(curve, aes(alpha,power,color=method)) +
      geom_line(size=1) +
      geom_abline(slope=1,intercept=0,linetype="dashed") +
      labs(title=paste(s,dtag),
           x=expression(alpha),
           y="Rejection rate") +
      theme_minimal()
    
    ggsave(
      file.path(fig_dir,paste0("size_power_",s,"_",dtag,".png")),
      p,width=5,height=4
    )
  }
}

cat("All size-power plots created.\n")