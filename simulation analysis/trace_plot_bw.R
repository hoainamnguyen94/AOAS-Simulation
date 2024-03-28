library(ggplot2)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

id <- 15000:30000

title.name <- c(expression(paste(beta[1]^G)),
                expression(paste(beta[1]^D[1])),
                expression(paste(beta[1]^D[2])),
                expression(paste(beta[2]^G)),
                expression(paste(beta[2]^D[1])), 
                expression(paste(beta[2]^D[2])))

true.beta <- c(4,2,3,3,3,2)

trace <- function(id, beta, true.beta, title) {
  x <- 1:length(beta)
  data <- data.frame(x = x, beta = beta)
  p <- ggplot(data, aes(x = x, y = beta)) +
    geom_line(color = 'gray', linewidth = 0.4) +
    geom_hline(yintercept = median(beta[id]), linewidth = 0.4, color = 'black') +
    geom_hline(yintercept = true.beta, linewidth = 0.4, color = 'black', linetype = 'dotted') +
    geom_hline(yintercept = quantile(beta[id],0.025), linewidth = 0.4, color = 'black', linetype = 'dashed') +
    geom_hline(yintercept = quantile(beta[id],0.975), linewidth = 0.4, color = 'black', linetype = 'dashed') +
    xlab('Iteration') +
    scale_y_continuous(name = 'Estimate', limits = c(true.beta-1,max(beta,true.beta)+0.25)) +
    ggtitle(title) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
  return(p)
}

hist <- function(id, beta, title) {
  data <- data.frame(beta = beta[id])
  p <- ggplot(data, aes(x = beta)) +
    geom_histogram(aes(y = after_stat(density)), color = 'black', fill = 'gray', size = 0.1) +
    stat_function(fun = dnorm, args = list(mean = mean(data$beta), sd = sd(data$beta)), 
                  color = 'black', size = 0.5) +
    xlab('Estimate') +
    ylab('Density') +
    ggtitle(title) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
  return(p)
}

###########################################################################

post.beta <- read.csv("Fam-simulation-300-Miss-50/post_beta_seed1.csv") # path to the folder

p.trace <- p.hist <- as.list(1:6)

for (j in 1:6) {
  beta <- post.beta[,j]
  title <- title.name[j]
  p.trace[[j]] <- trace(id, beta, true.beta[j], title)
  p.hist[[j]] <- hist(id, beta, title)
}

ggarrange(p.trace[[1]], p.hist[[1]], p.trace[[2]], p.hist[[2]], p.trace[[3]], p.hist[[3]],
          p.trace[[4]], p.hist[[4]], p.trace[[5]], p.hist[[5]], p.trace[[6]], p.hist[[6]], 
          ncol = 4, nrow = 3)
ggsave("figures/sim_beta_bw.pdf", width = 10, height = 6)


