#' This makes a heatmap with intergrated significance values (one line for p<0.05 and two lines for p<0.005)
#' 
#' It's my first attempt at a package, don't expect to much from it! 2021-11-15 and has very little easy control for now, just as a utility for specific needs of the lab
#' 
#' Call the function timsBasicHeatmap(means, pvalues)
#' 
#' means is a matrix of means with appropriate column and row names that will be used in the heatmap graphic
#' 
#' pvalue is a matching matrix with the pvalue that accompanies the mean, The size and layout of the matrices must be identical
#' 
#' Columns should be phenotypes, rows should be groups

#' @export

timsBasicHeatmap <- function(means, pvalues){ #just takes means and pvalues of graphs, makes a reasonable sized heatmap
  
  savePDF <- "no"
  if(savePDF == "yes"){pdf(paste0(directory,"/pdf/",name,"-Heatmap.pdf"), width = 16, height = 8)}
  #means <- values
  
  X <- ncol(means)
  Y <- nrow(means)
  
  for(y in 1:Y) {
    for(x in 1:X) {
      INT <- 1
      
      ## Create a viewport, rectangle, and text
      vp<-viewport(x=x/50+0.05,y=.95-y/25,width=0.019,height = 0.038)
      pushViewport(vp)
  
      values <- means[y,x]
      if (values > 1) {values <- 1}
      if (values < -1) {values <- -1}
      RED = -0.398*(values^2) + 0.2826*values + 0.9637
      GREEN = -0.6374*(values^2) - 0.1367*values + 0.9623
      BLUE = 0.2225*(values^3) - 0.3484*(values^2) - 0.4966*(values) + 0.7778
      if (RED > 1) {RED <- 1}
      if (GREEN > 1) {GREEN <- 1}
      if (BLUE > 1) {BLUE <- 1}
      
      grid.rect(x=0.5, y=0.5, height=1,width=1,gp=gpar(fill=rgb(RED,GREEN,BLUE,INT),col=rgb(0,0,0,1)))
      
      sig <- pvalues[y,x]
        
      if (sig<.005) {
        grid.rect(x=.5,y=.5,width=1,height = .15, gp=gpar(fill=rgb(0,0,0,0), col=rgb(1,1,1,1)))
      } else if (sig<.05) {
        grid.rect(x=.5,y=.5,width=1,height = 0.001, gp=gpar(fill=rgb(0,0,0,0), col=rgb(1,1,1,1)))
      }
      grid.rect(x=.5,y=.5,width=1,height = 1, gp=gpar(fill=rgb(0,0,0,0), col=rgb(0,0,0,1)))
      vp <- popViewport(1)
      sig <- 1
      
    }
  }
  
  for(x in 1:X) {
    grid.text(colnames(means)[x], just=0,x=x/50+0.05, y=.95-(Y+.7)/25, rot=-90,
              gp=gpar(fontsize=10, col="black"))
  }
  
  for(y in 1:Y) {
    grid.text(rownames(means)[y], just=0,x=(X+.7)/50+0.05, y=.95-y/25, rot=,
              gp=gpar(fontsize=12, col="black"))
  }
    
  if(savePDF == "yes"){dev.off()}
  
}

