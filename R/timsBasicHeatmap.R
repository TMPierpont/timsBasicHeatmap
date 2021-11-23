#' This makes a heatmap with integrated significance values (one line for p<0.05 and two lines for p<0.005)
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
PVHeatmap <- function(
  values, #These are the values that generate the main colors on the boxes
  crunched_data = FALSE, #if this is being passed from one of Tim's crunched data sets, let it know!
  pvalues = NULL, #These are the pvalues, ignored if left blank
  sem_high = NULL, #These add high SEM values, should be values+SEM and then normalized appropriately.
  sem_low = NULL, #These add low SEM values, should be values+SEM and then normalized appropriately.

  transpose = FALSE, #Not yet implemented, don't use
  
  title = NULL, #If you want to give the graph a title
  font_size = NULL, #Font size, will set base size of fonts, title will be this * 1.3
  font_size_x = NULL, #Specific font size for x axis
  font_size_y = NULL, #Specific font size for y axis
  font_size_title = NULL, #Speicify font size for title.
  
  #These will let you push title, map, or either axis around a bit
  nudge_title_verticle = 0,
  nudge_title_horizontal = 0,
  nudge_x_axis_verticle = 0,
  nudge_x_axis_horizontal = 0,
  nudge_y_axis_verticle = 0,
  nudge_y_axis_horizontal = 0,
  nudge_map_verticle = 0,
  nudge_map_horizontal = 0,
  gap = .1, #Fraction of boxes that are given between boxes
  box_size = NULL, #fraction of screen size that a box size will be
  
  justify_y = 0, #will take -1 to 1
  sort_y = NULL #If this exists then it'll reorganize the data to match
  ){ #just takes means and pvalues of graphs, makes a reasonable sized heatmap
  
  
  #If this is crunched data, then sort that out properly first.
  if (crunched_data) {
    crunched_data <- values
    
    if (!is.null(sort_y)) {
      reorder_y <- match(sort_y,rownames(crunched_data))
      if (any(is.na(reorder_y))) {
        print("Not all y axis names given match actual row names given")
        print(sort(sort_y))
        print(sort(rownames(crunched_data)))
        }
      crunched_data <- crunched_data[reorder_y,,,]
    }
    
    #this is a hacky way to get this info, it must be easier, but this will work...
    num_of_matrix <- length(crunched_data[,,,])/length(crunched_data[,,,1])
    
    #if there are multiple subgroups or just one
    if (num_of_matrix == 1) {
      values <- crunched_data[,,1,1]
      sem_high <- crunched_data[,,2,1]
      sem_low <- crunched_data[,,3,1]
      pvalues <- crunched_data[,,4,1]
    } else {
      pvalues <- sem_low <- sem_high <- values <- array(dim=c(nrow(values),ncol(values),num_of_matrix))
      for(a in 1:num_of_matrix) {
        values[,,a] <- crunched_data[,,1,a]
        sem_high[,,a] <- crunched_data[,,2,a]
        sem_low[,,a] <- crunched_data[,,3,a]
        pvalues[,,a] <- crunched_data[,,4,a]
      }
      rownames(values) <- rownames(sem_high) <- rownames(sem_low) <- rownames(pvalues) <- rownames(crunched_data)
      colnames(values) <- colnames(sem_high) <- colnames(sem_low) <- colnames(pvalues) <- colnames(crunched_data)
    }
  }
  layers <- length(dim(values)) #Check if this is a simple heatmap, or one with several layers
  graphics.off()#clear plot
  ratio <- dev.size("px")[1] / dev.size("px")[2] #Get screen ratio
  
  if (transpose) {
    Y <- ncol(values)
    X <- nrow(values) 
  } else {
    X <- ncol(values)
    Y <- nrow(values)
  }
  
  
  if (is.null(box_size)) { box_size <- .8/(max(X/ratio,Y)) } #Establishes box size based on window
  gap <- 1-gap

  #If there is a font size, we'll want to adjust that here
  if (is.null(font_size)) {
    hm_font_y = 12
    hm_font_x = 12
    hm_font_title = 18
  } else {
    hm_font_y = font_size
    hm_font_x = font_size
    hm_font_title = as.numeric(font_size*1.3)
  }
  if (!is.null(font_size_y)) {hm_font_y=font_size_y}
  if (!is.null(font_size_x)) {hm_font_x=font_size_x}
  if (!is.null(font_size_title)) {hm_font_title=font_size_title}
    
  #Process all of the nudges and justify, sets them relative to the box size
  nudge_title_verticle = nudge_title_verticle*box_size
  nudge_title_horizontal = nudge_title_horizontal*box_size
  nudge_x_axis_verticle = nudge_x_axis_verticle*box_size
  nudge_x_axis_horizontal = nudge_x_axis_horizontal*box_size
  nudge_y_axis_verticle = nudge_y_axis_verticle*box_size
  nudge_y_axis_horizontal = nudge_y_axis_horizontal*box_size
  nudge_map_verticle = nudge_map_verticle*box_size
  if (!is.null(title)) {nudge_map_verticle = nudge_map_verticle-2*box_size}
  nudge_map_horizontal = nudge_map_horizontal*box_size
  

  
  grid.text(title, just=0,x=.5+nudge_title_horizontal, y=.95+nudge_title_verticle, rot=0,
            gp=gpar(fontsize=hm_font_title, col="black"))
  SEMsize = .1
  
  #For single box experiments
  if (layers == 2) {
    
    #Make sure the order is what the person wants
    if (!is.null(sort_y)) {
      reorder_y <- match(sort_y,rownames(values))
      if (any(is.na(reorder_y))) {
        print("Not all y axis names given match actual row names given")
        print(sort(sort_y))
        print(sort(rownames(crunched_data)))
      }
      values <- values[reorder_y,]
      if (!is.null(sem_high)) {sem_high <- sem_high[reorder_y,]}
      if (!is.null(sem_low)) {sem_low <- sem_low[reorder_y,]}
      if (!is.null(pvalues)) {pvalues <- pvalues[reorder_y,]}
    }
    
    for(y in 1:Y) {
      for(x in 1:X) {
        ## Create a viewport, rectangle
        vp<-viewport(x=x*box_size/ratio+nudge_map_horizontal,y=1+nudge_map_verticle-y*box_size,width=box_size/ratio*gap,height=box_size*gap)
        pushViewport(vp)
        
        #Add base color
        color <- PVHeatmapColor(values[y,x])
        grid.rect(x=0.5, y=0.5, height=1,width=1,gp=gpar(fill=color,col=rgb(0,0,0,1)))
        
        #Was SEM provided? then draw the SEMs
        if (!is.null(sem_high)) {
          color <- PVHeatmapColor(sem_high[y,x])
          grid.rect(x=0.5, y=SEMsize/2, height=SEMsize,width=1,gp=gpar(fill=color,col=rgb(0,0,0,0)))
        }
        if (!is.null(sem_low)) {
          color <- PVHeatmapColor(sem_low[y,x])
          grid.rect(x=0.5, y=1-SEMsize/2, height=SEMsize,width=1,gp=gpar(fill=color,col=rgb(0,0,0,0)))
        }      
        
        #were pvalues given? add them if so
        if (!is.null(pvalues)) {
          PVHeatmapDrawP(pvalues[y,x])
        }
        
        grid.rect(x=.5,y=.5,width=1,height = 1, gp=gpar(fill=rgb(0,0,0,0), col=rgb(0,0,0,1))) #Draws border around box
        vp <- popViewport(0) #Pop the port
      }
    }
    
  } else if (layers > 2) { #For drawing several boxes inside a single box
    exps <- dim(values)[3] #then get the number of matrices, this will be how many times we loop
    
    #Make sure the order is what the person wants
    if (!is.null(sort_y)) {
      reorder_y <- match(sort_y,rownames(values))
      if (any(is.na(reorder_y))) {
        print("Not all y axis names given match actual row names given")
        print(sort(sort_y))
        print(sort(rownames(crunched_data)))
      }
      values <- values[reorder_y,,]
      if (!is.null(sem_high)) {sem_high <- sem_high[reorder_y,,]}
      if (!is.null(sem_low)) {sem_low <- sem_low[reorder_y,,]}
      if (!is.null(pvalues)) {pvalues <- pvalues[reorder_y,,]}
    }
    
      for(y in 1:Y) {
        for(x in 1:X) {
          ## Create a viewport, rectangle
          vp<-viewport(x=x*box_size/ratio+nudge_map_horizontal,y=1+nudge_map_verticle-y*box_size,width=box_size/ratio*gap,height=box_size*gap)
          pushViewport(vp)
          
          for(z in 1:exps) {
            vp<-viewport(x=z*1/exps-1/exps/2,y=.5,width=1/exps,height=1)
            pushViewport(vp)
            
            #Add base color
            color <- PVHeatmapColor(values[y,x,z])
            grid.rect(x=0.5, y=0.5, height=1,width=1,gp=gpar(fill=color,col=rgb(0,0,0,1)))
            
            #Was SEM provided? then draw the SEMs
            if (!is.null(sem_high)) {
              color <- PVHeatmapColor(sem_high[y,x,z])
              grid.rect(x=0.5, y=SEMsize/2, height=SEMsize,width=1,gp=gpar(fill=color,col=rgb(0,0,0,0)))
            }
            if (!is.null(sem_low)) {
              color <- PVHeatmapColor(sem_low[y,x,z])
              grid.rect(x=0.5, y=1-SEMsize/2, height=SEMsize,width=1,gp=gpar(fill=color,col=rgb(0,0,0,0)))
            }      
            
            #were pvalues given? add them if so
            if (!is.null(pvalues)) {
              PVHeatmapDrawP(pvalues[y,x,z])
            }
            grid.rect(x=.5,y=.5,width=1,height = 1, gp=gpar(fill=rgb(0,0,0,0), col=rgb(0,0,0,.5))) #Draws border around box
            vp <- popViewport(1)
          }
          grid.rect(x=.5,y=.5,width=1,height = 1, gp=gpar(fill=rgb(0,0,0,0), col=rgb(0,0,0,1))) #Draws border around box
          vp <- popViewport(0) #Pop the port
        }
      }
    
  }
  
  for(x in 1:X) {
    grid.text(colnames(values)[x], just=0,x=x*box_size/ratio+nudge_map_horizontal+nudge_x_axis_horizontal, y=(1+nudge_map_verticle+nudge_x_axis_verticle)-(Y+1)*box_size, rot=-90,
              gp=gpar(fontsize=hm_font_x, col="black"))
  }
  
  for(y in 1:Y) {
    grid.text(rownames(values)[y], just=justify_y,x=(X+1)*box_size/ratio+nudge_y_axis_horizontal+nudge_map_horizontal, y=1+nudge_map_verticle+nudge_y_axis_verticle-y*box_size, rot=,
              gp=gpar(fontsize=hm_font_y, col="black"))
  }
  
}

#' @export
PVHeatmapColor <- function(num){
  
  if (num > 1) {num <- 1}
  if (num < -1) {num <- -1}
  RED = -0.398*(num^2) + 0.2826*num + 0.9637
  GREEN = -0.6374*(num^2) - 0.1367*num + 0.9623
  BLUE = 0.2225*(num^3) - 0.3484*(num^2) - 0.4966*(num) + 0.7778
  if (RED > 1) {RED <- 1}
  if (GREEN > 1) {GREEN <- 1}
  if (BLUE > 1) {BLUE <- 1}
  color <- rgb(RED,GREEN,BLUE,1)
  return(color) 
  
}

#' @export
PVHeatmapDrawP <- function(sig) {
  
  if (sig<.005) {
    grid.rect(x=.5,y=.5,width=1,height = .15, gp=gpar(fill=rgb(0,0,0,0), col=rgb(1,1,1,1)))
  } else if (sig<.05) {
    grid.rect(x=.5,y=.5,width=1,height = 0.001, gp=gpar(fill=rgb(0,0,0,0), col=rgb(1,1,1,1)))
  }
}
