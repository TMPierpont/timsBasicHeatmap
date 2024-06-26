#' This takes a matrix of data with certain columns and converts them into a 3 or 4 D matrix that can be read by Tim's pvheatmap package

#' @export

flowDataCruncher <- function(
  my_data, #This is your data, in the form of a matrix, needs to have column names, but rownames are ignored.
  selected_type = NULL, #If you have multiple types, this is the type that you want to use right now
  type = NULL, #If you have multiple types, this identifies the column
  group, #These are the groups that get compared, give the column name
  subgroup = NULL, #If you have sub groups identified, this column will split them into however many sub boxes are needed.
  ID=NULL, #Column that gives the ID of the samples
  ignore=NULL, #The columns you want to ignore/ not be graphed or used
  combine_groups = NULL, #Set this to combine groups, takes a list c(name_of_new_column, list_of_columns_to_combine)
  control_group = NULL #The group that is the control group for the Dunnett test, #HOWEVER THIS NEEDS TO BE ADDED TO THE SCRIPT!
) {
  
  
  #If the group is numeric, we need to fix that
  if (apply(my_data[group], 2, is.numeric)) {
    max_sample <- nchar(toString(max(my_data[group]))) #get the biggest number, convert it to string, and then count the characters
    temp_list <- unlist(as.list(my_data[group])) #yes... I know I made it a list and unlisted, but somehow that's the only way I could find to get the desired type
    my_data[group] <- formatC(temp_list, width=max_sample, flag="0") #then convert them to 00x format as characters
  }

  
  #remove columns I won't need
  if (!is.null(ignore)){
    for (a in 1:length(ignore)) {
      if (length(which(colnames(my_data) == ignore[a])) == 1) {#Is there a type column?
        my_data <- my_data[-which(colnames(my_data) == ignore[a])] #Then delete it.
      }
    }
  }
  
  #remove samples for now as well
  if (!is.null(ID)){
    if (length(which(colnames(my_data) == ID)) == 1) {#Is there a type column?
      my_data <- my_data[-which(colnames(my_data) == ID)] #Then delete it.
    }
  }
  
  #Are there multiple types?
  if (length(which(colnames(my_data) == type)) == 1) {#Is there a type column?
    type_count <- nrow(unique(my_data[which(colnames(my_data) == type)]))
    if (nrow(unique(my_data[which(colnames(my_data) == type)])) == 1) {#Then there's only one unique name, no need to use this, delete it.
      my_data <- my_data[-which(colnames(my_data) == type)] 
    } else {
      #We only want to use one type, so select only that type if there's more than one.
      colnames(my_data)[which(colnames(my_data) == type)] <- "qqqaaa"
      my_data <- dplyr::filter(my_data, my_data$qqqaaa == selected_type)
      my_data <- my_data[-which(colnames(my_data) == type)]
    }
  }
  #Are there multiple subgroups?
  if (length(which(colnames(my_data) == subgroup)) == 1) {#Is there a type column?
    type_subgroup <- nrow(unique(my_data[which(colnames(my_data) == subgroup)]))
    if (nrow(unique(my_data[which(colnames(my_data) == subgroup)])) == 1) {#Then there's only one unique name, no need to use this, delete it.
      my_data <- my_data[-which(colnames(my_data) == subgroup)]
      subgroup_names <- NULL
    } else {
      subgroup_names <- unique(my_data[which(colnames(my_data) == subgroup)])
      colnames(my_data)[which(colnames(my_data) == subgroup)] <- "qqqwww"
    }
  } else { 
    type_subgroup = 1
    subgroup_names = NULL} 

  #Are there multiple groups? There kind of needs to be
  if (length(which(colnames(my_data) == group)) == 1) {#Is there a groups column?
    group_count <- nrow(unique(my_data[which(colnames(my_data) == group)]))
    if (nrow(unique(my_data[which(colnames(my_data) == group)])) == 1) {#Then there's only one unique name
      print("not enough groups")
    } else {
      colnames(my_data)[which(colnames(my_data) == group)] <- "group"
    }
  }

  #Are there group names to combine? If so, combine them.
  if (!is.null(combine_groups)) {
    new_name <- combine_groups[1]
    old_names <- combine_groups[-1]
    for (a in 1:length(old_names)) {#Loop and get rid of the old group names
      my_data[which(my_data$group == old_names[a]), 'group'] <- new_name
    }
  }

  #This needs to happen to make sure the control group is always the first group... Only because I can't figure out how to get the Dunnett test to work properly otherwise.
    my_data <- as.data.frame(my_data)
    my_data$group[which(my_data$group == control_group)] <- "Control"
    my_data$group <- paste0("zz-zz", my_data$group)
    my_data$group[which(my_data$group == "zz-zzControl")] <- "Control"  
  
  
  if (!is.null(subgroup)) {colnames(my_data)[which(colnames(my_data) == "qqqwww")] <- "subgroup"}
  
  for(i in 1:type_subgroup) { #if there's more than one subgroup, we'll want to build a multidimensional arrays
    for(ii in 2:(ncol(my_data))) {
      #print(paste("sub_box",ii,"is",subgroup_names))
      if (colnames(my_data)[ii] == "group" || colnames(my_data)[ii] == "subgroup") {next} #don't graph these guys
      
      graphDataMain <- my_data[which(colnames(my_data) == "group")]

      graphDataMain[2] <- my_data[ii]
      population <- colnames(my_data)[ii]
      colnames(graphDataMain)[2] <- "value"
      if (!is.null(subgroup_names)) {
        graphDataMain[3] <- my_data[which(colnames(my_data) == "subgroup")]
        }

      if (!is.null(subgroup_names)) {#if there's a subgroup name, we'll pull down the current samples for that subgroup, otherwise, use them all
        graphData <- dplyr::filter(graphDataMain, graphDataMain$subgroup == toString(subgroup_names[i,1]))
      } else { #If there isn't any subgroups, then just pull the whole thing over
        graphData <- graphDataMain
      }

            #rename control group to "control"
      graphData <- dplyr::filter(graphData, !is.na(graphData$value)) #Just delete any rows without numbers
      graphData$value <- as.numeric(graphData$value) #The behavior changed and these need to be explicately converted to numbers now
      means <- aggregate(graphData['value'], list(graphData$group), mean)
      conMeans <- dplyr::filter(means, means$Group.1 == "Control")[1,2]
      normed <- log2(means[2] / as.numeric(conMeans))
      pvalues <- SEM_high <- SEM_low <- SEM <- means #lazy way to recreate the new matrix, meh, probably a better way but this setup works
      for(a in 1:nrow(SEM)) {
        SEM[a,2] = std.error(graphData[which(graphData$group==SEM[a,1]),'value'])
      }
            
      SEM_high[2] <- log2((means[2] + SEM[2]) / as.numeric(conMeans))
      SEM_low[2] <- log2((means[2] - SEM[2]) / as.numeric(conMeans))
      
      #Get the stats
      ANOVA <- aov(value ~ group, data = graphData)
      ANOVA <- summary(ANOVA)[[1]][["Pr(>F)"]][[1]]  
      Dunnett <- DescTools::DunnettTest(x = graphData$value, control = "Control", g = as.factor(graphData$group))
      pvalues[1,2] <- 1
      for(a in 2:nrow(SEM)) {
        pvalues[a,2] <- Dunnett$Control[[a-1,4]]
      }

      #That's all for the hard stuff, now just package it up...
      
      ############
      #I'll want to add the bit about the ANOVA testing and stuff, but for now this will test things
      ############


      if (exists('Nmeans_2d')) {
        Nmeans_2d[ncol(Nmeans_2d)+1] <- normed
        colnames(Nmeans_2d)[ncol(Nmeans_2d)] = population
      } else {
        Nmeans_2d <- normed
        rownames(Nmeans_2d) <- means[,1]
        colnames(Nmeans_2d)[1] = population
      }
      
      if (exists('sem_high_2d')) {
        sem_high_2d[ncol(sem_high_2d)+1] <- SEM_high[2]
        colnames(sem_high_2d)[ncol(sem_high_2d)] = population
      } else {
        sem_high_2d <- SEM_high[2]
        rownames(sem_high_2d) <- means[,1]
        colnames(sem_high_2d)[1] = population
      }
      
      if (exists('sem_low_2d')) {
        sem_low_2d[ncol(sem_low_2d)+1] <- SEM_low[2]
        colnames(sem_low_2d)[ncol(sem_low_2d)] = population
      } else {
        sem_low_2d <- SEM_low[2]
        rownames(sem_low_2d) <- means[,1]
        colnames(sem_low_2d)[1] = population
      }
      
      if (exists('pvalues_2d')) {
        pvalues_2d[ncol(pvalues_2d)+1] <- pvalues[2]
        colnames(pvalues_2d)[ncol(pvalues_2d)] = population
      } else {
        pvalues_2d <- pvalues[2]
        rownames(pvalues_2d) <- means[,1]
        colnames(pvalues_2d)[1] = population
      }
 
      if (exists('means_2d')) {
        means_2d[ncol(means_2d)+1] <- means[,2]
        colnames(means_2d)[ncol(means_2d)] = population
      } else {
        means_2d <- means[2]
        rownames(means_2d) <- means[,1]
        colnames(means_2d)[1] = population
      }
           
    }
    
    if (exists('full_data')) {
      full_data[,,1,i] <- as.matrix(Nmeans_2d)
      full_data[,,2,i] <- as.matrix(sem_high_2d)
      full_data[,,3,i] <- as.matrix(sem_low_2d)
      full_data[,,4,i] <- as.matrix(pvalues_2d)
      full_data[,,5,i] <- as.matrix(means_2d)
      
    } else {
      matrix_names <- c("normalized means", "SEM high", "SEM low", "pvalues", "means")
      if (!is.null(subgroup_names)){test <- c(as.matrix(subgroup_names))
      full_data <- array(c(as.matrix(Nmeans_2d), as.matrix(sem_high_2d), as.matrix(sem_low_2d), as.matrix(pvalues_2d), as.matrix(means_2d)), 
        dim = c(nrow(Nmeans_2d),ncol(Nmeans_2d),5,type_subgroup),
        dimnames = list(rownames(Nmeans_2d), colnames(Nmeans_2d), matrix_names, test))
      } else {
        full_data <- array(c(as.matrix(Nmeans_2d), as.matrix(sem_high_2d), as.matrix(sem_low_2d), as.matrix(pvalues_2d), as.matrix(means_2d)), 
                           dim = c(nrow(Nmeans_2d),ncol(Nmeans_2d),5,type_subgroup),
                           dimnames = list(rownames(Nmeans_2d), colnames(Nmeans_2d), matrix_names))
        
        
      }
    }
    
    rm(Nmeans_2d)
    rm(sem_high_2d)
    rm(sem_low_2d)
    rm(pvalues_2d)
    rm(means_2d)    
  }
  
  #change back the control group, these two changes were to get the control group to always sort as the top even if the unique subgroups are integers
  rownames(full_data)[rownames(full_data) == "Control"] = control_group
  rownames(full_data) <- str_remove(rownames(full_data), "zz-zz")
  
  return(full_data)
}

