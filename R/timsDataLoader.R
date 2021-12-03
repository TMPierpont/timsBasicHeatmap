#' This grabs flowjo files and extracts data
#' I could make this combine the files, but it'll be easier to use and more versitle if we just combine them outside this function
#' @export
flowJoDataLoader <- function(
sample_list, #sample info reference list
file, #flowjo File
fcs_files, #directory with FCS files
sample_info, #The columns that are not counts or samples that are to be transferred from the sample info sheet
main_event = NULL, #The population you want to treat as the main population when calculating totals.
debug=FALSE
) {
  if (debug) {cat("\nopening wsp file")}
  ws <- open_flowjo_xml(file)

  if (debug) {cat("... good!\nLoading fcs files")}
  gs <- flowjo_to_gatingset(ws, name = 1, path = fcs_files)
  nodelist <- gs_get_pop_paths(gs, path = "full")
  shortnodeList  <- gs_get_pop_paths(gs, path = "auto")

  ######################
  #match the sample file to the myData DF
  ######################
  myData <- as.data.frame(gs_pop_get_stats(gs, nodes = nodelist[1], xml = TRUE, type="percent")[,1])
  colnames(myData)[1] <- "sample"

  if (debug) {print("cleaning up wsp file IDs")}  
  #first clean up the fcs file names
  for (i in 1:nrow(myData)) {
    myData[i,'sample'] <- str_split(myData[i,'sample'],"\\.",simplify = TRUE)[1]
    sample_list[i,'sample'] <- str_split(sample_list[i,'sample'],"\\.",simplify = TRUE)[1]
  }
  
  if (debug) {cat("... good!\nremoving any '(1)' from sample list")}  
  sample_list$sample <- str_remove(sample_list$sample, "\\(1\\)") #sometimes if a file is saved twice, the fcs file gets a (1) and it's a pain in the...
  myData$sample <- str_remove(myData$sample, "\\(1\\)") #sometimes if a file is saved twice, the fcs file gets a (1) and it's a pain in the...

  if (debug) {cat("... good!\nPulling sample list info from into new dataframe and combining it with cleaned IDs from the WSP file\n")}   
  
    for (i in 1:nrow(myData)) {
      if (debug) {cat(paste("\r    Hunting for:",myData[i,'sample']))}
    replaceNum <- stringr::str_which(sample_list$sample, paste0(myData[i,'sample'],"$"))
    myData[i,'counts'] <- sample_list[replaceNum,'counts']
      for (a in 1:length(sample_info)) {
        myData[i, sample_info[a]] <- sample_list[replaceNum, sample_info[a]]
      }
    }

  #########################
  #Now load all the numbers
  #########################
  
  if (debug) {cat("... good!\nCHecking to if we should use another pop instead of root")}   
  #If someone wants to use something other than root as their main population to normalize to:
  if (is.null(main_event)) {root <- nodelist[1]} else { 
    if (debug) {cat(paste("\nUsing", main_event, "instead"))}   
    root <- main_event
    if (!any(nodelist == main_event)) { 
      print(paste("failed to find population", main_event, "only found:", nodelist))
      return(NULL)
    }
  }

  if (debug) {cat("... good!\nGenerating Totals and applying totals")}   
  #Totals first
  for (i in 2:length(nodelist)) {
    myData['new_col'] <- gs_pop_get_stats(gs, nodes = nodelist[i])[,3]

    #Get % total by dividing by total events ("Root")
    myData['new_col'] <- myData$new_col / gs_pop_get_stats(gs, nodes = root)[,3]
    #Multiply by total counts
    myData['new_col'] <- myData$new_col * as.numeric(myData$counts)
    colnames(myData)[ncol(myData)] = toString(paste("Total", shortnodeList[i]))
  }
  if (debug) {cat("... good!\nApplying percents")}   
  #Then Percentages
  for (i in 2:length(nodelist)) {
    myData['new_col'] <- gs_pop_get_stats(gs, nodes = nodelist[i], type="percent")[,3]
    colnames(myData)[ncol(myData)] = toString(paste("Percent", shortnodeList[i]))
  }
  if (debug) {cat("... good!\nDone! returning data!")}
  return(myData)
}
