saneMarkers <- function(markers){
  #' @title Sanitizing Marker Names
  #' @author Julian Spagnuolo
  #'
  #' @param markers Character vector, taken from first row of each FACS data file imported into FACkit GUI
  #'
  #' @importFrom stringr str_split
  #'
  #' @export
  #'
  #'

# Split flowJo marker names on "::" (if present), or remove other punctuation from name
x <- str_split(markers, pattern="[:punct:]", simplify = F)

for(i in 1:length(x))
{
  markers[i] <- x[[i]][length(x[[i]])]
}

# remove white space
markers <- gsub("[^a-zA-Z0-9]", "", markers)


# Check for numeric prefix and add X to start if detected
x <- str_split(markers, pattern="", simplify = F)

for(i in 1:length(x))
{
  if(x[[i]][1] %in% as.character(c(0,1,2,3,4,5,6,7,8,9)))
  {
    x[[i]] <- append(x[[i]], "X", after=0)
    markers[i] <- paste(x[[i]], collapse = "")
  }else{
    markers[i] <- paste(x[[i]], collapse = "")
  }
}

return(markers)
}
