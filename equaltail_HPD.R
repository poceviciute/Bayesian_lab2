#Equal tail credible interval
#perc is the percentile that should be excluded from both ends. If 95% equal tail credible interval
#then perc=0.025
#values is a dataframe or matrix
#function finds lower and upper limits for the interval by row (if dataframe/matrix contains more than 
#one column)
#returns dataframe with lowertails and uppertails

equal_tail <- function(perc, values){
  values <- as.matrix(values)
  #if values has one column find interval like this:
  if(ncol(values)==1){
    nr <- nrow(values)
    #the index corresponding to 100*perc%
    ind <- perc*nr
    #order the values and choose the (ind+1)th value as lowertail
    low <- values[order(values, decreasing=FALSE)[ind+1]]
    #order the values and choose the (nr-(ind+1))th value as uppertail
    upp <- values[order(values, decreasing=FALSE)[nr-ind-1]]
  }
  #otherwise find interval like this:
  else{
    nr <- ncol(values)
    #the index corresponding to 100*perc%
    ind <- perc*nr
    #for each row order the values in the columns and choose the (ind+1)th value as lowertail
    low <- apply(values,1,function(x){x[order(x, decreasing = FALSE)[ind+1]]})
    #for each row order the values in the columns and choose the (nr-(ind+1))th value as uppertail
    upp <- apply(values,1,function(x){x[order(x, decreasing = FALSE)[nr-ind-1]]}) 
  }
  #return dataframe with the lowertail and uppertail values
  data.frame(lowertail=low, uppertail=upp)
}

#highest posterior density interval
#int_length is the length of the interval. Eg. if we want 95% HPD interval, then int_length=0.95
#df is a dataframe with x-values and y-values (need to be named as x and y)
#returns the x-values that belong in the HDP interval

HPD <- function(int_length, df){
  #new dataframe
  dataframe <- data.frame(y=df$y,x=df$x) 
  #order the y values
  dataframe <- dataframe[order(dataframe$y,decreasing = TRUE),] 
  #calculate density of the y values
  dataframe$dens <- cumsum(dataframe$y)/sum(dataframe$y) 
  #check which densities are smaller than int_length (ex. 0.95 if 95% HPD interval)
  which_x<-dataframe$dens<int_length
  #return the corresponding x-values
  dataframe$x[which_x]
}