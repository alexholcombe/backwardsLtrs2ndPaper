#Turn dataframe with array columns into list columns, so dplyr can work

#d<- data.frame(x=c(1,2,3,4), y= matrix(c(1,2,3, 11,12,13, 20,21,22, 23,24,25) ,nrow=4,ncol=3) )
#Above gives different results, variables auto-turned into y.1, y.2, etc
d<- data.frame(x=c(1,2,3,4))
y= matrix(c(1,2,3, 11,12,13, 20,21,22, 23,24,25) ,nrow=4,ncol=3)
d$y = y
d$z= y
#d now looks like my dataframe

#non-vectorized
d$l <- NaN
for (i in 1:nrow(d)) {
  d[i,]$l <- list(d$y[i,])
}

#vectorized https://stackoverflow.com/questions/51622017/change-all-array-columns-of-a-data-frame-into-lists
library(purrr)
e<- map_if(d,is.matrix,~split(.,seq(nrow(.)))) %>% as_tibble
#map_if(d,is.matrix,~split(.,seq(nrow(.)))) %>% as_tibble %>% print.data.frame

f<- e %>% mutate(a=x)  #works

###
d<- data.frame(x=c(1,2,3,4))
y= matrix(c(1,2,3, 11,12,13, 20,21,22, 23,24,25) ,nrow=4,ncol=3)
d$y = y
d$z= y
e<- map_if(d,is.matrix,~split(.,seq(nrow(.)))) 


#turn the weird array things into a list column
#As explained here, tibbles can have list-columns https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html
# Here is how to build one:
tt<-tibble(x=1:3)
tt$y=list(1:5,1:10,1:20)
tt$z = NaN
tt$z[1] <- list((letterSeq[1,]))