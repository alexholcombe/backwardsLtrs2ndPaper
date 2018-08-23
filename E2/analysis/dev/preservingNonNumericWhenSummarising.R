#Original idea from
#https://stackoverflow.com/a/45802176/302378

#My question
#https://stackoverflow.com/questions/51978519/summarise-numeric-columns-return-last-value-of-non-numeric/51980086#51980086

# dataset('iris')
# iris %>% group_by(Sepal.Length) %>% summarise_all(funs(if_else(is.numeric(.), mean(.), first(.))))
# iris %>% select(-Species) %>% group_by(Sepal.Length) %>% 
#   summarise_all(funs(if_else(is.numeric(.), mean(.), first(.))))


set.seed(1234)
category <- (c('A','A','E','E','B','B','C'))
date <- seq(as.Date("2017-01-01"), by = "month", length.out = 7)
value1 <- sample(seq(from = 91, to = 97, by = 1))
dt <- data.frame(category, date, value1)
dt<- as_tibble(dt)
#works
dt2<- dt %>% select(-marsupial) %>%
  group_by(category) %>%
  summarise_all(funs(if_else(is.numeric(.), mean(.), last(.))))
print(dt2)

marsupial <- c("quoll","phascogale",'triok','opossum','antechinus','bandicoot','Fat-tailed dunnart')
dt$marsupial <- marsupial
dt3<- dt %>% #doesn't work
  group_by(category) %>%
  summarise_all(funs(if_else(is.numeric(.), mean(.), last(.))))

#WORKS except fucks up date, turns it into double. Also doesn't preserve columns that are lists
dt4<- dt %>% group_by(category) %>%
  summarise_all(function(x){ifelse(is.numeric(x),mean(x),last(x))}) 

#Below works and the author said works on list columns, but doesn't, gives errors like
# Error in summarise_impl(.data, dots) : 
# Column `targetStreamOrientations` must be length 1 (a summary value), not 16
dt5<- dt %>% group_by(category) %>%
  summarise_all( function(x){ if(is.numeric(x)){ return(mean(x)) } else{ nth(x,-1) } })
#Below gives error related to date
dt5<- dt %>% group_by(category) %>%
  summarise_all(function(x){ if (is.numeric(.)) mean(.) else last(.) })

dt6<- dt %>% group_by(category) %>%
  summarise_if(is.numeric,mean)

y = if (3==5) 4 else 5
