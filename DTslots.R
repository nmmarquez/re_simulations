rm(list=ls())
library(data.table)

# create new class with a single slot that is a data table
objDT <- setClass("objDT", slots=c(x="data.table"))

# assign new variable
test <- objDT(x=as.data.table(sleep))
test@x[,V1:=rnorm(nrow(test@x))]
head(test@x)

# save to a temp file
tf1 <- tempfile()
save(test, file = tf1)

# remove everything except temp file reference
rm(list=setdiff(ls(), "tf1"))

# load and try adding a column to the data table
load(tf1)
test@x[,V2:=rnorm(nrow(test@x))]
# no v2 here
head(test@x)

# try again using DT's fix
rm(list=setdiff(ls(), "tf1"))
load(tf1)
truelength(test@x)
setDT(test@x)
truelength(test@x) == 0
# still no luck
test@x[,V2:=rnorm(nrow(test@x))]
head(test@x)

# this works
test@x$V2 <- rnorm(nrow(test@x))
head(test@x)

# so does this
test@x <- as.data.table(test@x)
test@x[,V3:=rnorm(nrow(test@x))]
head(test@x)
