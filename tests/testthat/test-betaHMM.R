##-----------------------------------------------------------------------------
context("betaHMM function tests")

set.seed(12345)
obj <- matrix(round(runif(50, 1, 100)), ncol=8, nrow=10)
obj<-as.data.frame(obj)
x=seq(1,10)
IlmnID<-paste("cg",x,sep = "")
obj1<-data.frame(IlmnID=IlmnID,obj)
obj2<-data.frame(IlmnID=IlmnID,CHR=1,MAPINFO=x)

test_that("M, N and R are scalar positive integers", {
  expect_error(betaHMM(obj1, obj2, M=-3,N =-3,R=-3))
  expect_error(betaHMM(obj1, obj2, M=matrix(rnorm(6), nrow=3),
                       N=matrix(rnorm(6), nrow=3),R=matrix(rnorm(6), nrow=3)))
  expect_error(betaHMM(obj1, obj2, M=c(-3, -2, -1, 0, 1, 2, 3),
                       N=c(-3, -2, -1, 0, 1, 2, 3),
                       R=c(-3, -2, -1, 0, 1, 2, 3),))
  expect_error(betaHMM(obj1, obj2, M=c(0, 1, 2, 3),
                       N=c(0, 1, 2, 3),R=c(0, 1, 2, 3)))
  expect_error(betaHMM(obj1, obj2, M=c(1.5, 3.3, 5.2),  N=c(1.5, 3.3, 5.2),
                       R=c(1.5, 3.3, 5.2)))
})
