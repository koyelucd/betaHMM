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
                       R=c(-3, -2, -1, 0, 1, 2, 3)))
  expect_error(betaHMM(obj1, obj2, M=c(0, 1, 2, 3),
                       N=c(0, 1, 2, 3),R=c(0, 1, 2, 3)))
  expect_error(betaHMM(obj1, obj2, M=c(1.5, 3.3, 5.2),  N=c(1.5, 3.3, 5.2),
                       R=c(1.5, 3.3, 5.2)))
})

obj<-new("betaHMMResults")
test_that("AUC_threshold and uncertainty are scalar positive floating point
          number and within 0 and 1", {
    expect_error(dmc_identification(obj, AUC_threshold = 1.2,
                                    uncertainty=0.8))
              expect_error(dmc_identification(obj, AUC_threshold = 0.2,
                                              uncertainty=1.6)
                           )
              expect_error(dmc_identification(obj, AUC_threshold = -0.2,
                                              uncertainty=-0.8))

})

obj3<-new("dmcResults")
test_that("DMC_count is a scalar positive integer greater than 1.", {
              expect_error(dmr_identification(obj3, DMC_count = 0))
              expect_error(dmr_identification(obj3, DMC_count=0.2)
              )
              expect_error(dmr_identification(obj3, DMC_count=-2))
              expect_error(dmr_identification(obj3, DMC_count=1))
          })

test_that("M and N are scalar positive integers", {
    expect_error(threshold_identification(obj1[,1:5], M=-3,N =-3))
    expect_error(threshold_identification(obj1[,1:5],M=matrix(rnorm(6),nrow=3),
                         N=matrix(rnorm(6), nrow=3)))
    expect_error(threshold_identification(obj1[,1:5], M=c(-3, -2,-1, 0, 1,2,3),
                         N=c(-3, -2, -1, 0, 1, 2, 3)))
    expect_error(threshold_identification(obj1[,1:5], M=c(0, 1, 2, 3),
                         N=c(0, 1, 2, 3)))
    expect_error(threshold_identification(obj1[,1:5], M=c(1.5, 3.3, 5.2),
                                          N=c(1.5, 3.3, 5.2)))
})

test_that("Uncertainty is a positive floating point number lying between
          0 and 1",{
    expect_error(plot(obj,chromosome="1",what="fitted density",
                      uncertainty_threshold=-2))
    expect_error(plot(obj,chromosome="1",what="fitted density",
                      uncertainty_threshold=2))
    expect_error(plot(obj,chromosome="1",what="fitted density",
                      uncertainty_threshold=c(-3,-2)))
    expect_error(plot(obj,chromosome="1",what="fitted density",
                      uncertainty_threshold=c(1.5,2)))
    expect_error(plot(obj,chromosome="1",what="fitted density",
                      uncertainty_threshold=matrix(rnorm(6),nrow=3)))

})
test_that("what is a scalar string.",
          {expect_error(plot(obj,chromosome="1",what=1))
              expect_error(plot(obj,chromosome="1",what=c("fitted",
                                                                "kernel")))
              expect_error(plot(obj,chromosome="1",what=matrix(rnorm(6),
                                                                     nrow=3)))
          })

test_that("start_CpG is scalar and of type character and
          end_CpG is either scalar numeric or scalar and character type",{
              expect_error(plot(obj3,start_CpG=c("cg1","cg2"),
                                end_CpG=2))
              expect_error(plot(obj3,start_CpG="cg1",
                                end_CpG=c(1,10)))
              expect_error(plot(obj3,start_CpG="cg1",
                                end_CpG=NULL))
              expect_error(plot(obj3,start_CpG=NULL,
                                end_CpG=2))
              expect_error(plot(obj3,start_CpG=1,
                                end_CpG=2))

          })

obj4<-new("threshold_Results")

test_that("plot_threshold is carrying TRUE/FALSE value.",
          {expect_error(plot(obj4,plot_threshold="ABC"))
          expect_error(plot(obj4,plot_threshold=1))
          expect_error(plot(obj4,plot_threshold=matrix(rnorm(6),nrow=3)))
          expect_error(plot(obj4,plot_threshold=c(TRUE,FALSE,TRUE)))})
test_that("what is a scalar string.",
          {expect_error(plot(obj4,plot_threshold=FALSE,what=1))
              expect_error(plot(obj4,plot_threshold=TRUE,what=c("fitted",
                                                                "kernel")))
              expect_error(plot(obj4,plot_threshold=TRUE,what=matrix(rnorm(6),
                                                                     nrow=3)))
              })
