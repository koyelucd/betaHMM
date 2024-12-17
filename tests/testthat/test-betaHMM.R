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

test_that("methylation and annotation data are in dataframe/matrix/
          RangedSummarizedExperiment/ GRanges object", {
    expect_error(betaHMM(methylation_data=list(a=1,b=2),
                         annotation_file=list(a=1,b=2),
                         M=3,N =4,R=2))
      expect_error(betaHMM(methylation_data=seq(1,10),
                           annotation_file=seq(1,10),
                           M=3,N =4,R=2))
      expect_error(betaHMM(methylation_data=rep("a",10),
                           annotation_file=rep("a",10),
                           M=3,N =4,R=2))
})

test_that("treatment_group is a vector of string values",{
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         treatment_group=c(1,1)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         treatment_group=1))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         treatment_group=matrix(rnorm(10),ncol=5)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         treatment_group=list(s=1,j=2)))
})

test_that("parallel_process takes a TRUE/FALSE value",{
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         parallel_process=c(1,2)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         parallel_process=1))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         parallel_process=list(a=1,b=1)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         parallel_process=matrix(rnorm(6),ncol=2)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         parallel_process="abc"))
})


test_that("seed takes a numerical scalar value",{
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         seed="abc"))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         seed=c(1,10)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         seed=list(a=1,b=1)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         seed=matrix(rnorm(6),ncol=2)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         seed=TRUE))
})

test_that("iterations takes a numerical scalar value",{
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         iterations="abc"))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         iterations=c(1,10)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         iterations=list(a=1,b=1)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         iterations=matrix(rnorm(6),ncol=2)))
    expect_error(betaHMM(obj1,obj2,M=3,N=4,R=2,
                         iterations=TRUE))
})

test_that("M, N and R are scalar positive integers", {
    expect_error(betaHMMrun(obj1, obj2, M=-3,N =-3,R=-3,
                 treatment_group = NULL, parallel_process = FALSE,
                 seed = NULL,iterations=100))
    expect_error(betaHMMrun(obj1, obj2, M=matrix(rnorm(6), nrow=3),
                         N=matrix(rnorm(6), nrow=3),R=matrix(rnorm(6), nrow=3),
                         treatment_group = NULL, parallel_process = FALSE,
                         seed = NULL,iterations=100))
    expect_error(betaHMMrun(obj1, obj2, M=c(-3, -2, -1, 0, 1, 2, 3),
                         N=c(-3, -2, -1, 0, 1, 2, 3),
                         R=c(-3, -2, -1, 0, 1, 2, 3),
                         treatment_group = NULL, parallel_process = FALSE,
                         seed = NULL,iterations=100))
    expect_error(betaHMMrun(obj1, obj2, M=c(0, 1, 2, 3),
                         N=c(0, 1, 2, 3),R=c(0, 1, 2, 3),
                         treatment_group = NULL, parallel_process = FALSE,
                         seed = NULL,iterations=100))
    expect_error(betaHMMrun(obj1, obj2, M=c(1.5, 3.3, 5.2),N=c(1.5, 3.3, 5.2),
                         R=c(1.5, 3.3, 5.2),
                         treatment_group = NULL, parallel_process = FALSE,
                         seed = NULL,iterations=100))
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

test_that("AUC_threshold and uncertainty are scalar positive floating point
          number and within 0 and 1", {
              expect_error(dmc_identification_run(obj, AUC_threshold = 1.2,
                                              uncertainty=0.8))
              expect_error(dmc_identification_run(obj, AUC_threshold = 0.2,
                                              uncertainty=1.6)
              )
              expect_error(dmc_identification_run(obj, AUC_threshold = -0.2,
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

obj6<-matrix(rnorm(100),ncol=10)
test_that("DMC_count is a scalar positive integer greater than 1.", {
    expect_error(dmr_identification_run(obj6, DMC_count = 0))
    expect_error(dmr_identification_run(obj6, DMC_count=0.2)
    )
    expect_error(dmr_identification_run(obj6, DMC_count=-2))
    expect_error(dmr_identification_run(obj6, DMC_count=1))
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

test_that("parameter_estimation_only has a TRUE/FALSE input",
          {
          expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                parameter_estimation_only=
                                                    "abc"))
          expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                parameter_estimation_only=
                                                    10))
          expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                parameter_estimation_only=
                                                    c(10,10)))
          })

test_that("package_workflow has a TRUE/FALSE input",
          {
              expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                    package_workflow = "abc"))
              expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                    package_workflow = 10))
              expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                    package_workflow =
                                                        c(10,10)))
          })

test_that("seed has a positive scalar input",
          {
              expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                    seed  = "abc"))
              expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                    seed = -10))
              expect_error(threshold_identification(obj1[,1:5], M=3,N =4,
                                                    seed = c(10,10)))
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

test_that("start_CpG is scalar and of type character",{
              expect_error(plot(obj3,start_CpG=c("cg1","cg2"),
                                end_CpG=2))
              expect_error(plot(obj3,start_CpG=NULL,
                                end_CpG=2))
              expect_error(plot(obj3,start_CpG=TRUE,
                                end_CpG=2))
              expect_error(plot(obj3,start_CpG=1,
                                end_CpG=2))

          })

test_that("end_CpG is either scalar numeric or scalar and character type",{
              expect_error(plot(obj3,start_CpG="cg1",
                                end_CpG=TRUE))
              expect_error(plot(obj3,start_CpG="cg1",
                                end_CpG=c(1,10)))
              expect_error(plot(obj3,start_CpG="cg1",
                                end_CpG=NULL))
              expect_error(plot(obj3,start_CpG="cg1",
                                end_CpG=matrix(rnorm(6),nrow=2)))

          })

test_that("Extra input parameters are not passed when plotting DMR
          detection plot",{
    expect_error(plot(obj3,start_CpG="cg1",
                      end_CpG=2,what="fitted density"))
    expect_error(plot(obj3,start_CpG="cg1",
                      end_CpG=2,plot_threshold=TRUE))

})

test_that("Title is a string",{
    expect_error(plot(obj3,start_CpG="cg1",
                      end_CpG=2,title=1))
    expect_error(plot(obj3,start_CpG="cg1",
                      end_CpG=2,title=TRUE))

})

test_that("dmcResults object passed when accessing plotDMR function",{
              expect_error(plot(obj1,start_CpG="cg1",
                                end_CpG=2))
              expect_error(plot(matrix(rnorm(6),nrow=2),start_CpG="cg1",
                                end_CpG=c2))
              expect_error(plot(c(1,1),start_CpG="cg1",
                                end_CpG=2))

          })

test_that("title is correct in plotDMR function",{
    expect_error(plot(obj3,start_CpG="cg1",
                      end_CpG=2,title=1))
    expect_error(plot(obj3,start_CpG="cg1",
                      end_CpG=c2,title=c(1,1)))

})
obj4<-new("threshold_Results")

test_that("plot_threshold is carrying TRUE/FALSE value.",
          {expect_error(plot(obj4,plot_threshold="ABC"))
          expect_error(plot(obj4,plot_threshold=1))
          expect_error(plot(obj4,plot_threshold=matrix(rnorm(6),nrow=3)))
          expect_error(plot(obj4,plot_threshold=c(TRUE,FALSE,TRUE)))})

test_that("For plotting threshold the correct object is passed as input",
          {expect_error(plot(obj1,plot_threshold=TRUE))
              expect_error(plot(obj3,plot_threshold=TRUE))
             })

test_that("what is a scalar string.",
          {expect_error(plot(obj4,plot_threshold=FALSE,what=1))
              expect_error(plot(obj4,plot_threshold=TRUE,what=c("fitted",
                                                                "kernel")))
              expect_error(plot(obj4,plot_threshold=TRUE,what=matrix(rnorm(6),
                                                                     nrow=3)))
              })

test_that("AUC_DM_analysis has positive values as input for M, N, R and K",
          {expect_error(AUC_DM_analysis(M=-3,N=-3, R=-3, K=-3,
                                        tau=c(0.1,0.3,0.6), A=matrix(rnorm(81),
                                                                     nrow=9),
                                        phi=list(sp_1=c(1,1,1),
                                                 sp_2=c(1,1,1))))})

test_that("AUC_DM_analysis has positive values as input for M, N, R and K",
          {expect_error(AUC_DM_analysis(M=-3,N=-3, R=-3, K=-3,
                                        tau=c(0.1,0.3,0.6), A=matrix(rnorm(81),
                                                                     nrow=9),
                                        phi=list(sp_1=c(1,1,1),
                                                 sp_2=c(1,1,1))))})

test_that("AUC_DM_analysis has scalar values for M",
          {expect_error(AUC_DM_analysis(M=c(1,1,1),N=c(4,4), R=2,
                                        K=9,
                                        tau=c(0.1,0.3,0.6), A=matrix(rnorm(81),
                                                                     nrow=9),
                                        phi=list(sp_1=c(1,1,1),
                                                 sp_2=c(1,1,1))))
              expect_error(AUC_DM_analysis(M=matrix(rnorm(4),nrow=2),N=c(4,4),
                                           R=2,
                                           K=9,
                                           tau=c(0.1,0.3,0.6),
                                           A=matrix(rnorm(81),nrow=9),
                                           phi=list(sp_1=c(1,1,1),
                                                    sp_2=c(1,1,1))))})

test_that("AUC_DM_analysis has scalar values for R",
          {expect_error(AUC_DM_analysis(M=3,N=c(4,4), R=c(1,1,1),
                                        K=9,
                                        tau=c(0.1,0.3,0.6), A=matrix(rnorm(81),
                                                                     nrow=9),
                                        phi=list(sp_1=c(1,1,1),
                                                 sp_2=c(1,1,1))))
              expect_error(AUC_DM_analysis(M=3,N=c(4,4),
                                           R=matrix(rnorm(4),nrow=2),
                                           K=9,
                                           tau=c(0.1,0.3,0.6),
                                           A=matrix(rnorm(81),nrow=9),
                                           phi=list(sp_1=c(1,1,1),
                                                    sp_2=c(1,1,1))))})
test_that("AUC_DM_analysis has scalar values for K",
          {expect_error(AUC_DM_analysis(M=3,N=c(4,4), R=2,
                                        K=c(1,1,1),
                                        tau=c(0.1,0.3,0.6), A=matrix(rnorm(81),
                                                                     nrow=9),
                                        phi=list(sp_1=c(1,1,1),
                                                 sp_2=c(1,1,1))))
              expect_error(AUC_DM_analysis(M=3,N=c(4,4),
                                           R=2,
                                           K=matrix(rnorm(4),nrow=2),
                                           tau=c(0.1,0.3,0.6),
                                           A=matrix(rnorm(81),nrow=9),
                                           phi=list(sp_1=c(1,1,1),
                                                    sp_2=c(1,1,1))))})
test_that("AUC_DM_analysis has a vector input for tau",
          {expect_error(AUC_DM_analysis(M=3,N=c(4,4), R=2,
                                        K=9,
                                        tau=0.5, A=matrix(rnorm(81),
                                                                     nrow=9),
                                        phi=list(sp_1=c(1,1,1),
                                                 sp_2=c(1,1,1))))
              expect_error(AUC_DM_analysis(M=3,N=c(4,4),
                                           R=2,
                                           K=9,
                                           tau=matrix(rnorm(4),nrow=2),
                                           A=matrix(rnorm(81),nrow=9),
                                           phi=list(sp_1=c(1,1,1),
                                                    sp_2=c(1,1,1))))})
test_that("AUC_DM_analysis has a matrix input for A",
          {expect_error(AUC_DM_analysis(M=3,N=c(4,4), R=2,
                                        K=9,
                                        tau=c(0.1,0.2,0.7), A=c(1,2,3),
                                        phi=list(sp_1=c(1,1,1),
                                                 sp_2=c(1,1,1))))
              expect_error(AUC_DM_analysis(M=3,N=c(4,4),
                                           R=2,
                                           K=9,
                                           tau=c(0.1,0.2,0.7),
                                           A=3,
                                           phi=list(sp_1=c(1,1,1),
                                                    sp_2=c(1,1,1))))})
test_that("AUC_DM_analysis has a list input for phi",
          {expect_error(AUC_DM_analysis(M=3,N=c(4,4), R=2,
                                        K=9,
                                        tau=c(0.1,0.2,0.7), A=matrix(rnorm(81),
                                                                     nrow = 9),
                                        phi=c(1,1)))
              expect_error(AUC_DM_analysis(M=3,N=c(4,4),
                                           R=2,
                                           K=9,
                                           tau=c(0.1,0.2,0.7),
                                           A=matrix(rnorm(81),nrow = 9),
                                           phi=3))})

test_that("backward function has a matrix input for probabilities",
          {expect_error(backward(probabilities=c(1,1,1),
                                trained_params = list(A = matrix(rnorm(81),
                                                                 nrow=9),
                                                      tau = c(0.1,0.2,0.7),
                                                      phi = list(sp_1=c(1,1,1),
                                                                 sp_2=c(1,1,1)
                                                                 ))))})
test_that("backward function has a list input for trained_params",
          {expect_error(backward(probabilities=matrix(rnorm(81),ncol=9),
                                 trained_params = c(1,1,1) ))
              expect_error(backward(probabilities=matrix(rnorm(81),ncol=9),
                                    trained_params=matrix(rnorm(81),ncol=9)))})

test_that("forward function has a matrix input for probabilities",
          {expect_error(forward(probabilities=c(1,1,1),
                                 trained_params = list(A = matrix(rnorm(81),
                                                                  nrow=9),
                                                       tau = c(0.1,0.2,0.7),
                                                       phi=list(sp_1=c(1,1,1),
                                                                sp_2=c(1,1,1)
                                                       ))))})
test_that("forward function has a list input for trained_params",
          {expect_error(forward(probabilities=matrix(rnorm(81),ncol=9),
                                 trained_params = c(1,1,1) ))
              expect_error(forward(probabilities=matrix(rnorm(81),ncol=9),
                                    trained_params=matrix(rnorm(81),ncol=9)))})

test_that("initialise_parameters has a positive scalar input for M parameter",
          {expect_error(initialise_parameters(data=obj1, M=-3, N=c(4,4), R=2,
                                              seed = NULL))
              expect_error(initialise_parameters(data=obj1, M=c(1,1), N=c(4,4),
                                                 R=2,
                                                 seed = NULL))})

test_that("initialise_parameters has a positive scalar input for R parameter",
          {expect_error(initialise_parameters(data=obj1, M=3, N=c(4,4), R=-2,
                                              seed = NULL))
              expect_error(initialise_parameters(data=obj1, M=3, N=c(4,4),
                                                 R=c(1,1),
                                                 seed = NULL))})

test_that("initialise_parameters has a positive vector input for N parameter",
          {expect_error(initialise_parameters(data=obj1, M=3, N=c(-4,-4), R=2,
                                              seed = NULL))
              expect_error(initialise_parameters(data=obj1, M=3, N=c("a","a"),
                                                 R=2,
                                                 seed = NULL))})

obj1<-as.data.frame(obj1)
test_that("initialise_parameters has a dataframe input for data parameter",
          {expect_error(initialise_parameters(data=c(1,1), M=3, N=c(4,4), R=2,
                                              seed = NULL))
              expect_error(initialise_parameters(data=1, M=3, N=c("a","a"),
                                                 R=2,
                                                 seed = NULL))
              expect_error(initialise_parameters(data=as.matrix(obj1), M=3,
                                                 N=c(4,4),
                                                 R=2,
                                                 seed = NULL))})


test_that("initialise_parameters has a positive scalar input for seeed
          parameter",
          {expect_error(initialise_parameters(data=obj1, M=3, N=c(4,4), R=2,
                                              seed = -100))
              expect_error(initialise_parameters(data=obj1, M=3, N=c("a","a"),
                                                 R=2,
                                                 seed = "abc"))
              expect_error(initialise_parameters(data=obj1, M=3, N=c("a","a"),
                                                 R=2,
                                                 seed = c(1,1)))})

test_that("Viterbi has positive scalar values for M, R and K",
          {expect_error(Viterbi(data=obj1[,-1], M=-3, N=c(4,4), R=-2,
                                tau=c(0.1,0.2,0.3),
                                A=matrix(rnorm(81),nrow=9),
                                phi=list(sp_1=c(1,1,1),
                                         sp_2=c(1,1,1)),K=-9))
              expect_error(Viterbi(data=obj1[,-1], M=c(1,1), N=c(4,4),R=c(1,1),
                                   tau=c(0.1,0.2,0.3),
                                   A=matrix(rnorm(81),nrow=9),
                                   phi=list(sp_1=c(1,1,1),
                                            sp_2=c(1,1,1)),K=c(1,1)))})

test_that("Viterbi has positive vector input for tau",
          {expect_error(Viterbi(data=obj1[,-1], M=3, N=c(4,4), R=2,
                                tau=c(-3,-3),
                                A=matrix(rnorm(81),nrow=9),
                                phi=list(sp_1=c(1,1,1),
                                         sp_2=c(1,1,1)),K=9))
              expect_error(Viterbi(data=obj1[,-1], M=3, N=c(4,4), R=2,
                                   tau=-3,
                                   A=matrix(rnorm(81),nrow=9),
                                   phi=list(sp_1=c(1,1,1),
                                            sp_2=c(1,1,1)),K=9))
              expect_error(Viterbi(data=obj1[,-1], M=3, N=c(4,4), R=2,
                                   tau=list(),
                                   A=matrix(rnorm(81),nrow=9),
                                   phi=list(sp_1=c(1,1,1),
                                            sp_2=c(1,1,1)),K=9))})

test_that("The phi has a list input for Viterbi function",
          {expect_error(Viterbi(data=obj1[,-1], M=3, N=c(4,4), R=2,
                                tau=c(0.1,0.2,0.7),
                                A=matrix(rnorm(81),nrow=9),
                                phi=c(1,1),K=9))
              expect_error(Viterbi(data=obj1[,-1], M=3, N=c(4,4), R=2,
                                   tau=c(0.1,0.2,0.7),
                                   A=matrix(rnorm(81),nrow=9),
                                   phi=1,K=9))
              expect_error(Viterbi(data=obj1[,-1], M=3, N=c(4,4), R=2,
                                   tau=c(0.1,0.2,0.7),
                                   A=matrix(rnorm(81),nrow=9),
                                   phi=matrix(rnorm(6),ncol=2),K=9))})

test_that("BaumWelch function has positive scalar values for M,K and R",
          {expect_error(BaumWelch(data=obj, trained_params=NULL,K=-3,
                                  M=-3, N=c(4,4), R=-3,
            seed = NULL,iterations= 100))
              expect_error(BaumWelch(data=obj, trained_params=NULL,K=c(1,1),
                                     M=c(1,1), N=c(4,4), R=c(1,1),
                                     seed = NULL,iterations= 100))
          })

test_that("BaumWelch function has positive scalar or vector
          integer input for N",
          {expect_error(BaumWelch(data=obj, trained_params=NULL,K=9,
                                  M=3, N=c("a","b"), R=2,
                                  seed = NULL,iterations= 100))
              expect_error(BaumWelch(data=obj, trained_params=NULL,K=9,
                                     M=3, N=c(0.4,0.2), R=2,
                                     seed = NULL,iterations= 100))
              expect_error((BaumWelch(data=obj, trained_params=NULL,K=9,
                                      M=3, N=0.2, R=2,
                                      seed = NULL,iterations= 100)))
              expect_error((BaumWelch(data=obj, trained_params=NULL,K=9,
                                      M=3,
                                      N=matrix(rnorm(6),ncol=2), R=2,
                                      seed = NULL,iterations= 100)))
              expect_error((BaumWelch(data=obj, trained_params=NULL,K=9,
                                      M=3, N=list(a=1,b=1), R=2,
                                      seed = NULL,iterations= 100)))
          })

test_that("BaumWelch function has a dataframe as input for data variable",
          {expect_error(BaumWelch(data=3, trained_params=NULL,K=9,
                                  M=3, N=c(4,4), R=2,
                                  seed = NULL,iterations= 100))
              expect_error(BaumWelch(data=c(1,1), trained_params=NULL,K=9,
                                     M=3, N=c(4,4), R=2,
                                     seed = NULL,iterations= 100))
          })
test_that("BaumWelch function has a scalar positive value as input for seed
            and iterations",
          {expect_error(BaumWelch(data=obj, trained_params=NULL,K=9,
                                  M=3, N=c(4,4), R=2,
                                  seed = c(1,1) ,iterations= c(1,1)))
              expect_error(BaumWelch(data=obj, trained_params=NULL,K=9,
                                     M=3, N=c(4,4), R=2,
                                     seed = -3,iterations= -3))
          })

test_that("initialise_parameters_th has a positive scalar input for
          M parameter",
          {expect_error(initialise_parameters_th(data=obj1, M=-3, N=c(4,4),
                                                 R=2,
                                              seed = NULL))
              expect_error(initialise_parameters_th(data=obj1, M=c(1,1),
                                                    N=c(4,4),
                                                 R=2,
                                                 seed = NULL))})

test_that("initialise_parameters_th has a positive scalar input for
          R parameter",
          {expect_error(initialise_parameters_th(data=obj1, M=3, N=c(4,4),
                                                 R=-2,
                                              seed = NULL))
              expect_error(initialise_parameters_th(data=obj1, M=3, N=c(4,4),
                                                 R=c(1,1),
                                                 seed = NULL))})

test_that("initialise_parameters_th has a positive vector input for
          N parameter",
          {expect_error(initialise_parameters_th(data=obj1, M=3, N=c(-4,-4),
                                                 R=2,
                                              seed = NULL))
              expect_error(initialise_parameters_th(data=obj1, M=3,
                                                    N=c("a","a"),
                                                 R=2,
                                                 seed = NULL))})

test_that("initialise_parameters_th has a positive scalar input for seeed
          parameter",
          {expect_error(initialise_parameters_th(data=obj1, M=3, N=c(4,4),
                                                 R=2,
                                              seed = -100))
              expect_error(initialise_parameters_th(data=obj1, M=3,
                                                    N=c("a","a"),
                                                 R=2,
                                                 seed = "abc"))
              expect_error(initialise_parameters_th(data=obj1, M=3,
                                                    N=c("a","a"),
                                                 R=2,
                                                 seed = c(1,1)))})


test_that("threshold function has a numerical vector input for tau",
          {expect_error(threshold_values(data=obj1,tau=1,
                                         phi=list(sp_1=c(1,1,1),sp_2=c(1,1,1)
                                                  )))
              expect_error(threshold_values(data=obj1,
                                            tau=matrix(rnorm(6),ncol=2),
                                            phi=list(sp_1=c(1,1,1),
                                                     sp_2=c(1,1,1)
                                            )))
              expect_error(threshold_values(data=obj1,
                                            tau=list(tau1=0.5,tau2=0.5),
                                            phi=list(sp_1=c(1,1,1),
                                                     sp_2=c(1,1,1)
                                            )))
              expect_error(threshold_values(data=obj1,
                                            tau=c("a","b"),
                                            phi=list(sp_1=c(1,1,1),
                                                     sp_2=c(1,1,1)
                                            )))
              })

test_that("threshold function has a list input for phi",
          {expect_error(threshold_values(data=obj1,tau=c(0.1,0.2,0.7),
                                         phi=1
                                         ))
              expect_error(threshold_values(data=obj1,tau=c(0.1,0.2,0.7),
                                            phi=c(1,2)
              ))
              expect_error(threshold_values(data=obj1,tau=c(0.1,0.2,0.7),
                                            phi=matrix(rnorm(6),ncol=2)
              ))
              expect_error(threshold_values(data=obj1,tau=c(0.1,0.2,0.7),
                                            phi=as.data.frame(
                                                matrix(rnorm(6),ncol=2))
              ))
          })



test_that("kernel density plot has correct input",{
    expect_error(kernel_density_plot(data=obj1, hidden_states=1,
                        auc=matrix(rnorm(6),ncol=3), C=c(1,10),
                         treatment_group=c("B","T"), title_text=NULL))
    expect_error(kernel_density_plot(data=obj1,
                                     hidden_states=matrix(rnorm(6),ncol=2),
                                     auc=matrix(rnorm(6),ncol=3), C=c(1,10),
                                     treatment_group=c("B","T"),
                                     title_text=NULL))
    expect_error(kernel_density_plot(data=obj1,
                                     hidden_states=list(a="1",b="2"),
                                     auc=matrix(rnorm(6),ncol=3), C=c(1,10),
                                     treatment_group=c("B","T"),
                                     title_text=NULL))
})

test_that("fitted density plot has correct input",{
    expect_error(fitted_density_plot(phi=c(1,10), hidden_states=1,
                                     auc=matrix(rnorm(6),ncol=3),R=c(2,2),
                                     K=c(2,2), C=c(1,10),
                                     treatment_group=c("B","T"),
                                     title_text=NULL))
    expect_error(fitted_density_plot(phi=1,
                                     hidden_states=matrix(rnorm(6),ncol=2),
                                     auc=matrix(rnorm(6),ncol=3),R=c(2,2),
                                     K=c(2,2), C=c(1,10),
                                     treatment_group=c("B","T"),
                                     title_text=NULL))
    expect_error(fitted_density_plot(phi=matrix(rnorm(6),ncol=2),
                                     hidden_states=1,
                                     auc=matrix(rnorm(6),ncol=3),
                                     R=matrix(rnorm(6),ncol=2),
                                     K=matrix(rnorm(6),ncol=2),
                                     C=matrix(rnorm(6),ncol=2),
                                     treatment_group=c("B","T"),
                                     title_text=NULL))
})


test_that("Uncertainty plot has correct form of input",{
    expect_error(uncertainty_plots(data_comp=obj2,
                                   z=seq(1,10), hidden_states=c(1,10),C=c(1,2),
                                   uncertainty_threshold=
                                       matrix(rnorm(6),ncol=2),
                                   title_text=NULL))
})

test_that("individual function for data processing has correct input",{
    expect_error(beta_value_numeric(data=1,anno_file=1))
    expect_error(beta_value_numeric(data=list(a=1,b=1),
                                    anno_file = list(a=1,b=1)))
})
