#' @title Generate the raw ROC data after training a machine learning model using rpart model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     rpart = ROC_rpart(task,data)
#'     rpart
#'     }
ROC_rpart = function(task,data,
                     split_ratio = 0.7,
                     inner_resamples = 5,
                     outer_resamples =5,
                     tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner_rpart = lrn("classif.rpart",
                      predict_type = "prob",
                      id ="Classification_Tree")
  search_space_rpart = ps(cp = p_dbl(lower = 0.001, upper = 0.1))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at = at_create(tuner,learner_rpart,search_space_rpart,inner_resample,measure)
  rr = resample(task, at,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  at$train(task = task,split$train)
  pred = at$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}


#' @title Generate the raw ROC data after training a machine learning model using knn model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     knn = ROC_knn(task,data)
#'     knn
#'     }
ROC_knn = function(task,data,
                   split_ratio = 0.7,
                   inner_resamples = 5,
                   outer_resamples =5,
                   tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner = lrn("classif.kknn",
                predict_type = "prob",
                id="KNN")
  search_space = ps(k = p_int(lower = 1, upper = 10))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at = at_create(tuner,learner,search_space,inner_resample,measure)
  rr = resample(task, at,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  at$train(task = task,split$train)
  pred = at$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}


#' @title Generate the raw ROC data after training a machine learning model using svm model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     svm = ROC_svm(task,data)
#'     svm
#'     }
ROC_svm = function(task,data,
                   split_ratio = 0.7,
                   inner_resamples = 5,
                   outer_resamples =5,
                   tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner = lrn("classif.svm",
                type = "C-classification",
                kernel = "radial",
                predict_type = "prob",
                id = "SVM")
  search_space = ps(cost = p_dbl(log(0.1), log(10),
                                 trafo = function(x) exp(x)),
                    gamma = p_dbl(log(0.1), log(10),
                                  trafo = function(x) exp(x)))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at = at_create(tuner,learner,search_space,inner_resample,measure)
  rr = resample(task, at,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  at$train(task = task,split$train)
  pred = at$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}


#' @title Generate the raw ROC data after training a machine learning model using xgboost model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     xgboost = ROC_xgboost(task,data)
#'     xgboost
#'     }
ROC_xgboost = function(task,data,
                       split_ratio = 0.7,
                       inner_resamples = 5,
                       outer_resamples =5,
                       tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner = lrn("classif.xgboost",
                predict_type = "prob",
                id ="xgboost")
  search_space = ps(
    eta = p_dbl(lower = 0.1, upper = 1),
    max_depth = p_int(lower = 1, upper = 10),
    nrounds = p_int(lower = 1, upper = 16))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at = at_create(tuner,learner,search_space,inner_resample,measure)
  rr = resample(task, at,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  at$train(task = task,split$train)
  pred = at$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}


#' @title Generate the raw ROC data after training a machine learning model using random_forest model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     rr = ROC_rr(task,data)
#'     rr
#'     }
ROC_rr = function(task,data,
                  split_ratio = 0.7,
                  inner_resamples = 5,
                  outer_resamples =5,
                  tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner = lrn("classif.ranger",
                predict_type = "prob",
                id ="Random_Forest")
  search_space = ps(
    mtry = p_int(lower = 1, upper = task$ncol-1),
    min.node.size = p_int(lower = 1, upper = 10),
    num.trees = p_int(lower = 1, upper = 500))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at = at_create(tuner,learner,search_space,inner_resample,measure)
  rr = resample(task, at,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  at$train(task = task,split$train)
  pred = at$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}

#' @title Generate the raw ROC data after training a machine learning model using glmnet model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3verse
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     glmnrt = ROC_glmnet(task,data)
#'     glmnrt
#'     }
ROC_glmnet = function(task,data,
                      split_ratio = 0.7,
                      inner_resamples = 5,
                      outer_resamples =5,
                      tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner = lrn("classif.glmnet",
                predict_type = "prob",
                id ="Glmnet")
  search_space = ps(
    lambda = p_dbl(lower = 0.001, upper = 1, logscale = TRUE),
    alpha = p_dbl(lower = 0, upper = 1))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at = at_create(tuner,learner,search_space,inner_resample,measure)
  rr = resample(task, at,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  at$train(task = task,split$train)
  pred = at$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}


#' @title Generate the raw ROC data after training a machine learning model using LDA model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     lda = ROC_lda(task,data)
#'     lda
#'     }
ROC_lda = function(task,data,
                   split_ratio = 0.7,
                   inner_resamples = 5,
                   outer_resamples =5,
                   tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner = lrn("classif.lda", predict_type = "prob")

  rr  = resample(task, learner, outer_resample)
  learner$train(task = task,split$train)
  pred = learner$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}


#' @title Generate the raw ROC data after training a machine learning model using Logreg model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     logreg = ROC_logreg(task,data)
#'     logreg
#'     }
ROC_logreg = function(task,data,
                      split_ratio = 0.7,
                      inner_resamples = 5,
                      outer_resamples =5,
                      tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  #定义学习器
  learner = lrn("classif.log_reg", predict_type = "prob")
  # 进行交叉验证
  rr = resample(task, learner, outer_resample)
  learner$train(task = task,split$train)
  pred = learner$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}


#' @title Generate the raw ROC data after training a machine learning model using Naive_Bayes model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     naive_bayes =  ROC_naive_bayes(task,data)
#'     naive_bayes
#'     }
ROC_naive_bayes = function(task,data,
                           split_ratio = 0.7,
                           inner_resamples = 5,
                           outer_resamples =5,
                           tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner= lrn("classif.naive_bayes", predict_type = "prob")
  search_space = ps(
    laplace = p_dbl(lower = 0.1, upper = 1),
    threshold = p_int(lower = 1, upper = 10))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at = at_create(tuner,learner,search_space,inner_resample,measure)
  rr = resample(task, at,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  at$train(task = task,split$train)
  pred = at$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}

#' @title Generate the raw ROC data after training a machine learning model using Nnet model
#'      and making predictions.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param split_ratio A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resamples \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return If the data contained any NA value and na.rm=FALSE, NA is returned.
#'      Otherwise, if smooth=FALSE, a list of class “roc” with the following fields:\code{\link[pROC]{roc}}
#' @export
#' @import mlr3
#' @import mlr3tuning
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import pROC
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     nnet = ROC_nnet(task,data)
#'     nnet
#'     }
ROC_nnet = function(task,data,
                    split_ratio = 0.7,
                    inner_resamples = 5,
                    outer_resamples =5,
                    tuner_resolution =5){
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  split = partition(task, ratio = split_ratio)

  learner = lrn("classif.nnet",predict_type = "prob")
  search_space = ps(
    size = p_int(lower = 1, upper = 10),
    decay = p_dbl(lower = 0.1, upper = 0.9),
    MaxNWts= p_int(lower = 10000, upper = 10000))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at = at_create(tuner,learner,search_space,inner_resample,measure)
  rr = resample(task, at,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  at$train(task = task,split$train)
  pred = at$predict(task,split$test)
  roc <- roc(pred$truth %>% as.numeric(), pred$response %>% as.numeric())

  roc
}


#' @title Plot the ROC curves using the ROC results.
#'
#' @param roc \code{\link[pROC]{roc}} The results of the functions, ROC_rpart, ROC_knn, ROC_svm, ROC_xgboost, ROC_rr,
#'     ROC_glmnet, ROC_lda, ROC_logreg, ROC_naive_bayes, and ROC_nnet, with class "roc".
#'
#' @return ROC curve
#' @export
#' @import pROC
#' @import ggsci
#' @import ggplot2
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test[c(1:400),c(1:9)])
#'     task = task_creat(data)
#'     rpart = ROC_rpart(task,data)
#'     roc_plot(rpart)
#'     }
roc_plot <- function(roc) {
  plot =  ggroc(roc,legacy.axes = TRUE,if(class(roc) == "list"){aes = c("color")}else{NULL}) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="red",
                 linetype=6)+  theme_bw()+
    if(class(roc) == "list"){scale_color_npg()}else{NULL}
  plot
}


#' @title ROC curves across 10 different machine learning models.
#'
#' @param gene_exp For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param multi_threaded A logic variable to enable multi-threading during execution.
#' @param split_ratios A numeric variable representing the ratio used to split the data into a testing set and a training set;
#'      it indicates the proportion of the total data that will be allocated to the testing set.
#' @param inner_resampless \code{\link[mlr3]{rsmp}}
#' @param outer_resampless \code{\link[mlr3]{rsmp}}
#' @param tuner_resolutions \code{\link[mlr3tuning]{tnr}}
#'
#' @return ROC curve
#' @export
#' @import parallel
#' @import parallelly
#' @import doParallel
#' @import foreach
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     ROC = COMP_ROC(test[c(1:400),c(1:9)])
#'     ROC
#'     }
COMP_ROC <- function(gene_exp,multi_threaded = T,
                     split_ratios = 0.7,
                     inner_resampless = 5,
                     outer_resampless =5,
                     tuner_resolutions =5) {
  set.seed(520)
  data = data_handle(gene_exp)
  task = task_creat(data)
  if (multi_threaded == T) {
    num_cores <- availableCores() %>% as.numeric()
    cl <- makeCluster(num_cores-1) # 创建并行计算集群
    registerDoParallel(cl) # 注册并行计算集群

    models <- c(ROC_rpart,ROC_knn,ROC_svm,ROC_xgboost,ROC_rr,
                ROC_glmnet,ROC_lda,ROC_logreg,ROC_naive_bayes,ROC_nnet)
    result <- foreach(model = models,
                      .packages = c("mlr3verse", "tidyverse","pROC")) %dopar% {
                        mod_dat <- model(task, data,
                                         split_ratio = split_ratios,
                                         outer_resamples = outer_resampless,
                                         inner_resamples = inner_resampless,
                                         tuner_resolution = tuner_resolutions)
                        return(mod_dat)
                      }
    stopCluster(cl)

  }

  if (multi_threaded != T) {
    rpart = ROC_rpart(task, data,
                      split_ratio = split_ratios,
                      outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                      tuner_resolution = tuner_resolutions)
    knn = ROC_knn(task, data,
                  split_ratio = split_ratios,
                  outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                  tuner_resolution = tuner_resolutions)
    svm = ROC_svm(task, data,
                  split_ratio = split_ratios,
                  outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                  tuner_resolution = tuner_resolutions)
    xgboost = ROC_xgboost(task, data,
                          split_ratio = split_ratios,
                          outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                          tuner_resolution = tuner_resolutions)
    random_forest = ROC_rr(task, data,
                           split_ratio = split_ratios,
                           outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                           tuner_resolution = tuner_resolutions)
    glmnet = ROC_glmnet(task, data,
                        split_ratio = split_ratios,
                        outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                        tuner_resolution = tuner_resolutions)
    lda = ROC_lda(task, data,
                  split_ratio = split_ratios,
                  outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                  tuner_resolution = tuner_resolutions)
    logreg = ROC_logreg(task, data,
                        split_ratio = split_ratios,
                        outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                        tuner_resolution = tuner_resolutions)
    naive_bayes = ROC_naive_bayes(task, data,
                                  split_ratio = split_ratios,
                                  outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                                  tuner_resolution = tuner_resolutions)
    nnet = ROC_nnet(task, data,
                    split_ratio = split_ratios,
                    outer_resamples = outer_resampless,inner_resamples = inner_resampless,
                    tuner_resolution = tuner_resolutions)
    result = list(rpart,knn,svm,xgboost,random_forest,glmnet,lda,logreg,naive_bayes,nnet)
  }
  models_name = c("rpart","knn","svm","xgboost","random_forest",
                  "glmnet","lda","logreg","naive_bayes","nnet")
  for (i in 1:length(result)) {
    auc = result[[i]]$auc %>% as.numeric() %>% round(digits = 2)
    models_name[[i]] = str_c(models_name[[i]],auc,sep = "  AUC:")
  }
  names(result) = models_name
  plot = roc_plot(result)
  plot
}













