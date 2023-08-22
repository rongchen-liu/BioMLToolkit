#' @title Performance evaluation of rpart model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = rpart_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
rpart_model <- function(task,data,
                        inner_resamples = 5,
                        outer_resamples =5,
                        tuner_resolution =5,
                        finial_score_msr = c("classif.auc",
                                             "classif.ce",
                                             "classif.acc",
                                             "classif.precision",
                                             "classif.recall",
                                             "classif.sensitivity",
                                             "classif.specificity")
) {

  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_rpart = lrn("classif.rpart",
                      predict_type = "prob",
                      id ="Classification_Tree")
  search_space_rpart = ps(cp = p_dbl(lower = 0.001, upper = 0.1))
  at_rpart = at_create(tuner,learner_rpart,search_space_rpart,inner_resample,measure)
  rr_rpart = resample(task, at_rpart,outer_resample, store_models = TRUE)  # 自动调参过程(部分)

  res_classif = rr_rpart$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}


#' @title Performance evaluation of KNN model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = knn_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
knn_model <- function(task,data,
                      inner_resamples = 5,
                      outer_resamples =5,
                      tuner_resolution =5,
                      finial_score_msr = c("classif.auc",
                                           "classif.ce",
                                           "classif.acc",
                                           "classif.precision",
                                           "classif.recall",
                                           "classif.sensitivity",
                                           "classif.specificity")
) {
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_kknn = lrn("classif.kknn",
                     predict_type = "prob",
                     id="KNN")
  search_space_kknn = ps(k = p_int(lower = 1, upper = 10))
  at_kknn = at_create(tuner,learner_kknn,search_space_kknn,inner_resample,measure)
  rr_kknn = resample(task, at_kknn,outer_resample, store_models = TRUE)  # 自动调参过程(部分)

  res_classif = rr_kknn$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}


#' @title Performance evaluation of SVM model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = svm_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
svm_model<- function(task,data,
                     inner_resamples = 5,
                     outer_resamples =5,
                     tuner_resolution =5,
                     finial_score_msr = c("classif.auc",
                                          "classif.ce",
                                          "classif.acc",
                                          "classif.precision",
                                          "classif.recall",
                                          "classif.sensitivity",
                                          "classif.specificity")
) {

  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_svm = lrn("classif.svm",
                    type = "C-classification",
                    kernel = "radial",
                    predict_type = "prob",
                    id = "SVM")
  search_space_svm = ps(cost = p_dbl(log(0.1), log(10),
                                     trafo = function(x) exp(x)),
                        gamma = p_dbl(log(0.1), log(10),
                                      trafo = function(x) exp(x)))
  at_svm = at_create(tuner,learner_svm,search_space_svm,inner_resample,measure)
  rr_svm = resample(task, at_svm,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  res_classif = rr_svm$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif

}


#' @title Performance evaluation of Xgboost model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = xgboost_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
xgboost_model<- function(task,data,
                         inner_resamples = 5,
                         outer_resamples =5,
                         tuner_resolution =5,
                         finial_score_msr = c("classif.auc",
                                              "classif.ce",
                                              "classif.acc",
                                              "classif.precision",
                                              "classif.recall",
                                              "classif.sensitivity",
                                              "classif.specificity")
) {
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_xgboost = lrn("classif.xgboost",
                        predict_type = "prob",
                        id ="xgboost")
  search_space_xgboost = ps(
    eta = p_dbl(lower = 0.1, upper = 1),
    max_depth = p_int(lower = 1, upper = 10),
    nrounds = p_int(lower = 1, upper = 16))
  at_xgboost = at_create(tuner,learner_xgboost,search_space_xgboost,inner_resample,measure)
  rr_xgboost = resample(task, at_xgboost,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  res_classif = rr_xgboost$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}


#' @title Performance evaluation of Random_Forest model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = rr_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
rr_model <- function(task,data,
                     inner_resamples = 5,
                     outer_resamples =5,
                     tuner_resolution =5,
                     finial_score_msr = c("classif.auc",
                                          "classif.ce",
                                          "classif.acc",
                                          "classif.precision",
                                          "classif.recall",
                                          "classif.sensitivity",
                                          "classif.specificity")
) {
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_rr = lrn("classif.ranger",
                   predict_type = "prob",
                   id ="Random_Forest")
  search_space_rr = ps(
    mtry = p_int(lower = 1, upper = task$ncol-1),
    min.node.size = p_int(lower = 1, upper = 10),
    num.trees = p_int(lower = 1, upper = 500))

  at_rr = at_create(tuner,learner_rr,search_space_rr,inner_resample,measure)
  rr_rr = resample(task, at_rr,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  res_classif = rr_rr$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}


#' @title Performance evaluation of glmnet model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = glmnet_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
glmnet_model <- function(task,data,
                         inner_resamples = 5,
                         outer_resamples =5,
                         tuner_resolution =5,
                         finial_score_msr = c("classif.auc",
                                              "classif.ce",
                                              "classif.acc",
                                              "classif.precision",
                                              "classif.recall",
                                              "classif.sensitivity",
                                              "classif.specificity")
) {
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_glmnet = lrn("classif.glmnet",
                       predict_type = "prob",
                       id ="Glmnet")
  search_space_glmnet = ps(
    lambda = p_dbl(lower = 0.001, upper = 1, logscale = TRUE),
    alpha = p_dbl(lower = 0, upper = 1))

  at_glmnet = at_create(tuner,learner_glmnet,search_space_glmnet,inner_resample,measure)
  rr_glmnet = resample(task, at_glmnet,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  res_classif = rr_glmnet$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}

#' @title Performance evaluation of LDA model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = lda_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
lda_model <-function(task,data,
                     inner_resamples = 5,
                     outer_resamples =5,
                     tuner_resolution =5,
                     finial_score_msr = c("classif.auc",
                                          "classif.ce",
                                          "classif.acc",
                                          "classif.precision",
                                          "classif.recall",
                                          "classif.sensitivity",
                                          "classif.specificity")
) {
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_LDA = lrn("classif.lda", predict_type = "prob")

  rr_LDA  = resample(task, learner_LDA, outer_resample)

  res_classif = rr_LDA$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}


#' @title Performance evaluation of logreg_model model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = logreg_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
logreg_model <-function(task,data,
                        inner_resamples = 5,
                        outer_resamples =5,
                        tuner_resolution =5,
                        finial_score_msr = c("classif.auc",
                                             "classif.ce",
                                             "classif.acc",
                                             "classif.precision",
                                             "classif.recall",
                                             "classif.sensitivity",
                                             "classif.specificity")
) {
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  #定义学习器
  learner_logreg = lrn("classif.log_reg", predict_type = "prob")
  # 进行交叉验证
  rr_logreg = resample(task, learner_logreg, outer_resample)

  res_classif = rr_logreg$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}

#' @title Performance evaluation of naive_bayes model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = naive_bayes_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
naive_bayes_model <-function(task,data,
                             inner_resamples = 5,
                             outer_resamples =5,
                             tuner_resolution =5,
                             finial_score_msr = c("classif.auc",
                                                  "classif.ce",
                                                  "classif.acc",
                                                  "classif.precision",
                                                  "classif.recall",
                                                  "classif.sensitivity",
                                                  "classif.specificity")
) {
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_naive_bayes = lrn("classif.naive_bayes", predict_type = "prob")
  search_space_naive_bayes = ps(
    laplace = p_dbl(lower = 0.1, upper = 1),
    threshold = p_int(lower = 1, upper = 10))
  at_naive_bayes = at_create(tuner,learner_naive_bayes,search_space_naive_bayes,inner_resample,measure)
  rr_naive_bayes = resample(task, at_naive_bayes,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  res_classif = rr_naive_bayes$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}


#' @title Performance evaluation of nnet model
#' @param task A task of mlr3 \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @return A data.frame containing the values of specific performance metrics
#'      for each iteration of model training.
#' @export
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = nnet_model(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5,
#'         finial_score_msr = c("classif.auc","classif.ce","classif.acc",
#'         "classif.precision","classif.recall","classif.sensitivity",
#'          "classif.specificity")
#'         )
#'     result[c(1:5),]
#'      }
nnet_model <- function(task,data,
                       inner_resamples = 5,
                       outer_resamples =5,
                       tuner_resolution =5,
                       finial_score_msr = c("classif.auc",
                                            "classif.ce",
                                            "classif.acc",
                                            "classif.precision",
                                            "classif.recall",
                                            "classif.sensitivity",
                                            "classif.specificity")
) {
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)#内层重采样策略
  outer_resample <- rsmp("cv", folds = outer_resamples)#外层重采样策略
  tuner <- tnr("grid_search",resolution = tuner_resolution)#超参数搜索策略
  final_score <- msrs(finial_score_msr)

  learner_nnet = lrn("classif.nnet",predict_type = "prob")
  search_space_nnet = ps(
    size = p_int(lower = 1, upper = 10),
    decay = p_dbl(lower = 0.1, upper = 0.9),
    MaxNWts= p_int(lower = 10000, upper = 10000))

  at_nnet = at_create(tuner,learner_nnet,search_space_nnet,inner_resample,measure)
  rr_nnet = resample(task, at_nnet,outer_resample, store_models = TRUE)  # 自动调参过程(部分)
  res_classif = rr_nnet$score(final_score) %>%
    as.data.table() %>%
    as_tibble() %>%
    select(task_id,learner_id,iteration,contains("classif")) %>%
    mutate(mean_auc = classif.auc %>% mean()) %>%
    mutate(mean_ce = classif.ce %>% mean) %>%
    mutate(mean_acc = classif.acc %>% mean()) %>%
    mutate(mean_precision = classif.precision %>% mean()) %>%
    mutate(mean_recall = classif.recall %>% mean()) %>%
    mutate(mean_sensitivity = classif.sensitivity %>% mean()) %>%
    mutate(mean_specificit = classif.specificity %>% mean())
  res_classif
}



#' @title Integrate the model performance evaluation across 10 different machine learning models.
#'
#' @param gene_exp For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param result_type For a character variable, the parameter "result_type" can be one of the following:
#'     "data," "plot," or "all."
#'     If "result_type" is set to "data,"it will return a data frame that records all the data from training 10 different models.
#'     If "result_type" is set to "plot," it will return a plotted image.
#'     If the parameter "all" is used, it will return a list containing the results for the previous two options.
#' @param inner_resampless For a character variable,\code{\link[mlr3]{rsmp}}
#' @param outer_resampless For a character variable,\code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @param finial_score_msr A character vector representing performance metrics for the model evaluation,
#'      including "classif.auc," "classif.ce," "classif.acc," "classif.precision,"
#'       "classif.recall," "classif.sensitivity," and "classif.specificity."
#' @param multi_threaded A logic variable to enable multi-threading during execution.
#'
#' @return Return a list containing four elements: "data","plot" In the "data" element,
#'      you will find the specific results of the 10 models across multiple resampling iterations.
#'      The "plot" element will contain the model performance plot.
#' @export
#' @import doParallel
#' @import parallelly
#' @import foreach
#' @import mlr3
#' @importFrom mlr3tuning tnr auto_tuner
#' @import paradox
#' @import dplyr
#' @import stringr
#' @import ggplot2
#' @import tidyr
#' @examples \donttest{
#'     data(test)
#'     res = efficiency_mod(test[,c(1:10)])
#'     res$plot
#'     res$data[c(1:10),]
#'      }
efficiency_mod <- function(gene_exp,
                           result_type = "all",
                           inner_resampless = 5,
                           outer_resampless =5,
                           tuner_resolutions =5,
                           finial_score_msrs = c("classif.auc",
                                                "classif.ce",
                                                "classif.acc",
                                                "classif.precision",
                                                "classif.recall",
                                                "classif.sensitivity",
                                                "classif.specificity"),
                           multi_threaded = T) {
  data = data_handle(gene_exp)
  task = task_creat(data)
  if (multi_threaded == T) {
    num_cores <- availableCores() %>% as.numeric()
    cl <- makeCluster(num_cores-1) # 创建并行计算集群
    registerDoParallel(cl) # 注册并行计算集群

    models <- c(rpart_model,knn_model,svm_model,xgboost_model,rr_model,
                glmnet_model,lda_model,logreg_model,naive_bayes_model,nnet_model)
    # 并行计算
    result_classif <- foreach(model = models, .combine = rbind,
                              .packages = c("tidyverse","mlr3verse")) %dopar% {
      mod_dat = model(task,data,
                      inner_resamples = inner_resampless,
                      outer_resamples = outer_resampless,
                      finial_score_msr = finial_score_msrs)
      return(mod_dat)
    }
    # 关闭并行计算集群
    stopCluster(cl)
  }
  if (multi_threaded != T) {

    models <- c(rpart_model,knn_model,svm_model,xgboost_model,rr_model,
                glmnet_model,lda_model,logreg_model,naive_bayes_model,nnet_model)

    rpart = rpart_model(task,data,
                inner_resamples = inner_resampless,
                outer_resamples = outer_resampless,
                finial_score_msr = finial_score_msrs)
    knn = knn_model(task,data,
                    inner_resamples = inner_resampless,
                    outer_resamples = outer_resampless,
                    finial_score_msr = finial_score_msrs)
    svm = svm_model(task,data,
                    inner_resamples = inner_resampless,
                    outer_resamples = outer_resampless,
                    finial_score_msr = finial_score_msrs)
    xgboost = xgboost_model(task,data,
                            inner_resamples = inner_resampless,
                            outer_resamples = outer_resampless,
                            finial_score_msr = finial_score_msrs)
    rr = rr_model(task,data,
                  inner_resamples = inner_resampless,
                  outer_resamples = outer_resampless,
                  finial_score_msr = finial_score_msrs)
    glmnet = glmnet_model(task,data,
                         inner_resamples = inner_resampless,
                         outer_resamples = outer_resampless,
                         finial_score_msr = finial_score_msrs)
    lda = lda_model(task,data,
                    inner_resamples = inner_resampless,
                    outer_resamples = outer_resampless,
                    finial_score_msr = finial_score_msrs)
    logreg = logreg_model(task,data,
                          inner_resamples = inner_resampless,
                          outer_resamples = outer_resampless,
                          finial_score_msr = finial_score_msrs)
    naive_bayes = naive_bayes_model(task,data,
                                    inner_resamples = inner_resampless,
                                    outer_resamples = outer_resampless,
                                    finial_score_msr = finial_score_msrs)
    nnet = nnet_model(task,data,
                      inner_resamples = inner_resampless,
                      outer_resamples = outer_resampless,
                      finial_score_msr = finial_score_msrs)
    result_classif = rpart %>%
      bind_rows(knn) %>%
      bind_rows(svm) %>%
      bind_rows(xgboost) %>%
      bind_rows(rr) %>%
      bind_rows(glmnet) %>%
      bind_rows(lda) %>%
      bind_rows(logreg) %>%
      bind_rows(naive_bayes) %>%
      bind_rows(nnet)
  }
  result_classif = result_classif %>% mutate(
    model = learner_id %>% str_replace_all(".tuned",""))
  result_classif$model = result_classif$model %>% str_replace_all("classif.","")
  colnames(result_classif) = colnames(result_classif) %>% str_replace_all("classif.","")

  hp_data = result_classif %>% dplyr::select(model,contains("mean")) %>% unique()%>%
    pivot_longer(-model,
                 names_to = "measures",
                 values_to = "value")

  hp_data$value = hp_data$value %>% round(digits = 2)
  plot_hp = ggplot(hp_data, aes(x= measures,y = model))+
    geom_tile(aes(fill = value),size = 1)+
    scale_fill_gradient2(low = "yellow",
                         high = "#F39B7F")+
    geom_text(aes(label=value),col ="black")+
    theme_bw()+ theme(axis.text = element_text(size = 10,
                                               colour = "gray7"),
                      axis.text.x = element_text(angle = 90)) +labs(x = NULL, y = NULL)
  res_EM = list(data_EM = result_classif,heatmap_EM = plot_hp)

  if (result_type == "data") {return(result_classif)}
  if (result_type == "plot") {return(plot_hp)}
  if (result_type == "all") {return(res_EM)}
}
