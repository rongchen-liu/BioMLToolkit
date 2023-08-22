#' @title Standardize Data
#' @description Change "-" in column names to "_",
#'     Convert the type to a factor variable.
#' @param data For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @return Return the processed data.
#' @export
#' @importFrom stringr str_detect str_replace_all %>%
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)}
data_handle  <- function(data) {
  if (
    str_detect(colnames(data),"-")  %>% unique() %>% length() != 1
  ) {
    colnames(data) = colnames(data) %>% str_replace_all("-","_")
  }
  data$type = data$type %>% as.factor()
  data
}

#' @title Create Task
#'
#' @param data Result of data_handle().
#' @param target Name of the target column.Default is "target"
#'
#' @return TaskClassif,a task base on mlr3
#' @export
#' @importFrom mlr3 as_task_classif
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)}
task_creat <- function(data,target = "type") {
  task = as_task_classif(data,target = "type")
  task
}

#' @title Create an automatic tuning machine for mlr3.
#'
#' @param tuner \code{\link[mlr3tuning]{auto_tuner}}
#' @param learner \code{\link[mlr3tuning]{auto_tuner}}
#' @param search_space \code{\link[mlr3tuning]{auto_tuner}}
#' @param inner_resample \code{\link[mlr3tuning]{auto_tuner}}
#' @param measure \code{\link[mlr3tuning]{auto_tuner}}
#'
#' @return The AutoTuner wraps a mlr3::Learner and augments it with an automatic tuning process for a given set of hyperparameters.
#' @export
#' @importFrom mlr3tuning auto_tuner
#' @examples \donttest{
#'     library(mlr3verse)
#'     measure <- msr("classif.auc")
#'     inner_resample <- rsmp("cv", folds = 5)
#'     outer_resample <- rsmp("cv", folds = 5)
#'     tuner <- tnr("grid_search",resolution = 5)
#'     learner_rpart = lrn("classif.rpart",predict_type = "prob",id ="Classification_Tree")
#'     search_space_rpart = ps(cp = p_dbl(lower = 0.001, upper = 0.1))
#'     at_rpart = at_creat(tuner,learner_rpart,search_space_rpart,inner_resample,measure)
#'     }
at_creat <- function(tuner,learner,search_space,inner_resample,measure) {
  at = auto_tuner(
    tuner = tuner,
    learner = learner,
    search_space = search_space,
    resampling = inner_resample,
    measure = measure)
  at
}

#' @title Obtain variable importance using a rpart(Classification_Tree) model.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     rpart = rpart_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5 )}
rpart_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

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

  at_rpart = at_create(tuner,learner_rpart,search_space_rpart,inner_resample,measure)
  rr_rpart = resample(task, at_rpart,outer_resample, store_models = TRUE)

  at_rpart$train(task = task)
  rpart_exp = DALEXtra::explain_mlr3(at_rpart,
                                     data = data,
                                     y = data$type %>% as.numeric(),
                                     label = "rpart",
                                     colorize = FALSE)
  importance_rpart<-variable_importance(
    rpart_exp,
    loss_function = loss_root_mean_square
  )

  imp_data = importance_data(importance_rpart)
  imp_data
}

#' @title Obtain variable importance using a KNN model.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     KNN = knn_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5 )}
knn_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_kknn = lrn("classif.kknn",
                     predict_type = "prob",
                     id="KNN")
  search_space_kknn = ps(k = p_int(lower = 1, upper = 10))
  at_create <- function(tuner,learner,search_space,inner_resample,measure) {
    at = auto_tuner(
      tuner = tuner,
      learner = learner,
      search_space = search_space,
      resampling = inner_resample,
      measure = measure)

    at
  }
  at_kknn = at_create(tuner,learner_kknn,search_space_kknn,inner_resample,measure)
  rr_kknn = resample(task, at_kknn,outer_resample, store_models = TRUE)
  at_kknn$train(task)
  kknn_exp = DALEXtra::explain_mlr3(at_kknn,
                                    data = data,
                                    y = data$type %>% as.numeric(),
                                    label = "KNN",
                                    colorize = FALSE)
  importance_kknn<-variable_importance(
    kknn_exp,
    loss_function = loss_root_mean_square
  )
  imp_data = importance_data(importance_kknn)
  imp_data
}

#' @title Obtain variable importance using a SVM model.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     SVM = svm_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     SVM[c(1:10),]
#'      }
svm_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_svm = lrn("classif.svm",
                    type = "C-classification",
                    kernel = "radial",
                    predict_type = "prob",
                    id = "SVM")
  search_space_svm = ps(cost = p_dbl(log(0.1), log(10),
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
  at_svm = at_create(tuner,learner_svm,search_space_svm,inner_resample,measure)
  rr_svm = resample(task, at_svm,outer_resample, store_models = TRUE)

  at_svm$train(task)
  svm_exp = DALEXtra::explain_mlr3(at_svm,
                                   data = data,
                                   y = data$type %>% as.numeric(),
                                   label = "SVM",
                                   colorize = FALSE)
  importance_svm <- variable_importance(
    svm_exp,
    loss_function = loss_root_mean_square)
  imp_res = importance_data(importance_svm)

  imp_res
}

#' @title Obtain variable importance using a Xgboost model.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = xgboost_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     result[c(1:10),]
#'      }
xgboost_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_xgboost = lrn("classif.xgboost",
                        predict_type = "prob",
                        id ="xgboost")
  search_space_xgboost = ps(
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
  at_xgboost = at_create(tuner,learner_xgboost,search_space_xgboost,inner_resample,measure)
  rr_xgboost = resample(task, at_xgboost,outer_resample, store_models = TRUE)

  at_xgboost$train(task)
  xgboost_exp = DALEXtra::explain_mlr3(at_xgboost,
                                       data = data,
                                       y = data$type %>% as.numeric(),
                                       label = "xgboost",
                                       colorize = FALSE)
  imp_res <-variable_importance(
    xgboost_exp,
    loss_function = loss_root_mean_square
  )
  imp_data = importance_data(imp_res)
  imp_data
}


#' @title Obtain variable importance using a Random_Forest model.
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = rr_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     result[c(1:10),]
#'      }
rr_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_rr = lrn("classif.ranger",
                   predict_type = "prob",
                   id ="Random_Forest")
  search_space_rr = ps(
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
  at_rr = at_create(tuner,learner_rr,search_space_rr,inner_resample,measure)
  rr_rr = resample(task, at_rr,outer_resample, store_models = TRUE)
  at_rr$train(task)
  rr_exp = DALEXtra::explain_mlr3(at_rr,
                                  data = data,
                                  y = data$type %>% as.numeric(),
                                  label = "Random_Forest",
                                  colorize = FALSE)
  imp_res <-variable_importance(
    rr_exp,
    loss_function = loss_root_mean_square
  )
  imp_data = importance_data(imp_res)
  imp_data
}

#' @title Obtain variable importance using a SVM model.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     SVM = svm_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     SVM[c(1:10),]
#'      }
svm_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_svm = lrn("classif.svm",
                    type = "C-classification",
                    kernel = "radial",
                    predict_type = "prob",
                    id = "SVM")
  search_space_svm = ps(cost = p_dbl(log(0.1), log(10),
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
  at_svm = at_create(tuner,learner_svm,search_space_svm,inner_resample,measure)
  rr_svm = resample(task, at_svm,outer_resample, store_models = TRUE)
  at_svm$train(task)
  svm_exp = DALEXtra::explain_mlr3(at_svm,
                                   data = data,
                                   y = data$type %>% as.numeric(),
                                   label = "SVM",
                                   colorize = FALSE)
  importance_svm <- variable_importance(
    svm_exp,
    loss_function = loss_root_mean_square)
  imp_res = importance_data(importance_svm)

  imp_res
}

#' @title Obtain variable importance using a Xgboost model.
#'
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#'
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = xgboost_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     result[c(1:10),]
#'      }
xgboost_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_xgboost = lrn("classif.xgboost",
                        predict_type = "prob",
                        id ="xgboost")
  search_space_xgboost = ps(
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
  at_xgboost = at_create(tuner,learner_xgboost,search_space_xgboost,inner_resample,measure)
  rr_xgboost = resample(task, at_xgboost,outer_resample, store_models = TRUE)

  at_xgboost$train(task)
  xgboost_exp = DALEXtra::explain_mlr3(at_xgboost,
                                       data = data,
                                       y = data$type %>% as.numeric(),
                                       label = "xgboost",
                                       colorize = FALSE)
  imp_res <-variable_importance(
    xgboost_exp,
    loss_function = loss_root_mean_square
  )
  imp_data = importance_data(imp_res)
  imp_data
}


#' @title Obtain variable importance using a glmnet model.
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = glmnet_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     result[c(1:10),]
#'      }
glmnet_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_glmnet = lrn("classif.glmnet",
                       predict_type = "prob",
                       id ="Glmnet")
  search_space_glmnet = ps(
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
  at_glmnet = at_create(tuner,learner_glmnet,search_space_glmnet,inner_resample,measure)
  rr_glmnet = resample(task, at_glmnet,outer_resample, store_models = TRUE)

  at_glmnet$train(task)
  glmnet_exp = DALEXtra::explain_mlr3(at_glmnet,
                                      data = data,
                                      y = data$type %>% as.numeric(),
                                      label = "glmnet",
                                      colorize = FALSE)
  imp_res<-variable_importance(
    glmnet_exp,
    loss_function = loss_root_mean_square)
  imp_data = importance_data(imp_res)
  imp_data
}



#' @title Obtain variable importance using a LDA model.
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = lda_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     result[c(1:10),]
#'      }
lda_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")#模型训练评价指标
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_LDA = lrn("classif.lda", predict_type = "prob")

  rr_LDA  = resample(task, learner_LDA, outer_resample)

  learner_LDA$train(task)
  LDA_exp = DALEXtra::explain_mlr3(learner_LDA,
                                   data = data,
                                   y = data$type %>% as.numeric(),
                                   label = "lda",
                                   colorize = FALSE)
  imp_res<-variable_importance(
    LDA_exp,
    loss_function = loss_root_mean_square
  )
  imp_data = importance_data(imp_res)
  imp_data
}


#' @title Obtain variable importance using a logreg model.
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = logreg_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     result[c(1:10),]
#'      }
logreg_imp <-function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_logreg = lrn("classif.log_reg", predict_type = "prob")
  rr_logreg = resample(task, learner_logreg, outer_resample)

  learner_logreg$train(task)
  logreg_exp = DALEXtra::explain_mlr3(learner_logreg,
                                      data = data,
                                      y = data$type %>% as.numeric(),
                                      label = "logreg",
                                      colorize = FALSE)
  imp_res<-variable_importance(
    logreg_exp,
    loss_function = loss_root_mean_square
  )
  imp_data = importance_data(imp_res)
  imp_data
}

#' @title Obtain variable importance using a naive_bayes model.
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = naive_bayes_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     result[c(1:10),]
#'      }
naive_bayes_imp <-function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_naive_bayes = lrn("classif.naive_bayes", predict_type = "prob")
  search_space_naive_bayes = ps(
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
  at_naive_bayes = at_create(tuner,learner_naive_bayes,search_space_naive_bayes,inner_resample,measure)
  rr_naive_bayes = resample(task, at_naive_bayes,outer_resample, store_models = TRUE)

  at_naive_bayes$train(task)
  naive_bayes_exp = DALEXtra::explain_mlr3(at_naive_bayes,
                                           data = data,
                                           y = data$type %>% as.numeric(),
                                           label = "Naive_bayes",
                                           colorize = FALSE)
  imp_res<-variable_importance(
    naive_bayes_exp,
    loss_function = loss_root_mean_square)
  imp_data = importance_data(imp_res)
  imp_data
}

#' @title Obtain variable importance using a nnet model.
#' @param task \code{\link[mlr3]{resample}}
#' @param data Result of data_handle().For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param inner_resamples \code{\link[mlr3]{rsmp}}
#' @param outer_resample \code{\link[mlr3]{rsmp}}
#' @param tuner_resolution \code{\link[mlr3tuning]{tnr}}
#' @return A data.frame records the dropout_loss of a model for each sampling operation,
#'     which represents the degree of impact of removing a variable on the model.
#'     The evaluation is done using RMSE (Root Mean Square Error)
#' @export
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows left_join
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     data = data_handle(test)
#'     task = task_creat(data)
#'     result = nnet_imp(task,data,
#'         inner_resamples = 5,
#'         outer_resample = 5,
#'         tuner_resolution = 5)
#'     result[c(1:10),]
#'      }
nnet_imp <- function(task,data,inner_resamples = 5,outer_resample =5,tuner_resolution =5) {
  importance_data = function(imp_model){
    importance_handle <- function(importance_result) {
      imp_data = importance_result %>% as.data.frame()
      variable_name = imp_data[,"variable"] %>% unique()
      model_name = imp_data$label %>% unique()
      result = data.frame()
      for (i in 1:length(variable_name)) {
        vb_name = variable_name[i]
        db_int = imp_data %>% dplyr::filter(variable == vb_name)
        mean_dropout_loss = db_int$dropout_loss %>% mean()
        res = data.frame(variable = vb_name,mean_dropout_loss,model = model_name)
        result = result %>% bind_rows(res)
      }
      result
    }
    imp_data = imp_model %>% as.data.frame()%>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_") %>% filter(variable != "type")
    imp_res = importance_handle(imp_model)
    imp_res = imp_res %>% filter(variable != "_full_model_") %>%
      filter(variable != "_baseline_")
    imp_data = imp_data %>% left_join(imp_res[,c("variable","mean_dropout_loss")],
                                      by = c("variable"))
    imp_data
  }
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = inner_resamples)
  outer_resample <- rsmp("cv", folds = outer_resample)
  tuner <- tnr("grid_search",resolution = tuner_resolution)

  learner_nnet = lrn("classif.nnet",predict_type = "prob")
  search_space_nnet = ps(
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
  at_nnet = at_create(tuner,learner_nnet,search_space_nnet,inner_resample,measure)
  rr_nnet = resample(task, at_nnet,outer_resample, store_models = TRUE)
  at_nnet$train(task)
  nnet_exp = DALEXtra::explain_mlr3(at_nnet,
                                    data = data,
                                    y = data$type %>% as.numeric(),
                                    label = "nnet",
                                    colorize = FALSE)
  imp_res <-variable_importance(
    nnet_exp,
    loss_function = loss_root_mean_square)
  imp_data = importance_data(imp_res)
  imp_data
}


#' @title Integrate the variable (gene) importance evaluation across 10 different machine learning models.
#'
#' @param gene_exp For a gene expression matrix, it is required to have sample names as row names and gene names as column names.
#'     Additionally, there should be an additional column called "type" consisting of 0 and 1, representing whether a sample has been treated or not.
#'     Here, 0 represents the control group (untreated), and 1 represents the treated group.
#' @param result_type For a character variable, the parameter "result_type" can be one of the following:
#'     "table," "plot," "gene," "gene_model," or "all."
#'     If "result_type" is set to "table,"it will return a data frame that records all the data from training 10 different models.
#'     If "result_type" is set to "plot," it will return a plotted image.
#'     If "result_type" is set to "gene," it will return the filtered gene names.
#'     If "result_type" is set to "gene_model," it will return the filtered gene names followed by the model.
#'     If the parameter "all" is used, it will return a list containing the results for the previous four options.
#' @param num For a numeric variable，Specify the number of top-ranked variables to select for filtering.
#' @param inner_resampless For a character variable,\code{\link[mlr3]{rsmp}}
#' @param outer_resampless For a character variable,\code{\link[mlr3]{rsmp}}
#' @param tuner_resolutions For a character variable\code{\link[mlr3tuning]{tnr}}
#' @param multi_threaded A logic variable to enable multi-threading during execution.
#'
#' @return Return a list containing four elements: "table," "plot," "gene," and "gene_model." In the "table" element,
#'      you will find the specific results of the 10 models across multiple resampling iterations.
#'      The "plot" element will contain the gene importance plot. The "gene" element will include the filtered gene names.
#'       And the "gene_model" element will provide the filtered gene names followed by the model name.
#' @export
#' @import future
#' @import parallel
#' @import parallelly
#' @import doParallel
#' @import foreach
#' @importFrom mlr3 msr rsmp lrn resample
#' @importFrom mlr3tuning tnr auto_tuner
#' @importFrom paradox ps
#' @import ggplot2
#' @import tidyverse
#' @importFrom ggsci scale_fill_npg
#' @import dplyr
#' @import stringr
#' @importFrom DALEXtra explain_mlr3
#' @importFrom DALEX variable_importance loss_root_mean_square
#' @examples \donttest{
#'     data(test)
#'     res = Comp_imp(test[c(1:400),c(1:9)])
#'     res$plot
#'      }
Comp_imp = function(gene_exp,
                    result_type = "all",
                    num = 5,
                    inner_resampless = 5,
                    outer_resampless = 5,
                    tuner_resolutions = 5,
                    multi_threaded = T){

  data = data_handle(gene_exp)
  task = task_creat(data)
  plot_data = function(imp_data,num){
    select_gene = imp_data[order(imp_data$mean_dropout_loss,decreasing = T),] %>%
      pull(var = "variable") %>% unique() %>% .[1:num]
    data_plot = imp_data %>% filter(variable %in% select_gene) %>%
      arrange(mean_dropout_loss) %>%
      mutate(variable = factor(variable,levels = variable %>% unique()))
    data_plot
  }
  importance_barplot = function(data_plot){
    plot = ggplot(data = data_plot,aes(x = dropout_loss,y = variable))+
      geom_violin(aes(fill = label))+
      geom_boxplot(aes(fill = label))+
      scale_fill_npg()+
      theme_bw()+
      labs(title = "Feature importance", x = "Root mean square error (RMSE) loss after permutations ")

    plot
  }
  if (multi_threaded == T) {
    set.seed(520)
    num_cores <- availableCores() %>% as.numeric()

    cl <- makeCluster(num_cores-1)
    registerDoParallel(cl)
    models = c(rpart_imp,knn_imp,svm_imp,xgboost_imp,rr_imp,
               glmnet_imp,lda_imp,logreg_imp,naive_bayes_imp,nnet_imp)

    models <- foreach(model = models,
                      .packages = c("mlr3verse", "tidyverse","DALEX","DALEXtra"),
                      .export = c("inner_resampless","outer_resampless","tuner_resolutions")) %dopar% {
                        mod <- model(task,data,
                                     inner_resamples = inner_resampless,
                                     outer_resample = outer_resampless,
                                     tuner_resolution =tuner_resolutions )
                        return(mod)
                      }
    stopCluster(cl)
  }else{
    rpart = rpart_imp(task,data,
                      inner_resamples = inner_resampless,
                      outer_resample = outer_resampless,
                      tuner_resolution =tuner_resolutions )
    knn = knn_imp(task,data,
                  inner_resamples = inner_resampless,
                  outer_resample = outer_resampless,
                  tuner_resolution =tuner_resolutions )
    svm = svm_imp(task,data,
                  inner_resamples = inner_resampless,
                  outer_resample = outer_resampless,
                  tuner_resolution =tuner_resolutions )
    xgboost = xgboost_imp(task,data,
                          inner_resamples = inner_resampless,
                          outer_resample = outer_resampless,
                          tuner_resolution =tuner_resolutions )
    rr = rr_imp(task,data,
                inner_resamples = inner_resampless,
                outer_resample = outer_resampless,
                tuner_resolution =tuner_resolutions )
    glmnet = glmnet_imp(task,data,
                        inner_resamples = inner_resampless,
                        outer_resample = outer_resampless,
                        tuner_resolution =tuner_resolutions )
    lda = lda_imp(task,data,
                  inner_resamples = inner_resampless,
                  outer_resample = outer_resampless,
                  tuner_resolution =tuner_resolutions )
    logreg = logreg_imp(task,data,
                        inner_resamples = inner_resampless,
                        outer_resample = outer_resampless,
                        tuner_resolution =tuner_resolutions )
    naive_bayes = naive_bayes_imp(task,data,
                                  inner_resamples = inner_resampless,
                                  outer_resample = outer_resampless,
                                  tuner_resolution =tuner_resolutions )
    nnet = nnet_imp(task,data,
                    inner_resamples = inner_resampless,
                    outer_resample = outer_resampless,
                    tuner_resolution =tuner_resolutions )


    models = list(rpart,knn,svm,xgboost,rr,
                  glmnet,lda,logreg,naive_bayes,nnet)
  }

  models_name = c("rpart","knn","svm","xgboost","random_forest",
                  "glmnet","lda","logreg","naive_bayes","nnet")

  names(models) = models_name
  data_plot_all = data.frame()
  for (i in 1:length(models)) {
    model = models[[i]]
    data_plot = plot_data(model,num)
    data_plot_all = data_plot_all %>% bind_rows(data_plot)
  }
  gene = pull(data_plot_all,var = "variable") %>% unique() %>% as.character()
  data_plot_all = data_plot_all %>% mutate(
    variable = str_c(variable,label,sep = "_")
  )

  gene_model = pull(data_plot_all,var = "variable") %>% unique() %>% as.character()
  data_plot_all = data_plot_all %>%  arrange(label,mean_dropout_loss) %>%
    mutate(variable = factor(variable,levels = variable %>% unique()))

  data = models$rpart %>%
    bind_rows(models$knn) %>%
    bind_rows(models$svm) %>%
    bind_rows(models$xgboost) %>%
    bind_rows(models$random_forest) %>%
    bind_rows(models$glmnet) %>%
    bind_rows(models$lda) %>%
    bind_rows(models$logreg) %>%
    bind_rows(models$naive_bayes) %>%
    bind_rows(models$nnet)

  plot = importance_barplot(data_plot_all)
  res = list(plot = plot ,data = data,gene = gene,gene_model = gene_model)
  if (result_type == "table") {return(data)}
  if (result_type == "plot") {return(plot)}
  if (result_type == "gene") {return(gene)}
  if (result_type == "gene_model") {return(gene_model)}
  if (result_type == "all") {return(res)}
}
