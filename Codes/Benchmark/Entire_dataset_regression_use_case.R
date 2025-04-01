# - Benchmark example usage a Regression Problem -
#                   No Feature Selection
#-------------------------------------------------------------------------------
# Usage with Sentinel-2 spectral bands and vegetation indices in predicting Leaf Chlorophyll Content (LCC)
# ------------------------------------------------------------------------------
# Library
library(ggplot2)
library(ggpubr)
library(grid)
library(ggExtra)
library(RColorBrewer)
library(grDevices)
library(GGally)
library(gridExtra)
library(corrplot)
library(psych)
# ------------------------Prediction algorithms---------------------------------
# This code implements various regression algorithms in MLR package. 
# Some portions of the code are generic and must not be changed, unless you are familiar with the MLR package. 
# Please run one algorithm at a time. The parameters override previous ones. 
#------------------------------------------------------------------
library(mlr)
# -------------------Algorithm Settings----------------------------
#------------------------------------------------------------------
as.data.frame(mlr::listLearners())
# -------------------Generic Settings------------------------------
# Set Cross-validation strategy
rdesc = makeResampleDesc("CV", stratify = F, iters =10)
# ------------------- (1) Random Forest Parameters-----------------
# Do not other algorithms after running this section. Jump to parameter tuning. 
# Learner selection and parameter values
# getLearnerPackages("rregr.randomForest")
lrn <- makeLearner("regr.randomForest", predict.type = "response")
# set parameter space
# getParamSet("regr.randomForest")
params <- makeParamSet(makeIntegerParam("ntree",lower = 100L, upper = 1000L),
                       makeIntegerParam("mtry",lower = 1L,upper = 10L))
# Grid-search strategy
ctrl <- makeTuneControlGrid()
# ------------------- (2) eXtreme Gradient Boosting Parameters------------------
# Do not other algorithms after running this section. Jump to parameter tuning.
# Learner selection and parameter values
# getLearnerPackages("regr.xgboost")
lrn <- makeLearner("regr.xgboost", predict.type = "response")

lrn$par.vals <- list(objective="reg:squarederror", eval_metric="rmse", nrounds=1000L, eta=0.1)

# set parameter space
# getParamSet("regr.xgboost")
params <- makeParamSet(makeDiscreteParam("booster",values = c("gbtree","gblinear")), makeIntegerParam("max_depth",lower = 3L,upper = 10L), makeNumericParam("min_child_weight",lower = 1L,upper = 10L), makeNumericParam("subsample",lower = 0.5,upper = 1), makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))

# Random search strategy
ctrl <- makeTuneControlRandom(maxit = 100L)

# ------------------- (3) Support Vector Machines Parameters--------------------
# Do not other algorithms after running this section. Jump to parameter tuning.
# set parameter space
# getParamSet("regr.svm")
lrn <- makeLearner("regr.svm", predict.type = "response")

lrn$par.vals <- list(objective="reg:squarederror", eval_metric="rmse", nrounds=1000L, eta=0.1)

# set parameter space
# getParamSet("regr.xgboost")
params <- makeParamSet(makeDiscreteParam("type", values = c("eps-regression","nu-regression")), 
                       makeDiscreteParam("kernel", values = c("radial","linear","polynomial", "sigmoid")), 
                       makeNumericParam("gamma",lower = 0L,upper = 10L), makeNumericParam("cost",lower = 0,upper = 1), 
                       makeNumericParam("degree",lower = 0, upper = 1))

# Random search strategy
ctrl <- makeTuneControlRandom(maxit = 100L)
# ------------------- (4) RF-Conditional Inference Trees Parameters-------------
# Do not other algorithms after running this section. Jump to parameter tuning.
# Random Forest Based on Conditional Inference Trees
# set parameter space
# getParamSet("regr.cforest")
lrn <- makeLearner("regr.cforest", predict.type = "response")

# set parameter space
# getParamSet("regr.cforest"")
params <- makeParamSet(makeIntegerParam("ntree",lower = 100L, upper = 1000L),
                       makeIntegerParam("mtry",lower = 1L,upper = 10L))
# Random search strategy
ctrl <- makeTuneControlRandom(maxit = 100L)

# ---------(4) Multivariate Adaptive Regression Splines Parameters--------------
# Do not other algorithms after running this section. Jump to parameter tuning.
# https://uc-r.github.io/mars
# set parameter space
# getParamSet("regr.earth")
lrn <- makeLearner("regr.earth", predict.type = "response")

# set parameter space
# getParamSet("regr.cforest"")
params <- makeParamSet(makeIntegerParam("degree",lower = 1L, upper = 5L),
                       makeIntegerParam("nprune",lower = 2L,upper = 100L))

# Random search strategy
ctrl <- makeTuneControlRandom(maxit = 100L)

# ------(5) Evolutionary learning of globally optimal trees Parameters----------
# Do not other algorithms after running this section. Jump to parameter tuning.
# install.packages("evtree")
# set parameter space
# getParamSet("regr.evtree")
lrn <- makeLearner("regr.evtree", predict.type = "response")

# set parameter space
# getParamSet("regr.evtree")
params <- makeParamSet(makeIntegerParam("ntrees",lower = 100L, upper = 1000L),
                       makeIntegerParam("minbucket",lower = 1L,upper = 10L),
                       makeIntegerParam("minsplit",lower = 1L,upper = 40L),
                       makeIntegerParam("maxdepth",lower = 1L,upper = 20L),
                       makeIntegerParam("psplit",lower = 1L,upper = 40L),
                       makeIntegerParam("pprune",lower = 1L,upper = 40L))
# Random search strategy
ctrl <- makeTuneControlRandom(maxit = 100L)

# ----------(6) Gaussian Process Regression Parameters--------------------------
# Do not other algorithms after running this section. Jump to parameter tuning.
# set parameter space
# getParamSet("regr.gausspr")
lrn <- makeLearner("regr.gausspr", predict.type = "response")

# set parameter space
# getParamSet("regr.gausspr")
params <- makeParamSet(makeDiscreteParam("kernel",values = c("rbfdot")),
                       makeIntegerParam("sigma",lower = 1e-5,upper = 1e+5))#,
#makeIntegerParam("degree",lower = 1,upper = 5))
# Random search strategy
ctrl <- makeTuneControlRandom(maxit = 100L)
# ----------(7) Cubist Parameters--------------------------
# Do not other algorithms after running this section. Jump to parameter tuning.
# set parameter space
getParamSet("regr.cubist")
lrn <- makeLearner("regr.cubist", predict.type = "response")

# set parameter space
# getParamSet("regr.cubist")
params <- makeParamSet(makeIntegerParam("committees",lower = 1,upper = 100),
                       makeIntegerParam("neighbors",lower = 1,upper = 9),
                       makeIntegerParam("rules",lower = 10,upper = 20))

# Random search strategy
ctrl <- makeTuneControlRandom(maxit = 100L)

# -------------------Computing settings----------------------------
#----------------------------------------------------------
# Load the data -----------------------------------------------------------
LCC_dataset <- read.csv("./Data/Regression-LCC/Bothaville_LCC_2021.csv", head = TRUE, as.is=F)
head(LCC_dataset)

# Pre-process the data -----------------------------------------------------------
# Remove unnecessary columns
LCC_dataset <- subset(LCC_dataset, select = -c(LAI, CCC))
head(LCC_dataset)

# Remove rows with NAs
LCC_dataset <- na.omit(LCC_dataset) 
head(LCC_dataset)

sum(is.na(LCC_dataset))
sum(is.nan(as.matrix(LCC_dataset)))

str(LCC_dataset)
dim(LCC_dataset)
head(LCC_dataset)

# Prepare dataset
S2_bands <- LCC_dataset[,c(3:24)]
head(S2_bands)
dim(S2_bands)
y <- LCC_dataset$LCC
S2_bands2 <- cbind(S2_bands,y)
#---------------------------------------------------------
#create a regression task
#----------------------------------------------------
# Regression analysis using Entire dataset (benchmark)
#--------------------------------
# No Feature selection / Entire dataset task
optispect_best_task <- makeRegrTask(data = S2_bands2[c(1:23)], target = "y", id="Cv all features")
optispect_best_task
#---------------------------------------------------------
# set parallel back-end
require(parallel)
require(parallelMap) 
parallelStartSocket(cpus = detectCores())
#------------------------------------------------------------------
# -------------------Parameter tuning------------------------------
set.seed(165)
system.time(tune.model <- tuneParams(learner = lrn, 
                                     task = optispect_best_task, 
                                     resampling = rdesc, 
                                     measures = list(mlr::mape, mlr::rmse, mlr::rsq), 
                                     par.set = params, 
                                     control = ctrl, 
                                     show.info = T))
# Print the optimal hyper-parameters
print(tune.model$x)
# Print CV metrics of hyper-parameters
print(tune.model$y)

# Set the optimal hyper-parameters
lrn_tune.model <- setHyperPars(lrn, par.vals = tune.model$x)
#------------------------------------------------------------------
# -------------------Model Train and Evaluation--------------------
# Train  the model 
# Returns an object of class Wrapped Model which encapsulates the fitted model, i.e., the output of the underlying R learning method.
model.split.nofeat <- resample(learner = lrn_tune.model,task = optispect_best_task, resampling = rdesc, measures = list(mlr::mape, mlr::rmse, mlr::rsq))

print(model.split.nofeat)

#------------------------------------------------------------------
#------------------------- Accuracy metrics -----------------------
## Accuracy metrics (from resampled data)
#Create a summary table from the aggregate values
SummaryTbl_test_metrics <- t(as.data.frame(model.split.nofeat$aggr))
head(SummaryTbl_test_metrics)

SummaryTbl_test_metrics <- round(subset(SummaryTbl_test_metrics, select= c(rsq.test.mean, rmse.test.rmse, mape.test.mean)),2)
head(SummaryTbl_test_metrics)

colnames(SummaryTbl_test_metrics) <- c("R2cv", "RMSEcv", "MAPEcv")
head(SummaryTbl_test_metrics)
#----------------------------------------------------------
#--------------------- Scatter plot ----------------------
# Set theme to allow for plotmath expressions
#tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tt <- ttheme_minimal(base_size = 10, base_colour = "black", base_fill="darkolivegreen1",base_font="bold",base_family = "serif",core = list(fg_params=list(cex = 3.0), fontface="bold"), colhead = list(fg_params=list(cex = 2.0), fontface="bold", parse=TRUE), rowhead = list(fg_params=list(cex = 3.0, fontface="plain"), fontface="bold"))
#----------------------------------------------------------
tbl <- tableGrob(SummaryTbl_test_metrics, rows = "", theme = tt)
separators <- replicate(ncol(tbl) - 2,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)

tbl <- gtable::gtable_add_grob(tbl, grobs = separators,
                               t = 1, b = nrow(tbl), l = seq_len(ncol(tbl)-2)+2)

# Plot table
grid.arrange(tbl,
             nrow=2,
             as.table=TRUE,
             heights=c(2,1))

#----------------------------------------------------------
# Scatter Plot settings
regr_result <-model.split$pred$data
#----------------------------------------------------------
# Plot (approx. 1250*997 as Tiff)
p <- ggplot(regr_result, aes(x=truth, y=response)) + geom_point(size=5, colour="#728a25") + xlim(0, 70) +  ylim(0, 70) +
  scale_colour_hue(l=90) + # Use a slightly darker palette than normal
  geom_smooth(method=lm,   # Add linear regression lines
              se=F,    # Don't add shaded confidence region
              fullrange=TRUE, level=0.95, size=2, color="#0c6cc6")+ theme_bw()+ labs(x= bquote("Observed LCC "~(mu*g~cm^-2)), y=bquote("Observed LCC "~ (mu*g~cm^-2))) + theme(plot.title=element_text(size=12, face="bold"), axis.text.x=element_text(size=30, color = "black"), axis.text.y=element_text(size=30, color = "black"), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30),legend.title = element_text(size=30, face="bold"), legend.text = element_text(size=30), legend.justification=c(1,0), legend.position=c(1,0), text = element_text(family = "serif", face="bold")) + annotation_custom(tbl, xmin = 10, ymin = 55,xmax = 50)

p + geom_abline(intercept = 0, slope = 1, color = "#c6290c", size=2)
#----------------------------------------------------------

