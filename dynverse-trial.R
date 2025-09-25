####trajectory analysis
install.packages("devtools")
devtools::install_github("dynverse/dyno")
library(dyno)
library(dynwrap)
library(dynmethods)
library(Seurat)
library(reticulate)

# Tell reticulate to use the virtualenv
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)

# Check that Scanpy is available
scanpy <- import("scanpy")
scanpy$pp$neighbors  # should return the function object

load("2w-stroma-subset_seurat.RData")
View(subset_seurat)



# Extract counts and normalized expression

dataset <- wrap_expression(
  counts = GetAssayData(subset_seurat, slot = "counts"),
  expression = GetAssayData(subset_seurat, slot = "data")
)

# Define method: PAGA
model <- ti_paga()


guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected

# Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = TRUE, 
  expect_topology = NULL, 
  expected_topology = NULL, 
  n_cells = 16613, 
  n_features = 18783, 
  memory = "100GB", 
  prior_information = c("start_id", "end_id", "end_n", "start_n", "leaves_n", "groups_n", "features_id", "dimred"), 
  docker = FALSE
)
guidelines <- dynguidelines::guidelines(answers = answers) 