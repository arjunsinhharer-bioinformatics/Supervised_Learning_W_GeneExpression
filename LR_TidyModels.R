## 1.) DOWNLOAD DATA & PACKAGES

library(GEOquery)
library(dplyr)
library(tibble)
library(tidymodels)

# Load the GEO dataset
gset <- getGEO("GSE11223", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Extract expression data
ex <- exprs(gset)

# Extract sample metadata
meta <- pData(gset)

# Extract platform annotation
platform <- annotation(gset); gpl <- getGEO(platform); platform_data <- Table(gpl)

#USE THE PLATFORM DATA TO ADD GENE_SYMBOLS AS ROW NAMES TO EXPRESSION DATA
# Convert 'ex' to a data frame and add row names as a column
ex_df <- ex %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID")  # Keep row names as 'ProbeID'

# Ensure both ID columns are character
ex_df <- ex_df %>%
  mutate(ID = as.character(ID))


platform_subset <- platform_data %>%
  dplyr::select(`ID`, `GENE_NAME`) %>%
  mutate(`ID` = as.character(`ID`))  # Convert ID to character


# Now perform the left join
annotated_ex <- ex_df %>%
  left_join(platform_subset, by = "ID")

# # Move 'GENE_NAME' to the first column
# annotated_ex <- annotated_ex %>%
#   select(`GENE_NAME`, everything())

annotated_ex <- annotated_ex %>%
  dplyr::select(GENE_NAME, everything())  # Move GENE_NAME to first position


annotated_ex <- annotated_ex %>% filter(!is.na(GENE_NAME))

annotated_ex <- annotated_ex %>%
  mutate(GENE_NAME = trimws(GENE_NAME)) %>%  # Remove extra spaces
  filter(GENE_NAME != "")  # Remove empty values

annotated_ex <- annotated_ex %>% drop_na()

annotated_ex_trim <- annotated_ex %>%
  dplyr::select(-c(ID, GENE_NAME))  # Removes ID and GENE_NAME, keeping only expression data



# Ensure that annotated_ex_trim and annotated_ex have matching rows
annotated_ex_trim <- annotated_ex_trim %>%
  mutate(GENE_NAME = annotated_ex$GENE_NAME[match(rownames(annotated_ex_trim), rownames(annotated_ex))])

# Remove any rows where GENE_NAME is missing
annotated_ex_trim <- annotated_ex_trim %>%
  filter(!is.na(GENE_NAME))

# Aggregate duplicate gene names by averaging their expression values
annotated_ex_trim <- annotated_ex_trim %>%
  group_by(GENE_NAME) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()

# Convert to a standard dataframe (not a tibble)
annotated_ex_trim <- as.data.frame(annotated_ex_trim)

# Assign GENE_NAME as row names
rownames(annotated_ex_trim) <- annotated_ex_trim$GENE_NAME

# Remove the GENE_NAME column as it's now the row names
annotated_ex_trim <- annotated_ex_trim[, !colnames(annotated_ex_trim) %in% "GENE_NAME"]

## 2.) FIND HIGHLY CORRELATED GENES (HCGs)

# Compute variance for each gene
gene_variances <- apply(annotated_ex_trim, 1, var, na.rm = TRUE)  # Variance per gene

# Select the top 1000 most variable genes
top_1000_genes <- names(sort(gene_variances, decreasing = TRUE)[1:1000])

# Filter the expression matrix to include only these genes
annotated_ex_trim_top1000 <- annotated_ex_trim[top_1000_genes, ]

# Transpose the matrix (samples as rows, genes as columns)
annotated_ex_trim_top1000_t <- t(annotated_ex_trim_top1000)

# Compute Spearman correlation
cor_matrix <- cor(annotated_ex_trim_top1000_t, method = "spearman", use = "pairwise.complete.obs")


# Convert the correlation matrix to a long format for easy filtering
cor_long <- as.data.frame(as.table(cor_matrix))

# Rename columns
colnames(cor_long) <- c("Gene1", "Gene2", "Correlation")

# Remove self-correlations (where Gene1 == Gene2)
cor_long <- cor_long %>%
  filter(Gene1 != Gene2)

# Find the most highly correlated genes (abs correlation >= 0.7)
highly_correlated_genes <- cor_long %>%
  filter(abs(Correlation) >= 0.7) %>%
  arrange(desc(abs(Correlation)))

## 3.) NORMALIZATION CHECK

#Check for log transformation
hist(as.numeric(unlist(annotated_ex_trim)), breaks = 50, main = "Distribution of Data", col = "blue")
all(annotated_ex_trim > 0, na.rm = TRUE)  # Log-transformed data should have only positive values

boxplot(annotated_ex_trim, main = "Boxplot of Samples", col = "lightblue")

## 4.) PREPARE OUTCOME DATA FOR TIDYMODELS
outcome <- meta %>%
  dplyr::select(`inflammation_status:ch1`) %>%
  dplyr::rename(Status = `inflammation_status:ch1`) %>%
  dplyr::mutate(Binary = ifelse(Status == "Inflamed", 1, 0))

# Remove duplicate correlation values
highly_correlated_genes_trimmed <- highly_correlated_genes %>%
  distinct(Correlation, .keep_all = TRUE)

# Convert rownames to a column and select only Sample_ID and Binary
outcome_fixed <- outcome %>%
  rownames_to_column(var = "Sample_ID") %>%  # Convert rownames to a new column
  select(Sample_ID, Binary)  # Keep only relevant columns


## 5.) Find the Top Gene
library(limma)

# Create model matrix for the binary condition
design <- model.matrix(~ Binary, data = outcome_fixed)

# Run limma analysis
fit <- lmFit(annotated_ex_trim, design)
fit <- eBayes(fit)

# Get top differentially expressed gene
top_de_gene <- rownames(topTable(fit, coef = 2, number = 1, sort.by = "p"))

cat("Most differentially expressed gene:", top_de_gene, "\n")

## 6.) WORKSPACE ORGANIZATION
#CODE TO DE-CLUTTER THE WORKSPACE
# List all objects in the environment
all_objects <- ls()

# Specify the datasets to keep
datasets_to_keep <- c("highly_correlated_genes_trimmed", "annotated_ex_trim", "outcome_fixed", "top_de_gene")

# Remove all objects except for those specified in datasets_to_keep
rm(list = all_objects[!all_objects %in% datasets_to_keep])

## 7.) LOGISTIC REGRESSION IMPLEMENTATION

# Extract expression values for "dual oxidase 2"
top_gene_expression <- annotated_ex_trim %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  filter(Gene == "dual oxidase 2") %>%
  column_to_rownames(var = "Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID") %>%
  rename(DUOX2_Expression = "dual oxidase 2")  # Rename for clarity

# Merge with outcome data
logistic_data <- outcome_fixed %>%
  inner_join(top_gene_expression, by = "Sample_ID")

# Convert to tibble for tidymodels
logistic_data <- as_tibble(logistic_data)

# Print structure to confirm
print(logistic_data)

### Split the Data

set.seed(123)
data_split <- initial_split(logistic_data, prop = 0.8, strata = Binary)
train_data <- training(data_split)
test_data <- testing(data_split)

train_data <- train_data %>%
  mutate(Binary = as.factor(Binary))

test_data <- test_data %>%
  mutate(Binary = as.factor(Binary))


# Verify structure before proceeding
print(colnames(train_data))

logistic_recipe <- recipe(Binary ~ DUOX2_Expression, data = train_data) %>%
  step_normalize(all_predictors())  # Normalize DUOX2 expression

logistic_model <- logistic_reg(mode = "classification") %>%
  set_engine("glm")

logistic_workflow <- workflow() %>%
  add_model(logistic_model) %>%
  add_recipe(logistic_recipe)

# Define and Train the model
logistic_fit <- logistic_workflow %>%
  fit(data = train_data)

### Make predictions on the test set
predictions <- logistic_fit %>%
  predict(new_data = test_data) %>%
  bind_cols(test_data)

# Print some predictions
print(head(predictions))

### Evaluate Model Performance
# Get predicted probabilities instead of just class
predictions <- logistic_fit %>% 
  predict(new_data = test_data, type = "prob") %>% 
  bind_cols(test_data)  # Merge with test data to retain Binary column

# Verify column names
colnames(predictions)

# Ensure Binary is a factor
predictions <- predictions %>%
  mutate(Binary = as.factor(Binary))

# Compute AUC
auc <- roc_auc(predictions,
               truth = Binary, 
               .pred_0,  # Make sure this column exists, change if needed
               event_level = "first")$.estimate

# Print AUC
print(auc)

### Confusion Matrix
predictions <- predictions %>%
  mutate(.pred_class = ifelse(.pred_1 > 0.5, 1, 0)) %>%
  mutate(.pred_class = as.factor(.pred_class)) # Ensure it's a factor for classification

conf_mat(predictions, truth = Binary, estimate = .pred_class)

# Convert Binary to factor for classification
predictions <- predictions %>%
  mutate(Binary = as.factor(Binary))

# Define a set of classification metrics
classification_metrics <- metric_set(accuracy, f_meas, spec, sens, npv, ppv)

# Apply the metrics
classification_metrics(predictions, truth = Binary, estimate = .pred_class)

### Confusion Matrix
predictions <- predictions %>%
  mutate(.pred_class = ifelse(.pred_1 > 0.5, 1, 0)) %>%
  mutate(.pred_class = as.factor(.pred_class)) # Ensure it's a factor for classification

conf_mat(predictions, truth = Binary, estimate = .pred_class)

# Convert Binary to factor for classification
predictions <- predictions %>%
  mutate(Binary = as.factor(Binary))

# Define a set of classification metrics
classification_metrics <- metric_set(accuracy, f_meas, spec, sens, npv, ppv)

# Apply the metrics
classification_metrics(predictions, truth = Binary, estimate = .pred_class)

### Compute AUC and Generate ROC
# Compute the AUC (Area Under the Curve)
auc <- roc_auc(predictions, truth = Binary, .pred_0, event_level = "first")$.estimate

# Generate ROC Curve
g_roc <- predictions %>%
  roc_curve(truth = Binary, .pred_0, event_level = "first") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(color = "red") +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  annotate(geom = "text", x = 0.75, y = 0.1, label = paste0("AUC: ", round(auc, 3)), color = "red") +
  ggtitle("Logistic Regression ROC Curve")

# Print ROC plot
print(g_roc)


