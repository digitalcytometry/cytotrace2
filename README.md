<p align="center">
  <img width="500" src="images/logo.jpg"> 
</p>

<h2> <p align="center">
      CytoTRACE 2 - Version 1.1.0 Released
</p> </h2>

We are thrilled to introduce **CytoTRACE 2 version 1.1.0**, packed with significant performance enhancements to elevate your single-cell transcriptomic analyses. Here's what's new in this release:

### 🔍 Major Updates and Enhancements

- **Retrained CytoTRACE 2 Framework**  
  The CytoTRACE 2 model has been retrained, yielding additional performance gains in granular potency prediction and enhancing cross-platform robustness.

- **Expanded Ensemble Model**  
  The ensemble now comprises **19 models** instead of 17, improving the predictive power and stability of the framework.

- **Background Expression Matrix**  
  Introduced a background expression matrix generated during training for improved regularization.

- **Enhanced Data Representations**  
  Added Log2-adjusted representation of the input expression data to be used for prediction on top of ranked expression profiles, to capture detailed transcriptomic signals. *This changes the requirement for the input expression data to contain only raw or CPM/TPM normalized counts.*

- **Adaptive Nearest Neighbor Smoothing**  
  Modified the KNN smoothing step to employ an adaptive nearest neighbor smoothing strategy.

### 💻 Codebase and Distribution Updates

- **Codebase Updates**  
  - Updated both R and Python package codebases to reflect all the above changes.
  - Optimized for **time and memory efficiency**, ensuring faster computations and scalability.

- **Enhanced Python Package Distribution**  
  The Python version of CytoTRACE 2 is now available on [PyPI](https://pypi.org/project/cytotrace2-py/), making installation easier for Python users.

### 📚 Documentation and Guides

- Updated Vignettes to align with the new model features and usage instructions.
- Refreshed README with new information, detailed explanations, and FAQ items tailored to the new framework.

---

We deeply appreciate the contributions from our community that made this release possible. Thank you for your continued support! 🙏

<h2> <p align="center">
      Prediction of absolute developmental potential <br> using  single-cell expression data
</p> </h2>

**CytoTRACE 2** is a computational method for predicting cellular potency categories and absolute developmental potential from single-cell RNA-sequencing data. 

Potency categories in the context of CytoTRACE 2 classify cells based on their developmental potential, ranging from totipotent and pluripotent cells with broad differentiation potential to lineage-restricted oligopotent, multipotent and unipotent cells capable of producing varying numbers of downstream cell types, and finally, differentiated cells, ranging from mature to terminally differentiated phenotypes.

The predicted potency scores additionally provide a continuous measure of developmental potential, ranging from 0 (differentiated) to 1 (totipotent).

Underlying this method is a novel, interpretable deep learning framework trained and validated across 34 human and mouse scRNA-seq datasets encompassing 24 tissue types, collectively spanning the developmental spectrum. 

This framework learns multivariate gene expression programs for each potency category and calibrates outputs across the full range of cellular ontogeny, facilitating direct cross-dataset comparison of developmental potential in an absolute space. 

<p align="center">
    <img width="900" src="images/schematic.png">
</p>

This documentation page details the R package for applying CytoTRACE 2. <strong> For the python package, see <a href="/cytotrace2_python" target="_blank">CytoTRACE 2 Python</a>.</strong>


## Installation


  We recommend installing the CytoTRACE 2 package using the **devtools** package from the R console. If you do not have **devtools** installed, you can install it by running ```install.packages("devtools")``` in the R console.

```R
  devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") #installing
  library(CytoTRACE2) #loading
```

See alternative installation and package management methods (including an easy-to-use conda environment that precisely solves all dependencies) in the [__Advanced options__](#advanced-options) section below.

The installation of the CytoTRACE 2 package itself typically takes about one minute on a standard computer. Optional installation of the provided conda environment generally takes 5-10 minutes but can vary substantially, sometimes requiring up to an hour depending on system and conda version.

NOTE: We recommend using Seurat v4 or later for full compatibility with CytoTRACE 2 package. If you don't have Seurat installed, you can install it by running ```install.packages("Seurat")``` in the R console prior to installing CytoTRACE 2 or use the provided conda environment.

<details><summary>Dependencies</summary>

The following list includes the versions of packages used during the development of CytoTRACE 2. While CytoTRACE 2 is compatible with various (older and newer) versions of these packages, it's important to acknowledge that specific combinations of dependency versions can lead to conflicts. The only such conflict known at this time happens when using Seurat v4 in conjunction with Matrix v1.6. This issue can be resolved by either upgrading Seurat or downgrading Matrix.

``` bash
    R (4.2.3)
    data.table (1.14.8)
    doParallel (1.0.17)
    dplyr (1.1.3)
    ggplot2 (3.4.4)
    HiClimR (2.2.1)
    magrittr (2.0.3)
    Matrix (1.5-4.1)
    parallel (4.2.3)
    plyr (1.8.9)
    RANN (2.6.1)
    Rfast (2.0.8)
    RSpectra (0.16.1)
    Seurat (4.3.0.1)
    SeuratObject (4.1.3)
    stringr (1.5.1)
  ```

   </details>


## Running CytoTRACE 2 

Running ```CytoTRACE 2``` is easy and straightforward. After loading the library, simply execute the [cytotrace2()](#cyto-trace2) function, with one required input, [__expression data__](#input-files), to obtain [__potency score and potency category__](#cytotrace-2-outputs)  predictions. Subsequently, running [plotData](#plot-data) will generate informative visualizations based on the predicted values, and external annotations, if available. Below, find two vignettes showcasing the application on a mouse dataset and a human dataset.


<details open><summary><span style="font-size: 15px;"><strong>Vignette 1: Development of mouse pancreatic epithelial cells (data table input, ~2 minutes)</strong></span></summary>

To illustrate use of CytoTRACE 2 with a mouse dataset, we will use the dataset Pancreas_10x_downsampled.rds, originally from [Bastidas-Ponce et al., 2019](https://doi.org/10.1242/dev.173849), filtered to cells with known ground truth developmental potential and downsampled, available to download [here](https://drive.google.com/uc?export=download&id=1TYdQsMoDIJjoeuiTD5EO_kZgNJUyfRY2), containing 2 objects:
- expression_data: gene expression matrix for a scRNA-seq (10x Chromium) dataset encompassing 2850 cells from murine pancreatic epithelium
- annotation: phenotype annotations for the scRNA-seq dataset above. 

After downloading the .rds file, we apply CytoTRACE 2 to this dataset as follows:

```r
# load the CytoTRACE 2 package
library(CytoTRACE2) 

# download the .rds file (this will download the file to your working directory)
download.file("https://drive.google.com/uc?export=download&id=1TYdQsMoDIJjoeuiTD5EO_kZgNJUyfRY2", "Pancreas_10x_downsampled.rds")

# load rds
data <- readRDS("Pancreas_10x_downsampled.rds")

# extract expression data
expression_data <- data$expression_data

# running CytoTRACE 2 main function - cytotrace2 - with default parameters
cytotrace2_result <- cytotrace2(expression_data)

# extract annotation data
annotation <- data$annotation

# generate prediction and phenotype association plots with plotData function
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation,
                  expression_data = expression_data
                  )

```

Expected prediction output, dataframe ```cytotrace2_result``` looks as shown below (can be downloaded from [here](./Vignette1_CytoTRACE2_results.csv)):

<p align="center">
    <img width="600" src="images/Vignette1_predictions.png">
</p>



<br>

This dataset contains cells from 4 different embryonic stages of a murine pancreas, and has the following cell types present:
- Multipotent pancreatic progenitors
- Endocrine progenitors and precursors
- Immature endocrine cells
- Alpha, Beta, Delta, and Epsilon cells

<p align="center">
        <img width="600" src="images/Vignette1_phenotype_umap.png">
</p>

Each of these cell types is at a different stage of development, with progenitors and precursors having varying potential to differentiate into other cell types, and mature cells having no potential for further development. We use CytoTRACE 2 to predict the absolute developmental potential of each cell, which we term as "potency score", as a continuous value ranging from 0 (differentiated) to 1 (stem cells capable of generating an entire multicellular organism). The discrete potency categories that the potency scores cover are ```Differentiated```, ```Unipotent```, ```Oligopotent```, ```Multipotent```, ```Pluripotent```, and ```Totipotent```.

In this case, we would expect to see:
- close to 0 potency scores alpha, beta, delta, and epsilon cells as those are known to be differentiated, 
- scores in the higher mid-range for multipotent pancreatic progenitors as those are known to be multipotent, 
- for endocrine progenitors, precursors and immature cells, the ground truth is not unique, but is in the range for unipotent category. So we would expect to see scores in the lower range for these cells, closer to differentiated.

Visualizing the results we can directly compare the predicted potency scores with the known developmental stage of the cells, seeing how the predictions meticulously align with the known biology. Take a look!

- ***Potency score vs. Ground truth*** 
 <br> UMAP embedding of predicted absolute potency score, which is a continuous value ranging from 0 (differentiated) to 1 (totipotent), binned into discrete potency categories as described [__below__](#cytotrace-2-outputs), indicating the absolute developmental potential of each cell. <br>
  ```bash
  plots$CytoTRACE2_UMAP
  ```

<div align="center">
  <div style="display: flex;">
    <img width="400" src="images/Vignette1_potency_score_umap.png">
    <img width="400" src="images/Vignette1_ground_truth_umap_with_pheno.png">
  </div>
</div>

<br>



 - <details> <summary> <strong>Other output plots</strong> </summary>
    
    - ***Potency score distribution by phenotype***
    <br> A boxplot of predicted potency score separated by phenotype/group from the annotation file. Can be used to assess the distribution of predicted potency scores across different cell phenotypes. <br>
      ```bash
      plots$CytoTRACE2_Boxplot_byPheno
      ```

      <p align="center">
        <img width="600" src="images/Vignette1_potencyBoxplot_byPheno.png">
      </p>


    - ***Potency category***
    <br> The UMAP embedding plot of predicted potency category reflects the discrete classification of cells into potency categories, taking possible values of ```Differentiated```, ```Unipotent```, ```Oligopotent```, ```Multipotent```, ```Pluripotent```, and ```Totipotent```. <br>
      ```bash
      plots$CytoTRACE2_Potency_UMAP
      ```
      <p align="center">
        <img width="600" src="images/Vignette1_potency_category_umap.png">
      </p>

    - ***Relative order***
    <br> UMAP embedding of predicted relative order, which is based on absolute predicted potency scores normalized to the range 0 (more differentiated) to 1 (less differentiated). Provides the relative ordering of cells by developmental potential <br>
      ```bash
      plots$CytoTRACE2_Relative_UMAP
      ```
      <p align="center">
        <img width="600" src="images/Vignette1_rel_order_umap.png">
      </p>

    - ***Phenotypes***
    <br> UMAP colored by phenotype annotation. Used to assess the distribution of cell phenotypes across the UMAP space. <br>
      ```bash
      plots$Phenotype_UMAP
      ```
      <p align="center">
        <img width="600" src="images/Vignette1_phenotype_umap.png">
      </p>
</details>

</details>

***Additional specifications***

- By default, CytoTRACE 2 expects mouse data. To provide human data, users should specify ```species = "human"```
- When running on computers with less than 16GB memory, we recommend reducing ```ncores``` to 1 or 2 to avoid memory issues.
- To pass a loaded Seurat object or an .rds file path containing a Seurat object as gene expression input to the ```cytotrace2()``` function, users should specify ```is_seurat = TRUE``` and ```slot_type``` as the name of the assay slot containing the gene expression matrix to use for prediction (can be either ```counts``` or ```data```; default is ```counts```). This will return a Seurat object with metadata containing all CytoTRACE 2 cell potency predictions, which can be further passed to the ```plotData()``` function for visualization with ```is_seurat = TRUE```, as shown in the [__vignette__](#vignette-2:-human-cord-blood-mononuclear-cells-(seurat-object-input,-~2-minutes)) below. 
- If the passed Seurat object already contains PCA and UMAP embeddings, the ```plotData()``` function will use these embeddings for visualization. Otherwise, it will generate new PCA and UMAP embeddings internally.
- **NOTE**: To reproduce the results in the manuscript, use the following parameters: <br>
      ```parallelize_models = TRUE```  <br>
      ```parallelize_smoothing = TRUE```  <br>
      ```batch_size = 100000```  <br>
      ```smooth_batch_size = 10000```  <br>

More details on expected function input files and output objects can be found in [__Input Files__](#input-files) and [__CytoTRACE 2 outputs__](#cytotrace-2-outputs) sections below.

For full usage details with additional options, see the [__Extended usage details__](#extended-usage-details) section below.



<details><summary><span style="font-size: 15px;"><strong>Vignette 2: Human cord blood mononuclear cells (Seurat object input, ~2 minutes) </strong> </span></summary>

To illustrate use of CytoTRACE 2 with a human dataset stored in a Seurat object, we will use the Seurat object Cord_blood_CITE_seq_SeuratObj.rds, originally from [Stoeckius et al., 2017](https://doi.org/10.1038/nmeth.4380), filtered to cells with known ground truth developmental potential and available to download [here](https://drive.google.com/uc?export=download&id=1_OhZdz4y0R0MeB6gNlYAFC8P8rutgqzJ), containing:
- gene expression data for a scRNA-seq (CITE-seq) dataset comprising of 2308 cells. Stored in the RNA assay, 'counts' slot of the Seurat object.
- phenotype annotations for the scRNA-seq dataset above. Stored in the metadata of the Seurat object.

After downloading the .rds file, we apply CytoTRACE 2 to this dataset as follows:


```r
# load the CytoTRACE 2 package
library(CytoTRACE2) 

# download the .rds file (this will download the file to your working directory)
download.file("https://drive.google.com/uc?export=download&id=1_OhZdz4y0R0MeB6gNlYAFC8P8rutgqzJ", "Cord_blood_CITE_seq_downsampled_SeuratObj.rds")# load the Seurat object

seurat_obj <- readRDS("Cord_blood_CITE_seq_downsampled_SeuratObj.rds")

# running CytoTRACE 2 main function - cytotrace2
# setting is_seurat = TRUE as input is a Seurat object
# setting slot_type = "counts" as the gene expression data is stored in the "counts" slot of the Seurat object
# setting species = 'human' 
cytotrace2_result <- cytotrace2(seurat_obj, is_seurat = TRUE, slot_type = "counts", species = 'human')

# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = seurat_obj@meta.data$standardized_phenotype) %>% set_rownames(., colnames(seurat_obj))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = TRUE)

```
Expected prediction output ```cytotrace2_result``` is your input Seurat object with predictions added to the metadata.  The predictions in the metadata look as shown below (the resulting Seurat object can be downloaded from [here](./Vignette2_CytoTRACE2_results.rds)):

<p align="center">
    <img width="600" src="images/Vignette2_predictions.png">
</p>

Expected plotting outputs: 

- ***Potency score***
<br> UMAP embedding of predicted absolute potency score, which is a continuous value ranging from 0 (differentiated) to 1 (totipotent), indicating the absolute developmental potential of each cell. As expected, hematopoietic stem and progenitor cells have much higher potency compared to other phenotypes, which are mostly mature differentiated cells and are predicted as such.
    ```bash
      plots$CytoTRACE2_UMAP
    ```

  <p align="center">
    <img width="600" src="images/Vignette2_potency_score_umap.png">
  </p>

- ***Potency category*** 
<br> UMAP embedding of predicted potency category, reflecting the discrete classification of cells into potency categories, taking possible values of ```Differentiated```, ```Unipotent```, ```Oligopotent```, ```Multipotent```, ```Pluripotent```, and ```Totipotent```.

  ```bash
  plots$CytoTRACE2_Potency_UMAP
  ```
  <p align="center">
    <img width="600" src="images/Vignette2_potency_category_umap.png">
  </p>

- ***Relative order***
<br> UMAP embedding of predicted relative order, which is based on absolute predicted potency scores normalized to the range 0 (more differentiated) to 1 (less differentiated). Provides the relative ordering of cells by developmental potential.
  ```bash
  plots$CytoTRACE2_Relative_UMAP
  ```
  <p align="center">
    <img width="600" src="images/Vignette2_rel_order_umap.png">
  </p>

- ***Phenotypes***
<br> UMAP colored by phenotype annotation. Used to assess the distribution of cell phenotypes across the UMAP space.
  ```bash
  plots$Phenotype_UMAP
  ```
  <p align="center">
    <img width="600" src="images/Vignette2_phenotype_umap.png">
  </p>

- ***Potency score distribution by phenotype***
<br> A boxplot of predicted potency score separated by phenotype/group from the annotation file. Can be used to assess the distribution of predicted potency scores across different cell phenotypes.
  ```bash
  plots$CytoTRACE2_Boxplot_byPheno
  ```
  <p align="center">
    <img width="600" src="images/Vignette2_potencyBoxplot_byPheno.png">
  </p>

</details>


## Input files

<details><summary>Expand section</summary>
  
CytoTRACE 2 requires a single-cell RNA-sequencing gene expression object as input. This can include either raw counts or CPM/TPM normalized counts, and should not be log-transformed.
This input can be provided in any of the following formats:

1. **A data table**. This object of class data.frame or of another class which allows storing row and column names. This should have genes as rows and cells as columns, with row and column names set accordingly.
2. **A filepath to a tab-delimited file**. This file should contain a gene expression matrix with genes as rows and cells as columns. The first row must contain the cell IDs (header), and the first column must contain the gene names that can have a column name or an empty header.
3. **A Seurat object**. This object should contain gene expression values to be used, stored in `object[["RNA"]]@slot_type`, where `slot_type` is the name of the assay slot containing the gene expression matrix to use for prediction (can be either `counts` or `data`). 
4. **A filepath to an .rds file containing a Seurat object**. This file should contain a Seurat object as described in `option 3`.

Please make sure that the input does not contain duplicate gene names or cell IDs.

<p align="center">
    <img width="600" src="images/data.png">
</p>

CytoTRACE 2 also accepts cell phenotype annotations as an optional input. This is a non-required argument used for plotting outputs, and if provided, will be used to generate UMAP and predicted potency score box plots by phenotype. It should be passed as a loaded dataframe containing cell IDs as rownames and phenotype labels (string formatted without special characters) in the first column. Any additional columns will be ignored.

<p align="center">
    <img width="600" src="images/annotation.png">
</p>

</details>

## CytoTRACE 2 outputs

<details><summary>Expand section</summary>

### ```cytotrace2()```

```cytotrace2()``` function returns the CytoTRACE 2 cell potency predictions in a data.frame format. If the gene expression input was provided as a Seurat object or as a filepath to a Seurat object, the results will be returned in the form of a Seurat object with metadata elements containing all CytoTRACE 2 cell potency predictions.

#### CytoTRACE 2 cell potency predictions

For each cell, the CytoTRACE 2 predictions include:

1. *CytoTRACE2_Score*: The final predicted cellular potency score following postprocessing. Possible values are real numbers ranging from 0 (differentiated) to 1 (totipotent), which are binned into potency categories according to the following ranges:
    <div style="text-align: center;">
        <table style="margin-left: auto; margin-right: auto;">
            <tr>
                <td>Range</td>
                <td>Potency</td>
            </tr>
            <tr>
                <td>0 to 1/6</td>
                <td>Differentiated</td>
            </tr>
            <tr>
                <td>1/6 to 2/6</td>
                <td>Unipotent</td>
            </tr>
            <tr>
                <td>2/6 to 3/6</td>
                <td>Oligopotent</td>
            </tr>
            <tr>
                <td>3/6 to 4/6</td>
                <td>Multipotent</td>
            </tr>
            <tr>
                <td>4/6 to 5/6</td>
                <td>Pluripotent</td>
            </tr>
            <tr>
                <td>5/6 to 1</td>
                <td>Totipotent</td>
            </tr>
        </table>
    </div>



2. *CytoTRACE2_Potency*: The final predicted cellular potency category following postprocessing. Possible values are ```Differentiated```, ```Unipotent```, ```Oligopotent```, ```Multipotent```, ```Pluripotent```, and ```Totipotent```. 
3. *CytoTRACE2_Relative*: The predicted relative order of the cell, based on the absolute predicted potency scores, normalized to the range [0,1] (0 being most differentiated, 1 being least differentiated).
4. *preKNN_CytoTRACE2_Score*: The cellular potency score predicted by the CytoTRACE 2 model before KNN smoothing (See ‘Postprocessing’ in the Methods section of the manuscript).
5. *preKNN_CytoTRACE2_Potency*: The cellular potency category  predicted by the CytoTRACE 2 model before KNN smoothing (See ‘Postprocessing’ in the Methods section of the manuscript). Possible values are ```Differentiated```, ```Unipotent```, ```Oligopotent```, ```Multipotent```, ```Pluripotent```, and ```Totipotent```.


For more details about postprocessing, see the [__Under the hood__](#under-the-hood) section below.

### ```plotData()```

The results of ```plotData()``` are returned as a named list of plots. To access any of the plots simply run ```plots$name_of_the_plot``` (i.e. plots$CytoTRACE2_UMAP).

By default, ```plotData``` produces three plots depicting the UMAP embedding of the input single-cell gene expression data, each colored according to a CytoTRACE 2 output prediction type.

- **Potency category UMAP**: a UMAP colored by predicted potency category (named *CytoTRACE2_UMAP* in the output list)
- **Potency score UMAP**: a UMAP colored by predicted potency score (named *CytoTRACE2_Potency_UMAP* in the output list)
- **Relative order UMAP**: a UMAP colored by predicted relative order (named *CytoTRACE2_Relative_UMAP* in the output list)

If a phenotype annotation file is provided, ```plotData``` will return two more plots.

- **Phenotype UMAP**: a UMAP colored by phenotype annotation (named *Phenotype_UMAP* in the output list)
- **Phenotype potency box plot**: a boxplot of predicted potency score separated by phenotype/group the annotation file (named *CytoTRACE2_Boxplot_byPheno* in the output list)

If the input is a Seurat object containing predictions, with ```is_seurat = TRUE``` and ```slot_type``` argument properly specified ("counts" or "data"), if it contains PCA and UMAP embeddings (named "pca" and "umap" if `seurat_object@reductions``), the function will use these embeddings for visualization. Otherwise, it will generate new PCA and UMAP embeddings internally.

</details>



## Extended usage details

<details id="cyto-trace2"><summary>cytotrace2</summary>


Required input:

- *input*: Single-cell RNA-seq data, which can be an expression matrix
(rows as genes, columns as cells), a filepath to a tab-separated file containing the data table,
a Seurat object, or the filepath to an .rds file containing a Seurat object.

Optional arguments:

- *species*: String indicating the species name for the gene names in the input data
(options: **"human"** or **"mouse"**, default is **"mouse"**).
- *is_seurat*: Logical indicating whether the input is a Seurat object/filepath to a Seurat object (default is **FALSE**).
- *slot_type*: Character indicating the type of slot to access from "RNA" assay if provided is a Seurat object & is_seurat = TRUE (options: **counts** or **data**, default is  **counts**)
- *parallelize_smoothing*: Logical indicating whether to run the smoothing function
on subsamples in parallel on multiple threads (default is **TRUE**).
- *parallelize_models*: Logical indicating whether to run the prediction function on
models in parallel on multiple threads (default is **TRUE**).
- *ncores*: Integer indicating the number of cores to utilize when parallelize_models
and/or parallelize_smoothing are TRUE (default is **NULL**; the pipeline detects the number of
available cores and runs on half of them; for Windows, it will be set to 1; when running on computers with less than 16GB memory, we recommend reducing it to 1 or 2 to avoid memory issues).
- *batch_size*: Integer or NULL indicating the number of cells to process at once, including subsampling for KNN smoothing.
No subsampling if NULL (default is **10000**; recommended for input data size > 10K cells).
- *smooth_batch_size*: Integer or NULL indicating the number of cells to subsample further
within the batch_size for the smoothing by diffusion step of the pipeline. No subsampling if NULL
(default is **1000**; recommended for input data size > 1K cells).
- *seed*: Integer specifying the seed for reproducibility in random processes (default is **14**).

Information about these arguments is also available in the function's manual, which can be accessed by running ```help(cytotrace2)``` within an R session.

A typical snippet to run the function with full argument specification on a file path containing human data: 

```r
cytotrace2_result <- cytotrace2("path/to/input/expression_file.tsv",
                       species = "human",
                       is_seurat = FALSE,
                       batch_size = 10000,
                       smooth_batch_size = 1000,
                       parallelize_models = TRUE,
                       parallelize_smoothing = TRUE,
                       ncores = NULL,
                       seed = 14)               
```
A typical snippet to run the function with full argument specification on a loaded Seurat object of mouse data: 

```r
seurat_obj <- loadData("path/to/input/Seurat_object.rds")
cytotrace2_result <- cytotrace2(seurat_obj,   
                  species = "mouse",
                  is_seurat = TRUE,
                  slot_type = "counts",
                  batch_size = 10000,
                  smooth_batch_size = 1000,
                  parallelize_models = TRUE,
                  parallelize_smoothing = TRUE,
                  ncores = NULL,
                  seed = 14)  

```
</details>



<details id="plot-data"><summary>plotData</summary>

Required input:
- *cytotrace2_result*: Output of the cytotrace2 function as described in [_CytoTRACE 2 outputs_](#cytotrace-2-outputs).

Optional input:
- *annotation*: The annotation data (optional, default is ***NULL***, used to prepare phenotype UMAP and box plots)
- *expression_data*: The expression data to be used for plotting (default is ***NULL***, if cytotrace2 is a Seurat object containing expression data and is_seurat = TRUE, can be left NULL).
- *pc_dims*: The number of principal components to use for UMAP visualization (default is ***30***).
- *is_seurat*: Logical, indicating whether the input is a Seurat object (default is ***FALSE***). If `TRUE` and the input Seurat object contains PCA and UMAP embeddings, those will be automatically used for visualization, otherwise, PCA and UMAP will be calculated.

- *slot_type*: Character indicating the type of slot to access from "RNA" assay if provided is a Seurat object & is_seurat = TRUE (options: ***counts*** or ***data***, default is  ***counts***).
- *seed*: Integer specifying the seed for reproducibility in random processes (default is **14**).

A typical snippet to run the function with full argument specification following CytoTRACE 2 prediction: 

```r
annotation <- read.csv('path/to/input/annotation_file.tsv', header = TRUE, row.names = 1)
data <- read.csv('path/to/input/expression_file.tsv', sep = "\t", header = TRUE, row.names = 1)
cytotrace2_result <- cytotrace2(data)
plots <- plotData(cytotrace2_result,
                  annotation = annotation,
                  expression_data = data,
                  is_seurat = FALSE,
                  pc_dims = 30,
                  seed = 14)
```
</details>

## Under the hood 

<details><summary>Expand section</summary> 
  
Underlying CytoTRACE 2 is a novel deep learning framework designed to handle the complexities of single-cell potency assessment while achieving direct biological interpretability. The core of this framework is a set of Gene Set Binary Network (GSBN) modules, in which binary neural networks learn gene sets associated with each potency category. This network was trained over 19 datasets from 16 diverse human and mouse tissues, and the package here relies on an ensemble of these per-dataset trained models. 
<p align="center">
    <img width="700" src="images/BNN_schematic.png">
</p>
Following initial prediction by the core model, CytoTRACE 2 implements a postprocessing step to leverage the information across transcriptionally similar cells to smooth potency score and correct potency category outliers using a combination of Markov diffusion and k-nearest neighbor smoothing. 
<!-- For more details about the CytoTRACE 2 method, please see the [_associated publication_](#Citation). -->

</details>

## Advanced options

<details><summary>Manual download and installation</summary>

1. You can clone the package folder from the GitHub repository and do local installation using the following commands:
``` bash
  git clone https://github.com/digitalcytometry/cytotrace2.git cytotrace2 #run in terminal to clone the repository
``` 
2. Navigate to `cytotrace2` directory:
    ```bash
      cd cytotrace2
    ```
3. (optional) To ensure your environment contains all the necessary dependencies, you can use the provided conda environment file to create a new environment with all the necessary packages. 
    <details><summary>Expand section</summary>

    1. Install <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank">Miniconda</a> if not already available.

    2. (5-10 minutes) Create a conda environment with the required dependencies:
    ```bash
      conda env create -f environment_R.yml
    ```

    3. Activate the `cytotrace2` environment you just created:
    ```bash
      conda activate cytotrace2
    ``` 
   </details>
4. (Activate R from terminal) Install and load the package into your R environment using the following commands:
```R
  devtools::install_local("./cytotrace2_r") #run to install the package locally (done once)
  library(CytoTRACE2) #run to load the package into the current R session
```

or

```R
  devtools::load_all("./cytotrace2_r") #to load the package into the current R session
```

Make sure you specify the path to the subdirectory "cytotrace2_r" within the cloned repository, if you're running outside the cloned repository.

Now you can use the package as described in the [__Running CytoTRACE 2__](#running-cytotrace-2) section.

NOTE: If you are running on a M1/M2 Mac, you may get an error solving the conda environment triggered by issues installing some of the dependencies from conda. To fix that, please try the following:
```bash
    conda create -n cytotrace2
    conda activate cytotrace2
    conda config --env --set subdir osx-64
    conda env update --file environment_R.yml
```

</details>


<details><summary>Updating your local installation</summary>

If you have made local updates to your version of the CytoTRACE 2 source code, you should execute 
```bash
  devtools::install_local("./cytotrace2_r") #to re-install the package locally (done once)
  library(CytoTRACE2) #to load the updated package into the current R session
```
or 

```bash
  devtools::load_all("./cytotrace2_r") #to load the updated package into the current R session
``` 
in the package folder, before running. 

</details>

## Frequently asked questions

<details><summary>Expand section</summary>

1. **What are the CytoTRACE 2 potency categories?**
CytoTRACE 2 classifies cells into six potency categories:

  - **Totipotent**: Stem cells capable of generating an entire multicellular organism
  - **Pluripotent**: Stem cells with the capacity to differentiate into all adult cell types
  - **Multipotent**: Lineage-restricted multipotent cells capable of producing >3 downstream cell types
  - **Oligopotent**: Lineage-restricted immature cells capable of producing 2-3 downstream cell types
  - **Unipotent**: Lineage-restricted immature cells capable of producing a single downstream cell type
  - **Differentiated**: Mature cells, including cells with no developmental potential
  
2. **What organism can my data be from?**

    CytoTRACE 2 was developed over mouse and human data, and this package accepts data from either. If human data is provided (with ```species = 'human'``` specified), the algorithm will automatically perform an orthology mapping to convert human genes to mouse genes for the CytoTRACE 2 feature set. 

3. **Should I normalize the data before running the main function?**

    There is no need to normalize data prior to running CytoTRACE 2, provided there are no missing values and all values are non-negative. The input needs to be **raw counts or CPM/TPM normalized counts, and should not be log-transformed**.
    
    For the UMAP plots produced by ```plotData```, if the input is not a Seurat object already containing UMAP and PCA embeddings, the input slot's expression is log-normalized internally unless the maximum value in the provided slot is less than 20.

4. **What if I have multiple batches of data? Should I perform any integration?**

    No batch integration is required. Instead, we recommend running CytoTRACE 2 separately over each dataset. While raw predictions are made per cell without regard to the broader dataset, the postprocessing step to refine predictions  adjusts predictions using information from other cells in the dataset, and so may be impacted by batch effects. Note that CytoTRACE 2 outputs, except for CytoTRACE2_Relative, are calibrated to be comparable across datasets without further adjustment. Therefore, no integration is recommended over the predictions either. As for CytoTRACE2_Relative, it should not be compared across different runs, as it is based on scaling of CytoTRACE2_Scores within the context of the specific input dataset. 

5. **Do the R and Python packages produce equivalent output?**

    When run without batching (i.e., downsampling the input dataset into batches [or chunks] for parallel processing or to save memory), these packages produce equivalent output. When batching is performed, package outputs will vary, but remain highly correlated in practice.

6. **What strategies are recommended for managing very large datasets (>100K cells) with CytoTRACE 2?**

    For large datasets, subdividing the data into smaller segments, each containing up to 100,000 cells, is advisable. This division not only facilitates more efficient memory management and processing but also preserves the integrity of your analysis. Depending on the dataset’s characteristics, you can segment by experimental conditions (technical batches) or samples. Additionally, when choosing your subset size, please be mindful of the computational resources available to you--some systems may support ~100,000 cells while for others, a further reduced subset size may be preferable.

7. **What if my dataset includes rare cell types?**

    CytoTRACE 2 implements an adaptive nearest neighbor smoothing step as the final component of postprocessing. When analyzing phenotypes expected to have five or fewer cells, we recommend bypassing the KNN smoothing step so that predictions for these rare cells are not forced toward more abundant phenotypes. In practice, you can simply use the preKNN score output (preKNN_CytoTRACE2_Score) instead of the final KNN-smoothed value (CytoTRACE2_Score).

8. **Does CytoTRACE 2's performance depend on the number of UMI/gene counts per cell?**

    Although generally insensitive to variation in gene/UMI counts per cell, CytoTRACE 2 requires further optimization for cells that have exceedingly low gene expression levels, particularly those with fewer than 500 genes expressed per cell, as its performance can become less reliable in these cases. For best results, **a minimum gene count of 500-1000 per cell is recommended.**


9. **Why does the UMAP generated by CytoTRACE 2's `plotData` function differ from the cell embeddings saved in my input Seurat object, and how can I get them to match?**

    If the input Seurat object already contains cell embeddings, those need to be accessible as `"pca"` and `"umap"` from `seurat_object@reductions` layer, and the `plotData` function will automatically use these embeddings for visualization. Make sure to have `is_seurat` set to `TRUE`.

    Otherwise, the `plotData` function utilizes the `RunPCA` and `RunUMAP` functions to generate cell embeddings internally. `pc_dims` and `seed` parameters of the `plotData` function (corresponding to the `dims` and `seed.use` parameters in `RunUMAP`) can be tweaked to adjust the UMAP plot.
      
    ```R
        # Example of setting pc_dims and seed in plotData
        plotData(data_object, pc_dims = 15, seed = 42)
    ```

    If you have a different type/name of reduction (other than UMAP/umap) in your Seurat object that you would like to use for visualization, you can either rename it to "umap" or manually replace the embeddings on `plotData` output plots:

      ```R
      # Assuming 'data_object' contains the prefered cell embeddings in "tsne" reduction

      emb_1 <- data_object@reductions$tsne@cell.embeddings[,1]  # First dimension
      emb_2 <- data_object@reductions$tsne@cell.embeddings[,2]  # Second dimension

      # Replace CytoTRACE 2 generated UMAP coordinates in the plot output
      plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
      plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_2"] <- emb_2
      ```
      Alternatively, you can use predicted values from `cytotrace2()` function and plot those directly on your prefered original plot, without using the `plotData()` function.

</details>



## GitHub package content

<details><summary>cytotrace2_r</summary>

- the main folder of the R package
  - *R*: the subfolder of R scripts
    - *cytotrace2.R*: the main function that executes the pipeline from input to final predictions
    - *plotting.R*: plotting functions based on generated prediction and provided annotation file (if available)
    - *postprocessing.R*: functions used for postprocessing (smoothing, binning, adaptive kNN smoothing) raw predicted values
    - *prediction.R*: functions used in the prediction part of the pipeline
    - *preprocess.R*: functions used in loading and processing (orthology mapping, feature selection, feature ranking per cell, log2 CPM-transformation) of the input data
  - *inst/extdata*: the subfolder of files internally loaded and used in the pipeline, includes:
    - model parameter matrices (.rds objects containing learned parameters of all the ensemble models)
    - preselected feature set for the model
    - orthology mapping dictionary (human to mouse)
    - alias mapping dictionary (gene symbol to alias gene symbol)
  - *man*: the subfolder of function manuals
    - can be accessed by ```help(function_name)```
  - *NAMESPACE*: the file containing the package namespace
  - *DESCRIPTION*: the file containing the package description
</details>

<details><summary>cytotrace2_python</summary>

- the <a href="/cytotrace2_python" target="_blank">main folder of the Python scripts and files</a> 
  - *cytotrace2_py*: the subfolder of Python scripts
    - *cytotrace2_py.py*: the main function that executes the pipeline from input to final predictions
    - *resources*: the subfolder of files internally loaded and used in the pipeline, includes:
      - model parameter matrices (.pt objects containing learned parameters of all the ensemble models)
      - preselected feature set for the model
      - human and mouse orthology table from Ensembl
      - alias table from HGNC (gene symbol to alias gene symbol)
      - alias table from MGI (gene symbol to alias gene symbol)
    - *common*: the subfolder of utility functions used in the pipeline
  - *environment_py.yml*: conda environment file for creating a new environment with all the necessary packages.
  - *setup.py*: the file containing the package description

</details>

## CytoTRACE 2 web application

An interactive RShiny web application can be accessed at [cytotrace2.stanford.edu](https://cytotrace2.stanford.edu), allowing users to run analyses on their own data, browse results for 33 ground-truth–annotated datasets, explore potency-associated genes and gene-set enrichment across the single-cell potency atlas, download the atlas, and access Python vignettes for model training and custom GSBN architectures.

## Authors
CytoTRACE 2 was developed in the <a href="https://anlab.stanford.edu/" target="_blank">Newman Lab</a> by Minji Kang, Gunsagar Gulati, Erin Brown, Susanna Avagyan, Jose Juan Almagro Armenteros and Rachel Gleyzer.

## Contact
If you have any questions, please contact the CytoTRACE 2 team at cytotrace2team@gmail.com.

## License
Please see the <a href="LICENSE" target="_blank">LICENSE</a> file.

## Citation
If you use CytoTRACE 2, please cite: 

**Improved reconstruction of single-cell developmental potential with CytoTRACE 2.** *Nature Methods*, 2025. <br>
Minji Kang<sup>\*</sup>, Gunsagar S. Gulati<sup>\*</sup>, Erin L. Brown<sup>\*</sup>, Zhen Qi<sup>\*</sup>, Susanna Avagyan, Jose Juan Almagro Armenteros, Rachel Gleyzer, Wubing Zhang, Chloé B. Steen, Jeremy Philip D’Silva, Janella Schwab, Michael F. Clarke, Aadel A. Chaudhuri, and Aaron M. Newman<sup>†</sup>.
[doi.org/10.1038/s41592-025-02857-2](https://doi.org/10.1038/s41592-025-02857-2)


<!-- ## Preprint
Mapping single-cell developmental potential in health and disease with interpretable deep learning.

bioRxiv 2024.03.19.585637; <br> doi: https://doi.org/10.1101/2024.03.19.585637 -->

