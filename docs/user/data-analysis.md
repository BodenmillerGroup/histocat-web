# Single-cell data analysis

## Pipelines

When a proper dataset is uploaded and selected in the Datasets panel, then users can run single-cell data analysis on them.
**Pipeline** consists of several processing steps that reflect the output user wants to get as a result.
You can run the same pipeline on different datasets, different acquisitions as many time as you wish.
All processing outputs are stored as **Results** (see next documentation section).

At the beginning, user needs to design processing pipeline in visual pipeline editor. The processing flow is from top to bottom:

![Pipeline editor](../assets/pipeline-editor.png)

First, select a step you want to add from the list of available processing steps:

![](../assets/add-pipeline-step.png)

There are three pipeline step categories:

* **Preprocessing**: such actions as marker filtering, data transformation, scaling, etc.
* **Embeddings**: dimensionality reduction via PCA, t-SNE, UMAP
* **Clustering**: community detection by Louvain or Leiden algorithms

Each step has his own sets of parameters. By clicking on parameter field you can see its short description:

![t-SNE perplexity parameter description](../assets/step-parameter-description.png)

One can change step order in a pipeline by clicking **Up** / **Down** arrow buttons or delete it:  ![](../assets/common-step-buttons.png) 

!!! warning "Warning"
    Some steps require additional steps defined beforehand in the pipeline to run properly. For instance, UMAP embedding or any clustering step needs Neighborhood graph detection step defined upper in the pipeline. If you miss such steps, there will be a warning message preventing the submission to processing.

![Pipeline step order](../assets/pipeline-step-order.png)

As soon as pipeline is complete, you can submit it to back-end for a processing by background workers.
After clicking ![](../assets/process-button.png) button, a dialog window will appear where one can select acquisitions that have to be analysed by the pipeline:

![Acquisitions selection dialog](../assets/acquisition-selection-dialog.png)

Depending on the pipeline's complexity, processing can take some time. For example, t-SNE / UMAP or clustering detection are quite long-running processes.
When processing is complete, there will be a notification message and the analysis result will appear in the **Results** view (please see next section in documentation).

In case you want to re-use a pipeline in the future or to share it with other group members, it can be saved in the Pipeline tab by clicking **Save** button:

![](../assets/save-pipeline-button.png)

At any moment you can rename saved pipeline, give it some description or delete it altogether.
To load stored pipeline please click Load Pipeline button: ![](../assets/load-pipeline-button.png) in the `Pipelines` panel:

![](../assets/pipelines-list.png)

## Results

All results of individual pipeline runs are grouped in the lower half of **Datasets** tab.
By default, their name is a timestamp of the processing completion.

![Results list](../assets/results-list-view.png)

!!! info "Info"
    Please keep in mind that result will appear in the list only when processing on the server side is complete, which can take a lot of time (up to hours).

Because results of pipeline data processing are stored in _AnnData_ file format, it is easy to download them and continue further analysis on the client side.
To download result as a zip archive, please click **Download** button: ![](../assets/result-download-button.png) 

Other available commands are: **Load** result, **Edit** result and **Delete** result.
By selecting **Edit** command user can give a custom, meaningful name and description to the result entry instead of the automatically assigned name (which is a timestamp by default).

!!! info "Info"
    Result also stores the whole pipeline that has been used for its processing so you can easily check which parameters were used for analysis or adjust these parameters and re-run pipeline in order to get new result.

Depending on the pipeline design, after selecting a result user will see Scatter, PCA, t-SNE and UMAP plots in the corresponding panels:

![Result view](../assets/result-view.png)

!!! info "Info"
    Scatter plot allows user to select markers on X and Y axes and display correlation between them on a 2D-scatterplot.

If you ran community detection in your pipeline, then additional tab called Clustering will appear.
There you may see violin and dot plots for Leiden / Louvain clustering:

![](../assets/clustering-view.png)

## Gating

Bi-directional gating on both raster images and single-cell plots can be used to select subsets of cells for further visualization strategies.
Furthermore, users may assign cell type annotations to selected gates and run Random Forest classification similarly to QuPath implementation. 

![Gating](../assets/gating.png)
