# Interactivity and Gene Informativeness
Often times we find our selves with too many genes to realistically interpret. Fortunately we've included three tools to help you in this search for informative genes. Here we present 
1. **RandomForest**: 1 vs all, or all vs all
2. **PCA** (Principal Component Analysis) with interactive singular vector selection
3. **Linear regression**: for statistical significance of genes.

## Random Forest
Random forest provides you with a ranking of most informative genes. Instead of simply returning the highest expressed genes, we return the genes that are most informative for defining either,
  * A single cell class vs everything else
  * All cell classes vs all cell classes.

1. To get started make sure you have a large selection of genes to work with. Its best to press `g` to select a search that will return a lot of genes. This also works with a simple list of genes. After these genes have been loaded press **`f`** (for forest). 
    * The First question is 1 vs all or all vs all. I have best result with 1 vs all.
    *  The next question ask which label/labels to select. You can choose more than 1, but the statistical significane on the boxplot will no longer work.
    * If you have a lot of genes selected, the console will now ask you how many genes to return. If you would like statistical significance returned on the heatmap, please choose **<200** genes.
    * Once complete notice how the label you've selected is now red. Also, the heatmap has statistical significance listed for each gene.
    ![][image1]


## BiPlot
The biPlot we use in this visualization is created using singular vector decomposition. For a better understandinf of this please see Dr. Jeff Phillips textbook [*Mathematical Foundations for Data Analysis*](http://www.cs.utah.edu/~jeffp/M4D/M4D.html) chapter 17. Principal component analysis represents a high dimensional structure (the heatmap) in a 2 dimensional structure (the bi-plot). The biplot is composed of 2 main features
1. Principal components, represented by the cells (colored lables).
2. principal directions represented by the genes (gray gene names with spectral segments projecting towards it).

During the singular vector decomposition singular vectors are found. Often only the first and second singular vectors are use to produce the principal component analysis. Here, we give you the option to pick and choose which singular vectors you would like to use for the bi plot display.

1. Press **`v`** (for vectors). This will bring up a plot. This plot shows the total variance explained for each singular vector. As you will continue to notice singluar vector 1 often dominates, which can lead to an uninformative bi-plot display. Selecting the 2nd and 3rd, although they do not explain a lot of variance, can produce informative and more interesting bi-plots.
![][image2]

2. Performing this analysis we are able to observe the genes which are driving this separation of the data. Now, to interact with the bi-plot, press **`F2`**. This will now allow you to select genes on this plot.
    * Selected genes will become selected on the heatmap and the boxplot.
    * Only the top 200 genes will be displayed (genes producing the most variance).
    ![][image3]

## BoxPlot
The labels define categories of cell, and the box-plot provides the best representation of these different classes of neurons. Each box-plot's width is varied based on the number of cells within that category. Additionally, each data point is plotted ontop of the box-plot. If the cell is defined by a subclassification, then this is display rather than a **[*]** 

1. The boxplot is the most informative plot in the analysis. The y axis is the level of the gene expressed. Additional gene aliases, and other designations are added to the plot to provide more information about the selected gene.
2. If a specific label was selected for random forest (**`f`** button), then this designates the cell class to be compared against for linear regression. If other cell classes are found to be significantly different, the significance is added to the plot with red asteriks.

![][image4]

## Heatmap
The heatmap provides a representation of all selected gene. It also shows the label separation, by different color representation, and horizontal lines added to the plot. 

The heatmap is interactive and clicks on this heatmap updates all three windows. 
1. The boxPlot updates with the gene selected
2. The calcium imaging trace updates with the cell selected. 
3. The bi-plot changes the color of the selected gene to red, and make the font bold.

To make the heatmap clickable press **`F1`**. The heatmap is now active and accepting clicks.
![][image5]




[image1]: ../howToGifs/20_randomForest.gif
[image2]: ../howToGifs/21_biPlot.gif
[image3]: ../howToGifs/22_biPlot_clicks.gif
[image4]: ../howToGifs/23_boxPlot.gif
[image5]: ../howToGifs/24_heatMap.gif

