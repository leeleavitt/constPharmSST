# Visualization
The portion of the walk through will introduce you to choosing the label to wok with, subsetting data, creating different names for the cells.

Lets start from scratch. You have created a new profile. Now what (all keyboard shortcuts are case sensitive)?

1. Load up genes.
    1. press <kbd>g</kbd>for the search terms to be applied against gene names
    2. press <kbd>ctrl</kbd><kbd>g</kbd> for the search to be applied against the go-Terms. This is significantly slower than the first option if many searches are included. 

![][image1]

2. Remember to add new searches, or edit current searches navigate to `~\searches` within your profile folder,

![][image2]

3. It is now important to select the way we should label the cells. Cells are the rows in this heatmap. Starting off they are all red. Pressing `l` with allow you to choose the correct labeling. Experiment at this point to view all the different types of labels we can use.

![][image3] 

4. Normalization of the heatmap can undergo various types of normalization techniques. This can also affect the various other analysis, and should be used frequently. Press **`n`** will allow you to select the method for normalization. **NOTE** combinations are allowed, even if they do not make sense.

![][image7]

4. Although we have the cells represented by specific lables, we might want the names to be different. For example, I'd also want to know what the subclasses are within the, pressing `r` (for rename), will bring up all options to rename the neurons.

![][image4]

5. Selecting a gene subset can be a useful way to view a smaller portion of the data pressing `G` will allow you to select genes within your gene super list.

![][image5]

6. Choosing a cell subset can also be quite useful. To choose a subset of cells press `c`. 
    1. Notice though, at this point you can only select a subset of cell, based on the label you have selected. 
    2. The first question posed is asking if you would like to subset based on experiment. Answering cancel, will simply take you onto the next question of which label to display.

![][image6]


[image1]: ../howToGifs/5_selectingSearch.gif
[image2]: ../howToGifs/4_editingSearchTerms.gif
[image3]: ../howToGifs/7_labelSelection.gif
[image4]: ../howToGifs/17_rename_classes.gif
[image5]: ../howToGifs/18_geneSubset.gif
[image6]: ../howToGifs/19_cellSubset.gif
[image7]: ../howToGifs/30_normalization.gif
