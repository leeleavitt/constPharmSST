## Propriocepter subclass analysis
In this example I select proprioceptors/L1 only to analyze. Initally there doesn't seem to be any contrast. from there i press **`l`**. This then allows me to label the neurons by cell-type and also subclass. I then press **`f`** to perform random forest analysis. After this i click through significant genes and observe the boxplot. **NOTICE** before I started this analysis i pressed **`g`** to select a large/custom gene list to begin this analysis.
![][image1]

## IPSI cci analysis
In this example i've start off with looking at all cell-types. In this labeling example we are unable to define whether the animal has undergone a chronic constriction injury to provoke a pain state. By pressing **`l`** I change my label to be both `label_experiment` and `label_cellType`. This operation combines these labels to creates the labels. From there, I perfomr randomforest for R13_ipsi vs all others. This analysis shows `MRAP` to be a significant gene postively expressed only in the R13 cell class. This could be a potential drug target for pain targeting.
![][image2]

## Multiple label Selection
You can select more than one label for the random forest analysis. This might be useful if you have a drug that hits two or three cell classes specifically. In this example I select L2(a-delta-LTMRS), and N14(c-LTMRS). This now shows me all genes that are expressed or not expressed exlcusively in these cell classes.
![][image3]

## Total Gene Analysis
The software now supports analysis on the entirety of the genetic data. To perform this analysis.

  1. Add text file to the search in your profile which does not contain any terms
  2. Press <kbd>g</kbd> and select this empty search. This will return all gene discovered during the transcriptomics.
  3. Now that all gene are returned performing random forest analysis will be performed on all genes. <kbd>f</kbd>
  4. To work through all genes during it is useful to simply press <kbd>shift</kbd><kbd>g</kbd>. The genes are sorted based on the random forest. At this time the biplot become most useful for understanding the genes which play a positive role in defining the cell class (these genes will be closest to the cells of interest), vs the genes which play a negative role in defining the cell class(these gene will be furthest away and opposing the cell class).


[image1]: ../howToGifs/28_advancedCellLabeling.gif
[image2]: ../howToGifs/29_advancedCciAnalysis.gif
[image2]: ../howToGifs/34_InteractiveLabeling.gif
[image3]: ../howToGifs/31_advancedTwoClassAnalysis.gif

