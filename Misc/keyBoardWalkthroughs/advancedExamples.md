## Propriocepter subclass analysis
In this example I select proprioceptors/L1 only to analyze. Initally there doesn't seem to be any contrast. from there i press **`l`**. This then allows me to label the neurons by cell-type and also subclass. I then press **`f`** to perform random forest analysis. After this i click through significant genes and observe the boxplot. **NOTICE** before I started this analysis i pressed **`g`** to select a large/custom gene list to begin this analysis.
![][image1]

## IPSI cci analysis
In this example i've start off with looking at all cell-types. In this labeling example we are unable to define whether the animal has undergone a chronic constriction injury to provoke a pain state. By pressing **`l`** I change my label to be both `label_experiment` and `label_cellType`. This operation combines these labels to creates the labels. From there, I perfomr randomforest for R13_ipsi vs all others. This analysis shows `MRAP` to be a significant gene postively expressed only in the R13 cell class. This could be a potential drug target for pain targeting.
![][image2]


[image1]: ../howToGifs/28_advancedCellLabeling.gif
[image2]: ../howToGifs/29_advancedCciAnalysis.gif


