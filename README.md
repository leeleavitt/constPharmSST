# Transcriptome Surfer

**NOT COMPATIBLE WITH MAC**

This a package developed to combine calcium imaging data with single cell transcriptomic data. To install this software on your local computer please follow these [instructions](./Misc/keyBoardWalkthroughs/installation.md). It is also highly reccomended you learn [regular expressions](https://regex101.com/). Making a new profile will automatically save two different gene searches into the `./profiles/yourProfile/searches`. Look at it for examples of how to search. Also, feel free to delete these files.

1. To get started on using this software navigate to the directory and simply double click `clickMe.Rdata`. This will open an R console, and load in all software and data.
![][image1]
2. The first question that will be asked is what profile you would like to use. If you want to make a new profile simply click **New profile** or click **Cancel**. 
![][image2]
3. Answering `New profile` will next ask you the name of the new profile.
![][image3]
4. Starting from line three we are now working from scratch. The important thing to realize here is that we need to select search terms to find genes. 

  * First navigate into the folder named profiles
  * Find the profile you just created
  * Open the searches to find the first search terms we will use. Notice each search term is separated by new line. 
  * Spend some time exploring the folder [Lee Leavitt](./profiles/Lee_Leavitt/Searches) for a few examples of different types of searches you can do.
  ![][image4]

5. Now we will use the keyboard to guide us during the rest of the data analysis. Follow this.
    1. [Visualizations](./Misc/keyBoardWalkthroughs/visualizations.md):
       1. Genes, and gene-subset
       2. Cell labels
       3. Matrix normalization
    2. [Interactivity and gene Informativeness](./Misc/keyBoardWalkthroughs/labelsInformation.md) : Using labels to identify informative genes
       1. **RandomForest** against selected labels to rank genes in terms of importance
       2. Using linear regression to identify statistical significance of difference across all labels vs selected label.
    3. [Saving and reloading](./Misc/keyBoardWalkthroughs/saving.md)
    4. To gain insite for how you can potentially use this, check out the [advanced examples](./Misc/keyBoardWalkthroughs/advancedExamples.md) page.




[image1]: ./Misc/howToGifs/1_startup.gif
[image2]: ./Misc/howToGifs/2_profileSelection.gif
[image3]: ./Misc/howToGifs/3_makingNewProfile.gif
[image4]: ./Misc/howToGifs/4_editingSearchTerms.gif