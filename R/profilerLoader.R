profileLoader <- function(){
    # Now ask the user to select their profile
    bringToTop(-1)
    # Go into profiles and recursively list what is there
    profilesToChoosePretty <- list.files('./profiles', recursive = F)
    profilesToChoose <- c('New Profile', profilesToChoosePretty)

    cat("\nSelect the profile you would like to use\n")
    chosenProfile <- select.list(profilesToChoose, title="Profiles")

    if(chosenProfile == 'New Profile' | chosenProfile == ""){
        cat('\nEnter the name of your profile, buddy\n')
        chosenProfile <- scan(what = 'character', n=1, quiet = T, sep=">")
        
        # Now that we have the name of the new profile, lets make a new directory
        dir.create(paste0("./profiles/", chosenProfile))
        # Make the savedCsv folder
        dir.create(paste0("./profiles/", chosenProfile,'/savedCsv'))
        # Make the searches folder
        dir.create(paste0("./profiles/", chosenProfile,'/searches'))    
        
        # Give them a default searches to begin working with
        invisible(file.copy("./Misc/rawData/goTermSearch.txt", paste0("./profiles/",chosenProfile,'/searches')))
        invisible(file.copy("./Misc/rawData/leeTermSearch.txt", paste0("./profiles/",chosenProfile,'/searches')))
    }
    setwd(paste0("./profiles/", chosenProfile))

    settingLogic <- list.files(pattern = 'SETTINGS')
    if(length(settingLogic) == 0){
        SETTINGS <<- NULL
    }else{
        load("SETTINGS.Rdata", globalenv())
    }
    invisible(tsInteract(SETTINGS))
}