# Trichromatic advantage in wild primates

This repository contains all analysis code used for the capuchin color vision project.

The manuscript is titled:

**Trichromacy increases fruit intake rates of wild capuchins (*Cebus capucinus imitator*)**

All code consists of a single script:

* capuchin.R

The script requires four input files containing data from the Santa Rosa Capuchin Project's database. These files are located in:

* data/CapuchinForagingData.txt

	A table with rows consisting of feeding events and associated metadata on the behavioral state, the feeding bout, the animal subject, and the feeding patch (including its phenological conditions).

* data/FruitColoration.txt

	A table with rows consisting of fruit resources utilized by Santa Rosa capuchin monkeys and associated taxonomic and color metadata.

* data/2007_2008Dominance.xlsx

	An Excel spreadsheet with date, group, individualID, and ordinal dominance rank information

* data/FPVDetailwithVirtualCBH_Duration(Min)_May30_2017.xlsx

	An Excel spreadsheet with phenological information, mostly importantly a ripe fruit score (SimpleRipeFruitScore) and a SeasonID, the latter of which has values corresponding to the PhenologyID in other datasets.

* data/energy.txt

	Energetic information about fruits involved in this study. The most important columns are the Conspicuity, Diameter, DryEnergy, and FinalEnergyRate.

For questions regarding the data, contact [Amanda Melin](mailto:amanda.melin@ucalgary.ca).

For questions regarding the code, contact [Kenny Chiou](mailto:kchiou@uw.edu).