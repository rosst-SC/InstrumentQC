[![GPL-3.0][gpl3-shield]][gpl3]
[![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

# InstrumentQC 

Having reached this webpage, you likely discovered the [InstrumentQC](https://umgccfcss.github.io/InstrumentQC/) website I implemented at the University of Maryland, Baltimore to track daily QC for the Cytek and BD instruments. I am providing this documentation for those of you who would like to understand how it is implemented, whether for curiosity or because you'd like to modify the code to set up something similar for your own cytometers.

Having placed significant effort into this pet-project, I would love if others in the cytometery community could also benefit. Since assembling everything from scratch requires some familiarity with both R and Git, I am writing up this detailed explanation of the process behind setting up and customizing the dashboards to give the average cytometry aficionado a fighting-chance. 

To get started, navigate to the walk-through explanation website by [clicking here](https://davidrach.github.io/InstrumentQC_Install/)

If you get stuck, or parts of the documentation are unclear, please reach out by opening a [discussion](https://github.com/DavidRach/InstrumentQC_Install/discussions) and I will try my best to help you troubleshoot. What makes free-and-open source software great (and in my opinion,fun to work with) is seeing how the small collective contributions of users improve the projects over time. 

Best-
David

# High-level Overview

After daily QC is run on a Cytek instrument, the relevant information is stored as .fcs and .csv files in specific folders. Using functions incorporated in the [Luciernaga](https://github.com/DavidRach/Luciernaga) R package, the data associated with these newly generated files is processed using [R](https://www.r-project.org/) once a day at a user-designated time (implemented via Windows Task Sceduler). The version control software [Git](https://git-scm.com/) keeps track of changes to the processed data, passing updates to the online [GitHub](https://github.com/) repository for storage. This data is then referenced when building the dashboard website using [Quarto](https://quarto.org/), which is published as a GitHub page allowing for url access.

For multiple instruments, the above process is repeated on each instruments computer, and the Quarto webpages are modified to display data from each instrument. 

# Contributing

Trying to implement this on your end and encountering issues? Glad to help, please open a [new discussion](https://github.com/DavidRach/InstrumentQC_Install/discussions)

Encountered a bug and want to report it? Found a glaring typo in my pre-morning coffee English? I always appreciate the help, please open a [new issue](https://github.com/DavidRach/InstrumentQC_Install/issue) and I will work on fixing/editing when I get a chance.

[![GPL-3.0][gpl3-image]][gpl3]
[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]


This material is licensed under the [Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].
The code examples are licensed under the [GNU GENERAL PUBLIC LICENSE (GPL-3.0)]()



[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
[gpl3]: https://www.gnu.org/licenses/gpl-3.0.en.html
[gpl3-image]: https://www.gnu.org/graphics/gplv3-88x31.png
[gpl3-shield]: https://img.shields.io/badge/license-GPLv3-blue

