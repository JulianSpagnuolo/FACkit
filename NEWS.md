# FACkit 0.1.3

* Modified shiny app data input method to use rhandsontables instead of DT - much easier/quicker to input.
* Improved marker name sanitization

# FACkit 0.1.2

* Added Marker name sanitizing function to the shiny app
* Removed the `NEWS.md` file integration (it wasn't working)
* Explicit conversions between Rcpp::vector and std::vector in bindist.cpp to prevent overloaded "=".

# FACkit 0.1.1

* Increased App Upload Limit to 100 Mb
* Removed app loading of ggalt (for server purposes) and switched instances of geom_bkde to geom_density 

# FACkit 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added documentation to FACkit GUI.


