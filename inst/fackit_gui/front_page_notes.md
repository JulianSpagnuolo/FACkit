
<!-- front_page_notes.md is generated from front_page_notes.Rmd. Please edit that file -->

# FACkit GUI

Running FACkit, a FACS Analysis Tool Kit, version 0.1.2 developed by
Julian Spagnuolo.

## Introduction

This GUI was developed to enable intuitive and streamlined exploration
and analysis of high-dimension FACS data using tSNE.

Documentation and explanation of the tools and steps used are either
located beside each box or can be accessed by clicking “help” buttons.

In the event of crashes or unexpected behaviour please contact the
package developer; Julian Spagnuolo <julianspagnuolo@gmail.com>.

# VERSION NEWS

## FACkit 0.1.2

  - Added Marker name sanitizing function to the shiny app
  - Removed the `NEWS.md` file integration (it wasn’t working)

## FACkit 0.1.1

  - Increased App Upload Limit to 100 Mb
  - Removed app loading of ggalt (for server purposes) and switched
    instances of geom\_bkde to geom\_density

## FACkit 0.1.0

  - Added a `NEWS.md` file to track changes to the package.
  - Added documentation to FACkit GUI.
