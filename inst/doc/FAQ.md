# FAQs

This document intends to be the first reference for any doubts that you may have regarding s2dv. If you do not find the information you need, please open an issue for your problem.

## Index
1. **How to**
   1. [Global Map with non-standard longitudinal boundaries](#1-global-map-with-non-standard-longitudinal-boundaries)

2. **Something goes wrong...**  
   1. [CDORemap() returns errors or warnings with specific module versions](#1-cdoremap-returns-errors-or-warnings-with-specific-module-versions)


## 1. How to

### 1. Global Map with non-standard longitudinal boundaries

Usually, the global maps are displayed:
  a) from -180 to 180 degrees Est being the Grenwich meridian in the center of the map or
  b) from 0 to 360 degrees Est being the Grenwich meridian in the left marging of the map.

You can run the following code to test both cases using PlotEquiMap:

```
library(s2dv)
a <- 1:(180*360)
dim(a) <- c(lat = 180, lon = 360)
PlotEquiMap(a, lon = -179.5 : 179.5, lat = -89.5 : 89.5) # case a)
PlotEquiMap(a, lon = 1 : 360, lat = -89.5 : 89.5) # case b)
```

What if I want to use different boundaries of the region? For instance, if I want to display Atlantic, Indic and Pacific Oceans being centered. Then, you should do some extra steps:

```
library(ClimProjDiags)
a <- Subset(a, along = 'lon', indices = c(20 : 360, 1 : 19))
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(3, 1))
PlotEquiMap(a, lon = 21 : 380, lat = -89.5 : 89.5, drawleg = FALSE,
            coast_width = 0.0, filled.continents = FALSE)
map("world", wrap = c(20, 390), add = TRUE)
ColorBar(var_limits = c(1, max(a)))
```

Note: You can adjust many parameters to visualize the plot, here we are just showing how to move the boundaries.

If you want to add other information to the plot (e.g.: hatching, points, countours, ...), you can add it just before ColorBar() function.



## 2. Something goes wrong...  

### 1. CDORemap() returns errors or warnings with specific module versions  
CDORemap() uses cdo and ncdf4 inside, and the performance is impacted by those tools a lot. 
Some instances may work with a specific set of module combination but not with another. 
Since the incompatibility is not from the R code, it is hard to improve or prevent the failure. 
Here are some detected cases that specific versions need to be used.  
(1) The 'grid' parameter is a file  
- The workable version combination:  
CDO/1.9.8-foss-2015a  
R/3.6.1-foss-2015a-bare  
HDF5/1.8.14-foss-2015a  
- The unworkable version combination:
_It returns a warning about HDF5._  
CDO/1.6.3-foss-2015a  
R/3.6.1-foss-2015a-bare  
HDF5/1.10.5-foss-2015a  





