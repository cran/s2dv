# s2dv 1.0.0 (Release date: 2021-06-16)
- New functions:
ACC, Ano_CrossValid, BrierScore, CDORemap, Cluster, Consistent_Trend, EOF, EuroAtlanticTC, Filter, Histo2Hindcast, 
NAO, Plot2VarsVsLTime, PlotACC, PlotBoxWhisker, PlotVsLTime, ProbBins, ProjectField, RatioRMS, 
RatioSDRMS, REOF, Spectrum, Spread, StatSeasAtlHurr, UltimateBrier
- Season(): Accept one-dimension input.  
- Persistence(): Add parameters checks for 'start' and 'end'; correct the output 'AR.lowCI' and 'AR.highCI'.  
- Corr(): Add parameter 'member' and 'memb_dim'. They allow the existence of the member dimension
 which can have different length between exp and obs, and users can choose to do the ensemble mean 
first before correlation or calculate the correlation for individual member. 
- InsertDim(): Remove Apply() to improve the efficiency.  
- Reorder(): Improve efficiency.  
- Indices functions take the case without 'memb_dim' into consideration. The climatology calculation for the anomaly is member-dependent if member exists.  
- PlotStereoMap(): Add contour and arrow feature.  
- PlotAno(): Add parameter check for 'sdates'.  
- PlotEquiMap(): Add new arguments 'contour_draw_label', 'lake_color', 'lab_dist_x', 'lab_dist_y', and 'degree_sym'. Fix the border error; the border grids are fully plotted now. Add ocean mask feature.

# s2dv 0.1.1 (Release date: 2020-11-16)
- Change the lincense to Apache License 2.0.
 
# s2dv 0.1.0 (Release date: 2020-11-12)
- New functions: Ano(), Composite(), PlotAno(), Smoothing(), AMV(), GSAT(), SPOD(), TPI(), GMST(), Persistence().
- Change the default value of PlotClim() parameter 'fileout' to NULL.
- Change Regression() parameter 'time_dim' to 'reg_dim', and enable the inputs to be vectors.
- Change Trend() parameter 'time_dim' default value from 'sdate' to 'ftime'.
- Change the default of Season() parameter 'time_dim' from 'sdate' to 'ftime'.
- Bugfix for Regression() na.action. 'na.action = na.fail' was not functional before.
- Add p-value by ANOVA in Trend().
- Bugfix for Trend() slope, detrended, and p-value.
- Change MeanDims() na.rm default to FALSE to be in line with mean().
- Remove unecessary parameter checks in Clim().
- Change parameter 'memb_dim' to 'dat_dim', and the default value from 'member' to 'dat' in Corr(), RMS(), and RMSSS().
- Allow RMS() and RMSSS() to have vector data input.
- Bugfix for Load() when start date and first lead time is not consistent.
- Improve Season() performance by using apply() when 'ncores' is not bigger than 1

# s2dv 0.0.1 (Release date: 2020-02-07)
- The package is the advanced version of package 's2dverification', adopting the regime of package 'multiApply' for all the analytic functions. Most of the other functions for plotting and data retrieval in 's2dverification' are also preserved in this package.
- Because of the adoption of 'multiApply' regime, the functions work well with package 'startR'. 
- All the packages mentioned above are developed by BSC-CNS.

