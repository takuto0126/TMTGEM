#################################################
#                                               #
# Control file for COMCOT program (v1.7)        #
#                                               #
#################################################
#--+-----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
#===============================================:===============================
# General Parameters for Simulation             : Value Field                  |
#===============================================:===============================
#Job Description: NZ30sec bathymetry, Spherical Coordinates for code testing
 Total run time (Wall clock, seconds)           :  300.000
 Time interval to Save Data    ( unit: sec )    :    5.0
 Output Zmax & TS (0-Max Z;1-Timeseries;2-Both) :     0
 Start Type (0-Cold start; 1-Hot start)         :     0
 Resuming Time If hot start (Seconds)           :   800.00
 Specify Min WaterDepth offshore  (meter)       :    10.00
 Initial Cond. (0:FLT,1:File,2:WM,3:LS,4:FLT+LS):     0
 Specify BC  (0-Open;1-Sponge;2-Wall;3-FACTS)   :     0
 Specify Input Z filename (for BC=3, FACTS)     : mw94_n22_nz_ha.xyt
 Specify Input U filename (for BC=3, FACTS)     : mw94_n22_nz_ua.xyt
 Specify Input V filename (for BC=3, FACTS)     : mw94_n22_nz_va.xyt

#===============================================:===============================
# Parameters for Fault Model (Segment 01)       :Values                        |
#===============================================:===============================
 No. of FLT Planes (With fault_multi.ctl if >1) :   1
 Fault Rupture Time (seconds)                   :   60.0
 Faulting Option (0: Model; 1- Data;)           :   0
 Focal Depth                             (meter):   10000.000
 Length of source area                   (meter):   100000.000
 Width of source area                    (meter):   50000.000
 Dislocation of fault plate              (meter):    10.000
 Strike direction (theta)               (degree):   193.000
 Dip  angle       (delta)               (degree):    10.000
 Slip angle       (lamda)               (degree):    81.000
 Origin of Comp. Domain (Layer 01) (Lat, degree):    38.8
 Origin of Comp. Domain (Layer 01) (Lon, degree):   143.5
 Epicenter: Latitude                    (degree):    38.8
 Epicenter: Longitude                   (degree):   143.5
 File Name of Deformation Data                  : segment_parameter.dat
 Data Format Option (0-COMCOT; 1-MOST; 2-XYZ)   :     2

#===============================================:===============================
#  Parameters for Wave Maker                    :Values                        |
#===============================================:===============================
 Wave type  ( 1:Solit, 2:given, 3:focusing )    :     1
 FileName of Customized Input (for Type=2)      : fse.dat
 Incident direction( 1:top,2:bt,3:lf,4:rt,5:ob ):     2
 Characteristic Wave Amplitude        (meter)   :     0.500
 Typical Water depth                  (meter)   :  2000.000 
 
#===============================================:===============================
#  Parameters for Submarine LS/Transient Motion :ValUes                        |
#===============================================:===============================
 X Coord. of Left/West Edge of Landlide Area    :  177.00
 X Coord. of Right/East Edge of Landlide Area   :  179.00
 Y Coord. of Bottom/South Edge of Landlide Area :  -41.00
 Y Coord. of Top/North Edge of Landlide Area    :  -39.00
 File Name of landslide Data                    : landslide_test.dat
 Data Format Option (0-Old; 1-XYT; 2-Function)  :     2
 
#===============================================:===============================
# Configurations for all grids                  :Values                        |
#===============================================:===============================
# Parameters for 1st-level grid -- layer 01     :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     0
 Coordinate System    (0:spherical, 1:cartesian):     0
 Governing Equations  (0:linear,    1:nonlinear):     0
 Grid Size  (dx, sph:minute, cart:meter)        :     0.50
 Time step                            ( second ):     1.0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     0
 X_start                                        :  131.0
 X_end                                          :  154.0
 Y_Start                                        :   33.0
 Y_end                                          :   43.0
 File Name of Bathymetry Data                   :../topo/W130E155S33N45_1min.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    01
 Grid Level                                     :     1
 Parent Grid's ID Number                        :    -1

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 02    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     3
 X_start                                        : 178.011643
 X_end                                          : 178.501284
 Y_start                                        : -39.506430
 Y_end                                          : -39.025147
 FileName of Water depth data                   : depth125m_utm.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    02
 Grid Level                                     :     2
 Parent Grid's ID Number                        :     1

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 03    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     3
 X_start                                        :  178.25
 X_end                                          :  178.45
 Y_start                                        :  -39.45
 Y_end                                          :  -39.25
 FileName of Water depth data                   : depth_nz30sec_LL.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    03
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 04    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :    11.
 X_end                                          :    30.
 Y_start                                        :    11.
 Y_end                                          :    30.
 FileName of Water depth data                   : layer23.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    04
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 05    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :    61.
 X_end                                          :    80.
 Y_start                                        :    61.
 Y_end                                          :    80.
 FileName of Water depth data                   : layer24.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    05
 Grid Level                                     :     3
 Parent Grid's ID number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 06    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   115.
 X_end                                          :   233.
 Y_start                                        :   407.
 Y_end                                          :   573.
 FileName of Water depth data                   : layer31.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    06
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 07    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   140.
 X_end                                          :   233.
 Y_start                                        :   143.
 Y_end                                          :   310.
 FileName of Water depth data                   : layer32.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    07
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 08    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Use Bottom friction ?(only cart,nonlin,0:y,1:n):     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   274.
 X_end                                          :   329.
 Y_start                                        :   143.
 Y_end                                          :   235.
 FileName of Water depth data                   : layer33.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    08
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 09    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :    41.
 X_end                                          :    60.
 Y_start                                        :    41.
 Y_end                                          :    60.
 FileName of Water depth data                   : layer34.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    09
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 10    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   129.
 X_end                                          :   247.
 Y_start                                        :   471.
 Y_end                                          :   588.
 FileName of Water depth data                   : layer41.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    10
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 11    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   858.
 X_end                                          :   968.
 Y_start                                        :   246.
 Y_end                                          :   388.
 FileName of Water depth data                   : layer42.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    11
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 12    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :  2154.
 X_end                                          :  2258.
 Y_start                                        :   671.
 Y_end                                          :   825.
 FileName of Water depth data                   : layer43.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    12
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 13    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :    41.
 X_end                                          :    60.
 Y_start                                        :    41.
 Y_end                                          :    60.
 FileName of Water depth data                   : layer44.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0
 Grid Identification Number                     :    13
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01



