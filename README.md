odfinal.py contains code for orbital determination. Code takes in three observation data points, and returns velocity vector and using Method of Gauss.
  The final function is used to output user friendly data in the format of Position Vector, Velocity Vector, Range, Orbital Components, and Mean Anomaly.
  The finalmonte function outputs data used in montecarlo.py
montecarlo.py was the function used for generating Monte Carlo Data. 
  Code contains constants used in the determination of 2253 Espinette. 
  Outputs graph.



odfinal.py is a code to determine orbit given three observations. It contains all helper functions that are needed and takes in a txt file. 
Each line should contain an observation. If you input more than three lines, the program will ask for which of the observations you would like to select.
Each input line should contain date/time/RA/Dec/R_Vector(Sun Vector) in the form:
Year Month Day UTC Time RA Dec  R_x, R_y, R_z
ie. 2021 06 25 00:00:00.000 12:18:12.50 +10:11:28.0  -5.996994128131034E-02, 9.309517649933454E-01, 4.035390946349719E-01
To skip to the main function, you may skip to line 176.
The output contains position vector, velocity vector, range, orbital components(a,e,i,OMEGA,omega,M), and Mean Anomaly for July 22, 2025 @ 6:00:00 UTC in the ecliptic.
Since the function uses an iterative loop, a difference is also printed per cycle in case of an error in converging the difference to 10e-10.
A txt file will be included for an example for input. 
Note that Light Speed has been accounted for.
