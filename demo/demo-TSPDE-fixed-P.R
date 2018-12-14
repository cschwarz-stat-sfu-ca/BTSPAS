
# 2018-12-15 CJS demonstration of fixeding logitP for some entries (made up data)
# 2014-09-01 CJS removed all of the prompts
# 2010-11-25
# This is a demonstration of how to call the Time Stratified Petersen with Diagonal Entries (TSPDE) program
# It is based on the analysis of California Junction City 2003 Chinook data.
# It is slightly different from the data used in the Trinity River Project because I
# (arbitrarily) set some bad.n1 and bad.m2 values to test this feature. 
# If you want to run the Trinity example, remove the bad.m2 and bad.n1 assignment statements.
#
# In each julian week j, n1[j] are marked and released above the rotary screw trap.
# Of these, m2[j] are recaptured. All recaptures take place in the week of release, i.e. the matrix of
# releases and recoveries is diagonal.
# The n1[j] and m2[j] establish the capture efficiency of the trap.
#
# At the same time, u2[j] unmarked fish are captured at the screw trap.
# The simple stratified Petersen estimator would inflate the u2[j] buy 1/capture efficiency[j]
# giving U2[j] = total fish passing the trap in julian week [j] = u2[j] * n1[j]/ m2[j].
#
# The program assumes that the trap was operating all days of the week. The sampfrac[j] variable
# gives the proportion of days the trap was operating. For example, if the trap was operating for 3 of the
# 7 days in a week, then sampfrac[j]<- 3/7
#
#
# Notes:
#    - the number of recaptures in sample week 33 (julian week 41) is far too low. 
#      This leads to an estimate of almost 13 million fish from the simple stratified Petersen. 
#      Consequently, the recaptures for this
#      week are set to missing and the program will interpolate the number of fish for this week
#
#    - the number of days operating is 8 in sample weeks 2 (julian week 10) and 
#      6 in sample week 3 (julian week 11). The 8 days in sample week 2 is "real" as
#      the code used on the marked fish was used for 8 days. The program will automatically 
#      "reduce" the number of unmarked fish captured in this week to a "7" day week 
#      and will increase the number of unmarked fish captured in week 3 to "7" days as well. 
# 
#  The program tries to fit a single spline to the entire dataset. However, in julian weeks
#  23 and 40, hatchery released fish started to arrive at the trap resulting in sudden jump
#  in abundance. The jump.after vector gives the julian weeks just BEFORE the suddent jump,
#  i.e. the spline is allowed to jump AFTER the julian weeks in jump.after.
#
#  The vector bad.m2 indicates which julian weeks something went wrong. For example, the
#  number of recoveries in julian week 41 is far below expectations and leads to impossible
#  Petersen estimate for julian week 41.
#  
#  The vector bad.u2 indicates which julian weeks, the number of unmarked fish is suspect.
#  I arbitrarily chose the third julian week to demonstrate this feature.
# 
#  The prefix is used to identify the output files for this run.
#  The title  is used to title the output.

library("BTSPAS")

# Create a directory to store the results
if(file.access("demo-TSPDE-fixed-P")!=0){ dir.create("demo-TSPDE-fixed-P", showWarnings=TRUE)}  # Test and then create the directory
setwd("demo-TSPDE-fixed-P")

# Get the data. In many cases, this is stored in a *.csv file and read into the program
# using a read.csv() call. In this demo, the raw data is assigned directly as a vector.
#

demo.data.csv <- textConnection("
jweek   ,  n1   ,    m2   ,      u2   , sampfrac
 9   ,      0   ,     0   ,    4135   ,   3
10   ,   1465   ,    51   ,   10452   ,   8
11   ,   1106   ,   121   ,    2199   ,   6
12   ,    229   ,    25   ,     655   ,   7
13   ,     20   ,     0   ,     308   ,   7
14   ,    177   ,    17   ,     719   ,   7
15   ,    702   ,    74   ,     973   ,   7
16   ,    633   ,    94   ,     972   ,   7
17   ,   1370   ,    62   ,    2386   ,   7
18   ,    283   ,    10   ,     469   ,   7
19   ,    647   ,    32   ,     897   ,   7
20   ,    276   ,    11   ,     426   ,   7
21   ,    277   ,    13   ,     407   ,   7
22   ,    333   ,    15   ,     526   ,   7
23   ,   3981   ,   242   ,   39969   ,   7
24   ,   3988   ,    55   ,   17580   ,   7
25   ,   2889   ,   115   ,    7928   ,   7
26   ,   3119   ,   198   ,    6918   ,   7
27   ,   2478   ,    80   ,    3578   ,   7
28   ,   1292   ,    71   ,    1713   ,   7
29   ,   2326   ,   153   ,    4212   ,   6
30   ,   2528   ,   156   ,    5037   ,   7
31   ,   2338   ,   275   ,    3315   ,   7
32   ,   1012   ,   101   ,    1300   ,   7
33   ,    729   ,    66   ,     989   ,   7
34   ,    333   ,    44   ,     444   ,   7
35   ,    269   ,    33   ,     339   ,   7
36   ,     77   ,     7   ,     107   ,   7
37   ,     62   ,     9   ,      79   ,   7
38   ,     26   ,     3   ,      41   ,   7
39   ,     20   ,     1   ,      23   ,   7
40   ,   4757   ,   188   ,   35118   ,   7
41   ,   2876   ,     8   ,   34534   ,   7
42   ,   3989   ,    81   ,   14960   ,   7
43   ,   1755   ,    27   ,    3643   ,   7
44   ,   1527   ,    30   ,    1811   ,   7
45   ,    485   ,    14   ,     679   ,   7
46   ,    115   ,     4   ,     154   ,   5")

demo.data <- read.csv(demo.data.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
demo.data$sampfrac <- demo.data$sampfrac/ 7


# After which weeks is the spline allowed to jump?
demo.jump.after <- c(22,39)  # julian weeks after which jump occurs

# Which julian weeks have "bad" recapture values. These will be set to missing and estimated.
demo.bad.m2     <- c(41)   # list julian weeks with bad m2 values. This is used in the Trinity Example
demo.bad.u2     <- c(11)   # list julian weeks with bad u2 values. [This was arbitrary to demostrate the feature.]

demo.logitP.fixed <- c(37,38,46)
demo.logitP.fixed.values <- rep(-10, length(demo.logitP.fixed))
demo.data$u2[ demo.data$jweek %in% demo.logitP.fixed] <- 0   # no captures when capture probability is 0
demo.data$m2[ demo.data$jweek %in% demo.logitP.fixed] <- 0

# The prefix for the output files:
demo.prefix <- "demo-JC-2003-CH-TSPDE" 

# Title for the analysis
demo.title <- "Junction City 2003 Chinook "

cat("*** Starting ",demo.title, "\n\n")

# Make the call to fit the model and generate the output files
demo.tspde.logitP.fixed <- TimeStratPetersenDiagError_fit(
                  title  =demo.title,
                  prefix =demo.prefix,
                  time   =demo.data$jweek,
                  n1     =demo.data$n1, 
                  m2     =demo.data$m2, 
                  u2     =demo.data$u2,
                  sampfrac=demo.data$sampfrac,
                  jump.after=demo.jump.after,
                  bad.n1 =demo.bad.n1,
                  bad.m2 =demo.bad.m2,
                  bad.u2 =demo.bad.u2,
                  logitP.fixed=demo.logitP.fixed, logitP.fixed.values=demo.logitP.fixed.values,
                  InitialSeed=890110,
                  debug=TRUE,  # this generates only 10,000 iterations of the MCMC chain for checking.
                  save.output.to.files=TRUE)

# Rename files that were created.

file.copy("data.txt",       paste(demo.prefix,".data.txt",sep=""),      overwrite=TRUE)
file.copy("CODAindex.txt",  paste(demo.prefix,".CODAindex.txt",sep=""), overwrite=TRUE)
file.copy("CODAchain1.txt", paste(demo.prefix,".CODAchain1.txt",sep=""),overwrite=TRUE)
file.copy("CODAchain2.txt", paste(demo.prefix,".CODAchain2.txt",sep=""),overwrite=TRUE)
file.copy("CODAchain3.txt", paste(demo.prefix,".CODAchain3.txt",sep=""),overwrite=TRUE)
file.copy("inits1.txt",     paste(demo.prefix,".inits1.txt",sep=""),    overwrite=TRUE)
file.copy("inits2.txt",     paste(demo.prefix,".inits2.txt",sep=""),    overwrite=TRUE)
file.copy("inits3.txt",     paste(demo.prefix,".inits3.txt",sep=""),    overwrite=TRUE)

file.remove("data.txt"       )       
file.remove("CODAindex.txt"  )
file.remove("CODAchain1.txt" )
file.remove("CODAchain2.txt" )
file.remove("CODAchain3.txt" )
file.remove("inits1.txt"     )
file.remove("inits2.txt"     )
file.remove("inits3.txt"     )
 
# save the results in a data dump that can be read in later using the load() command.
# Contact Carl Schwarz (cschwarz@stat.sfu.ca) for details.
save(list=c("demo.tspde.logitP.fixed"), file="demo-tspde-logitP-fixed-saved.Rdata")  # save the results from this run
cat("\n\n\n ***** FILES and GRAPHS saved in \n    ", getwd(), "\n\n\n")
print(dir())

# move up the directory
setwd("..")

cat("\n\n\n ***** End of Demonstration *****\n\n\n")

