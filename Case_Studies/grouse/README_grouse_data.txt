The grouse_data.csv contains the survival data for chick and juvenile Columbian sharp-tailed grouse in north-central Colorado based on radio-collar data collected during three spring/summer seasons of 2014 through 2016.

Individuals with a mortality event have two rows within the data, while individuals that are right censored have a single row. 

Column names define these data as follows:

left.time = the entry period (week) of the study, where 1 represents the first week of the study
right.time = the exit period (week)
censor = 1/0, whether an event is a censored event (1), or a mortality event (0)
id = unique identifier for each collared individuals
year.cap1 = an indicator whether an individual was captured and collared during the first year of the study (2014)
year.cap2 = an indicator whether an individual was captured and collared during the second year of the study (2015)
year.cap3 = an indicator whether an individual was captured and collared during the third year of the study (2016)
left.age = the age (in weeks) when the individual enters the study at left.time period
right.age = the age (in weeks) when the individual exits the study, either through mortality or right censoring
brood = the id number of the brood that the chick or juvenile were hatched into

w
