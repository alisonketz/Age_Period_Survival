The wtd_data.csv contains the survival data for adult female white-tailed deer in south-central Wisconsin based on GPS collar data collected from January 7,2017 through May 4th, 2021.

Deer with a mortality event have two rows within the data, while deer that are right censored have a single row. 

Column names define these data as follows:

left.time = the entry period (week) of the study, where 1 represents the first week of the study
right.time = the exit period (week)
censor = 1/0, whether an event is a censored event (1), or a mortality event (0)
id = unique identifier for each collared deer
yearCapture = year that the deer was captured and GPS collared
ageCapture = the age (in weeks) when the individual enters the study at left.time period
ageRight = the age (in weeks) when the individual exits the study, either through mortality or right censoring

 


