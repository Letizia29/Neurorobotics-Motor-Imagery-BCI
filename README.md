# Neurorobotics-Motor-Imagery-BCI
Anlysis of data collected during a 3-day Motor Imagery (MI)  Brain-Computer Interface (BCI) experiment involving 8 healthy participants.

## Data description:
The data was recorded using a 16-channel EEG amplifier (g.USBamp, g.Tec) at a sampling rate of 512 Hz. Electrodes were positioned according to the 10-20 international system. Each participant completed at least two recording days:
● Day 1: 3 "offline" runs (calibration, without real feedback) and 2 "online" runs (with real feedback).
● Day 2 and Day 3: 2 "online" runs per day.

## The task and the visual paradigm:
Participants performed two motor imagery tasks—imagining movements of both hands or both feet—and a rest task. The color of the cue indicated which motor imagery task to perform (e.g., both hands, both feet, or rest).
● During the calibration runs, feedback associated with the cue automatically moved in the correct direction.
● During the online runs, feedback movement was determined by the output of the classifier.

## The analysis:
1. Grand average analyses on the whole population and on representative subjects
   a. Process the data using appropriate methods.
   b. Identify and extract the most relevant features.
   c. Report the achieved results
2. Analyses on BMI decoding on each subject
   a. Calibration phase:
    ▪ Consider only the offline runs
    ▪ Process the data, compute the features, select the most discriminant features
    ▪ Create a classifier based on those features
   b. Evaluation phase:
    ▪ Consider only the online runs
    ▪ Process the data, compute the features, and extract those already selected during the calibration phase
    ▪ Use this data to evaluate the classifier created during the calibration phase
    ▪ Implement and apply an evidence accumulation framework on the posterior probabilities
