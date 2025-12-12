# Climate Data

Climate Data

## Format

Scientists are interested in the Earth's temperature change since the
last glacial maximum, about 20,000 years ago. The first study to
estimate the temperature change was published in 1980, and estimated a
change of -1.5 degrees C, +/- 1.2 degrees C in tropical sea surface
temperatures. The negative value means that the Earth was colder then
than now. Since 1980 there have been many other studies. `climate` is a
dataset with 63 measurements on 5 variables:

- *sdev*:

  a standard deviation for the calculated *deltaT*;

- *proxy*:

  a number 1-8 reflecting which type of measurement system was used to
  derive deltaT. Some proxies can be used over land, others over water.
  The proxies are coded as  
  1 "Mg/Ca"  
  2 "alkenone"  
  3 "Faunal"  
  4 "Sr/Ca"  
  5 "del 180"  
  6 "Ice Core"  
  7 "Pollen"  
  8 "Noble Gas"  

- *T/M*:

  , an indicator of whether it was a terrestrial or marine study (T/M),
  which is coded as 0 for Terrestrial, 1 for Marine;

- *latitude*:

  the latitude where the data were collected.

## Source

Data provided originally by Michael Lavine
