Data dictionaries for data sets not created by setup_data.R

CMV_primary_infection_data.csv - infection viral loads for infants after primary CMV infection

|**Variable name** | **Description** | **Type** | **Categories (if applicable)** |
|------------------|-----------------|------------------|-----------------|
|**PatientID2** | Randomized patient ID | chr | -P denotes infants (all)| 
|**days** | time since first observation (days) | num | |
|**times** | observation date | date | |
|**Virus** | Virus ID | chr | ORL_CMV|
|**infantInfection** | infant infection in household | bool | |
|**infantInfDate** | date of infant infection in household | date | |
|**count** | log10 viral concentration | num | |


demographics_data.csv - some time series metadata for breastfeeding, saliva sharing, and food sharing during the week prior to the date.

|**Variable name** | **Description** | **Type** | **Categories (if applicable)** |
|------------------|-----------------|------------------|-----------------|
|**PatientID2** | Randomized patient ID | chr | -P denotes infants (all)| 
|**times** | observation date | date | |
|**breastfed** | breastfeeding in the previous week | num | 0 = no, 1 = yes|
|**saliva** | saliva sharing recorded between mom and infant previous week|num | 0 = no, 1 = yes|
|**chewfood** | mom recorded chewing food for infant in previous week | num | 0 = no, 1 = yes|


raw_data.csv - time series viral load data for all infants for all viruses during the study.

|**Variable name** | **Description** | **Type** | **Categories (if applicable)** |
|------------------|-----------------|------------------|-----------------|
|**PatientID2** | Randomized patient ID | chr | -P denotes infants (all)| 
|**days** | time since first observation (days) | num | |
|**times** | observation date | date | |
|**Virus** | Virus ID | chr | ORL_HSV, ORL_CMV, ORL_EBV, ORL_HHV6, ORL_HHV8|
|**infantInfection** | infant infection in household | bool | |
|**infantInfDate** | date of infant infection in household | date | |
|**count** | log10 viral concentration | num | |