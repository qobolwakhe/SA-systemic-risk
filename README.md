# SA-systemic-risk
This repository contains the code for http://systemicrisk.org.za, a systemic risk ranking of South Africa's financial institutions.

__INSTRUCTIONS__

 The data used by the model is contained in the spreasheet MES_data.xlsx. All information is retrieved from Bloomberg and can be updated by opening the spreadsheet in an instance of Excel with the Bloomberg plug-in installed. The attached spreadsheet already contains data for 108 firms. You can alter the sample being imported into the model, by changing ONLY the tickers in the third row of the 'Share Prices' worksheet. You can also add firms not present in the dataset in the 'All', taking note of the format and structure to be maintained. You will also need to update the 'Ranges' worksheet with the new information and edit the lookup ranges in the relevnt worksheets. 
 
 The current size of the sample to be imported has been set to 16. You can change this by:
 
 1. Updating the spreadsheet as per instructions above
 2. Updating the ranges in the `Reading in data` `(line73)` section of the main.m script
 
The script requires additional functions to run, these can be found in the 'MFE'and 'SysRiskMeasures' folders located in this respository. To add these folders to Matlab's path, you must update the file paths under the `Add functions to path` section `(line 64)` and point the directories to the location where the folders have been saved.

__OUTPUT__

This script will output all data into the 'OUT\' directory. Output will include a csv file for each institution which will contain time series' of the MES, SRISK, SRISK contribution and SRISK ranking. Institutions with no contribution to systemic risk on a given day assigned a rank of 0. 

You an also plot a graph of each institution's MES, SRISK and SRISK contribution over time by changing the value of `j (line253)` to that of the index of the desired intitution in the variable `Series`. 
