# SA-systemic-risk
This repository contains the code for http://systemicrisk.org.za, a systemic risk ranking of South Africa's financial institutions.

__INSTRUCTIONS__

 The data used by the model is contained in the spreasheet MES_data.xlsx. All information is retrieved from Bloomberg and can be updated by opening the spreadsheet in an instance of Excel with the Bloomberg plug-in installed. The attached spreadsheet already contains data for 108 firms. You can alter the sample being imported into the model, by changing ONLY the tickers in the third row of the 'Share Prices' worksheet. You can also add firms not present in the dataset in the 'All', taking note of the format and structure to be maintained. You will also need to update the 'Ranges' worksheet with the new information and edit the lookup ranges in the relevnt worksheets. 
 
 The current size of the sample to be imported has been set to 16. You can change this by:
 1) Updating the spreadsheet as per instructions above
 2) Updating the ranges in the 'Reading in data' (line73) section of the main.m script
 
