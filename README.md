# WPFG-database
Collaborative project to collate existing water plant functional group assignments to species or genus levels 
using the 7- or 10-category classification system.

## To update database
Updating the database to add a new set of WPFG classifications involves four steps: 

1. update the master species list. 
Follow the data format of "wpfg_lists_Ver1.csv" (in folder data.in) to add your records. 

This is a four column comma separated value file with the following columns:

original_name = the scientific name used in your database to identify the taxon. 
WPFG = the water plant functional group you have assigned to the taxon.
source = an identifier for your data set (e.g., manuscript first author and year  "Smith_2026")
chrono = an integer identifier for your work. Assign this as the next value after the highest in the list.

Add your species assignments and then save the file with a new version number (This will be the database version).

2. Update the Australian Plant Census database version.
The database automatically aligns taxon names to aligns names to the Australian Plant Census (APC) and 
Australian Plant Name Index (APNI). The database records the latest stable versions of these resources, 
used. To identify this resource, load the APCalign package and run: 

'default_version()'

This value is updated in the next step so it is assigned to the revised WPFG database.


3. Update the run file. 
Now open 'run_all.R' and update the version number (next integer number) and the apc database. 

Save this file with the same name.

4. Re-run the code. 
the file 'run_all.R' re-creates the entire workflow. Open this script, select all and run. 

The output files will overwrite in the data_out folder. The user database is in the file "WPFG_user_database_Ver1.xlsx"

