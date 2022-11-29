# PCR samples preparation
Python script destined to perform in a Opentrons robot (OT-2) a PCR mix preparation, i.e., to prepare well in which there is polymerase mix, primers, water and colony.

For that purpose the script also needs a csv file with several customizable attributes, an OT-2 robot with its app and this script. Besides that, it may need also the json description files in case we are using custom labware, but this is optional.




## Script Structure

Even if the script structure can vary from version to version the main structure of the script will be described here. The newest version of the script can be checked at the end of this page under the section URL GitHub repository.




*Definition of labware context*

We need to extract several times information from the labware, manly names of wells, number of columns, rows, etc.


*Definition of classes*

`setted_parameters` Class that will set the parameters in the csv to variables and besides it has a method to process some of the variables to convert them from strings to dictionaries.

`calculated_volumes` In this class we will calculate and store the values of the reactives per reaction and tube. This way, we can access it easily during the protocol.

`positions` Class to store the positions of reactives and labware, mainly used for transfering and printing the information parts of the script.


*Definition of functions*

`check_setted_parameters` There is an internal validation after settign the parameters. This is to validate that there are o errors about the variables in the csv as labwares that are not in the system, number of cells in a plate is higher than the wells of that labware, etc.
define_tiprack. Everytime that a pipette does not have enough tips and there is a need to place a deck slot with a tiprack, this function will be called and return the name of the tiprack that will be assocated with the pipette.

`check_tip_and_pick` Function called to check if there are tips, in case that there is not is going to check if we defined in the variables file if we want to replace the tiprack or not. In case that that variable is True, it will pause the protocol, the user needs to change that tiprack and wait until you resume the protocol. If th evariable is set as False, it is going to check if there is space in the deck and in affirmative case, will call define_tiprack, in negative case, it will raise an error. After all this, it will pick up a tip.

`number_tubes_needed` Given a volume per reaction of a reactive, the number of reaction and the maximum volume of the tube in which the reactive is going to be stored will return the number of tubes needed and the number or reactions per tube. This function is used several times when setting the labware. To calculate the number of tubes needed it divides the reactions in different number of tubes until it finds the lower number of tubes in which the volume for that quantity of reactions fits the tube volume. This function does not garantee the lower total number of tubes for the reactives because it divides by reactions and not by total volume needed but this way we ensure that the volume can be picked with at least 1 of the mounted pipettes.

`setting_number_plates` Given a number of labware plates and the name of it, it sets the list of empty deck slots and if there is space, it will load the labwares and if there is not enough space, it will raise an error. This function is used one or more times in the setting labware phase of the script

`setting_reactives_in_places` Given the number of tubes, the name of the reactive and the possible positions in the labware the function will return a dictionary with the positions of that reactive in the labware for further use.

`setting_labware` Given the setted_parameters class (called variables in the script) will return list of the source labware, final labware and the tuberacks where reactives are. For that, it sets these labwares with the setting_number_plate, it calculates the tubes needed for each reactive with number_tubes_needed and calculates how many tuberacks are needed. Besides of this, it updates the variables class with the number and positions of each of the previously named labwares.

`z_positions_mix` Given a volume, the function will return based on that volume 3 heights (z positions) that will serve as bottom, medium and max position for the mixing positions based always in a 1.5mL eppendorf (it is tunned for that kind of tubes, it will not work as efficient with other volume tubes). This function will be used in the function mixing_eppendorf_15.

`distribute_vol_tracking` Given the volume to distribute to each well, the position(s) of the reactive tube(s) and the positions to dispense it will transfer the liquid to that wells tracking the volume of each tube. In this case this process is easy because of the kind of separation that we did of the total volume in the tubes. We have the reactions per tube, so we only need to transfer rthat many reactions.

`set_labware` Generator used to give sequentially the positions to dispense from a list of final wells. Each time that this function is called will return the next position of the list.

`give_me_optimal_pipette_to_use` Given the pipettes that are mounted in the robot arm and the volume of transferring it return the pipette that will perform the less movements to transfer that volume. In case that only one of the pipettes can transfer the volume, the function will return that one.

`pds_labware_source` Given the name of the source labware it will return a dataframe with the same dimensions so it can be filled after with the names of the wells where the samples are placed.

`print_information_user` This function is used at the end of the script to print in the terminal information that is needed for the user as some general information, the needed labware and the reactives positions and volumes in the tuberack(s)

`map_tracking` Given the source and final wells, the final table and the map of the source plate, the function will process the names and insert them in the corresponding cell in the final map.

`delete_certain_wells_from_list` We are going to delete a list of wells from a labware, i.e., we have 3 lists, one of all the wells, another that we want to remove and the third one is all the labware that we are working with. This is use in this script to remove the controls and the wells that the user has stated that they do not want to take from.

`mixing_eppendorf_15` This function will mix the liquid that is in the tube mixing at different heights, stir a little bit and transfering liquid from different heights. This mixing it is tunned to be perform with less than 1500uL using the function z_positions_mix.

`transfer_volume_tracking` This function will transfer reactives from one placve to another taking in account the pipettes an dbeing able to change from one to another in case that is more favorable to use the other pipette. In this script this function is used to create the PCR mix in the eppendorfs.

*Body of the script*

1. Read variables from csv file
2. Create the setted_parameters class, calling it variables
3. Validate the labware
4. Processing some variables
5. Validating rest of the variables
6. Define labware
7. Define volumes needed of each reactives
8. Define pipettes
9. Transfer reactives (water, primers and then polymerase)
10. Mixing and then distribute set of mixes
11. Define the wells of source plates that we have to transfer and the order (controls at the end)
12. Create the maps
13. Distribute samples from source labware
14. Export map if that is settled in the csv of variables
15. Print user information

All of this steps are embedded in a try except structure to catch errors in the process and give them in a more user-friendly output.

Everything in the script is inside of the run function due to the requirement of opentrons, for more information visit the [Opentrons page](https://docs.opentrons.com/v2/writing.html) and its API manual version 2.




## Requirements (packages, software)


Python 3.7 or higher (https://www.python.org/downloads/)

**Python packages:** opentrons, pandas, logging, math, numpy and copy

For the installation of the opentrons package you can visit the [opentrons blog page](https://support.opentrons.com/s/article/Simulating-OT-2-protocols-on-your-computer)





## Input(s) - (file format/number)


The script needs only the *Variables-PCRs-OT.csv*.

This file can be obtained by filling *Template_Variables_PCR.xlsx* and converting it or exporting it as a csv with that name, it is important that the file has that exact name, otherwise the script will not work.

In case that you have a labware error, one of the reasons can be that the OT-App is not recognising the labware. In that case you need to create a directory with the api_name of the labware and a file with the labware description (a json file created by the [custom labware creator](https://labware.opentrons.com/create/)) called *1.json*

This file should be located in the robot system, specifically in the folloring directory */data/labware/v2/custom_definitions/custom_labware*



## How to run script

1. Fill template
2. Convert to csv
3. Pass csv to directory */data/user_storage* and script to OT-2 in which the protocol will be run
4. Run script in OT and storing the output in a file
5. Pass that file to your computer and the map(s) that will be stored in the directory */data/user_storage*
6. Load script into OT-App
7. Load labware as stated in the App and the instruction file
8. Run protocol

For additional information of how to proceed to run the protocol check the [protocols.io page of the protocol](dx.doi.org/10.17504/protocols.io.n92ldpyznl5b/v1)

Besides, you can use the files in the example folder of this GitHub directory to test how the script works.



## Expected output files (file format/number)


You will obtain the instruction file in case that you stored the output of the script in a file and a map of each source plate with a list of all the positions the samples are in the final PCR plates.


Besides that, you will obtain plates with the samples and the PCR mix with the respective primers for all of them, i.e., we will have all the PCR set mixes inoculated with the samples, one per each combination set-sample.





## Parameters of creation (i.e python version)


Python 3.9.7

Microsoft Windows 10 Home 64-bit
