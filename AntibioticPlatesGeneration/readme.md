# Plate creation with reactive and specific sample

Python script destined to perform in a Opentrons robot (OT-2) a dispensation of one or more reactives, for example antibiotics, and the inoculation of microorganisms into these plates.

For that purpose it needs a csv file with several customizable attributes, an OT-2 robot with its app and this script. Besides that, it may need also the json description files in case we are using custom labware.




## Script Structure

Even if the script structure can vary from version to version the main structure of the script will be described here. The newest version of the script can be checked at the end of this page under the section URL GitHub repository)




*Definition of labware context*

We need to extract several times information from the labware, manly names of wells, number of columns, rows, etc.

*Definition of classes*

setted_parameters. Class that will set the parameters in the csv to variables and besides it has a method to process some of the variables to convert them from strings to dictionaries.

*Definition of functions*

 `check_setted_parameters` There is an internal validation after settign the parameters. This is to validate that there are o errors about the variables in the csv as labwares that are not in the system, number of cells in a plate is higher than the wells of that labware, etc.
define_tiprack. Everytime that a pipette does not have enough tips and there is a need to place a deck slot with a tiprack, this function will be called and return the name of the tiprack that will be assocated with the pipette.

`check_tip_and_pick` Function called to check if there are tips, in case that there is not is going to check if we defined in the variables file if we want to replace the tiprack or not. In case that that variable is True, it will pause the protocol, the user needs to change that tiprack and wait until you resume the protocol. If th evariable is set as False, it is going to check if there is space in the deck and in affirmative case, will call define_tiprack, in negative case, it will raise an error. After all this, it will pick up a tip.

`number_tubes_needed` Given a volume per reaction of a reactive, the number of reaction and the maximum volume of the tube in which the reactive is going to be stored will return the number of tubes needed and the number or reactions per tube. This function is used several times when setting the labware. To calculate the number of tubes needed it divides the reactions in different number of tubes until it finds the lower number of tubes in which the volume for that quantity of reactions fits the tube volume. This function does not garantee the lower total number of tubes for the reactives because it divides by reactions and not by total volume needed but this way we ensure that the volume can be picked with at least 1 of the mounted pipettes.
setting_number_plates. Given a number of labware plates and the name of it, it sets the list of empty deck slots and if there is space, it will load the labwares and if there is not enough space, it will raise an error. This function is used one or more times in the setting labware phase of the script

`setting_labware` Given the setted_parameters class (called variables in the script) will return list of the source labware, final labware and the tuberacks where reactives are. For that, it sets these labwares with the setting_number_plate, it calculates the tubes needed for each reactive with number_tubes_needed and calculates how many tuberacks are needed. Besides of this, it updates the variables class with the number and positions of each of the previously named labwares.

`position_dispense_aspirate` The dispense and aspirate action from the 15mL falcon tubes needs to be at a certain  height so the pipette does not get wet. These heights are measured by hand and given a position of the tube and the volume of it the function will return the position with a certain z accordingly to the volume of the falcon. This function is used in the distribution of the reactives in the different final plates.

`distribute_z_tracking` Given the volume to distribute to each well, the position(s) of the reactive tube(s) and the positions to dispense it will transfer the liquid to that wells tracking the height of aspiration by comparing the height before and after the volume aspiration. The position is the result of the position_dispense_aspirate function.

`set_labware` Generator used to give sequentially the positions to dispense from a list of final wells. Each time that this function is called will return the next position of the list.
print_information_user. This function is used at the end of the script to print in the terminal information that is needed for the user as some general information, the needed labware and the reactives positions and volumes in the tuberack(s)

*Body of the script*

1. Read variables from csv file
2. Create the setted_parameters class, calling it variables
3. Process them
4. Fill it with other variables, all extracted from the original ones
5. Validate variables
6. Set pipettes
7. Set labware
8. Distribute reactives
9. Distribute samples from source labware

All of this steps are embedded in a try except structure to catch errors in the process and give them in a more user-friendly output.

Everything in the script is inside of the run function due to the requirement of opentrons, for more information visit the [Opentrons page](https://docs.opentrons.com/v2/writing.html) and its API manual version 2.




## Requirements (packages, software)

Python 3.7 or higher (https://www.python.org/downloads/)

**Python packages:** opentrons, pandas, logging, math and copy

For the installation of the opentrons package you can visit the [opentrons blog page](https://support.opentrons.com/s/article/Simulating-OT-2-protocols-on-your-computer)





## Input(s) - (file format/number)

The script needs only 1 mandatory file, the Variables-AntibioticPlatesCreation-OT.csv

This file can be obtained by filling Template_Variables_PlateGeneration.xlsx and converting it or exporting it as a csv with that name, it is important that the file has that exact name, otherwise the script will not work.

In case that you have a labware error, one of the reasons can be that the OT-App is not recognising the labware. In that case you need to create a directory with the api_name of the labware and a file with the labware description (a json file created by the custom labware creator) called 1.json

This file should be located in the robot system, specifically in the folloring directory /data/labware/v2/custom_definitions/custom_labware



## How to run script
 1. Fill template
 2. Convert to csv
 3. Pass csv (to directory /data/user_storage) and script to OT-2 in which the protocol will be run
 4. Run script in OT and storing the output in a file
 5. Pass that file to your computer
 6. Load script into OT-App
 7. Load labware as stated in the App and the instruction file
 8. Run protocol

For additional information of how to proceed to run the protocol check the [protocols.io entry for this protocol](dx.doi.org/10.17504/protocols.io.q26g7yb3kgwz/v1.)

You can do a test try with the files in the example folder inside this GitHub directory.

## Expected output files (file format/number)

You will only obtain the instruction file in case that you stored the output of the script in a file.

Besides that, you will obtain plates with a reactive an dinoculated with the samples from the source plate(s)




## Parameters of creation (i.e python version)

Python 3.9.7

Microsoft Windows 10 Home 64-bit
