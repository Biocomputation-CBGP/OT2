
"""
Este sera un script con el shaker y con el termociclador, con su perfil, con el shaker y con las variables que habra que dar
Este ya se ha probado y funciona biologicamente (se probo con 7 mixes)
Cosas que hay que arreglar y cambair lo siguiente:
 - lista de las posicines que esten bien ya que como lo pongo yo no coge que puede poner en el 9 pero no en el 7 y 8
 - hay que arreglar la posicion del coldblock tambien, ya qu eno se puede poner alrededor del shaker y del termociclador y si hace falta ponerlo
poner un wanring de que tiene squ eestar seguro de que tu "coldblock" se ajuste bine con los termocicladores
 - Hay que probar lo de los tipracks de 1000 con los modulos por lo que nos dijo sofia, para saber si te alerta o no
 - Hacer controles adecuados al shaker
 - Hacer el mapa correcto de las PCRs, que se parezca al de los demas scripts y no como lo estamos dando en este
 - Subrilo a github
Este script no va a ser el de opentrons, va a ser el de la version 2
"""




## Packages needed for the running of the protocol
import opentrons.execute
from opentrons import protocol_api
import pandas as pd
import numpy as np
import logging
import math
import copy
import os

## Structure needed for the running of th escript in the OT
metadata = {
'apiLevel':'2.13' #For the heater-shaker needs to be >=2.13
}

def run(protocol: protocol_api.ProtocolContext):
    # Before everything because we will need it for the labware setting
    labware_context = opentrons.protocol_api.labware
    
    # We are going to create a logger to control errors that can occure in the script that are not taked in account and some of ours that will quit the program
    logging.basicConfig(format="""-------------------------------------------------------------------------------
--> SUM ERROR:
%(message)s
-------------------------------------------------------------------------------
--> INFO ERROR:""",
                        level=logging.ERROR)
    logger = logging.getLogger()

# Classes definitions
# ----------------------------------
# ----------------------------------
    
    class setted_parameters():
        """
        Class that will contain the parameters setted in the variables csv and will process them to work easily in the rest of the protocol
        The coding of this function is dependant of the variables in the Template of the protocol and the names have to be consistent with the rest of the code
        """
        def __init__(self, modular_parameters):
            self.name_coldblock = modular_parameters["API name of ColdBlock"]["Value"] # Needs to have 12 positions (3 rows and 4 columns) and take eppednorfs of 1.5mL
            self.number_primers = int(modular_parameters["Number of primers per set"]["Value"])
            self.sets = int(modular_parameters["Number of primer sets"]["Value"])
            self.volume_reactives = float(modular_parameters["Final volume of reactives (without colony)"]["Value"])
            self.volume_primer = float(modular_parameters["Volume of each primer transfer (uL)"]["Value"])
            self.volume_polymerase = float(modular_parameters["Volume of polymerase transfer (uL)"]["Value"])
            self.name_source_plate = modular_parameters["Source Plate (Name)"]["Value"]
            self.name_final_plate = modular_parameters["Final Plate (Name)"]["Value"]
            self.volume_colonie = float(modular_parameters["Volume of colony Transfer (uL)"]["Value"])
            self.volume_water = self.volume_reactives - self.volume_polymerase - self.number_primers*self.volume_primer
            if modular_parameters["Optional Map"]["Value"].lower() == "false":
                self.map = False
                self.name_map = None
            elif modular_parameters["Optional Map"]["Value"].lower() == "true":
                self.map = True
                self.name_map = "/data/user_storage/"+modular_parameters["Name of map"]["Value"]
                # self.name_map = modular_parameters["Name of map"]["Value"]
            self.number_source_plates = int(modular_parameters["Number of source plates"]["Value"])
            self.number_samples_source_plates = modular_parameters["Number of samples in every source plate"]["Value"]
            self.index_start_source_plate = modular_parameters["Index of start cell in source plate"]["Value"]
            self.index_start_final_plate = modular_parameters["Index of start cell in final plate"]["Value"]
            self.number_colonies = int(modular_parameters["Number of colonies to perform PCR"]["Value"])
            if self.index_start_final_plate == "-":
                self.index_start_final_plate = 0
            else:
                self.index_start_final_plate = int(self.index_start_final_plate)
            self.positions_not_perform_pcr = modular_parameters["Position colonies not to do PCR"]["Value"]
            self.positions_control = modular_parameters["Position of controls"]["Value"]
            self.number_controls = int(modular_parameters["Number of controls"]["Value"])
            self.name_right_pip = modular_parameters["Pipette Right Mount"]["Value"]
            self.name_left_pip = modular_parameters["Pipette Left Mount"]["Value"]
            self.starting_tip_right_pip = modular_parameters["Starting tip Pipette Right"]["Value"]
            self.starting_tip_left_pip = modular_parameters["Starting tip Pipette Left"]["Value"]
            self.extra_pip_factor = float(modular_parameters["Extra Pipetting Adjust Factor"]["Value"])+1.0
            self.primer_volume_per_set = self.volume_primer*self.number_primers # volume of each primer * number of primers per set
            if modular_parameters["Replace Tiprack"]["Value"].lower() == "false":
                self.replace_tiprack = False
            elif modular_parameters["Replace Tiprack"]["Value"].lower() == "true":
                self.replace_tiprack = True
            else:
                self.replace_tiprack = False
            self.need_change_tiprack = False
            self.number_reactions = self.number_colonies+self.number_controls
            self.right_pipette = None
            self.left_pipette = None
            self.warnings = []
            self.errors = []
            
            # Variables for the new type of map
            self.maps_source_plates_names = modular_parameters["Name Source Plates Maps"]["Value"][1:-1].replace(" ","").split(",") # We have the list of the maps files
            self.maps_source_plate = {} # We are going to fill this dictionary with the dataframes corresponding to the files
            for index_map, map_source in enumerate(self.maps_source_plates_names):
                self.maps_source_plate[index_map] = pd.read_csv("/data/user_storage/"+map_source+".csv", index_col = 0)
                # self.maps_source_plate[index_map] = pd.read_csv(map_source+".csv", index_col = 0)

            # Estos son variables del shaker que son customizables
            if modular_parameters["Heater-Shaker"]["Value"].lower() == "true":
                self.shaker = True
            elif modular_parameters["Heater-Shaker"]["Value"].lower() == "false":
                self.shaker = False
            else:
                raise Exception("'Heater-Shaker' variable can only be True or False")
            
            if modular_parameters["RMP Shaker"]["Value"] == "-" or modular_parameters["RMP Shaker"]["Value"].lower() == "none":
                self.revolutions = None
            else:
                self.revolutions = int(modular_parameters["RMP Shaker"]["Value"])
            self.hs_mod = None # Initial value
            self.name_labware_shaking = modular_parameters["API Name Labware in Heater-Shaker"]["Value"]
            self.labware_shaking = None # Initial value
            
            # Estos son variables del termociclador que son customizables
            if modular_parameters["Thermocycler"]["Value"].lower() == "true":
                self.thermocycler = True
            elif modular_parameters["Thermocycler"]["Value"].lower() == "false":
                self.thermocycler = False
            else:
                raise Exception("Thermocycler variable can only be True or False")
            self.tc_mod = None # Initial value
            self.thermocycler_profile_file = "/data/user_storage/"+modular_parameters["Program Thermocycler"]["Value"]+".csv"
            if modular_parameters["Final Open Lid"]["Value"].lower() == "true":
                self.final_open_lid = True
            elif modular_parameters["Final Open Lid"]["Value"].lower() == "false":
                self.final_open_lid = False
            elif modular_parameters["Final Open Lid"]["Value"].lower() == "-":
                self.final_open_lid = None
            else:
                raise Exception("'Final Open Lid' variable can only be True, False or '-'")
            if modular_parameters["Hold Block Temperature"]["Value"] == "-":
                self.final_block_temperature = None
            else:
                self.final_block_temperature = float(modular_parameters["Hold Block Temperature"]["Value"])
            self.lid_temperature = 0 # Initial value
            if modular_parameters["Pause before Temperature Program"]["Value"].lower() == "true":    
                self.pause_tc = True
            elif modular_parameters["Pause before Temperature Program"]["Value"].lower() == "false":
                self.pause_tc = False
            elif modular_parameters["Pause before Temperature Program"]["Value"] == "-":
                self.pause_tc = None
            else:
                raise Exception("'Pause before Temperature Program' variable can only be True, False or '-'")
            
            self.positions = {key: None for key in range(1,12)} # For controling if we can put more labware or not
            return
            
        def proccess_variables(self, name_parameter, value_default_empty, must_entries):
            """
            This function will process some text that is extracted from the variables and make it to a dictionary so we can work with it
            Some of the variables need processing to be usable in the script, usually is from a "(,),(,)" format to a dictionary
            
            The function need the name of the parameter in the class, default value of the variable and how many items the dictiony needs to have
            In case that the strcutre of the parameters are different, this function should be changed
            """
            
            value_parameter = getattr(self, name_parameter).replace(" ","")[1:-1].split('),(')
            
            # We separate the text value into the real values (they were separated by commas)
            for i, element in enumerate(value_parameter):
                    value_parameter[i] = element.split(",")
            
            dict_final = {}
            if len(value_parameter[0][0]) != 0: # This checks if the the variable has a default value (like 0 or -) or if it has other values
                                                # This works with the process that we have done in  the previous part
                for well, plate in value_parameter:
                    try:
                        dict_final[int(plate)] += [well]
                    except:
                        dict_final[int(plate)] = [well]
            else:
                pass
            
            # We make sure that we have the number of entries that are specified in must_entries and if there is one missing we fill it with the default value for the item
            # If the dictionary is empty (when value in template is the default one) we will cretae all the items with the default value of the items (for example can be the number of wells of a plate or the index)
            for i in range(1, must_entries+1):
                try:
                    dict_final[i]
                except:
                    dict_final[i] = [str(value_default_empty)]
            
            setattr(self, name_parameter, dict_final)
            return
        
        def check_setted_parameters(self):
            """
            Function that will check the variables of the Template and will raise errors that will crash the OT run
            It is a validation function of the variables checking errors or inconsistencies
    
            This function is dependant again with the variabels that we have, some checks are interchangable between protocols, but some of them are specific of the variables
            """
            
            # Volume of reactives is larger than the established one
            if (self.volume_primer*self.number_primers + self.volume_polymerase) > self.volume_reactives:
                self.errors.append("Volume of each reactive added is larger than the total volume of reactives")
            
            # Is there any space for water
            elif (self.volume_primer*self.number_primers + self.volume_polymerase) == self.volume_reactives:
                self.warnings.append("Volume of water in each reaction will be 0 but there will be one tube in the coldblock for it")
            
            # Number of controls equal to number of positions of controls
            number_positions_control = 0
            for list_controls_plate in self.positions_control.values():
                if "None" in list_controls_plate:
                    pass
                else:
                    number_positions_control += len(list_controls_plate)
            if number_positions_control == self.number_controls:
                pass
            else:
                self.errors.append("Number of controls and the number of positions given are not equal")
            
            # Check that the variables related to the source plate are correct: that there is as number of plates as the variable establish it and that it goes in the range(1, number_source_plates+1)
            if len(self.number_samples_source_plates.keys()) != self.number_source_plates or sorted(self.number_samples_source_plates.keys()) != list(range(1, self.number_source_plates + 1)):
                self.errors.append("There is a problem with the number of colonies per source plate. Remember that it is a relative position, it needs to go from 1 to the number of source plates")
            if len(self.index_start_source_plate.keys()) != self.number_source_plates or sorted(self.index_start_source_plate.keys()) != list(range(1, self.number_source_plates + 1)):
                self.errors.append("There is a problem with the start index per source plate. Remember that it is a relative position, it needs to go from 1 to the number of source plates")
            if len(self.positions_control.keys()) != self.number_source_plates or sorted(self.positions_control.keys()) != list(range(1, self.number_source_plates + 1)):
                self.errors.append("There is a problem with the positions of controls. Remember that it is a relative position, it needs to go from 1 to the number of source plates")
            if len(self.positions_not_perform_pcr.keys()) != self.number_source_plates or sorted(self.positions_not_perform_pcr.keys()) != list(range(1, self.number_source_plates + 1)):
                self.errors.append("There is a problem with the positions where the PCR should not be perform. Remember that it is a relative position, it needs to go from 1 to the number of source plates")
            
            # Control wells + Colonies pcr + Colonies not perform PCR = All colonies
            number_positions_not_pcr = 0
            number_all_colonies_plates = 0
            
            for list_notpcr_plate in self.positions_not_perform_pcr.values():
                if "None" in list_notpcr_plate:
                    pass
                else:
                    number_positions_not_pcr += len(list_notpcr_plate)
                    
            for list_samples_plate in self.number_samples_source_plates.values():
                number_all_colonies_plates += int(list_samples_plate[0])
                
            if self.number_colonies + self.number_controls + number_positions_not_pcr != number_all_colonies_plates:
                self.errors.append("The total number of samples is not equal to the addition of the controls, PCR samples and samples in which we should not perform the PCR")
            
            # We are going to check that the number of cells in each plate is not larger than the capacity of the source plates
            for number_cells_per_plate in self.number_samples_source_plates.values():
                if int(len(labware_context.get_labware_definition(self.name_source_plate)["wells"])) < int(number_cells_per_plate[0]):
                    self.errors.append("Number of cells is larger than the capacity of the source plate labware")
                else:
                    pass      
            
            # Index of any of the source plates is larger than the wells of that plate
            if self.index_start_final_plate > len(labware_context.get_labware_definition(self.name_final_plate)["wells"]):
                self.errors.append("Index of start for the final plate is larger than the number of wells of that labware")
            number_wells_source_plate = len(labware_context.get_labware_definition(self.name_source_plate)["wells"])
            for index_plate in self.index_start_source_plate.values():
                if int(index_plate[0]) > number_wells_source_plate:
                    self.errors.append("Index of start for one or more source plate(s) is larger than the number of wells of that labware")
                    break
        
            # We are goin to check that all the wells that are mentioned exist in the correspondent plates (i.e., we are looking for typos like AA1 or A111)
            all_values_wells_source_plate = [x for x in sum(list(self.positions_not_perform_pcr.values())+list(self.positions_control.values()),[]) if x != "None"]
            name_wells_source_plate = list(labware_context.get_labware_definition(self.name_source_plate)["wells"].keys())
            if all(item in name_wells_source_plate for item in all_values_wells_source_plate):
                pass
            else:
                self.errors.append("One or more wells from the 'not perform pcr' or 'controls' wells do not exist in the source labware")
            
            # Check if there is some typo in the the starting tip
            if self.name_right_pip not in ["None", "-"]:
                if not any(self.starting_tip_right_pip in columns for columns in labware_context.get_labware_definition(define_tiprack(self.name_right_pip))["ordering"]):
                    self.errors.append("Starting tip of right pipette is not valid, check for typos")
            if self.name_left_pip not in ["None", "-"]:
                if not any(self.starting_tip_left_pip in columns for columns in labware_context.get_labware_definition(define_tiprack(self.name_left_pip))["ordering"]):
                    self.errors.append("Starting tip of left pipette is not valid, check for typos")
    
    
    
            # Check if the RPM of the shaker are correct
            if self.shaker:
                if self.revolutions != "-" and (self.revolutions > 3000 or self.revolutions < 200):
                    self.errors.append("The revolutions set for the heater-shaker can not be produced by the module, the range is 200-3000rpm")
    
    
            # We are going to check if we need more than 1 final plate only if thermocycler is there
            if self.thermocycler:
                if self.index_start_final_plate+(self.number_controls+self.number_colonies)*self.sets > len(labware_context.get_labware_definition(self.name_final_plate)["wells"]):
                    self.errors.append("You need more than 1 final labware to run this protocol. There is only 1 thermocycler and we can use only 1 final labware.")
            
                # Check if the profile csv exists
                if os.path.isfile(self.thermocycler_profile_file) == False:
                    self.errors.append(f"The program file '{self.thermocycler_profile_file}' does not exist or it is not found!")
                    return
            
                # Lets check if the temperatures are inside the range of the thermocycler
                # First hold block temperature and lid temperature
                self.lid_temperature = max(pd.read_csv(self.thermocycler_profile_file)["Temperature (C)"].values)+10
                if self.lid_temperature > 110 or self.lid_temperature < 37:
                    self.errors.append("Lid temperature cannot be set with the thermocycler, the operative range of the thermocycler is 37-110C")
                    
                if self.final_block_temperature != None and (self.final_block_temperature > 99 or self.final_block_temperature < 4):
                    self.errors.append("Block temperature cannot be set with the thermocycler, the operative range of the thermocycler is 4-99C")
                
                # Now we check if the temperatures in the steps of the profile are in the range of temperature
                for row in pd.read_csv(self.thermocycler_profile_file).iterrows():
                    if row[1]["Temperature (C)"] > 110 or row[1]["Temperature (C)"] < 4:
                        self.errors.append("One step of the profile cannot be set with the thermocycler, the operative range of the thermocycler is 4-99C")
                    if row[1]["Cycle Status"].lower() not in ["start","end","-"]:
                        self.errors.append("One step of the profile has another value for 'Cycle Status' that is neither 'Start', 'End' nor '-'")
                
                if self.name_final_plate != "nest_96_wellplate_100ul_pcr_full_skirt":
                    self.warnings.append("The only labware oficially supported by OT-2 is 'nest_96_wellplate_100ul_pcr_full_skirt'.\n\t  Make sure that the thermocycler can work with yours.")
                    
                if self.pause_tc:
                    self.warnings.append("A pause will be effectuated before the thermocycler closes its lid.\n\t  Stay aware of the robot so you can resume the program.\n\t  If this step is not desired, set the variable 'Pause before Temperature Program' to False")
            
            return
    
    class calculated_volumes():
        """
        This class will store the calculated volumes that are needed for the mixes and other commands in the rest of the protocol
        
        All the volumes are calculated from the variables of the setted_parameters class
        """
        
        def __init__(self, variables, reactions_water, reactions_polymerase, reactions_primer, reactions_mix):
            # We are going to do the calculations with the extra pipetting factor
            self.water_per_reaction = variables.volume_water*variables.extra_pip_factor
            self.primer_per_reaction = variables.volume_primer*variables.extra_pip_factor
            self.polymerase_per_reaction = variables.volume_polymerase*variables.extra_pip_factor
            self.mix_per_reaction = variables.volume_reactives*variables.extra_pip_factor
            
            # To each tube there are going to be attached an amount of reactions so the reactives in that tube are the osurce for that ammount of reactions 
            self.reactions_per_tube_water = reactions_water
            self.reactions_per_tube_polymerase = reactions_polymerase
            self.reactions_per_tube_primer = reactions_primer
            self.reactions_per_tube_mix = reactions_mix
            
            # Volumes of each tube of reactives
            self.volume_tube_water = list(map(lambda item:math.ceil(item*self.water_per_reaction), self.reactions_per_tube_water))
            self.volume_tube_primer = list(map(lambda item:math.ceil(item*self.primer_per_reaction), self.reactions_per_tube_primer))
            self.volume_tube_polymerase = list(map(lambda item:math.ceil(item*self.polymerase_per_reaction), self.reactions_per_tube_polymerase))
            self.volume_tube_mix = list(map(lambda item:math.ceil(item*self.mix_per_reaction), self.reactions_per_tube_mix))
            
            return

    class positions():
        """
        Class to store the positions of reactives and labware
        
        Main function of this class is for transfering commands and for printing the information to the user
        """
        
        def __init__(self, positions_labwares_reactives):
            self.water = positions_labwares_reactives[0]
            self.polymerase = positions_labwares_reactives[1]
            self.primers = {key_primer:val_primer for dict_primer in positions_labwares_reactives[2] for key_primer,val_primer in dict_primer.items()}
            self.mixes = {key_mix:val_mix for dict_mixes in positions_labwares_reactives[3] for key_mix, val_mix in dict_mixes.items()}
            self.source_plates = positions_labwares_reactives[4]
            self.final_plates = positions_labwares_reactives[5]
            self.coldblocks = positions_labwares_reactives[6]
            
            return
    
# Functions
# ----------------------------------
# ----------------------------------

    def run_program_PCR(tc_mod, program, final_lid_state, final_block_state, volume_sample):
        """
        Function that will read the csv file with the steps of the program and will perform it
        in the thermocycler
        the program is already the pd.DataFrame, because it has already been read and establish that the lid temperature is okay,
        it exists, etc, i.e., control error of the file
        """
        
        # Set the lid temperature --> this is going to be estetic, but it could be dinamical if wanted
        lid_temperature = variables.lid_temperature
        
        # initialyze the state of the variable cycle that we will use to control if the step is a cycle or a step
        cycle = False
        
        # set the initail temperature of the lid
        tc_mod.set_lid_temperature(lid_temperature)
        for row in program.iterrows():
            # Check if it is a cycle or not, if it is a start of the end of it
            # This will work because we have already donde contorl of the values of this column which only can be -, Start or End
            if row[1]["Cycle Status"].lower() == "start":
                profile_termo =[{"temperature":float(row[1]["Temperature (C)"]),"hold_time_seconds":float(row[1]["Time (s)"])}]
                cycle = True
                continue
            elif row[1]["Cycle Status"].lower() == "end":
                profile_termo.append({"temperature":float(row[1]["Temperature (C)"]),"hold_time_seconds":float(row[1]["Time (s)"])})
                tc_mod.execute_profile(steps = profile_termo,
                                       repetitions = int(row[1]["Number of cycles"]),
                                       block_max_volume = volume_sample)
                cycle = False
                continue
            
            # Now we know if we have to add a step to the cycle or do the step directly
            if cycle == True:
                profile_termo.append({"temperature":float(row[1]["Temperature (C)"]),"hold_time_seconds":float(row[1]["Time (s)"])})
            elif cycle == False:
                tc_mod.set_block_temperature(row[1]["Temperature (C)"],
                                             hold_time_seconds = float(row[1]["Time (s)"]),
                                             block_max_volume = volume_sample)
            
        # Now we are going to put the block at one temeprature/open lid if it is establish like that
        tc_mod.deactivate_lid()
        
        if final_lid_state:
            tc_mod.open_lid()
            
        if final_block_state != None:
            tc_mod.set_block_temperature(final_block_state,
                                         block_max_volume = volume_sample)
        else:
            tc_mod.deactivate_block()
    
        return


    def define_tiprack (pipette):
        """
        Define the correspondent tiprack to the pipette
        This function should be updated as more pipettes are included in the possible ones
        
        Also this function only has the not filtered tips associated
        """
        
        if pipette == "p20_single_gen2" or pipette == "p20_multi_gen2":
            return "opentrons_96_tiprack_20ul"
        elif pipette == "p300_single_gen2" or pipette == "p300_multi_gen2":
            return "opentrons_96_tiprack_300ul"
        elif pipette == "p1000_single_gen2":
            return "opentrons_96_tiprack_1000ul"
        else:
            raise Exception("The pippette that you are using is not established.")
        return


    def check_tip_and_pick(pipette_used, variables):
        """
        This functions is used to pick tips and in case of need, replace the tiprack or add one to the labware
        This way we add labwares in the process of the simulation and we do not need to calculate it a priori
        
        In the OT-App it will appear directly in the deck but it has been added with this function
        """
        # One future improvemnt of this function is to check if the pipettes use the same tipracks and add them to both pipettes, that way we will need less tipracks if they
        # have the same tips associated, for exmaple, if we are using a 300 multi and single pipette
        
        try:
            pipette_used.pick_up_tip()
            # When there are no tips left in the tiprack OT will raise an error
        except:
            # We only enter in this part of the function when there are not more tips in the respective tiprack
            if len(pipette_used.tip_racks) == 0: # There are no tipracks associated to this pipette yet (it is the beginning of the protocol)
                
                try:
                    first_position_free = [position for position, labware in variables.positions.items() if labware == None][0] # Check which is the first position free in the deck
                except:
                    raise Exception("There is not enough space in the deck for this protocol, try less samples")
                
                # Add to the pipette the tiprack that we are going to set in the first free position in the deck 
                pipette_used.tip_racks += [protocol.load_labware(define_tiprack(pipette_used.name), first_position_free)]
                
                # Update the positions
                variables.positions[first_position_free] = define_tiprack(pipette_used.name)
                
                # We establish now the starting tip, it will only be with the first addition, the rest will be establish that the first tip is in A1 directly
                if pipette_used.mount == "right":
                    pipette_used.starting_tip = pipette_used.tip_racks[0][variables.starting_tip_right_pip]
                elif pipette_used.mount == "left":
                    pipette_used.starting_tip = pipette_used.tip_racks[0][variables.starting_tip_left_pip]
                    
            else:
                if variables.replace_tiprack == False:
                    
                    try:
                        first_position_free = [position for position, labware in variables.positions.items() if labware == None][0]
                    except:
                        raise Exception("There is not enough space in the deck for this protocol, try less samples")
                        
                    pipette_used.tip_racks += [protocol.load_labware(define_tiprack(pipette_used.name), first_position_free)]
                    
                    # Update the positions
                    variables.positions[first_position_free] = define_tiprack(pipette_used.name)
                    
                else:
                    #Careful with this part if you are traspassing this script into jupyter because this will crash your jupyter (will wait until resume and it does not exist)
                    variables.need_change_tiprack = True
                    protocol.pause("Replace Empty Tiprack With A Full One And Press Resume In OT-APP")
                    pipette_used.reset_tipracks()
            
            #Finally, we pick up the needed tip        
            pipette_used.pick_up_tip()
        return
    
    def number_tubes_needed (vol_reactive_per_reaction_factor, number_reactions, vol_max_tube):
        """
        Given a maximum volume of the tube (vol_max_tube), the volume of that reactive/reaction (vol_reactive_per_reaction_factor) and the total number of reactions (number_reactions)
        this function will return the number of tubes needed for this reactive and how many reactions are filled in every tube
        
        This function does not garantee the lower number of tubes but it assures that everything can be picked with the pipettes that we have
        
        This function will be used mainly, sometimes exclusively, in setting_number_plates
        """
        
        number_tubes = 1 # Initizialice the number of tubes
        reactions_per_tube = [number_reactions] # Initialice the reactions per tube
        volumes_tubes = [vol_reactive_per_reaction_factor*number_reactions]*number_tubes # Initialice the number of tubes
        
        while any(volume > vol_max_tube for volume in volumes_tubes): # If there is some volume that is greater than the max volume we are goin to enter in the loop
            number_tubes += 1 # We add on emore tube so th evolume can fit in the tubes
            
            # Now we redistribute the reactions (and correspondant volume) to the tubes so it will be the most homogeneus way
            reactions_per_tube = [int(number_reactions/number_tubes)]*number_tubes
            tubes_to_add_reaction = number_reactions%number_tubes
            for i in range(tubes_to_add_reaction):
                reactions_per_tube[i] += 1
            
            # Calculate the new volumes
            volumes_tubes = [vol_reactive_per_reaction_factor*number_reactions_tube for number_reactions_tube in reactions_per_tube]
        
        # When the volume can fit every tube (exit from th ewhile loop) we return the number of tubes and the reactions that will fit in every tube
        return (number_tubes, reactions_per_tube)
    
    def setting_number_plates (number_plates, labware_name, coldblock = False):
        """
        In this function we will set how many labwares we need of every category (source labwares, final, coldblocks, falcon tube racks, etc)

		This function will only set the labwares in the different slots of the deck, with not calculate how many we need,
		this way we do not have to change this function and only change the setting_labware function from protocol to protocol
		"""

        position_plates = [position for position, labware in variables.positions.items() if labware == None] # We obtain the positions in which there are not labwares
        if coldblock:
            if variables.shaker and variables.thermocycler:
                position_not = [2, 4, 5]
                positions_plates_coldblock = [position for position in position_plates if position not in position_not]
            elif variables.shaker:
                position_not = [2, 4]
                positions_plates_coldblock = [position for position in position_plates if position not in position_not]
            elif variables.thermocycler:
                position_not = [4, 5]
                positions_plates_coldblock = [position for position in position_plates if position not in position_not]
                
            # Now we have to put some warning about putting the coldblocks next to the modules
            if (len(positions_plates_coldblock) < number_plates) and (len(position_plates) >= number_plates):
                variables.warnings.append("Due to the lack of space we need to put one or more of the coldblocks next to a module, make sure that it fits!\n\tIf it does not fit, try lowering the number of samples.")
                for needed_plates_from_not_optimal_positions in range(number_plates-len(positions_plates_coldblock)):
                    for position_deck_not_optimal in position_not:
                        if variables.positions[position_deck_not_optimal] == None:
                            positions_plates_coldblock.append(position_deck_not_optimal)
                        else:
                            next
                    
                # order the coldblock positions only to the map to be more intuitive
                positions_plates_coldblock.sort()
                position_plates = positions_plates_coldblock
            else:
                position_plates = positions_plates_coldblock
        all_plates = []
        
        if len(position_plates) < number_plates: # There were not enough items in the position_plate to fit the number of labwares that we needed
            raise Exception("There is not enough space in the deck for this protocol, try less samples")
        
        # There are enough deck space for the plates
        for index_number_need_plate in range(number_plates):
            for available_position_deck in position_plates:
                try:
                    plate = protocol.load_labware(labware_name, available_position_deck)
                except opentrons.protocols.geometry.deck_conflict.DeckConflictError: # There would be dimensions problems
                    continue
                except:
                    raise("There has been some problem establishing the coldblocks in the deck")
                
                all_plates.append(plate)
                
                # update the variables.positions
                variables.positions[available_position_deck] = labware_name
                position_plates.remove(available_position_deck)
                
                if len(position_plates) == 0 and index_number_need_plate != number_plates-1:
                    raise("Due to dimensions restrictions between the coldblocks and other labwares/modules the protolc cannot be done.\nTry less examples or taking out a module.")
                else:
                    break
        
        return all_plates

    
    def setting_reactives_in_places(number_tubes_reactives, name_reactive, possible_positions):
        """
        Function that will set tubes of reactives within a type of labware, in the case of PCR sample preparation is setting eppendorfs
        within one or more coldblocks
        
        Because we have already established how many labwares we needed, we can go through this function without having a possible error
        
        This function will be mainly or exclusivelly used in the setting_labware function
        
        We are going to obtain a dictionary with the name of the reactive and its positions
        """
        
        dict_position_reactive = {}
        for tube in range(int(number_tubes_reactives)):
            try:
                # It is not the first time the reactive is being assigned to a well so we are adding it to the list and not assigning it
                dict_position_reactive[name_reactive] = dict_position_reactive[name_reactive] + [possible_positions[tube]]
            except:
                dict_position_reactive[name_reactive] = [possible_positions[tube]]
        
        return (dict_position_reactive)
    
    def setting_labware(variables):
        """
        In this function we will calculate how many labwares of each type do we need to perform the protocol, and then we are going to call
        the function setting_number_plates that will set the labware in the respective positions
        
        In this function we are going to calculate and set the labwares in the deck according to the csv provided, with the exception of the tipracks
        
        This function should be changed depending of how many types of labware we need to set and how many different reactives do we need
        """
        # Setting thermocycler
        if variables.thermocycler:
            variables.tc_module = protocol.load_module("thermocycler")
            # Set the labware in the thermocycler
            labware_pcr = [variables.tc_module.load_labware(variables.name_final_plate)]
            for position_thermocycler in [7,8,10,11]:
                variables.positions[position_thermocycler] = variables.tc_module
            
        # Setting heater-shaker
        if variables.shaker:
            variables.hs_mod = protocol.load_module('heaterShakerModuleV1',1)
            variables.hs_mod.close_labware_latch()
            variables.labware_shaking = variables.hs_mod.load_labware(variables.name_labware_shaking)
            variables.positions[1] = variables.hs_mod
        
        # We start settign the source labware which number has been provided
        labware_source = setting_number_plates(variables.number_source_plates, variables.name_source_plate)
        
        # Set the final plates (PCR plates)
        # To do that we need to know how many wells the labware have and how many we will need
        if not variables.thermocycler:
            number_wells_labware_final = len(labware_context.get_labware_definition(variables.name_final_plate)["wells"])
            number_plates_pcr = math.ceil(((variables.number_reactions*variables.sets)+variables.index_start_final_plate)/number_wells_labware_final)
            labware_pcr = setting_number_plates(number_plates_pcr, variables.name_final_plate)
        
        # Set the coldblocks
        # For that we need to know the maximum volume of the tubes and how many tubes of the reactives we need in total
    
        max_volume_well = labware_context.get_labware_definition(variables.name_coldblock)["wells"]["A1"]["totalLiquidVolume"]
        
        # We are going to extract 5ul from the water, poly and primer max volumes so we make sure that we dont run out of reactives (that way we can always add 5 ul minimum)
        # Set number tubes needed for water
        tubes_water, reactions_per_tube_water = number_tubes_needed(variables.volume_water*variables.extra_pip_factor, variables.number_reactions*variables.sets, max_volume_well-5)
        # Set number tubes of the polymerase
        tubes_phyre, reactions_per_tube_phyre = number_tubes_needed(variables.volume_polymerase*variables.extra_pip_factor, variables.number_reactions*variables.sets, max_volume_well-5)
        # Set number of tubes for primers
        tubes_primer_each_primer, reactions_per_tube_primer = number_tubes_needed(variables.volume_primer*variables.extra_pip_factor, variables.number_reactions, max_volume_well-5)
        tubes_primer_total = tubes_primer_each_primer*variables.sets*variables.number_primers
        # Set number tubes for mixes
        tubes_mixes_each_set, reactions_per_tube_mix = number_tubes_needed(variables.volume_reactives*variables.extra_pip_factor, variables.number_reactions, max_volume_well)
        tubes_mixes_total = tubes_mixes_each_set*variables.sets
        
        # Set how many colblocks we need the tubes in the coldblock after we know how many tubes we need for everything
        number_wells_coldblock = len(labware_context.get_labware_definition(variables.name_coldblock)["wells"])
        # coldblocks_needed = math.ceil((tubes_water+tubes_phyre+tubes_primer_total+tubes_mixes_total)/number_wells_coldblock)
        coldblocks_needed = math.ceil((tubes_water+tubes_phyre+tubes_primer_total)/number_wells_coldblock)
        # labware_coldblock = setting_number_plates(coldblocks_needed, variables.name_coldblock) Hay que arreglar esto porque no se puede poner al lado del shaker
        # Let's find in which positions the coldblock can be put
        
        labware_coldblock = setting_number_plates(coldblocks_needed, variables.name_coldblock, coldblock = True)
               
        
        
        if tubes_mixes_total > len(labware_context.get_labware_definition(variables.name_labware_shaking)["wells"]):
            raise Exception(f"Tubes of mixes doesnt fit in shaker. It is needed {tubes_mixes_total} tubes of mixes and labware only fits {len(labware_context.get_labware_definition(variables.name_labware_shaking)['wells'])}")
        
        # Now we are going to set the reactives in the coldblock positions, we need to keep track of these positions for liquid movement        
        # Get the possible positions 
        positions_coldblock = []
        for labware in labware_coldblock:
            positions_coldblock += labware.wells()
        
        # Set the positions of the reactives and delete from the possible possitions the ones occupied by reactives
        # We set the water, polymerase, primers and mixes tubes respectivelly
        positions_water = positions_coldblock[:tubes_water]
        del positions_coldblock[:tubes_water]
        
        positions_phyre = positions_coldblock[:tubes_phyre]
        del positions_coldblock[:tubes_phyre]
        
        position_primers = []
        position_mixes = []
        for primer in range(1, (variables.sets*variables.number_primers)+1):
            name_reactive = "primer_"+str(primer)
            position_primers.append(setting_reactives_in_places(tubes_primer_each_primer, name_reactive, positions_coldblock))
            positions_coldblock = positions_coldblock[int(tubes_primer_each_primer):]
        
        positions_labware_shaker = variables.labware_shaking.wells()
        for mix in range(1,variables.sets+1):
            name_reactive = "setPrimers_"+str(mix)
            position_mixes.append(setting_reactives_in_places(tubes_mixes_each_set, name_reactive, positions_labware_shaker))
            positions_labware_shaker = positions_labware_shaker[tubes_mixes_each_set:]
            
        # Now we return the positions of the reactives and the labware we have setted with this function
        return (positions_water, positions_phyre, position_primers, position_mixes, labware_source, labware_pcr, labware_coldblock, reactions_per_tube_water, reactions_per_tube_phyre, reactions_per_tube_primer, reactions_per_tube_mix)
    
    def z_positions_mix(vol_mixing):
        """
        Function that will define the positions of mixing according to the volume of each tube of primer set
        
        These heights have been manually measured for 1.5mL eppendorfs to attach z to aproximatelly the volume associated
        
        We will have 3 mixing heights at the end, but not neccessarilly different within each other
        """
        
        position_bottom = 1
        position_100 = 6
        position_100_250 = 9
        position_250 = 11
        position_500 = 16
        position_750 = 20
        position_1000 = 25
        position_1250 = 30
        
        #Assigned to the volume the 3 positions [min, center, max] that we are going to use in the mixing process
        if vol_mixing <= 100: # The values of comparing are volumes (in uL)
            return [position_bottom, position_bottom, position_bottom]
        elif vol_mixing > 100 and vol_mixing <= 250:
            return [position_bottom, position_100, position_100_250]
        elif vol_mixing > 250 and vol_mixing <= 500:
            return [position_bottom, position_100, position_250]
        elif vol_mixing > 500 and vol_mixing <= 750:
            return [position_100, position_250, position_500]
        elif vol_mixing > 750 and vol_mixing <= 1000:
            return [position_100, position_250, position_750]
        elif vol_mixing > 1000 and vol_mixing <= 1250:
            return [position_100, position_500, position_1000]
        elif vol_mixing > 1250:
            return [position_100, position_500, position_1250]
        else:
            pass
        
        return    
        
    def distribute_vol_tracking(volume_distribute, positions_tube_source, reactions_tubes_source, positions_final, current_pipette_used, variables):
        """
        This function is destined to transfer small amounts (in comparasion to transfer_vol_tracking) from a source to different wells
        
        Specifically, here we used it to distribute the set mix to the wells and try to do it in the lower ammount of movements
        
        This function works because we have assigned previously the reactions to the tube, if we separate the reactives another way we
        should use other function
        
        This function of distribution is different from the other protocols because in this case the height is not an issue and the vol tracking is done in the volume calculating section
        """
        
        # We are not going to pick a tip because the selection of the tip already happened in the script
        
        # Lets distribute the reactive
        for index_tube, reactions_tube in enumerate(reactions_tubes_source):
            current_step = "Mixing tube(s) of primer set"
            variables.hs_mod.set_and_wait_for_shake_speed(variables.revolutions)
            protocol.delay(seconds=10)
            variables.hs_mod.deactivate_shaker()
            
            current_step = "Distributing primer set to final plates"
            current_pipette_used.distribute(volume_distribute, positions_tube_source[index_tube], positions_final[:reactions_tube], disposal_volume = 0, new_tip = "never")
            del positions_final[:reactions_tube]
        
        # The drop of the tip will be done outside of this function as well
        
        return    
        
    def set_labware(labware_wells_name):
        """
        Generator of the positions to transfer, it will give you the next element of a list every time it is called
        """
        for well in labware_wells_name:
            yield well

    def give_me_optimal_pipette_to_use(aVolume, pipette_r, pipette_l):
        """
        Function that will return the optimal pipette to use for the volume that we want to handle.
        
        In case that it is a great volume (higher than the maximal volume of both pipettes) will return the pipette that will give the minimal quantity of movements
        
        If none of the pipettes attached can pick the volume (because it is too small) the function will raise an error
        
        For the correct functioning of this function at least 1 pipette should be attached to the OT, otherwise, the function will raise an error
        """

        if pipette_r == None and pipette_l == None: #no pipettes attached
            raise Exception("There is not a pippette attached")
        
        # First we look if one of them is the only option
        elif pipette_r == None or pipette_l == None: # One mount is free, only need that the volume is more than the min of the pipette
            if pipette_r == None and aVolume >= pipette_l.min_volume:
                return pipette_l
            elif pipette_l == None and aVolume >= pipette_r.min_volume:
                return pipette_r
            else: # One of them does not exist and the other one is not valid
                raise Exception("Volume cannot be picked with pipette(s) associated. Try another pipette combination")
                
        else: # Both of them are in the OT
            # Define which has a bigger min volume so it can do the fewer moves to take the reactives or even distribute with fewer moves
            if pipette_l.min_volume > pipette_r.min_volume:
                max_pipette = pipette_l
                min_pipette = pipette_r
            else:
                max_pipette = pipette_r
                min_pipette = pipette_l
            
            if aVolume >= pipette_l.min_volume and aVolume >= pipette_r.min_volume:
                if pipette_l == max_pipette:
                    return pipette_l
                else:
                    return pipette_r
            elif aVolume >= pipette_l.min_volume and pipette_l == min_pipette:
                return pipette_l
            elif aVolume >= pipette_r.min_volume and pipette_r == min_pipette:
                return pipette_r
            else: # None of the pipettes can hold that volume
                raise Exception("Volume cannot be picked with pipette(s) associated. try another pipette combination")
                # Error control
        return

    def print_information_user(variables, output_settinglabware, volume_reactives):
        """
        This function will print the information that the user should know as the volumes of each tube and where to place it.
        
        This function should be customized to every protocol
        
        This will be printed when done the simulation because there is info needed to know before running the protocol in the OT
        """
        
        print("""--------------------------------------------------------------\n--------------------------------------------------------------\nINFORMATION FILE
    Data the users need to set the OT and perform the protocol
    It is the output of the python file for the protocol when executed with the opentrons_simulate option\n""")
        #Some general information that the user has setted
        print("--------------------------------------------------------------\nGENERAL INFORMATION\nThis details are set by the user, this is only a remainder of some variables")
        print("\t- Number of samples (colonies + controls): "+str(variables.number_colonies+variables.number_controls))
        print("\t- Number of primer sets: "+str(variables.sets))
        print("\t- Number of total reactions: "+str((variables.number_colonies+variables.number_controls)*variables.sets))
        print("\t- Final volume of output labware (uL): "+str(variables.volume_reactives+variables.volume_colonie))
        if variables.replace_tiprack == True and variables.need_change_tiprack == True:
            print("\t- You should stay near the OT to change the tiprack(s)\n")
        elif variables.replace_tiprack == True and variables.need_change_tiprack == False:
            print("\t- There will be no need to change tipracks\n")
        elif variables.replace_tiprack == False:
            print("\t- You have the replace tiprack option as False so there is no need of replacing none of the tipracks")
        
        if variables.map == True:
            print("\t- The map(s) is going to be  in "+variables.name_map+"_positionInDeck.csv")
        else:
            pass
        
        # Deck information of the labware
        print("\n--------------------------------------------------------------\nDECK LABWARE POSITIONS\nThis information will be provided by the OT App as well")
        for labware_position in variables.positions.values():
            print("\t"+str(labware_position))
        
        # Coldblock positions for the reactives
        coldblock_positions = output_settinglabware[6]
        water_positions = output_settinglabware[0]
        poly_positions = output_settinglabware[1]
        primer_positions = output_settinglabware[2]
        mixes_positions = output_settinglabware[3]
    
        # Printing of the reactives and their volume as it will bein the coldblock(s)
        print("\n--------------------------------------------------------------\nCOLDBLOCK REACTIVES POSITIONS\nReactives in each position of the coldblock(s) and their respective volume (minimum)")
        for coldblock in coldblock_positions:
            coldblock_table = pd.DataFrame(index=["A","B","C"],columns=["1","2","3","4"])
            print("\n--> Coldblock in Slot "+str(coldblock).split(" ")[-1])
            for index_water, water_tube in enumerate(water_positions):
                water_well = str(water_tube).split(" ")[0]
                if str(water_tube).split(" ")[-1] == str(coldblock).split(" ")[-1]:
                    coldblock_table[water_well[1:]][water_well[0]] = "Water - "+str(volume_reactives.volume_tube_water[index_water])+str("uL")
            for index_poly, poly_tube in enumerate(poly_positions):
                poly_well = str(poly_tube).split(" ")[0]
                if str(poly_tube).split(" ")[-1] == str(coldblock).split(" ")[-1]:
                    coldblock_table[poly_well[1:]][poly_well[0]] = "Polymerase - "+str(volume_reactives.volume_tube_polymerase[index_poly])+str("uL")
            for index_primer, primer_tubes in enumerate(primer_positions):
                for index_position, position in enumerate(list(primer_tubes.values())[0]):
                    primer_well = str(position).split(" ")[0]
                    if str(position).split(" ")[-1] == str(coldblock).split(" ")[-1]:
                        coldblock_table[primer_well[1:]][primer_well[0]] = "Primer "+str(index_primer+1)+" - "+str(volume_reactives.volume_tube_primer[index_position])+str("uL")
            # for index_mix, mix_tubes in enumerate(mixes_positions):
            #     for position in list(mix_tubes.values())[0]:
            #         mix_well = str(position).split(" ")[0]
            #         if str(position).split(" ")[-1] == str(coldblock).split(" ")[-1]:
            #             coldblock_table[mix_well[1:]][mix_well[0]] = "Set Reactives "+str(index_mix+1)
            print(coldblock_table.fillna("-")) # for esthetic purposes
            print("\n--------------------------------------------------------------")
    
        if variables.shaker == True:
            print("""\n--------------------------------------------------------------\nSHAKER REACTIVES POSITIONS\nReactives in each position of the heater-shaker
                  
--> Heater-Shaker in Slot 1""")
            labware_shaker_table = pd.DataFrame(index = ["A","B","C","D"], columns= ["1","2","3","4","5","6"])
            for index_mix, mix_tubes in enumerate(mixes_positions):
                for position in list(mix_tubes.values())[0]:
                    mix_well = str(position).split(" ")[0]
                    labware_shaker_table[mix_well[1:]][mix_well[0]] = "Set Reactives "+str(index_mix+1)
            print(labware_shaker_table.fillna("-")) # for esthetic purposes
            print("\n--------------------------------------------------------------")
        return
    
    ## The following functions are specific to the PCR sample preparation, they can be used in future protocols but for now they are only used in this one
  
    
    
    def pds_labware_sources(labware_source):
        """
        Function destined to create a DataFrame with the dimensions of the labware_source.
        
        Its main use is to creat emaps that will be filled during the running of the protocol to track where is every sample going
        
        It is destined to create a Dataframe in which the columns are numbers and the rows letters (only 1 letter) (as most of the 96 tisue plates)
        Future improvements will take in account than the rows could have more than one letter, for example, a well that will be AA5
        """
        
        # We obtain the names of the columns and the rows from the labware definition
        labware_wells = labware_context.get_labware_definition(labware_source[0].name)["groups"][0]["wells"]
        name_rows = sorted(list(set([well[0] for well in labware_wells])))
        number_rows = len(name_rows)
        number_columns = len(set([well[1:] for well in labware_wells]))
        
        # Initialize the dictiony of dataframes
        dataframes_source_with_index = {}
        
        # Create the dataframes for every labware
        for labware in labware_source:
            dataframes_source_with_index[int(str(labware).split(" ")[-1])] = pd.DataFrame(np.full((number_rows,number_columns),None),columns=list(range(1,number_columns+1)),index=name_rows)
            dataframes_source_with_index[int(str(labware).split(" ")[-1])].index.name = "Row/Column"
    
        return dataframes_source_with_index

    def map_tracking(pos_source, pos_final, dataframes_with_index, map_source_plate, primer_set):
        """
        Function that will track which sample from the source plate (pos_source in map_source_plate) is being dispense in the pos_final
        
        This way we will obtain a map of the identifiers of the plate + their position in the deck and we can backtrack the samples in case that we need it
        """
        
        # Process the pos_source and pos_final to get the well and the labware position in the deck to fill correctly the dataframe from dataframes_with_index 
        pos_source_deck = str(pos_source).split(" ")[-1]
        pos_source_labware = str(pos_source).split(" ")[0]
        
        column_pos_cource = str(pos_source_labware[1:])
        row_pos_source = str(pos_source_labware[0])
        
        pos_final_deck = str(pos_final).split(" ")[-1]
        pos_final_labware = str(pos_final).split(" ")[0]
    
        identifier_sample = map_source_plate[int(pos_source_deck)][column_pos_cource][row_pos_source]
        dataframes_with_index[int(pos_final_deck)][int(pos_final_labware[1:])][pos_final_labware[0]] = str(identifier_sample+" on "+pos_source_deck+" with "+str(primer_set+1))
        # We are going to still put the position on the deck because 2 samples can have the same identifier in different plates
        
        return
    
    def change_dict_key(d, old_key, new_key):
        """
        We are going to use this function to preserve the map_tracking of other scripts but in this case we dont have the source plates in the first positions
        Only we are going to use it if we have the shaker, because it is in the first 
        """
        d[new_key] = d.pop(old_key, d[old_key])
        return
    
    
    
    
    def delete_certain_wells_from_list(all_wells, delete_wells, plates_in_action):
        """
        This function is going to remove from a list of elements (wells) specific elements (delete_wells) from specific deck positions (plates in action)
        
        In the PCR sample preparation is going to be used to delete certain wells from the PCR final plates, that way if there are several samples that have been discarded
        we do not have to perform following PCRs with them
        """
        
        for index_plate, plate in enumerate(plates_in_action):
            try:
                # Take the values of the dictionary that should be erased
                delete_wells_in_this_plate = delete_wells[index_plate+1]
            except: # There is no well to eliminate from this source plate
                continue # Go to the next plate
            
            for well in delete_wells_in_this_plate: # For every element of that dictionary value
                if well == str(None): # Double check
                    next
                elif plate[well] in all_wells:
                    all_wells.remove(plate[well])

        return all_wells
    
    def mixing_eppendorf_15(locations_of_mixing, volumes_tubes_mixing, variables, warning):
        """
        This is a function to perfrom an extensive mixing of every eppendorf which should be done before distributing the reactives along the final plates
        
        Mixing is one of the most crucial parts of this workflow and that is why theer is a function only for it
        
        This function will perform the mixing and will return warnings (in case that is needed) and the final pipette that has been used
        """
        # This function is going to perform the mixing of the tubes in the coldblock (it is setted for the 1500ul eppendorf because the positions are done manually)
        
        original_pipette = variables.right_pipette #this is only to initialize as one, it will change if neccessary
        for index_location, location_mixing in enumerate(locations_of_mixing):
            mastermix_mixing_volume_theory = volumes_tubes_mixing[index_location] / 3
            try:
                aSuitablePippet_mixing = give_me_optimal_pipette_to_use(mastermix_mixing_volume_theory, variables.right_pipette, variables.left_pipette)
                max_vol_pipette = aSuitablePippet_mixing.max_volume
                if max_vol_pipette < mastermix_mixing_volume_theory:
                    volume_mixing = max_vol_pipette
                    warning.append("Volume of tube in "+str(location_mixing)+" is "+str(volumes_tubes_mixing[index_location])+"uL and mixing volume should be "+str(mastermix_mixing_volume_theory)+"uL, this volume cannot be picked and final mixing volume is going to be "+str(volume_mixing)+"uL.")
                else:
                    volume_mixing = mastermix_mixing_volume_theory
                    pass
            except: # If this happens it means that the the volume is too low to any of the pipettes
                if variables.right_pipette.min_volume < variables.left_pipette.min_volume:
                    volume_mixing = variables.right_pipette.min_volume
                    aSuitablePippet_mixing = variables.right_pipette
                else:
                    volume_mixing = variables.left_pipette.min_volume
                    aSuitablePippet_mixing = variables.left_pipette
                warning.append("Volume of tube in "+str(location_mixing)+" is "+str(volumes_tubes_mixing[index_location])+"uL and mixing volume should be "+str(mastermix_mixing_volume_theory)+"uL, this volume cannot be picked and final mixing volume is going to be "+str(volume_mixing)+"uL.")
            
            if index_location == 0:
                check_tip_and_pick(aSuitablePippet_mixing, variables)
                original_pipette = aSuitablePippet_mixing
            else:
                if aSuitablePippet_mixing == original_pipette: # We do not change the tip because it is the same mix (same reactives)
                    pass
                else:
                    original_pipette.drop_tip()
                    check_tip_and_pick(aSuitablePippet_mixing, variables)
                    original_pipette = aSuitablePippet_mixing
            
            # After calculating the mixing volume, choosing a pipette and picking up a tip we perform the mix
            positions_mixing = z_positions_mix(volume_mixing) # This is the part that is customized for the 1500uL eppendorfs
            
            # We are going to mix 7 times at different heighs of the tube
            aSuitablePippet_mixing.mix(7, volume_mixing, location_mixing.bottom(z=positions_mixing[1])) 
            aSuitablePippet_mixing.mix(7, volume_mixing, location_mixing.bottom(z=positions_mixing[0])) 
            aSuitablePippet_mixing.mix(7, volume_mixing, location_mixing.bottom(z=positions_mixing[2])) 
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -20, radius=0.7, speed=30)
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -20, radius=0.7, speed=30)
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -20, radius=0.7, speed=30)
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -20, radius=0.5, speed=30)
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -20, radius=0.5, speed=30)
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -20, radius=0.5, speed=30)
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -27, radius=0.3, speed=30)
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -27, radius=0.3, speed=30)
            aSuitablePippet_mixing.touch_tip(location_mixing,v_offset = -27, radius=0.3, speed=30)
            # Now we are going to aspirate and dispense 3 times at different heights to mix a little bit more the content of the tube
            for i in range(2):
                aSuitablePippet_mixing.aspirate(volume_mixing, location_mixing.bottom(z=positions_mixing[0]))
                aSuitablePippet_mixing.dispense(volume_mixing, location_mixing.bottom(z=positions_mixing[2]))
            for i in range(2):
                aSuitablePippet_mixing.aspirate(volume_mixing, location_mixing.bottom(z=positions_mixing[2]))
                aSuitablePippet_mixing.dispense(volume_mixing, location_mixing.bottom(z=positions_mixing[0]))
            aSuitablePippet_mixing.blow_out(location_mixing.center())
        
        return warning, aSuitablePippet_mixing
    
    def transfer_volume_tracking(vol_distribute, position_source, position_final, volumes_source, variables):
        """
        Function used to transfer bigger volumes than the distributing with volume tracking but this does not take in account the heights and the distributing takes in account
        that all the final locations (wells) will have the same quantity (distributing along the plate, for example) while this is used mainly to create the mix that will be distributed
        
        In a future transfer and distributing function will be unified
        """
        # This function is used to transfer bigger volumes (in comparasion to distribute)
        # The distribute and transfering functions are very similar, sometimes they can be interchanged (this is not the case for versions < 11 of the script)
        
        # Initialize different variables
        current_position_source = position_source[0]
        current_volume = volumes_source[0]
        i_position = 0
        pipette_used = give_me_optimal_pipette_to_use(vol_distribute[0], variables.right_pipette, variables.left_pipette)
        
        # Just in case, better to waste 1 tip than to contaminate reactives
        if variables.right_pipette != None and variables.right_pipette.has_tip == True:
            variables.right_pipette.drop_tip() 
        if variables.left_pipette != None and variables.left_pipette.has_tip == True:
            variables.left_pipette.drop_tip()
        
        # Now we loop over the different sets
        # Every set of mixes is going to have the same amount of liquid, the only difference is the primers that are going to have
        for index_position, position_final_current in enumerate(position_final):
            vol_distribute_position = vol_distribute[index_position]
            while vol_distribute_position > 0:
                if current_volume >= vol_distribute_position: #We use it only one time, usually at the beginning or the end of the loop
                    provisional_pipette = give_me_optimal_pipette_to_use(vol_distribute_position, variables.right_pipette, variables.left_pipette)
                    if provisional_pipette == pipette_used:
                        if pipette_used.has_tip == False: # This can be possible in the first movement of transfering the volume, it is only a double check
                            check_tip_and_pick(pipette_used, variables)
                        else:
                            pass
                    elif provisional_pipette != pipette_used and pipette_used.has_tip == True:
                        pipette_used.drop_tip()
                        pipette_used = provisional_pipette
                        check_tip_and_pick(pipette_used, variables)
                    else: # We do not ned to drop the tip because it does not have one
                        pipette_used = provisional_pipette
                        check_tip_and_pick(pipette_used, variables)
                    
                    pipette_used.transfer(vol_distribute_position, current_position_source, position_final_current, new_tip="never")
                    current_volume -= vol_distribute_position
                    vol_distribute_position -= vol_distribute_position
                elif current_volume < 1: #this will be the round of math.ceil volume that it is used only for math problems due to all decimals from the float format in python
                    i_position += 1
                    current_position_source = position_source[i_position]
                    current_volume = volumes_source[i_position]
                else: # This will mean that the tube doe snot have enough volume for all wells and we will need to aspirate from another tube also
                    provisional_pipette = give_me_optimal_pipette_to_use(current_volume, variables.right_pipette, variables.left_pipette)
                    if provisional_pipette == pipette_used:
                        if pipette_used.has_tip == False: #just in case
                            check_tip_and_pick(pipette_used, variables)
                        else:
                            pass
                    else:
                        pipette_used.drop_tip()
                        pipette_used = provisional_pipette
                        check_tip_and_pick(pipette_used, variables)
                    
                    pipette_used.transfer(current_volume, current_position_source, position_final_current, new_tip="never")
                    vol_distribute_position -= current_volume
                    #change the position of the reactive to the next tube of reactives
                    i_position += 1
                    current_position_source = position_source[i_position]
                    current_volume = volumes_source[i_position]
        
        #Need to return the last pipette used because it could not be the same as the first one and usually after transfering we need to drop the tip
        return pipette_used

# Body of the script
#--------------------------------
#--------------------------------

# We use a try and except so we can give some more user-friendly infromation about the errors that we can receive
    try:
        # Loading of the csv parameters and using the first column (the name of the variables) as index
        current_step = "Reading csv and transforming them to parameters/variables"
        variables_df = pd.read_csv("/data/user_storage/Variables-PCRs-OT-Shaker.csv", index_col = 0)
        # variables_df = pd.read_csv("Variables-PCRs-OT-Shaker.csv", index_col = 0)
        
        # We are going to convert these parameters into arguments of the class variables and we ar egoing to process some of them so they can be usable (they are going to be dictionaries)
        variables = setted_parameters(variables_df.to_dict(orient = "index"))
        
        # Due to the fact that even in the parameters checking we use the labwares definition, we are going to check that they exist outside that function
        try:
            labware_context.get_labware_definition(variables.name_source_plate)
            labware_context.get_labware_definition(variables.name_final_plate)
            labware_context.get_labware_definition(variables.name_coldblock)
        except:
            raise Exception("One or more of the introduced labwares are not in the custom labware directory of the opentrons. Check for any typo of the API labware name.")
        
        
        # We need the number of wells to set the default value of number_sampleS_source_plates
        number_wells_source_plate = len(labware_context.get_labware_definition(variables.name_source_plate)["wells"])
        setted_parameters.proccess_variables(variables, "number_samples_source_plates", number_wells_source_plate, variables.number_source_plates)
        setted_parameters.proccess_variables(variables, "positions_not_perform_pcr", None, variables.number_source_plates)
        setted_parameters.proccess_variables(variables, "index_start_source_plate", 0, variables.number_source_plates)
        setted_parameters.proccess_variables(variables, "positions_control", None, variables.number_source_plates)
        
        # Now we are going to submit the variables to some previous checking to prevent errors in the protocol
        # Not every error is going to be prevented
        variables.check_setted_parameters()
        if len(variables.errors) > 0:
            print("------------------------------------------------\nERROR(S):")
            for error in variables.errors:
                print("\t- "+error)
            print("------------------------------------------------\nBefore proceeding this error(s) should be fixed")
            quit()
        
        #----------------------------------------------------------------------------------------------------------
        
        # Set the labwares into the deck and get the positions of the reactive tubes
        current_step = "Setting labwares into deck (not the tipracks)"
        
        output_setting_labware = setting_labware(variables)
        positions_reactives_labware = positions(output_setting_labware[:7])
        # volume_tubes_water, volumes_tubes_polymerase, volumes_tubes_primer(each primer), volumes_tubes_mix (per set) --> input of the positions class
        
        source_cell_plates_locations = positions_reactives_labware.source_plates
        pcr_final_plates_locations = positions_reactives_labware.final_plates
        
        volumes_reactives = calculated_volumes(variables, output_setting_labware[7], output_setting_labware[8], output_setting_labware[9], output_setting_labware[10])
        # variables, volume_tubes_water, volumes_tubes_polymerase, volumes_tubes_primer(each primer), volumes_tubes_mix (per set) --> input of calculated_volumes class
        #----------------------------------------------------------------------------------------------------------
        current_step = "Changing keys of maps to their positions in the deck"
        
        # print("++++++")
        # print(positions_reactives_labware.source_plates)
        # print("++++++")
        for index, position in enumerate(positions_reactives_labware.source_plates):    
            change_dict_key(variables.maps_source_plate, index, int(str(position).split(" ")[-1]))
    
        # Create the dataframes for the mapping
        current_step = "Creating maps of final labware"
        # First we create the maps for filling, we are creating it even if we dont want it
        maps_final_plates = pds_labware_sources(output_setting_labware[5])
        
        # Here we set pipettes and tips
        current_step = "Stating pipettes and associated tip racks"
        if variables.name_right_pip == "None" or variables.name_right_pip == "-":
            variables.right_pipette = None
        else:
            variables.right_pipette = protocol.load_instrument(variables.name_right_pip, "right", tip_racks=[])
            
        if variables.name_left_pip == "None" or variables.name_left_pip == "-":
            variables.left_pipette = None
        else:
            variables.left_pipette = protocol.load_instrument(variables.name_left_pip, "left", tip_racks=[])

        #----------------------------------------------------------------------------------------------------------
        
        # Now we are going to transfer the liquids to the tubes of mixes
        current_step = "Creating mixes (set of primers + water + polymerase)"
        
        number_tubes_mixes_per_set = len(volumes_reactives.reactions_per_tube_mix)
        
        #----------------------------------------------------------------------------------------------------------
        
        # Water
        current_step = "Transfering water"
        
        # Calculate the water that we have to transfer to every tube of mix
        volume_water_per_tube = [reactions_tube*volumes_reactives.water_per_reaction for reactions_tube in volumes_reactives.reactions_per_tube_mix]
        
        # Transfer with volume tracking
        last_pipette_used = transfer_volume_tracking(volume_water_per_tube*variables.sets, positions_reactives_labware.water, sum(positions_reactives_labware.mixes.values(),[]), volumes_reactives.volume_tube_water, variables)
        last_pipette_used.drop_tip()
        
        #----------------------------------------------------------------------------------------------------------
        
        # Distribute the primers
        current_step = "Transfering primers"
        
        # Calculate the primer amount that we have to transfer to each tube of mix
        volume_distribute_primer = [reactions_tube*volumes_reactives.primer_per_reaction for reactions_tube in volumes_reactives.reactions_per_tube_mix]
        
        # Transfer the primers to the correspondant tubes
        for number_subset, mix_current_location in enumerate(positions_reactives_labware.mixes.values()):
            for index_primer in range(number_subset*variables.number_primers,(number_subset*variables.number_primers)+variables.number_primers):
            #Create a list of the indexes of primers that go to each mix
                primer = list(positions_reactives_labware.primers.values())[index_primer]
                last_pipette_used = transfer_volume_tracking(volume_distribute_primer, primer, mix_current_location, volumes_reactives.volume_tube_primer, variables)
                last_pipette_used.drop_tip()
        
        #----------------------------------------------------------------------------------------------------------
        
        # Phire
        current_step = "Transfering polymerase"
        
        # Calculate the polymerase amount that we have to transfer to each tube of mix
        volume_phyre_per_tube = [reactions_tube*volumes_reactives.polymerase_per_reaction for reactions_tube in volumes_reactives.reactions_per_tube_mix]
        
        # We change aspiration/dispensing speed to move liquid slower (for the Green Phyre or in general the polymerase)
        # We change it because usually the polymerases are more viscous and need a lower aspiration and dispense rate
        if variables.right_pipette != None:
            variables.right_pipette.flow_rate.aspirate = 10
            variables.right_pipette.flow_rate.dispense = 10
        if variables.left_pipette != None:
            variables.left_pipette.flow_rate.aspirate = 10
            variables.left_pipette.flow_rate.dispense = 10
        
        # TranSfer the polymerase to the tubes with volume tracking of the polymerase tubes
        last_pipette_used = transfer_volume_tracking(volume_phyre_per_tube*variables.sets, positions_reactives_labware.polymerase, sum(positions_reactives_labware.mixes.values(),[]), volumes_reactives.volume_tube_polymerase, variables)
        last_pipette_used.drop_tip()
        
        # We RESTORE aspiration/dispensing speed to defaults
        if variables.right_pipette != None:
            variables.right_pipette.flow_rate.aspirate = 92.86
            variables.right_pipette.flow_rate.dispense = 92.86
        if variables.left_pipette != None:
            variables.left_pipette.flow_rate.aspirate = 92.86
            variables.left_pipette.flow_rate.dispense = 92.86
        
        #----------------------------------------------------------------------------------------------------------
        if variables.thermocycler:
            variables.tc_module.open_lid()
        
        # Now we distribute each mix in the final plate(s)
        # Let's use a pippet that can transfer such small volume
        current_step = "Picking pipette for distributing reactives to final plate"
        aSuitablePippet = give_me_optimal_pipette_to_use(variables.volume_reactives, variables.right_pipette, variables.left_pipette)
        # Because it does not need to be that low we are going to change the dispense height, this is only for precaution, it is optional
        aSuitablePippet.well_bottom_clearance
        
        # We create a list of all the possible positions in which we can distribute our mixes
        positions_pcr_reactives = []
        for pcr_labware in pcr_final_plates_locations:
            positions_pcr_reactives += pcr_labware.wells()
        positions_pcr_reactives = positions_pcr_reactives[variables.index_start_final_plate:]
        
        # For every set we are going to mix all the tubes and then distribute them
        for index_set_primer, set_primer_to_distribute in enumerate(positions_reactives_labware.mixes.values()):
            initial_volume_tubes_mix = copy.deepcopy(volumes_reactives.volume_tube_mix)
            check_tip_and_pick(aSuitablePippet,variables)
            distribute_vol_tracking(variables.volume_reactives, set_primer_to_distribute, volumes_reactives.reactions_per_tube_mix, positions_pcr_reactives[variables.number_reactions*index_set_primer:variables.number_reactions*(index_set_primer+1)], aSuitablePippet, variables)
            aSuitablePippet.drop_tip()
        
        # We restore the height of dispense
        aSuitablePippet.well_bottom_clearance.dispense = 1
        
        #----------------------------------------------------------------------------------------------------------
        
        # And now we distribute the cells
        current_step = "Picking pipette for distributing cells"
        
        # Creation of the list of poistions to distribute the cells
        aSuitablePippet = give_me_optimal_pipette_to_use(variables.volume_colonie, variables.right_pipette, variables.left_pipette)
        position_source_cells = []
        
        #----------------------------------------------------------------------------------------------------------
        
        current_step = "Creating list of cell destinations"
        
        for index_labware, index_start_well in variables.index_start_source_plate.items():
            position_source_cells += source_cell_plates_locations[index_labware-1].wells()[int(index_start_well[0]):int(index_start_well[0])+int(variables.number_samples_source_plates[index_labware][0])]
        
        #----------------------------------------------------------------------------------------------------------
        
        # We are going to delete the positions that are not going to the thermocycler and put the control positions at the end so they can be at the end of each primer subset
        current_step = "Deleting from the picking wells the positions that should not be picked"
        
        # delete the not wanted positions and the controls of the positions to take
        positions_source_cells = delete_certain_wells_from_list(position_source_cells, variables.positions_not_perform_pcr, source_cell_plates_locations)
        positions_source_cells = delete_certain_wells_from_list(position_source_cells, variables.positions_control, source_cell_plates_locations)
        #----------------------------------------------------------------------------------------------------------
        
        current_step = "Stating the control positions (they are going to be distributed at the end of each primer set)"
        for plate, wells_plate in variables.positions_control.items():
            # Add the controls at the final of the list
            positions_source_cells += [source_cell_plates_locations[plate-1][well] for well in wells_plate if well != str(None)]
        
        # Define the generator of positions
        current_position_cell = set_labware(position_source_cells)
        current_position_pcr = set_labware(positions_pcr_reactives)
        #hay que quitar esta variable porque no se usa, si lo hacemos por indices
        #----------------------------------------------------------------------------------------------------------
        
        # And now lets distribute cells to the primer pairs
        current_step = "Distributing samples to final plates"
        
        # First we create the maps for filling, we are creating it even if we dont want it
        labwares_source_indexes = pds_labware_sources(positions_reactives_labware.source_plates)
        
        
        for colony_index in range(variables.number_reactions):
            source_cell = next(current_position_cell)
            for type_set_primer in range(variables.sets):
                check_tip_and_pick(aSuitablePippet, variables)
                aSuitablePippet.transfer(variables.volume_colonie, source_cell, positions_pcr_reactives[(colony_index+(type_set_primer*variables.number_reactions))], new_tip="never")
                # We put that transfer into the maps
                map_tracking(source_cell, positions_pcr_reactives[(colony_index+(type_set_primer*variables.number_reactions))], maps_final_plates, variables.maps_source_plate, type_set_primer)
                # map_tracking(source_cell, positions_pcr_reactives[(colony_index+(type_set_primer*variables.number_reactions))], labwares_source_indexes)
                aSuitablePippet.drop_tip()
        
        #----------------------------------------------------------------------------------------------------------
        
        # Now we create the map in case the user wants it
        current_step = "Importing the map of cells and primers for each final plate"
        if variables.map == True:
            for map_source in maps_final_plates:
                # print(map_source)
                maps_final_plates[map_source] = maps_final_plates[map_source].fillna(value="-")
                maps_final_plates[map_source].to_csv(str(variables.name_map)+"_"+str(map_source)+".csv")
        
        
        if variables.thermocycler:
	    if variables.pause_tc:
	    	protocol.pause("Protocol is pause so plate in thermocyler can be mix or user can put caps on it")
            variables.tc_module.close_lid()
            run_program_PCR(variables.tc_module, pd.read_csv(variables.thermocycler_profile_file), variables.final_open_lid, variables.final_block_temperature, variables.volume_reactives+variables.volume_colonie)
        
        
        
        
        #----------------------------------------------------------------------------------------------------------
        
        current_step = "Final homing of the robot arm"
        protocol.home()
        
        #----------------------------------------------------------------------------------------------------------
        
        current_step = "Printing of warnings and user information as position and volumes of reactives"
        if len(variables.warnings) > 0:
            print("--------------------------------------------------------------\nWARNING(S):")
            for warning in variables.warnings:
                print("\t- "+warning)
            print("--------------------------------------------------------------\nThis warning(s) will not stop the proccess but should be inspected\n")
        print_information_user(variables, output_setting_labware, volumes_reactives)
        
    except Exception as e:
        print("-------------------------------------------------------------------------------")
        print("--> ERROR STEP:\n" + current_step)
        logger.critical(e, exc_info=True)
        print("-------------------------------------------------------------------------------")
