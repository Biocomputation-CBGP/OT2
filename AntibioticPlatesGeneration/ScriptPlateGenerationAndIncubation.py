"""
Python script destined to OT-2
This script performs a creation of reactive plates mixed with the respective colonies of the source plate
This scirpt needs a csv attached to perform the running and will give an output file (txt) with some instructions
For more info go to https://github.com/Biocomputation-CBGP/OT2/tree/main/AntibioticPlatesGeneration and/or ### pagina del protocols.io del protocolo
"""

## Packages needed for the running of the protocol
import opentrons.execute
from opentrons import protocol_api
import pandas as pd
import logging
import math
import copy

## Structure needed for the running of th escript in the OT
metadata = {
'apiLevel':'2.12'
}

def run(protocol: protocol_api.ProtocolContext):
# Some neccesary statements
# ----------------------------------
# ----------------------------------

	## We need to assign to a variable some labware functions because we will need info of different labwares in the process of the script
	labware_context = opentrons.protocol_api.labware

	## We are going to create a logger to control errors that can occure in the script that are not taked in account and some of ours that will quit the program
	logging.basicConfig(format="""-------------------------------------------------------------------------------
--> SUM ERROR:
%(message)s
-------------------------------------------------------------------------------
--> INFO ERROR:""", level=logging.ERROR)
	logger = logging.getLogger()

# Classes definitions
# ----------------------------------
# ----------------------------------

	class setted_parameters():
		"""
		Class that will contain the parameters setted in the variables csv and will process them to work easily in the rest of the protocol
		The coding of this function is dependant of the variables in the Template of the protocol and the names have to be consistent with the rest of the code
		"""
		vol_max_tube = 15000 # The same always because the transfer is setted with the heights of the 15mL falcon
							 # One of the few labware restrictions of this protocol
		
		def __init__(self, variables_csv):
			"""
			This function will take the pandas dataframe that will be the table of the csv variables
			"""
			self.number_source_plates = int(variables_csv._get_value("Number of Source Plates","Value"))
			self.samples_per_plate = variables_csv._get_value("Samples Per Plate","Value") # This will need some processing
			self.number_samples = None # We need to process first samples_per_plate and then state the totla number of samples
			self.number_antibiotics = len(variables_csv._get_value("Name Antibiotics","Value").split(","))
			self.plates_final = 0 # Initial value
			self.volume_antibiotic = float(variables_csv._get_value("Volume of Antibiotics to Transfer (uL)","Value"))
			self.volume_colonies = float(variables_csv._get_value("Volume of Sample to Transfer (uL)","Value"))
			self.name_right_pipette = variables_csv._get_value("Name Right Pipette (Multichannel)","Value")
			self.name_left_pipette = variables_csv._get_value("Name Left Pipette (Singlechannel)","Value")
			self.right_pipette = None # First we are going to check the parameters (for errors) and then we are going to settle the pipettes
			self.left_pipette = None # First we are going to check the parameters (for errors) and then we are going to settle the pipettes
			self.starting_tip_right_pip = variables_csv._get_value("Initial Tip Right Pipette","Value")
			self.starting_tip_left_pip = variables_csv._get_value("Initial Tip Left Pipette","Value")
			self.name_source_plate = variables_csv._get_value("Name Source Plate","Value")
			self.name_final_plate = variables_csv._get_value("Name Final Plate","Value")
			self.name_rack_falcons = variables_csv._get_value("Name 15mL Tuberack","Value")
			if variables_csv._get_value("Replace Tipracks","Value").lower() == "true":
				self.replace_tiprack = True
			elif variables_csv._get_value("Replace Tipracks","Value").lower() == "false":
				self.replace_tiprack = False
			self.antibiotics = {} # Initialize
			for number_antibiotic in range(self.number_antibiotics):
				self.antibiotics[number_antibiotic] = {"name":variables_csv._get_value("Name Antibiotics","Value").strip().split(",")[number_antibiotic]}
			self.data_source_plates = {} # Will be filled in the course of the script
			self.data_final_plates = {} # Will be filled in the course of the script
			return
		
		def proccess_variables(self, name_parameter, value_default_empty, must_entries):
			"""
			This function will process some text that is extracted from the variables and make it to a dictionary so we can work with it
			Some of the variables need processing to be usable in the script, usually is from a "(,),(,)" format to a dictionary
			
			The function need the name of the parameter in the class, default value of the variable and how many items the dictiony needs to have
			In case that the structure of the parameters are different, this function should be changed
			"""
			
			value_parameter = getattr(self, name_parameter).replace(" ","")[1:-1].split('),(')
			
			# We separate the text value into the real values (they were separated by commas)
			for i, element in enumerate(value_parameter):
					value_parameter[i] = element.split(",")
			
			dict_final = {}
			if len(value_parameter[0][0]) != 0: # This checks if the the variable has a default value (like 0 or -) or if it has other values
												# This works with the process that we have done in  the previous part
				for well, plate in value_parameter:
						dict_final[int(plate)] = int(well)
			else:
				pass
			
			# We make sure that we have the number of entries that are specified in must_entries and if there is one missing we fill it with the default value for the item
			# If the dictionary is empty (when value in template is the default one) we will cretae all the items with the default value of the items (for example can be the number of wells of a plate or the index)
			for i in range(1, must_entries+1):
				try:
					dict_final[i]
				except:
					dict_final[i] = value_default_empty
			
			# Set the attribute after processing it
			setattr(self, name_parameter, dict_final)
			return

# Definition of the functions that are going to be used in the script
# ----------------------------------
# ----------------------------------

	def check_setted_parameters(variables):
		"""
		Function that will check the variables of the Template and will raise errors that will crash the OT run
		It is a validation function of the variables checking errors or inconsistencies
		
		This function is dependant again with the variabels that we have, some checks are interchangable between protocols, but some of them are specific of the variables
		"""
		# Initialize the errors list in which we are going to store the messages that we need to display for the user
		errors = []
		
		# We are going to check that the number of cells in each plate is not larger than the capacity of the source plates
		for number_cells_per_plate in variables.samples_per_plate.values():
			if int(len(labware_context.get_labware_definition(variables.name_source_plate)["wells"])) < int(number_cells_per_plate):
				errors.append("Number of cells is larger than the capacity of the source plate labware")
			else:
				pass      
		
		# Because we need both pipettes, a multichannel and a singlechannel we are going to check that both are provided
		# We are going to check if they have the appropiate channels after defining the pipettes
		if variables.name_left_pipette.lower() == "none" or variables.name_right_pipette.lower == "none":
			errors.append("We need 2 pipettes: a multichannel in the right and a singlechannel in the left")
		else:
			pass
		
		# We are going to check that the colonies + antibiotic is not more than the max volume of the wells in the final plates
		max_volume_well = float(list(labware_context.get_labware_definition(variables.name_final_plate)["wells"].values())[0]['totalLiquidVolume'])
		if variables.volume_antibiotic + variables.volume_colonies > max_volume_well:
			errors.append("The sum of sample and antibiotic volumes exceeds the max volume of final plate wells")
		else:
			pass
		
		# We are going to check, if for some reason, the dict of sampels per plate is major than the number of plates
		if len(variables.samples_per_plate) > variables.number_source_plates:
			errors.append("Check if the number of source plates and the samples per plate are correct")
		else:
			pass
		
		# We are going to check that the source and final plates have same dimensions (rows and columns)
		rows_source = len(labware_context.get_labware_definition(variables.name_source_plate)["ordering"][0])
		rows_final = len(labware_context.get_labware_definition(variables.name_final_plate)["ordering"][0])
		columns_source = len(labware_context.get_labware_definition(variables.name_source_plate)["ordering"])
		columns_final = len(labware_context.get_labware_definition(variables.name_final_plate)["ordering"])
		
		if rows_source != rows_final or columns_source != columns_final:
			errors.append("Source and final plates have not same dimensions (rows and columns)")
		else:
			pass
			
			return errors
	
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
					first_position_free = [position for position, labware in protocol.deck.items() if labware == None][0] # Check which is the first position free in the deck
				except:
					raise Exception("There is not enough space in the deck for this protocol, try less samples")
				
				# Add to the pipette the tiprack that we are going to set in the first free position in the deck 
				pipette_used.tip_racks += [protocol.load_labware(define_tiprack(pipette_used.name), first_position_free)]
				
				# We establish now the starting tip, it will only be with the first addition, the rest will be establish that the first tip is in A1 directly
				if pipette_used.mount == "right":
					pipette_used.starting_tip = pipette_used.tip_racks[0][variables.starting_tip_right_pip]
				elif pipette_used.mount == "left":
					pipette_used.starting_tip = pipette_used.tip_racks[0][variables.starting_tip_left_pip]
					
			else:
				if variables.replace_tiprack == False:
					
					try:
						first_position_free = [position for position, labware in protocol.deck.items() if labware == None][0]
					except:
						raise Exception("There is not enough space in the deck for this protocol, try less samples")
						
					pipette_used.tip_racks += [protocol.load_labware(define_tiprack(pipette_used.name), first_position_free)]
					
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
		
		while any(volume > vol_max_tube for volume in volumes_tubes): # If there is some volume that is greater than the max volume we are going to enter in the loop
			number_tubes += 1 # We add one tube so the volume can fit in the tubes
			
			# Now we redistribute the reactions (and correspondant volume) to the tubes so it will be the most homogeneus way
			reactions_per_tube = [int(number_reactions/number_tubes)]*number_tubes
			tubes_to_add_reaction = number_reactions%number_tubes
			for i in range(tubes_to_add_reaction):
				reactions_per_tube[i] += 1
			
			# Calculate the new volumes
			volumes_tubes = [vol_reactive_per_reaction_factor*number_reactions_tube for number_reactions_tube in reactions_per_tube]
		
		# When the volume can fit every tube (exit from th ewhile loop) we return the number of tubes and the reactions that will fit in every tube
		return (number_tubes, reactions_per_tube)
	
	def setting_number_plates (number_plates, labware_name):
		"""
		In this function we will set how many labwares we need of every category (source labwares, final, coldblocks, falcon tube racks, etc)
		
		This function will only set the labwares in the different slots of the deck, with not calculate how many we need,
		this way we do not have to change this function and only change the setting_labware function from protocol to protocol
		"""
		
		try:
			position_plates = [position for position, labware in protocol.deck.items() if labware == None] # We obtain the positions in which there are not labwares
			all_plates = []
			for i in range (number_plates):
				plate = protocol.load_labware(labware_name, position_plates[i])
				all_plates.append(plate)
			return all_plates
		
		except: # There were not enough items in the position_plate to fit the number of labwares that we needed
			raise Exception("There is not enough space in the deck for this protocol, try less samples")
	
	def setting_labware(variables):
		"""
		In this function we will calculate how many labwares of each type do we need to perform the protocol, and then we are going to call
		the function setting_number_plates that will set the labware in the respective positions
		
		In this function we are going to calculate and set the labwares in the deck according to the csv provided, with the exception of the tipracks
		
		This function should be changed depending of how many types of labware we need to set and how many different reactives do we need
		"""
		
		# We start settign the source labware which number has been provided
		labware_source = setting_number_plates(variables.number_source_plates, variables.name_source_plate)
		
		# Set the final plates which number has been previously calculated and it is in the setted_parameters class
		labware_final = setting_number_plates(variables.plates_final, variables.name_final_plate)
		
		list_labware_final = copy.deepcopy(labware_final)
		
		#Now we are going to assign to which final plates the samples from the source plates should go
		for data_source_plate in variables.data_source_plates.values():
			data_source_plate["final_plates"] = list_labware_final[:len(data_source_plate["antibiotics"])]
			del list_labware_final[:len(data_source_plate["antibiotics"])]
		
		
		# Set the antibiotic falcons
		# For that we need to know the maximum volume of the tubes and how many tubes of the reactives we need in total
		max_volume_well = variables.vol_max_tube
		
		tubes_antibiotic_total = 0
		for antibiotic in variables.antibiotics.keys():
			number_reactions_per_antibiotic = 0
			for plate in variables.data_source_plates.values():
				if antibiotic+1 in plate["antibiotics"]:
					# We need this antibiotic +1 because antibiotic indexes are assigned in the script so index starts with 0 and the antibiotics per plate are given by the user so the start index is 1
					number_reactions_per_antibiotic += plate["samples"]
			tubes_antibiotic, reactions_per_tube_antibiotic = number_tubes_needed(variables.volume_antibiotic, number_reactions_per_antibiotic, max_volume_well-10)
			tubes_antibiotic_total += tubes_antibiotic
			variables.antibiotics[antibiotic]["reactions_per_tube"] = reactions_per_tube_antibiotic
			variables.antibiotics[antibiotic]["positions_tubes"] = [] # Fill later, when we set the positions
		
		# Set how many tuberacks now that we now how many tubes of antibiotic we need
		number_wells_tuberack = len(labware_context.get_labware_definition(variables.name_rack_falcons)["wells"])
		tuberacks_needed = math.ceil(tubes_antibiotic_total/number_wells_tuberack)
		labware_falcons = setting_number_plates(tuberacks_needed, variables.name_rack_falcons)
		
		# Now we are going to set the reactives in the coldblock positions, we need to keep track of these positions for liquid movement
		# Get the possible positions merging all the labwares from the tuberacks
		positions_tuberack = []
		for labware in labware_falcons:
			positions_tuberack += labware.wells()
		
		# Assign to each antibiotic the positions of the falcons
		for antibiotic in variables.antibiotics.keys():
			number_tubes = len(list(variables.antibiotics[antibiotic]["reactions_per_tube"]))
			variables.antibiotics[antibiotic]["positions_tubes"] = positions_tuberack[:number_tubes]
			del positions_tuberack[:number_tubes]

		# Finally, we return the positions of the reactives and the labware we have setted with this function
		# Even if we do not return it, some variables from the class setted_parameters have been changed
		return (labware_source, labware_final, labware_falcons)
	
	def position_dispense_aspirate (vol_falcon, theory_position):
		"""
		This function will return the height in which the pipette should aspirate the volume
		
		It is manually measured, meaning that if you change the tubes you should test if this work or redo the heights
		"""
		if vol_falcon <= 100: # The values of comparing are volumes (in uL)
			final_position = theory_position.bottom(z=0.7)
		elif vol_falcon > 100 and vol_falcon <= 3000:
			final_position = theory_position.bottom(z=1)
		elif vol_falcon > 3000 and vol_falcon <= 6000:
			final_position = theory_position.bottom(z = 25)
		elif vol_falcon > 6000 and vol_falcon <= 9000:
			final_position = theory_position.bottom(z = 45)
		elif vol_falcon > 9000:
			final_position = theory_position.bottom(z = 65)
		return final_position

	def distribute_z_tracking(pipette_used, vol_source, vol_distribute_well, pos_source, pos_final):
		"""
		This function will distribute from a pos_source to pos_final (list) taking in account the height of aspiration of the 15mL falcon tube,
		this way the pipette will not get wet during the transfering of liquids
		
		This function is mainly design to distribute reactives to a plate making less pipette movements with the same pipette because the volume to distribute is the same
		"""
		
		while len(pos_final) != 0: # This will go on until there are no elements in the list pos_final (there are no more positions to transfer reactives to)
			# We are going to compare the height before and after aspirating
			if position_dispense_aspirate(vol_source, pos_source).point == position_dispense_aspirate((vol_source-(len(pos_final)*vol_distribute_well)), pos_source).point:
				# Heights are the same, so we are going to take the whole volume and distribute it
				pipette_used.distribute(vol_distribute_well, position_dispense_aspirate(vol_source, pos_source), pos_final, new_tip = "never", disposal_volume = 0)
				vol_source -= vol_distribute_well*len(pos_final)
				pos_final = []
			
			else:
				# Heights are not the same so we are going to take the maximum, distribute it, change the height and then distribute the rest
				pos_final_original = pos_final
				for number_positions in range(1, len(pos_final_original)+1):
					vol_needed = vol_distribute_well*len(pos_final[:number_positions])
					if position_dispense_aspirate(vol_source, pos_source) == position_dispense_aspirate((vol_source-vol_needed), pos_source):
						next
					else:
						pipette_used.distribute(vol_distribute_well, position_dispense_aspirate(vol_source, pos_source), pos_final[:number_positions], new_tip = "never", disposal_volume = 0)
						# If you change the disposal volume you need to take it in account when you are calculating the number of tubes
						vol_source -= vol_distribute_well*len(pos_final[:number_positions])
						pos_final = pos_final[number_positions:]
						break
		return
	
	def set_labware(labware_wells_name):
			"""
			Generator of the positions to transfer, it will give you the next element of a list every time it is called
			"""
			for well in labware_wells_name:
				yield well

	def print_information_user(variables, tuberack_positions):
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
		print("\t- Number total of samples (all source plates): "+str(sum(variables.samples_per_plate.values())))
		print("\t- Number of different antibiotics: "+str(variables.number_antibiotics))
		print("\t- Final volume of output labware (uL): "+str(variables.volume_antibiotic+variables.volume_colonies))
		if variables.replace_tiprack == True and variables.need_change_tiprack == True:
			print("\t- You should stay near the OT to change the tiprack(s)\n")
		elif variables.replace_tiprack == True and variables.need_change_tiprack == False:
			print("\t- There will be no need to change tipracks\n")
		elif variables.replace_tiprack == False:
			print("\t- You have the replace tiprack option as False so there is no need of replacing none of the tipracks\n")
		
		#Deck information of the labware
		print("--------------------------------------------------------------\nDECK LABWARE POSITIONS\nThis information will be provided by the OT App as well")
		for labware_position in protocol.loaded_labwares.values():
			print("\t"+str(labware_position))
	
		# Printing of the reactives and their volume as it will bein the coldblock(s)
		print("\n--------------------------------------------------------------\n15 ML FALCON TUBERACK REACTIVES POSITIONS\nReactives in each position of the coldblock(s) and their respective volume (minimum)")
		for tube_rack in tuberack_positions:
			print("Tuberack "+str(tube_rack))
			tuberack_table = pd.DataFrame(index=["A","B","C"], columns=["1","2","3","4","5"])
			print("\n--> Tuberack in Slot "+str(tube_rack).split(" ")[-1])
			
			for data_antibiotic in variables.antibiotics.values():
				for index_position, tube_antibiotic in enumerate(data_antibiotic['positions_tubes']):
					tube_well = str(tube_antibiotic).split(" ")[0]
					if str(tube_antibiotic).split(" ")[-1] == str(tube_rack).split(" ")[-1]:
						volume_antibiotic_tube = int(math.ceil(variables.volume_antibiotic*data_antibiotic["reactions_per_tube"][index_position]+10))
						tuberack_table[tube_well[1:]][tube_well[0]] = data_antibiotic["name"]+" - "+str(volume_antibiotic_tube)+"uL"
			print(tuberack_table.fillna("-")) # for esthetic purposes
			print("\n--------------------------------------------------------------")
		return

#Body of the script
#--------------------------------
#--------------------------------
	try:
		current_step = "Reading csv and transforming them to parameters/variables"
		# Setting variables and calculating others
		variables_csv = pd.read_csv("/data/user_storage/Variables-AntibioticPlatesCreation-OT.csv", index_col = 0)
		# We are going to convert these parameters into arguments of the class variables and we are going to process some of them so they can be usable (they are going to be dictionaries in their majority)
		variables = setted_parameters(variables_csv)
		
		# Due to the fact that even in the parameters checking we use the labwares definition, we are going to check that they exist outside that function
		try:
			labware_context.get_labware_definition(variables.name_source_plate)
			labware_context.get_labware_definition(variables.name_final_plate)
			labware_context.get_labware_definition(variables.name_rack_falcons)
		except:
			raise Exception("One or more of the introduced labwares are not in the custom labware directory of the opentrons. Check for any typo of the API labware name.")
		
		
		wells_source_plate = len(labware_context.get_labware_definition(variables.name_source_plate)["wells"])
		setted_parameters.proccess_variables(variables, "samples_per_plate", wells_source_plate, variables.number_source_plates)
		
		# Fill the dictionary for data_source_plates (samples and antibiotic), final plates will be filled after
		for number_plate in range(variables.number_source_plates):
				list_antibiotics = variables_csv._get_value("Antibiotics per plate","Value").strip().split(",(")[number_plate].replace("(","").replace(")","").split(",")
				list_antibiotics = [eval(i) for i in list_antibiotics]
				variables.data_source_plates[number_plate] = {"samples":variables.samples_per_plate[number_plate+1], "antibiotics":list_antibiotics,"final_plates":[]} # Final plates will be filled after
				# We need that +1 variables.samples_per_plate[number_plate+1] becaus ethis is a variable setted by the user and their index starts with 1, not 0
		
		# Set some variables of the final plates that we need
		# Number of final plates
		for data_plate in variables.data_source_plates.values():
			variables.plates_final += len(data_plate["antibiotics"])
		
		# Information of the final plates 
		for data_plate in variables.data_source_plates.values():
			for antibiotic in data_plate["antibiotics"]:
				try:
					index_plate = list(variables.data_final_plates.keys())[-1] + 1
				except:
					index_plate = 0
				variables.data_final_plates[index_plate] = {"samples":data_plate["samples"], "antibiotic":antibiotic}
		
		# Now we are going to submit the variables to some previous checking to prevent errors in the protocol
		# Not every error is going to be prevented but most of them can be checked at the beginning
		current_step = "Checking variables for errors"
		errors_variables = check_setted_parameters(variables)
		if len(errors_variables) > 0:
			print("------------------------------------------------\nERROR(S):")
			for error in errors_variables:
				print("\t- "+error)
			print("------------------------------------------------\nBefore proceeding this error(s) should be fixed")
			quit()
		
		current_step = "Defining pipettes and checking them"
		# We can define them because in the parameters check we have seen if there is a pipette that it is not in the mount
		variables.right_pipette = protocol.load_instrument(variables.name_right_pipette, mount = "right")
		variables.left_pipette = protocol.load_instrument(variables.name_left_pipette, mount = "left")
		
		# We are going to perform some pipette checkings because they need to be setted into the instrument context to these controls to be done
		# First we are going to check if the multi and single channel are put correctly
		if variables.right_pipette.channels > 1 and variables.left_pipette.channels == 1:
			# Now we are going to check if the pipettes can take the volumes (because the OT will not give an error about this)
			# In case it is not possible, we are going to raise an exception
			if variables.left_pipette.min_volume > variables.volume_antibiotic:
				raise Exception("Antibiotic volume to transfer can not be picked with "+variables.name_left_pipette)
			if variables.right_pipette.min_volume > variables.volume_colonies:
				raise Exception("Sample volume to transfer can not be picked with "+variables.name_right_pipette)
		else:
			raise Exception("Either left pipette is not single channel, right pipette is nor multichannel or both")
		
		current_step = "Setting labwares into deck (not the tipracks)"
		
		# Set the labware that we need
		labware_source, labware_final, labware_falcons = setting_labware(variables)
		
		# Distribute the antibiotics to the respective plates
		current_step = "Distributing antibiotic(s) to plate(s)"

		# Let's create the positions to distribute for every antibiotic and then distribute it
		for index_antibiotic, data_antibiotic in variables.antibiotics.items():
			positions_antibiotic = []
			for index_plate, data_plate in variables.data_final_plates.items():
				if index_antibiotic == data_plate["antibiotic"]:
					positions_antibiotic += labware_final[index_plate].wells()[:data_plate["samples"]]
				else:
					pass

			check_tip_and_pick(variables.left_pipette, variables)

			# We have the wells that we need to distribute the antibiotic by name
			for index_tube, reactions in enumerate(data_antibiotic["reactions_per_tube"]):
				current_tube = data_antibiotic["positions_tubes"][index_tube]
				wells_to_distribute = positions_antibiotic[:reactions]
				
				distribute_z_tracking(variables.left_pipette, len(wells_to_distribute)*variables.volume_antibiotic+10, variables.volume_antibiotic, current_tube, wells_to_distribute)
				
				del positions_antibiotic[:reactions]
			variables.left_pipette.drop_tip()
		
		# All the antibiotics are in the plates (in the corresponding wells)
		
		# Now we distribute the cells
		current_step = "Distributing samples to the plates with antibiotics"
		number_rows_source_plate = len(labware_context.get_labware_definition(variables.name_source_plate)["ordering"][0])
		for index_plate, data_plate in variables.data_source_plates.items():
			columns_plate = math.ceil(data_plate["samples"]/number_rows_source_plate)
			for final_plate_position in data_plate["final_plates"]:
				for column in range(columns_plate):
					check_tip_and_pick(variables.right_pipette, variables)
					variables.right_pipette.transfer(variables.volume_colonies, labware_source[index_plate].columns(column), final_plate_position.columns(column), new_tip = "never")
					variables.right_pipette.drop_tip()
		
		current_step ="Homing robot"
		protocol.home()
		
		current_step = "Printing User Information"
		print_information_user(variables, labware_falcons) #positions of the antibiotics and the reactions per tube of each one are in the varaibles class
		
	except Exception as e:
		print("-------------------------------------------------------------------------------")
		print("--> ERROR STEP:\n" + current_step)
		logger.critical(e, exc_info = True)
		print("-------------------------------------------------------------------------------")
