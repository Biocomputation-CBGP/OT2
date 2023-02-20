"""
Python script destined to OT-2
This script performs a transfer of a defined by the user number of samples from different plates to one or more final(s) plates.
One example of usage is to merge several libraries in 1 plate to test them more easily.
This scirpt needs number of source plates +1 csv attached to perform the running (variables and as much maps of identities as source plates)
and will give an output file (txt) with some instructions and a map of the final plate samples positions
For more info go to https://github.com/Biocomputation-CBGP/OT2/tree/main/ColonySelection and/or
www.protocols.io/view/ot-2-protocol-to-transfer-volume-from-several-plat-6qpvr4o62gmk/v1
"""

## Packages needed for the running of the protocol
import opentrons.execute
from opentrons import protocol_api
import pandas as pd
import logging
import math
import copy
import numpy as np
import random

## Structure needed for the running of th escript in the OT
metadata = {
'apiLevel':'2.12'
}

def run(protocol: protocol_api.ProtocolContext):
	# Before everything because we will need it for the labware setting
	labware_context = opentrons.protocol_api.labware

	# We are going to create a logger to control errors that can occure in the script that are not taked in account and some of ours that will quit the program
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
		
		def __init__(self, varibles_csv):
			self.number_source_plates = int(varibles_csv._get_value("Number of Source Plates","Value"))
			self.samples_per_plate = varibles_csv._get_value("Samples Per Plate","Value")
			self.name_right_pipette = varibles_csv._get_value("Name Right Pipette","Value")
			self.name_left_pipette = varibles_csv._get_value("Name Left Pipette","Value")
			# First we are going to check the parameters and then we are going to settle the pipettes
			self.right_pipette = None
			self.left_pipette = None
			self.starting_tip_right_pip = varibles_csv._get_value("Initial Tip Right Pipette","Value")
			self.starting_tip_left_pip = varibles_csv._get_value("Initial Tip Left Pipette","Value")
			if variables_csv._get_value("Replace Tipracks","Value").lower() == "true":
				self.replace_tiprack = True
			elif variables_csv._get_value("Replace Tipracks","Value").lower() == "false":
				self.replace_tiprack = False
			else:
				raise Exception("Replace Tiprack value can only be True or False")
			# self.name_map = "/data/user_storage/"+variables_csv._get_value("Name maps","Value")
			self.name_map = variables_csv._get_value("Name maps","Value")
			self.volume_transfer_sample = float(variables_csv._get_value("Volume Colony transfer (uL)","Value"))
			self.volume_transfer_water = float(variables_csv._get_value("Volume Water transfer (uL)","Value"))
			self.water = {"position_tubes":[], "reactions_per_tube":[]}
			self.name_rack_falcon = variables_csv._get_value("Name 15mL Falcon Rack","Value")
			self.name_source_plate = variables_csv._get_value("Name Source Plates","Value")
			self.name_final_plate = variables_csv._get_value("Name Final Plates","Value")
			self.maps_source_plates_names = variables_csv._get_value("Name Source Plates Maps","Value")[1:-1].replace(" ","").split(",") # We have the list of the maps files
			self.maps_source_plate = {} # We are going to fill this dictionary with the dataframes corresponding to the files
			for index_map, map_source in enumerate(self.maps_source_plates_names):
				# self.maps_source_plate[index_map] = pd.read_csv("/data/user_storage/"+map_source+".csv", index_col = 0)
				self.maps_source_plate[index_map] = pd.read_csv(map_source+".csv", index_col = 0)
			self.type_selection = variables_csv._get_value("Type of Sample Selection", "Value").lower()
			self.number_samples_source_plates = variables_csv._get_value("Number of samples in every source plate", "Value")
			self.index_start_source_plate = variables_csv._get_value("Index of start cell in source plate", "Value")
			self.need_change_tiprack = False # initialize with a value, it will be changed in case that there is a need to replace it
			
			return
		
		def proccess_variables(self, name_parameter, value_default_empty, must_entries):
			"""
			This function will process some text that is extracted from the variables and make it to a dictionary so we can work with it
			Some of the variables need processing to be usable in the script, usually is from a "(,),(,)" format to a dictionary
			
			The function need the name of the parameter in the class, default value of the variable and how many items the dictiony needs to have
			In case that the strcutre of the parameters are different, this function should be changed
			
			In this case we have 2 types of variables structures that should be treated differently
			"""
			dict_final = {}

			if name_parameter in ["number_samples_source_plates", "index_start_source_plate"]:
				value_parameter = getattr(self, name_parameter).replace(" ","")[1:-1].split('),(')
				
				# We separate the text value into the real values (they were separated by commas)
				for i, element in enumerate(value_parameter):
					value_parameter[i] = element.split(",")

				if len(value_parameter[0][0]) != 0:# This checks if the the variable has a default value (like 0 or -) or if it has other values
												   # This works with the process that we have done in  the previous part
					for well, plate in value_parameter:
						dict_final[int(plate)-1] = int(well)
				else:
						pass
			else:
				# We separate the text value into the real values (they were separated by commas) from the beginning, that is the difference between this category and the previous one
				value_parameter = getattr(self, name_parameter).replace(" ","")[1:-1].split(',')
				
				if len(value_parameter) != 0:# This checks if the the variable has a default value (like 0 or -) or if it has other values
											 # This works with the process that we have done in  the previous part
					for plate, well in enumerate(value_parameter):
							dict_final[int(plate)] = int(well)
				else:
						pass
			
			# We make sure that we have the number of entries that are specified in must_entries and if there is one missing we fill it with the default value for the item
			# If the dictionary is empty (when value in template is the default one) we will create all the items with the default value of the items (for example can be the number of wells of a plate or the index)
			for i in range(must_entries):
				try:
					dict_final[i]
				except:
					dict_final[i] = value_default_empty
			
			setattr(self, name_parameter, dict_final)
			return

# Functions
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
		
		# Index of any of the source plates is larger than the wells of that plate
		if any(index_labware > len(labware_context.get_labware_definition(variables.name_source_plate)["wells"]) for index_labware in variables.index_start_source_plate.values()):
			errors.append("One or more of the start indexes is higher than the number of wells in the source plate")
			
		# Number of samples of any of the source plates is larger than the wells of that plate
		if any(number_samples > len(labware_context.get_labware_definition(variables.name_source_plate)["wells"]) for number_samples in variables.number_samples_source_plates.values()):
			errors.append("One or more of the samples per source plate is higher than the number of wells in the source plate")
			
		# Index+number of samples of any of the source plates is larger than the wells of that plate
		if any(index+number > len(labware_context.get_labware_definition(variables.name_final_plate)["wells"]) for index, number in list(zip(dict(sorted(variables.index_start_source_plate.items())).values(),dict(sorted(variables.number_samples_source_plates.items())).values()))):
			errors.append("One or more of the samples per plate + starting index of the plate is higher than the numer of wells in the source plate")
		# Check if the samples that we are going to take are lower or equal to the number of samples in each of the source plate
		sorted_values_by_keys = [variables.number_samples_source_plates[key] for key in sorted(variables.number_samples_source_plates.keys())]
		# This sorting is needed because the list of samples per plate is ordered from the first source plate to the last one and we are going to "merge" both lists
		if any(samples_take > samples_plate for samples_take, samples_plate in list(zip(variables.samples_per_plate.values(), sorted_values_by_keys))):
			errors.append("In one or more of the source plates you are trying to take more samples to transfer than the number of samples that there are in that plate")

		# Check if we have given all the maps to the source plates
		if len(variables.maps_source_plates_names) != variables.number_source_plates:
			errors.append("The number of source plates and the number of maps do not match")
		
		# Check if there is any typo in the starting tip of both pipettes
		if variables.name_right_pipette not in ["None","none","-"]:
			if not any(variables.starting_tip_right_pip in columns for columns in labware_context.get_labware_definition(define_tiprack(variables.name_right_pipette))["ordering"]):
				errors.append("Starting tip of right pipette is not valid, check for typos")
		if variables.name_left_pipette not in ["None","none","-"]:
			if not any(variables.starting_tip_left_pip in columns for columns in labware_context.get_labware_definition(define_tiprack(variables.name_left_pipette))["ordering"]):
				errors.append("Starting tip of left pipette is not valid, check for typos")
		
		# Check if the type of selection variable is one of the established ones
		if variables.type_selection not in ["random","first","last"]:
			errors.append("'"+variables.type_selection+"' not recognised as a 'type of sample selection'. Options are 'random', 'first' and 'last'")
			
		# Check if the number of elements in samples per plate is the same as number of sourc eplates, because if we are not going to take from it, it doesnt make sense to have it in the deck
		if any(number_samples == 0 for number_samples in variables.samples_per_plate.values()):
			errors.append("You are not taking any samples from one of the source plates")
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
					# Careful with this part if you are traspassing this script into jupyter because this will crash your jupyter (will wait until resume and it does not exist)
					variables.need_change_tiprack = True
					protocol.pause("Replace Empty Tiprack With A Full One And Press Resume In OT-APP")
					pipette_used.reset_tipracks()
			
			# Finally, we pick up the needed tip        
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

		return

	def setting_labware(variables):
		"""
		In this function we will calculate how many labwares of each type do we need to perform the protocol, and then we are going to call
		the function setting_number_plates that will set the labware in the respective positions
		
		In this function we are going to calculate and set the labwares in the deck according to the csv provided, with the exception of the tipracks
		
		This function should be changed depending of how many types of labware we need to set and how many different reactives do we need
		"""
		
		# We start settign the source labware which number has been provided
		labware_source = setting_number_plates(variables.number_source_plates, variables.name_source_plate)
		
		# Set the final plates
		# We need to calculate how many plates do we need by sum (cells by plate) / wells per final plate
		wells_we_need = sum(variables.samples_per_plate.values())
		wells_final_plate = len(list(labware_context.get_labware_definition(variables.name_final_plate)["wells"]))
		labware_final = setting_number_plates(math.ceil(wells_we_need/wells_final_plate), variables.name_final_plate)
		
		list_labware_final = copy.deepcopy(labware_final)
			
		# Set the reactive falcons
		# For that we need to know the maximum volume of the tubes and how many tubes of the reactives we need in total
		max_volume_well = variables.vol_max_tube
		
		# We are going to put this in a condition because maybe we only want to perform a transfering and not a media+samples
		if variables.volume_transfer_water > 0:
			tubes_water, variables.water["reactions_per_tube"] = number_tubes_needed(variables.volume_transfer_water, wells_we_need, max_volume_well-10)
			# Set how many tuberacks
			number_wells_tuberack = len(labware_context.get_labware_definition(variables.name_rack_falcon)["wells"])
			tuberacks_needed = math.ceil(tubes_water/number_wells_tuberack)
			labware_falcons = setting_number_plates(tuberacks_needed, variables.name_rack_falcon)
			
			# Now we are going to set the reactives in the coldblock positions, we need to keep track of these positions for liquid movement
			# Get the possible positions 
			positions_tuberack = []
			for labware in labware_falcons:
				positions_tuberack += labware.wells()
			
			variables.water["position_tubes"] = positions_tuberack[:tubes_water]
		else: #We do not want to do media+sample
			labware_falcons = []
		
		# Finally, we return the positions of the reactives and the labware we have setted with this function
		return (labware_source, labware_final, labware_falcons)

	def position_dispense_aspirate (vol_falcon, theory_position):
		"""
		This function will return the height in which the pipette should aspirate the volume
		
		It is manually measured, meaning that if you change the tubes you should test if this work or redo the heights
		"""
		if vol_falcon <= 100:
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
				# If you change the disposal volume you need to take it in account when you are calculating the number of tubes
				vol_source -= vol_distribute_well*len(pos_final)
				pos_final = []
			
			else:
				# Heights are not the same so we are going to take the maximum, distribute it, change the height and then distribute the rest
				pos_final_original = pos_final
				for number_positions in range(1, len(pos_final_original)+1): #vamos mirando hasta cuanto puede distribuir sin tener que cambair la altura
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
				
	def give_me_optimal_pipette_to_use(aVolume, pipette_r, pipette_l):
		"""
		Function that will return the optimal pipette to use for the volume that we want to handle.
		
		In case that it is a great volume (higher than the maximal volume of both pipettes) will return the pipette that will give the minimal quantity of movements
		
		If none of the pipettes attached can pick the volume (because it is too small) the function will raise an error
		
		For the correct functioning of this function at least 1 pipette should be attached to the OT, otherwise, the function will raise an error
		"""
		
		if pipette_r == None and pipette_l == None: # No pipettes attached
			raise Exception("There is not a pippette attached")
		
		# First we look if one of them is the only option
		elif pipette_r == None or pipette_l == None: # One mount is free, only need that the volume is more than the min of the pipette
			if pipette_r == None and aVolume >= pipette_l.min_volume:
				return pipette_l
			elif pipette_l == None and aVolume >= pipette_r.min_volume:
				return pipette_r
			else: # One of them does not exist and the other one is not valid
				raise Exception("Volume cannot be picked with pipette(s) associated. Try another pipette combination")
				# Error control
				
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
		return
	
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

	def map_tracking(pos_source, pos_final, dataframes_with_index, map_source_plate):
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
		
		
		identifier_sample = map_source_plate[int(pos_source_deck)-1][column_pos_cource][row_pos_source]
	
		
		dataframes_with_index[int(pos_final_deck)][int(pos_final_labware[1:])][pos_final_labware[0]] = str(identifier_sample+" on "+pos_source_deck)
		# We are going to still put the position on the deck because 2 samples can have the same identifier in different plates
		
		return
	
	def wells_transfer(source_plates, final_plates, variables):
		"""
		Function that establish the final wells (as many as the sum of the samples from every source plate)
		and the list of wells which we are going to take from the source plates
		
		The list of wells from the osurce plate are going to be variable depending on the type of picking that is setted in variables
		"""
		# We state which are going to be the final wells, the ones that we are going to transfer to
	
		final_wells = sum([labware.wells() for labware in final_plates],[])[:sum(variables.samples_per_plate.values())]

		# Now lets take the source wells that we are going to aspirate from
		source_well_distribute = []
		for index_labware, samples_labware in variables.samples_per_plate.items():
			possible_wells_plate = source_plates[index_labware].wells()[variables.index_start_source_plate[index_labware]:variables.index_start_source_plate[index_labware]+variables.number_samples_source_plates[index_labware]]
			if variables.type_selection == "first":
				source_well_distribute += possible_wells_plate[:samples_labware]
			elif variables.type_selection == "last":
				source_well_distribute += list(reversed(possible_wells_plate))[:samples_labware]
			elif variables.type_selection == "random":
				source_well_distribute += random.sample(possible_wells_plate, samples_labware)
		
		return (source_well_distribute, final_wells)
		
	def print_information_user(variables, tuberack_positions):
		"""
		This function will print the information that the user should know as the volumes of each tube and where to place it.
		
		This function should be customized to every protocol
		
		This will be printed when done the simulation because there is info needed to know before running the protocol in the OT
		"""
		
		print("""--------------------------------------------------------------\n--------------------------------------------------------------\nINFORMATION FILE
	Data the users need to set the OT and perform the protocol
	It is the output of the python file for the protocol when executed with the opentrons_simulate option\n""")
		# Some general information that the user has setted
		print("--------------------------------------------------------------\nGENERAL INFORMATION\nThis details are set by the user, this is only a remainder of some variables")
		print("\t- Number total of samples (all source plates): "+str(sum(variables.samples_per_plate.values())))
		print("\t- Final volume of output labware (uL): "+str(variables.volume_transfer_water+variables.volume_transfer_sample))
		if variables.replace_tiprack == True and variables.need_change_tiprack == True:
			print("\t- You should stay near the OT to change the tiprack(s)\n")
		elif variables.replace_tiprack == True and variables.need_change_tiprack == False:
			print("\t- There will be no need to change tipracks\n")
		elif variables.replace_tiprack == False:
			print("\t- You have the replace tiprack option as False so there is no need of replacing none of the tipracks\n")
		
		# Deck information of the labware
		print("--------------------------------------------------------------\nDECK LABWARE POSITIONS\nThis information will be provided by the OT App as well")
		for labware_position in protocol.loaded_labwares.values():
			print("\t"+str(labware_position))
	
		# Printing of the reactives and their volume as it will be in the falcon tube racks(s) in case that there is one or more
		if variables.volume_transfer_water > 0:
			print("\n--------------------------------------------------------------\n15 ML FALCON TUBERACK REACTIVES POSITIONS\nReactives in each position of the coldblock(s) and their respective volume (minimum)")
			for tube_rack in tuberack_positions:
				print("Tuberack "+str(tube_rack))
				tuberack_table = pd.DataFrame(index=["A","B","C"], columns=["1","2","3","4","5"])
				print("\n--> Tuberack in Slot "+str(tube_rack).split(" ")[-1])
				for index_position, tube_water in enumerate(variables.water['position_tubes']):
					tube_well = str(tube_water).split(" ")[0]
					if str(tube_water).split(" ")[-1] == str(tube_rack).split(" ")[-1]:
						volume_water_tube = int(math.ceil(variables.volume_transfer_water*variables.water["reactions_per_tube"][index_position]+10))
						tuberack_table[tube_well[1:]][tube_well[0]] = "Water - "+str(volume_water_tube)+"uL"
				print(tuberack_table.fillna("-")) # for aesthetic purposes
				print("\n--------------------------------------------------------------")
		else:
			print("\n--------------------------------------------------------------")

		return
	
# Body of the script
#--------------------------------
#--------------------------------

# We use a try and except so we can give some more user-friendly infromation about the errors that we can receive
	try:
		current_step = "Reading csv and transforming them to parameters/variables"
		# Loading of the csv parameters and using the first column (the name of the variables) as index
		# variables_csv = pd.read_csv("/data/user_storage/Variables-SamplesMerging-OT.csv", index_col = 0)
		variables_csv = pd.read_csv("Variables-SamplesMerging-OT.csv", index_col = 0)

		# We are going to convert these parameters into arguments of the class variables and we ar egoing to process some of them so they can be usable (they are going to be dictionaries)
		variables = setted_parameters(variables_csv)
		
		# Due to the fact that even in the parameters checking we use the labwares definition, we are going to check that they exist outside that function
		try:
			labware_context.get_labware_definition(variables.name_source_plate)
			labware_context.get_labware_definition(variables.name_final_plate)
			labware_context.get_labware_definition(variables.name_rack_falcon)
		except:
			raise Exception("One or more of the introduced labwares are not in the custom labware directory of the opentrons. Check for any typo of the API labware name.")
		
		# We are going to process some of the variables so they can be usable (they are going to be dictionaries)
		setted_parameters.proccess_variables(variables, "samples_per_plate", 0, variables.number_source_plates) # We take 0 as default
		
		# We need to set the number of wells in the plates so we can use it as defaut number in some variables
		number_wells_source_plate = len(labware_context.get_labware_definition(variables.name_source_plate)["wells"])
		# Process some variables so they are turn usually from a character list to a dictionary that is going to be used in the rest of the protocol
		setted_parameters.proccess_variables(variables, "number_samples_source_plates", number_wells_source_plate, variables.number_source_plates)
		setted_parameters.proccess_variables(variables, "index_start_source_plate", 0, variables.number_source_plates)
		
		# Now that we have set the variables and process them, lets check if there is not an error than we can catch before running the rest of the script
		errors_variables = check_setted_parameters(variables)
		if len(errors_variables) > 0:
			print("------------------------------------------------\nERROR(S):")
			for error in errors_variables:
				print("\t- "+error)
			print("------------------------------------------------\nBefore proceeding this error(s) should be fixed")
			quit()
		
		# We define the pipettes
		current_step = "Stating pipettes and associated tip racks"
		if variables.name_right_pipette == "None" or variables.name_right_pipette == "-":
			variables.right_pipette = None
		else:
			variables.right_pipette = protocol.load_instrument(variables.name_right_pipette, "right", tip_racks=[])
			
		if variables.name_left_pipette == "None" or variables.name_left_pipette == "-":
			variables.left_pipette = None
		else:
			variables.left_pipette = protocol.load_instrument(variables.name_left_pipette, "left", tip_racks=[])
		
		current_step = "Setting labwares into deck (not the tipracks)"
		labware_source, labware_final, labware_falcons = setting_labware(variables)
		
		current_step = "Deciding the final wells and the samples that there are going to be taken"
		# We take the wells that we are going to need and the ones taht we are going to dispense to
		source_well_take, wells_dispense = wells_transfer(labware_source, labware_final, variables)
		
		# Create the dataframes for the mapping
		# First we create the maps for filling, we are creating it even if we dont want it
		maps_final_plates = pds_labware_sources(labware_final)
		
		#In case that we want to distribute some reactive to all wells, we are going to enter in that condition, otherwise, we are going to pass and distribute the samples only
		if variables.volume_transfer_water > 0:
			current_step = "Transfering water to final wells"
			
			wells_distribute_water = copy.deepcopy(wells_dispense)
			pipette_use = give_me_optimal_pipette_to_use(variables.volume_transfer_water, variables.right_pipette, variables.left_pipette)
			check_tip_and_pick(pipette_use, variables)
			
			# We are going to perform the distribution with height tracking for every tube, that is why we stored the reactions per tube
			for index_tube, reactions in enumerate(variables.water["reactions_per_tube"]):
					current_tube = variables.water["position_tubes"][index_tube]
					wells_to_distribute = wells_distribute_water[:reactions]
					
					distribute_z_tracking(pipette_use, (len(wells_to_distribute)*variables.volume_transfer_water)+10, variables.volume_transfer_water, current_tube, wells_to_distribute)
					
					del wells_distribute_water[:reactions]
			pipette_use.drop_tip()
		
		# Distribute and track in the dataframes corresponding to the final labwares
		current_step = "Transfering samples from source plate(s) to final merged plate(s)"
		
		# First we create the generator of wells from which we are going to take a sample from and the one sthat we are going to distribute them
		source_well_distribute = set_labware(source_well_take)
		final_well_distribute = set_labware(wells_dispense)
		used_pipette = give_me_optimal_pipette_to_use(variables.volume_transfer_sample, variables.right_pipette, variables.left_pipette)
		
		for i in range(len(wells_dispense)): # Iterate over all the wells that we are going to transfer
			# Define the wells
			well_source = next(source_well_distribute)
			well_final = next(final_well_distribute)
			
			# Transfer
			check_tip_and_pick(used_pipette, variables)
			used_pipette.transfer(variables.volume_transfer_sample, well_source, well_final, new_tip = "never")
			used_pipette.drop_tip()
			
			# Map
			map_tracking(well_source, well_final, maps_final_plates, variables.maps_source_plate)
		
		# Export the maps
		current_step = "Importing the map of cells and primers for each final plate"
		for map_source in maps_final_plates:
			maps_final_plates[map_source] = maps_final_plates[map_source].fillna(value="-")
			maps_final_plates[map_source].to_csv(str(variables.name_map)+"_"+str(map_source)+".csv")
		
		# Printing user's info
		current_step = "Printing of warnings and user information as position and volumes of reactives"
		print_information_user(variables, labware_falcons)
		
		# Homing
		current_step = "Final homing of the robot arm"
		protocol.home()
		
	except Exception as e:
		print("-------------------------------------------------------------------------------")
		print("--> ERROR STEP:\n" + current_step)
		logger.critical(e, exc_info = True)
		print("-------------------------------------------------------------------------------")
