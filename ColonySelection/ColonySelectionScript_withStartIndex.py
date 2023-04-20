"""
Python script destined to OT-2
This script performs a selection of colonies that have a certain value in 2 conditions, this can be OD in differnt antibiotics, OD in differnt time periods, GFP, etc
In one of the file has to have a lower value than the threshold and in other a higher value
This scirpt needs 3 csv attached to perform the running (variables and 2 values files) and will give an output file (txt) with some instructions
For more info go to https://github.com/Biocomputation-CBGP/OT2/tree/main/ColonySelection
This is a script that variates the https://www.protocols.io/view/ot-2-counter-selection-5qpvor5xdv4o/v1 adding 1 variable
"""

"""
Hay que hacer lo siguiente:
 - añadir la variante en la lectura del csv
 - control de que esta variable este bien y que se acepte el - y numeros nada mas
 - asegurarse de que se puede poner mas de 1 placa necesario por si el indice, y aun asi, yo eso lo podnria como un warning
 de la manera en la que si es mas de 1 se va a formar
 
Como cambiar lo del indice:
 1. averiguar cuantas hacen falta haciendo indice+samples_selected/numero_pocillos_labware para 1 replica, i.e, solo para 1 glycerol o solo para 1 pcr (es lo mismo)
 2. Sumar todas los plates necesarios para saber si caben
 3. Fusionar las listas de esos plates, la lista de los wells para cada uno de ellos, por una parte los de glycerol y por otro lso de PCR
 3. Seleccionar low wells a los que les ponemos con los indices:indices+samples
"""



# This script has an inmediate upgrade, making the reactives general, i.e., not only glycerol and water but general reactives


## Packages needed for the running of the protocol
import opentrons.execute
from opentrons import protocol_api
import pandas as pd
import numpy as np
import logging
import math
import copy

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
		number_source_plates = 1 # We only have 1 source plate that we are going to select from
		vol_max_tube = 15000 # The same always because the transfer is setted with the heights of the 15mL falcon
							 # One of the few labware restrictions of this protocol
		def __init__(self, variables_csv):
			"""
			This function will take the pandas dataframe that will be the table of the csv variables
			"""
			self.grow_OD = float(variables_csv._get_value("Grow OD","Value"))
			self.volume_transfer_colony = float(variables_csv._get_value("Vol Transfer Colony (uL)", "Value"))
			self.volume_transfer_glycerol = float(variables_csv._get_value("Vol Transfer Glycerol (uL)", "Value"))
			self.volume_transfer_water = float(variables_csv._get_value("Vol Transfer Water (uL)", "Value"))
			self.file_transform = "/data/user_storage/"+variables_csv._get_value("Antibiotic Transposition Genome File Name", "Value") + ".csv" # In this file the selection values are higher than the grow_OD
			self.file_integration = "/data/user_storage/"+variables_csv._get_value("Antibiotic Ampiciline Plasmid File Name", "Value") + ".csv" # In this file the selection values are lower than the grow_OD
			# self.file_transform = variables_csv._get_value("Antibiotic Transposition Genome File Name", "Value") + ".csv" # In this file the selection values are higher than the grow_OD
			# self.file_integration = variables_csv._get_value("Antibiotic Ampiciline Plasmid File Name", "Value") + ".csv" # In this file the selection values are lower than the grow_OD
			self.left_pipette_name = variables_csv._get_value("Name Pipette Left Mount", "Value")
			self.right_pipette_name = variables_csv._get_value("Name Pipette Right Mount", "Value")
			self.left_pipette = None
			self.right_pipette = None
			self.starting_tip_right_pip = variables_csv._get_value("Starting tip Right Pipette", "Value")
			self.starting_tip_left_pip = variables_csv._get_value("Starting tip Left Pipette", "Value")
			self.name_final_plate = variables_csv._get_value("Name Final Plate", "Value")
			self.name_source_plate = variables_csv._get_value("Name Source Plate", "Value")
			self.number_glycerol_plates = int(variables_csv._get_value("Number of Glycerol Plates", "Value"))
			self.number_pcr_plates = int(variables_csv._get_value("Number of PCR Plates", "Value"))
			self.name_rack_falcons = variables_csv._get_value("Name Rack Falcon 15mL", "Value")
			self.final_map_name = "/data/user_storage/"+variables_csv._get_value("Final Map Name", "Value")+".csv"
			# self.final_map_name = variables_csv._get_value("Final Map Name", "Value")
			self.reactives_info = {"glycerol":{"position_tubes":[],"reactions_per_tube":[],"destination_plates":[]},"water":{"position_tubes":[],"reactions_per_tube":[],"destination_plates":[]}}
			if variables_csv._get_value("Replace Tiprack","Value").lower() == "true":
				self.replace_tiprack = True
			elif variables_csv._get_value("Replace Tiprack","Value").lower() == "false":
				self.replace_tiprack = False
			else:
				raise Exception("Replace Tiprack value can only be True or False")
			self.index_start_final_plate = int(variables_csv._get_value("Index Start Final Plate", "Value"))
			self.maps_selected_cells = None
			
	def check_setted_parameters (variables):
		"""
		Function that will check the variables of the Template and will raise errors that will crash the OT run
		It is a validation function of the variables checking errors or inconsistencies
		
		This function is dependant again with the variabels that we have, some checks are interchangable between protocols, but some of them are specific of the variables
		"""
		# Initialize the errors list in which we are going to store the messages that we need to display for the user
		errors = []
		
		# We are going to check that the colonies + antibiotic is not more than the max volume of the wells in the final plates
		max_volume_well = float(list(labware_context.get_labware_definition(variables.name_final_plate)["wells"].values())[0]['totalLiquidVolume'])
		if (variables.volume_transfer_colony + variables.volume_transfer_water > max_volume_well) or (variables.volume_transfer_colony + variables.volume_transfer_glycerol > max_volume_well):
			errors.append("The sum of sample and reactive volumes exceeds the max volume of final plate wells")
		else:
			pass
		
		# We are going to check if the OD files exist in the /data/user_storage directory
		try:
			open(variables.file_transform,"r")
			open(variables.file_integration,"r")
		except:
			errors.append("One of the input csv files is not in the /data/user_storage/ directory")
		
		# Check if there is any typo in the starting tip of both pipettes
		try:
			if variables.starting_tip_right_pip not in labware_context.get_labware_definition(define_tiprack(variables.right_pipette_name))["groups"][0]["wells"]:
				errors.append("Starting tip of right pipette is not valid, check for typos")
			if variables.starting_tip_left_pip not in labware_context.get_labware_definition(define_tiprack(variables.left_pipette_name))["groups"][0]["wells"]:
				errors.append("Starting tip of left pipette is not valid, check for typos")
		except:
			errors.append("At least one of the pipettes is not established, check for typos in the name")
			
		# Check if the volume that we are trying to take is bigger than the one in the well
		vol_colonies_total = (variables.number_glycerol_plates + variables.number_pcr_plates)*variables.volume_transfer_colony
		max_volume_well_source_plates = float(list(labware_context.get_labware_definition(variables.name_source_plate)["wells"].values())[0]['totalLiquidVolume'])
		if vol_colonies_total > max_volume_well_source_plates:
			errors.append(f'To run the protocol {vol_colonies_total}ul of each colony is needed and the max volume of each well in source plate is {max_volume_well_source_plates}ul')
		
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

	def check_tip_and_pick (pipette_used, variables):
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
	
	def setting_labware(variables, total_number_colonies_pick):
		"""
		In this function we will calculate how many labwares of each type do we need to perform the protocol, and then we are going to call
		the function setting_number_plates that will set the labware in the respective positions
		
		In this function we are going to calculate and set the labwares in the deck according to the csv provided, with the exception of the tipracks
		
		This function should be changed depending of how many types of labware we need to set and how many different reactives do we need
		"""
		
		# We start setting the source labware which number has been provided
		labware_source = setting_number_plates(variables.number_source_plates, variables.name_source_plate)
		# Set the final plates
		# Find out how many plates do we need for 1 of glycerol or 1 PCR. We need to do this to set how many sets of glycerol plates we need in case that it doesnt fit with the index start + colonies to pick
		plate_per_samples = math.ceil((variables.index_start_final_plate+total_number_colonies_pick)/len(labware_context.get_labware_definition(variables.name_final_plate)["wells"]))
		variables.reactives_info["glycerol"]["destination_plates"] = []
		for i in range(variables.number_glycerol_plates):
			variables.reactives_info["glycerol"]["destination_plates"].append(setting_number_plates(plate_per_samples, variables.name_final_plate))
		variables.reactives_info["water"]["destination_plates"] = []
		for i in range(variables.number_pcr_plates):
			variables.reactives_info["water"]["destination_plates"].append(setting_number_plates(plate_per_samples, variables.name_final_plate))
		# We are going to end up havin a set for every set of glyceorl or pcr plates, for each "replica"
		
		# Set the antibiotic falcons
		# For that we need to know the maximum volume of the tubes and how many tubes of the reactives we need in total
		max_volume_well = variables.vol_max_tube
		
		# Calculate the neccessary volume of glycerol and water
		volume_needed_glycerol = variables.volume_transfer_glycerol*total_number_colonies_pick*variables.number_glycerol_plates
		volume_needed_water = variables.volume_transfer_water*total_number_colonies_pick*variables.number_pcr_plates
		
		# Set the number of tubes needed for both reactives
		tubes_glycerol, variables.reactives_info["glycerol"]["reactions_per_tube"] = number_tubes_needed(variables.volume_transfer_glycerol, total_number_colonies_pick*variables.number_glycerol_plates, max_volume_well-10)
		tubes_water, variables.reactives_info["water"]["reactions_per_tube"] = number_tubes_needed(variables.volume_transfer_water, total_number_colonies_pick*variables.number_pcr_plates, max_volume_well-10)
		
		# Set how many tuberacks
		number_wells_tuberack = len(labware_context.get_labware_definition(variables.name_rack_falcons)["wells"])
		tuberacks_needed = math.ceil((tubes_glycerol + tubes_water) / number_wells_tuberack)
		labware_falcons = setting_number_plates(tuberacks_needed, variables.name_rack_falcons)
		positions_labware_falcon = []
		for labware in labware_falcons:
			positions_labware_falcon += labware.wells()
		for reactive, data_reactive in variables.reactives_info.items():
			data_reactive["position_tubes"] = positions_labware_falcon[:len(data_reactive["reactions_per_tube"])] 
			del positions_labware_falcon[:len(data_reactive["reactions_per_tube"])]
		
		# Finally, we return the positions of the reactives and the labware we have setted with this function
		return (labware_source, labware_falcons)

	def position_dispense_aspirate (vol_falcon, theory_position):
		"""
		This function will return the height in which the pipette should aspirate the volume
		
		It is manually measured, meaning that if you change the tubes you should test if this work or redo the heights
		"""
		if vol_falcon <= 100: # The values of comparing are volumes (in uL)
			final_position = theory_position.bottom(z = 0.7)
		elif vol_falcon > 100 and vol_falcon <= 3000:
			final_position = theory_position.bottom(z = 1)
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
	

	def print_information_user(variables, tuberack_positions, samples_picked):
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
		print("\t- Number total of selected colonies (all source plates): "+str(len(samples_picked)))
		print("\t- Volume (uL) of each colony that is needed: "+str((variables.number_glycerol_plates + variables.number_pcr_plates)*variables.volume_transfer_colony))
		print("\t- Number of total final plates: "+str(variables.number_glycerol_plates+variables.number_pcr_plates))
		print("\t- Final volume of output glycerol labware (uL): "+str(variables.volume_transfer_colony+variables.volume_transfer_glycerol))
		print("\t- Final volume of output PCR labware (uL): "+str(variables.volume_transfer_colony+variables.volume_transfer_water))
		print("\t- Map of samples location with their respective identities is going to be stored at "+str(variables.final_map_name))
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
	
		# Printing of the reactives and their volume as it will bein the coldblock(s)
		print("\n--------------------------------------------------------------\n15 ML FALCON TUBERACK REACTIVES POSITIONS\nReactives in each position of the coldblock(s) and their respective volume (minimum)")
		for tube_rack in tuberack_positions:
			# We have a restriction of the positions, that is why the tube rack should have the well distribution of the original OT falcon tuberack
			print("Tuberack "+str(tube_rack))
			tuberack_table = pd.DataFrame(index=["A","B","C"], columns=["1","2","3","4","5"])
			print("\n--> Tuberack in Slot "+str(tube_rack).split(" ")[-1])
			
			for reactive, data_reactives in variables.reactives_info.items():
				for index_position, tube_reactive in enumerate(data_reactives['position_tubes']):
					tube_well = str(tube_reactive).split(" ")[0]
					if str(tube_reactive).split(" ")[-1] == str(tube_rack).split(" ")[-1]:
						if reactive == "glycerol":
							volume_reactive = math.ceil(variables.volume_transfer_glycerol*data_reactives["reactions_per_tube"][index_position]+10)
							tuberack_table[tube_well[1:]][tube_well[0]] = str(reactive)+" - "+str(volume_reactive)+"uL"
						elif reactive == "water":
							volume_reactive = math.ceil(variables.volume_transfer_water*data_reactives["reactions_per_tube"][index_position]+10)
							tuberack_table[tube_well[1:]][tube_well[0]] = str(reactive)+" - "+str(volume_reactive)+"uL"
			print(tuberack_table.fillna("-")) # for aesthetic purposes
			print("\n--------------------------------------------------------------")
		return
	
	def pds_labware_sources(labware_source):
			"""
			Function that will create a dataframe (pandas) with the dimension of the final plates so it can be filled with the positions where each sample is dispensed (list of positions)
			
			This function is setted to suit labware which has rows with alphabetical name (letters, only 1) and columns with numbers
			"""
			
			labware_wells = labware_context.get_labware_definition(labware_source[0].name)["groups"][0]["wells"]
			name_rows = sorted(list(set([well[0] for well in labware_wells])))
			number_rows = len(name_rows)
			number_columns = len(set([well[1:] for well in labware_wells]))
			
			dataframes_source_with_index = {}
			
			for labware in labware_source:
				dataframes_source_with_index[int(str(labware).split(" ")[-1])] = pd.DataFrame(np.full((number_rows,number_columns),None),columns=list(range(1,number_columns+1)),index=name_rows)
				dataframes_source_with_index[int(str(labware).split(" ")[-1])].index.name = "Row/Column"
	
			return dataframes_source_with_index

# Body of the script
#--------------------------------
#--------------------------------
	try:
		current_step = "Reading csv and transforming them to parameters/variables"
		variables_csv = pd.read_csv("/data/user_storage/Variables-ColonieScreening-OT-new.csv", index_col = 0)
		# variables_csv = pd.read_csv("Variables-ColonieScreening-OT-new.csv", index_col = 0)
		variables = setted_parameters(variables_csv)
		
		# Due to the fact that even in the parameters checking we use the labware definitions, we are going to check that they exist outside that function
		try:
			labware_context.get_labware_definition(variables.name_source_plate)
			labware_context.get_labware_definition(variables.name_final_plate)
			labware_context.get_labware_definition(variables.name_rack_falcons)
		except:
			raise Exception("One or more of the introduced labwares are not in the custom labware directory of the opentrons. Check for any typo of the api labware name.")
		
		current_step = "Variables verification"
		errors_variables = check_setted_parameters(variables)
		if len(errors_variables) > 0:
			print("------------------------------------------------\nERROR(S):")
			for error in errors_variables:
				print("\t- "+error)
			print("------------------------------------------------\nBefore proceeding this error(s) should be fixed")
			quit()
		
		current_step = "Defining pipettes"
		# We can define them because in the parameters check we have seen if there is a pipette that it is not in the mount
		if variables.right_pipette_name == "None" or variables.right_pipette_name == "-":
			variables.right_pipette = None
		else:
			variables.right_pipette = protocol.load_instrument(variables.right_pipette_name, "right", tip_racks=[])
			
		if variables.left_pipette_name == "None" or variables.left_pipette_name == "-":
			variables.left_pipette = None
		else:
			variables.left_pipette = protocol.load_instrument(variables.left_pipette_name, "left", tip_racks=[])
		
		current_step = "Load OD data"
		# Loading the csv from Amp and Ant antibiotics (for example)
		data_transform = pd.read_csv(variables.file_transform, index_col = 0)
		data_integration = pd.read_csv(variables.file_integration, index_col = 0)
		
		# Lets check if the files have the same dimensions as the source plate
		size_source_plate = (len(labware_context.get_labware_definition(variables.name_source_plate)["ordering"][0]), len(labware_context.get_labware_definition(variables.name_source_plate)["ordering"]))#(rows, columns)
		if data_transform.shape != size_source_plate or data_integration.shape != size_source_plate:
			raise Exception("One or both files given with data to compate does not have the same dimensions of the source plate")
		
		current_step = "Selecting colonies that fulfill the conditions"
		# Select the colonies that meet our conditions of growth (selective conditions)
		name_columns = list(data_transform.columns.values)
		name_index = list(data_transform.index.values)
		
		list_colonies_pip = []

		for i in name_columns:
			for j in name_index:
				if (float(data_transform.loc[j][i]) > variables.grow_OD) and (float(data_integration.loc[j][i]) < variables.grow_OD):
					# One possible update is having different values of selection for the different files
					list_colonies_pip.append(j+i)
				else:
					pass
		
		total_number_colonies_pick = len(list_colonies_pip)
		# Now we set all the labware with the information that we have calculated or obtained
		
		labware_source, labware_falcons = setting_labware(variables, total_number_colonies_pick)
		
		# We create a list of all the wells where we have to transfer glycerol and water
		wells_glycerol = []
		for set_labware_glycerol in variables.reactives_info["glycerol"]["destination_plates"]:
			all_wells_glycerol = []
			for labware_glycerol in set_labware_glycerol:
				all_wells_glycerol += labware_glycerol.wells()
			wells_glycerol += all_wells_glycerol[variables.index_start_final_plate:variables.index_start_final_plate+total_number_colonies_pick]
		
		wells_water = []

		for set_labware_pcr in variables.reactives_info["water"]["destination_plates"]:
			all_wells_pcr = []
			for labware_final_pcr in set_labware_pcr:
				all_wells_pcr += labware_final_pcr.wells()

			wells_water += all_wells_pcr[variables.index_start_final_plate:variables.index_start_final_plate+total_number_colonies_pick]

		# We are going to transfer the reactives
		# We can have more than one tube and that is why we are going to rack the volume of the tubes and also the height of aspirate
		
		# Transfer glycerol
		pipette_use = give_me_optimal_pipette_to_use(variables.volume_transfer_glycerol, variables.right_pipette, variables.left_pipette)

		check_tip_and_pick(pipette_use, variables)

		for index_tube, reactions in enumerate(variables.reactives_info["glycerol"]["reactions_per_tube"]):
				current_tube = variables.reactives_info["glycerol"]["position_tubes"][index_tube]
				wells_to_distribute = wells_glycerol[:reactions]
				
				distribute_z_tracking(pipette_use, (len(wells_to_distribute)*variables.volume_transfer_glycerol)+10, variables.volume_transfer_glycerol, current_tube, wells_to_distribute)
				
				del wells_to_distribute[:reactions]
		pipette_use.drop_tip()
		
		# Transfer water
		pipette_use = give_me_optimal_pipette_to_use(variables.volume_transfer_water, variables.right_pipette, variables.left_pipette)
		check_tip_and_pick(pipette_use, variables)
		for index_tube, reactions in enumerate(variables.reactives_info["water"]["reactions_per_tube"]):
				current_tube = variables.reactives_info["water"]["position_tubes"][index_tube]
				wells_to_distribute = wells_water[:reactions]
				
				distribute_z_tracking(pipette_use, (len(wells_to_distribute)*variables.volume_transfer_water)+10, variables.volume_transfer_water, current_tube, wells_to_distribute)
				
				del wells_to_distribute[:reactions]
		pipette_use.drop_tip()
		
		
		# Create the map for the colonies
		variables.maps_selected_cells = pds_labware_sources(variables.reactives_info["glycerol"]["destination_plates"][0])
		
		# Create the list in which the selected colonies have to go

		final_wells_each_selected_colony = []
		all_wells_with_reactives = wells_glycerol+wells_water
		for index_cell_selected in range(total_number_colonies_pick):
			list_append_final_wells = []
			for number_repetition_plate in range(int(len(all_wells_with_reactives)/total_number_colonies_pick)):
				list_append_final_wells.append(all_wells_with_reactives[index_cell_selected+(number_repetition_plate*total_number_colonies_pick)])
			final_wells_each_selected_colony.append(list_append_final_wells)
		# Lets move and fill the map
		# Set the position to move
		position_dispense = set_labware(final_wells_each_selected_colony)
		
		# We select the pipette that we are going to use
		suitablepip = give_me_optimal_pipette_to_use(variables.volume_transfer_colony, variables.right_pipette, variables.left_pipette)
		
		# Iterator to set the cell to which we are going to dispense in the different plates
		# Iterate over the colonies that we select
		for selected_colonie in list_colonies_pip:
			check_tip_and_pick(suitablepip, variables)
			
			# Final well position
			final_colonie_positions = next(position_dispense)
			# Transfer to all labwares (pcr plates and glycerol stocks)
			suitablepip.distribute(variables.volume_transfer_colony, labware_source[0][selected_colonie].bottom(z=0.7), final_colonie_positions, disposal_volume = 0, new_tip = "never")
			suitablepip.drop_tip()
			
			# Annotate in the map
			well_annotate_map_deck_position = str(final_colonie_positions[0]).split(" ")[-1]
			well_annotate_map_well_position = str(final_colonie_positions[0]).split(" ")[0]

			variables.maps_selected_cells[int(well_annotate_map_deck_position)][int(well_annotate_map_well_position[1:])][well_annotate_map_well_position[0]] = selected_colonie
			
		
		# We export the map (s)
		number_map = 1
		for labware_selected_colonies in variables.maps_selected_cells.values():
			labware_selected_colonies.index.name = "Row/Column"
			labware_selected_colonies.fillna("-",inplace=True)
			labware_selected_colonies.to_csv(f"{variables.final_map_name}_{number_map}.csv")
			number_map += 1
			
		# Home the robot
		current_step ="Homing robot"
		protocol.home()
		
		current_step = "Printing User Information"
		print_information_user(variables, labware_falcons, list_colonies_pip)
		
	except Exception as e:
		print("-------------------------------------------------------------------------------")
		print("--> ERROR STEP:\n" + current_step)
		logger.critical(e, exc_info = True)
		print("-------------------------------------------------------------------------------")
		
		
