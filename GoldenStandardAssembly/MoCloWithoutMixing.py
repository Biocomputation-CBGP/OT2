import opentrons.execute
from opentrons import protocol_api
import pandas as pd
import numpy as np
import logging
import copy
import csv
import os
import math

## Structure needed for the running of th escript in the OT
metadata = {
'apiLevel':'2.12'
}

def run(protocol: protocol_api.ProtocolContext):


#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
## FUNCIONES JUNTAS

	def generate_plate_maps(filename):
		"""
		Function that takes a file name, a csv file, and creates a dictionary with the
		name of the file and rows as list in which element is a cell of that csv row.
		
		Used in this script to create an easy way to access the map of DNA fragments when creating the MoClo fragments
		"""
		plate_maps = {}
		plate_map = []
		with open(filename, "r") as file:
			for row in csv.reader(file, dialect= "excel"):
				if len(row) == 0:
					continue
				if row[0]:
					if '\ufeff' in row[0]:
						row[0] = str(row[0].replace(u'\ufeff',''))
					if "ï»¿" in row[0]:
						row[0] = str(row[0].replace('ï»¿',''))
					plate_map.append(row)
		plate_name = os.path.splitext(os.path.basename(filename))[0]
		plate_maps[plate_name] = plate_map
	
		return plate_maps
	
	def generate_combinations(combinations_filename):
		"""
		Function that takes the name of a csv file in which each row is based in a name and the combination of fragments.
		
		In this script this is being use to create a dictionary with the combinations.
		
		It returns a dictionary with 2 elements, "name" (the first column of the csv) and "parts" with as elements as cells in the csv where filled
		"""
		combinations_to_make = []
		with open(combinations_filename, "r") as f:
			for row in csv.reader(f, dialect='excel'):
				if len(row) == 0:
					continue
				if row[0]:
					if '\ufeff' in row[0]:
						row[0] = str(row[0].replace(u'\ufeff',''))
					combinations_to_make.append({
												"name": row[0],
												"parts": [x for x in row[1:] if x]
												})
		return combinations_to_make
	
	def generate_and_save_output_plate_maps(combinations_to_make, output_filename, labware_name):
		"""
		This function digresses a little bit from the original one. We changed it so it can be more free to have any type of labware, not only agar plates
		
		What it takes is the dict of combinations and create as much plates as it is needed and record the places of each combination on them.
		After doing this, that map is exported as a csv
		"""
		
		# Split combinations_to_make into plate maps
		number_rows = len(labware_context.get_labware_definition(labware_name)["ordering"][0])
		number_cols = len(labware_context.get_labware_definition(labware_name)["ordering"])
		number_wells = number_rows*number_cols
		output_plate_map_flipped = []
		number_wells_in_current_plate = 0
		position_final_map = 1 # lo vamos a inicializar asi pero esto tambine lo podemos sacar de opentrons, lo vamos a hacer manual por ahora como prueba
							   # Hay que sacarlo de OT
		
		for i, combo in enumerate(combinations_to_make):
			name = combo["name"]
			if number_wells_in_current_plate % number_rows == 0 and number_wells_in_current_plate < number_wells: # In the current map but need another column
				# new column
				output_plate_map_flipped.append([name])
				number_wells_in_current_plate += 1
				
			elif number_wells_in_current_plate % number_rows != 0 and number_wells_in_current_plate < number_wells: # In the current map but we do not need another column
				output_plate_map_flipped[-1].append(name)
				number_wells_in_current_plate += 1
			
			else: # We need to Initialize another map
				# For that we just create this map and export it
				output_plate_map = []
				for i, row in enumerate(output_plate_map_flipped):
					for j, element in enumerate(row):
						if j >= len(output_plate_map):
							output_plate_map.append([element])
						else:
							output_plate_map[j].append(element)
			
				with open(output_filename+str(position_final_map)+".csv", 'w+', newline='') as f:
					writer = csv.writer(f)
					for row in output_plate_map:
						writer.writerow(row)
				f.close()
				
				position_final_map += 1
				number_wells_in_current_plate = 0
				output_plate_map_flipped = []
			
			# We export the map that is currently in place
			output_plate_map = []
			for i, row in enumerate(output_plate_map_flipped):
				for j, element in enumerate(row):
					if j >= len(output_plate_map):
						output_plate_map.append([element])
					else:
						output_plate_map[j].append(element)
		
			with open(output_filename+str(position_final_map)+".csv", 'w+', newline='') as f:
				writer = csv.writer(f)
				for row in output_plate_map:
					writer.writerow(row)
			f.close()
		
		return
	
	#--------
	
	def number_tubes_needed (vol_reactive_per_reaction_factor, number_reactions, vol_max_tube):
			"""
			Given a maximum volume of the tube (vol_max_tube), the volume of that reactive/reaction (vol_reactive_per_reaction_factor) and the total number of reactions (number_reactions)
			this function will return the number of tubes needed for this reactive and how many reactions are filled in every tube
			
			This function does not garantee the lower number of tubes but it assures that everything can be picked with the pipettes that we have
			
			This function will be used mainly, sometimes exclusively, in setting_number_plates
			"""
			# print("Number of tubes function")
			# print(vol_reactive_per_reaction_factor)
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
				# print(volumes_tubes)
			# When the volume can fit every tube (exit from th ewhile loop) we return the number of tubes and the reactions that will fit in every tube
			return (number_tubes, reactions_per_tube, volumes_tubes)
	
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
				if pipette_l == max_pipette and aVolume != pipette_r.max_volume:
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
	
	def transfer_volume_tracking(vol_distribute, position_source, position_final, volumes_source, variables):
		"""
		ESTA FUNCION HAY QUE CAMBIARLA, PORQUE EN REALIDAD DEBERIA ESTAR CON LAS COSAS DE NUMERO DE REACCIONES, NO CON VOLUMEN, PORQUE EL VOLUMEN YA LO HEMOS CALCULADO
		ANTES EN OTRA FUNCION, lo vamos a dejar asi por ahora, pero habra que cambiarlo, y tambien en la de PCR probablemente !!!!
		"""
		
		"""
		Function used to transfer bigger volumes than the distributing with volume tracking but this does not take in account the heights and the distributing takes in account
		that all the final locations (wells) will have the same quantity (distributing along the plate, for example) while this is used mainly to create the mix that will be distributed
		
		In a future transfer and distributing function will be unified
		"""
		# print("++++++++++++++")
		# print(position_source)
		# print(position_final)
		# This function is used to transfer bigger volumes (in comparasion to distribute)
		# The distribute and transfering functions are very similar, sometimes they can be interchanged (this is not the case for versions < 11 of the script)
		
		# Initialize different variables
		current_position_source = position_source[0]
		current_volume = volumes_source[0]
		i_position = 0
		pipette_used = give_me_optimal_pipette_to_use(vol_distribute[0], variables.right_pipette, variables.left_pipette)
		#pipette_used = give_me_optimal_pipette_to_use(vol_distribute, variables.right_pipette, variables.left_pipette)
		
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
					
					pipette_used.transfer(vol_distribute_position, current_position_source, position_final_current, new_tip = "never")
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
		
		# print("++++++++++++++")
		#Need to return the last pipette used because it could not be the same as the first one and usually after transfering we need to drop the tip
		return pipette_used

	class setted_parameters():
		"""
		Class that will contain the parameters setted in the variables csv and will process them to work easily in the rest of the protocol
		The coding of this function is dependant of the variables in the Template of the protocol and the names have to be consistent with the rest of the code
		"""
		
		def __init__(self, modular_parameters):
			self.name_file_source_labware = modular_parameters["Name File Source Reactives Maps"]["Value"]
			self.name_file_combinations = modular_parameters["Name Construct Combination File"]["Value"]
			self.name_final_maps = modular_parameters["Name Final Construct Maps"]["Value"]
			
			self.volume_acceptor = float(modular_parameters["Volume Acceptor Plasmid (uL)"]["Value"])
			self.volume_module_plasmid = float(modular_parameters["Volume Module Plasmid (uL)"]["Value"])
			self.volume_restriction_enzyme = float(modular_parameters["Volume Restriction Enzyme (uL)"]["Value"])
			self.volume_ligase = float(modular_parameters["Volume Ligase (uL)"]["Value"])
			self.volume_buffer = float(modular_parameters["Volume Buffer (uL)"]["Value"])
			self.volume_serum = float(modular_parameters["Volume ATP/Serum (uL)"]["Value"])
			self.volume_final_reaction = float(modular_parameters["Volume Final Each Reaction (uL)"]["Value"])
			
			self.name_coldblock = modular_parameters["Name Coldblock"]["Value"]
			self.name_source_plate = modular_parameters["Name Source Labware"]["Value"]
			self.name_final_plate = modular_parameters["Name Final Labware"]["Value"]
			self.index_start_final_plate = modular_parameters["Index of start cell in final labware"]["Value"]
			if self.index_start_final_plate == "-" or str(self.index_start_final_plate).lower() == "none":
				self.index_start_final_plate = 0
			else:
				self.index_start_final_plate = int(self.index_start_final_plate)
				
			self.name_right_pip = modular_parameters["Pipette Right Mount"]["Value"]
			self.name_left_pip = modular_parameters["Pipette Left Mount"]["Value"]
			self.starting_tip_right_pip = modular_parameters["Starting Tip Pipette Right"]["Value"]
			self.starting_tip_left_pip = modular_parameters["Starting Tip Pipette Left"]["Value"]
			self.extra_pip_factor = float(modular_parameters["Extra Pipetting Adjust Factor"]["Value"])+1.0
			if modular_parameters["Replace Tiprack"]["Value"].lower() == "false":
				self.replace_tiprack = False
			elif modular_parameters["Replace Tiprack"]["Value"].lower() == "true":
				self.replace_tiprack = True
			else:
				self.replace_tiprack = False
			
			# Inicializamos primero
			self.need_change_tiprack = False
			self.number_reactions = 0
			self.right_pipette = None
			self.left_pipette = None
			self.number_source_plates = 0
			self.number_final_plates = 0
			self.maps_reactives = {} # Each map is going to have his spot in the dictionary, maybe it will turn into a list, we will see because we dont go by names of the values
			self.mix_reactives = 0
			# Algunas variables necesarias para que se pueda printear la info y usarla para otras cosas del script que hay qu einicializar tambien
			self.positions_labware = {"final labware":[], "source labware":[], "cold block":[]}
			self.positions_reactives = {"restriction enzyme":{"positions":[],"reactions":[],"volume":[]},"ligase":{"positions":[],"reactions":[],"volume":[]},"buffer":{"positions":[],"reactions":[],"volume":[]},"serum":{"positions":[],"reactions":[],"volume":[]},"mix":{"positions":[],"reactions":[],"volume":[]}}
			self.map_dna_input = None
			self.water_tubes = {"volumes tubes":None, "volumes transfer from tube":None, "positions":[]}
			return

#------------------------------------------------
#------------------------------------------------

	def check_setted_parameters(variables):
		"""
		Function that will check the variables of the Template and will raise errors that will crash the OT run
		It is a validation function of the variables checking errors or inconsistencies
	
		This function is dependant again with the variabels that we have, some checks are interchangable between protocols, but some of them are specific of the variables
		"""
	
		errors = []
			
		# Check if some of the common reactives can be distributed with some pipette
		try:
			give_me_optimal_pipette_to_use(variables.volume_restriction_enzyme + variables.volume_ligase + variables.volume_buffer + variables.volume_serum, variables.right_pipette, variables.left_pipette)
		except:
			errors.append("Reactive mix volume cannot be picked by any of the set pipettes")
		
		try:
			give_me_optimal_pipette_to_use(variables.volume_module_plasmid, variables.right_pipette, variables.left_pipette)
			give_me_optimal_pipette_to_use(variables.volume_acceptor, variables.right_pipette, variables.left_pipette)
		except:
			errors.append("Either the volume of the acceptor or the volume of the module cannot be picked by set pipettes")
		
		# Index of any of the source plates is larger than the wells of that plate
		if variables.index_start_final_plate > len(labware_context.get_labware_definition(variables.name_final_plate)["wells"]):
			errors.append("Index of start for the final plate is larger than the number of wells of that labware")
		number_wells_source_plate = len(labware_context.get_labware_definition(variables.name_source_plate)["wells"])
		
		# Check if there is some typo in the the starting tip
		if variables.name_right_pip not in ["None", "-"]:
			if not any(variables.starting_tip_right_pip in columns for columns in labware_context.get_labware_definition(define_tiprack(variables.name_right_pip))["ordering"]):
				errors.append("Starting tip of right pipette is not valid, check for typos")
		if variables.name_left_pip not in ["None", "-"]:
			if not any(variables.starting_tip_left_pip in columns for columns in labware_context.get_labware_definition(define_tiprack(variables.name_left_pip))["ordering"]):
				errors.append("Starting tip of left pipette is not valid, check for typos")
		
		# We will check that the DNA map that the script has been provided to the script actually fits in the source labware
		labware_columns = len(labware_context.get_labware_definition(variables.name_source_plate)["ordering"])
		labware_rows = len(labware_context.get_labware_definition(variables.name_source_plate)["ordering"][0])
		for name, map_current in variables.map_dna_input.items():
			number_rows = len(map_current)
			number_columns = len(map_current[0])
			if labware_rows < number_rows or labware_columns < number_columns:
				errors.append("The map "+name+" does not fit in the labware with API name "+variables.name_source_plate)
			else:
				pass
			
		# We are going to check now if there is some repeated element in the map(s) of dna input
		# first we concatenate all elements
		all_dna_modules = []
		for name_map in variables.map_dna_input:
			for row in variables.map_dna_input[name_map]:
				all_dna_modules += row
		# second we remove the "", " " and "-" elements
		while "" in all_dna_modules: all_dna_modules.remove("")
		while " " in all_dna_modules: all_dna_modules.remove(" ")
		while "-" in all_dna_modules: all_dna_modules.remove("-")
		# Finally we see if there is some element that is repeated
		if len(all_dna_modules) != len(set(all_dna_modules)):
			errors.append("One or more elements of the DNA input map is repeated. For now this program only can take from one source :(")
		
		# We are going to check if a combination name is repeated
		all_combination_names = []
		for combination in variables.maps_reactives:
			all_combination_names.append(combination["name"])
		if len(all_combination_names) != len(set(all_combination_names)):
			errors.append("One or more names in the combinations file is repeated. Make sure that the combinations have different names so you can identify them in the final map")
		
		return errors
#------------------------------------------------
#------------------------------------------------

	#-----------
	
	def distribute_vol_tracking(volume_distribute, positions_tube_source, reactions_tubes_source, positions_final, variables):
		"""
		This function is destined to transfer small amounts (in comparasion to transfer_vol_tracking) from a source to different wells
		
		Specifically, here we used it to distribute the set mix to the wells and try to do it in the lower ammount of movements
		
		This function works because we have assigned previously the reactions to the tube, if we separate the reactives another way we
		should use other function
		
		This function of distribution is different from the other protocols because in this case the height is not an issue and the vol tracking is
		done in the volume calculating section
	
		"""
		# print("++++++++")
		# print(volume_distribute)
		# print(positions_tube_source)
		# print(reactions_tubes_source)
		# print(positions_final)
		# print("++++++++")		
		# print("entramos por primera vez al distribute")
		# Let's define the distributing pipette that is not going to change
		pipette_distributing = give_me_optimal_pipette_to_use(volume_distribute, variables.right_pipette, variables.left_pipette)
		# print("La pipeta para distribuir es "+str(pipette_distributing.name))
		
		# Lets distribute the reactive
		for index_tube, reactions_tube in enumerate(reactions_tubes_source):
			# print("este es el index tube donde estamos "+str(index_tube)+" con todas estas reacciones en el tubo "+str(reactions_tube))
			# Define the mixing pipette
			volume_mixing_tube = volume_distribute*reactions_tube/3
			pipette_mixing = give_me_optimal_pipette_to_use(volume_mixing_tube, variables.right_pipette, variables.left_pipette)
			# print("este es la pipeta de mixing "+str(pipette_mixing))
			if volume_mixing_tube > pipette_mixing.max_volume:
				volume_mixing_tube = pipette_mixing.max_volume
			else:
				pass
			
			if pipette_distributing == pipette_mixing:
				# print("las pipetas son las mismas")
				if pipette_mixing.has_tip == True:
					pass
				else:
					check_tip_and_pick(pipette_mixing, variables)
				# Mixing eppendorf ya cogera el tip y se supone que la pipeta sera la misma que la que hemos calculado previamente
				### mixing_eppendorf_15(positions_tube_source[index_tube], pipette_mixing.max_volume, pipette_mixing, variables)
				### pipette_mixing.blow_out(positions_tube_source[index_tube].bottom())
				protocol.pause("Stir mix tube so liquid is homogeneus!")
				for position_well_final in positions_final[:reactions_tube]:
					pipette_distributing.transfer(volume_distribute, positions_tube_source[index_tube], position_well_final, new_tip = "never")
					pipette_distributing.blow_out(position_well_final.bottom())
				del positions_final[:reactions_tube]
			else:
				# print("las pipetas no son iguales")
				if pipette_distributing.has_tip == True:
					# print("hay que tirar el tip de la pipeta de distribuir")
					pipette_distributing.drop_tip()
				# Mixing eppendorf ya cogera el tip y se supone que la pipeta sera la misma que la que hemos calculado previamente
				### check_tip_and_pick(pipette_mixing, variables)
				### mixing_eppendorf_15(positions_tube_source[index_tube], pipette_mixing.max_volume, pipette_mixing, variables)
				### pipette_mixing.blow_out(positions_tube_source[index_tube].bottom())
				# Sin embargo, no tira el tip, por lo que hay que tirarlo nosotros
				### pipette_mixing.drop_tip()
				protocol.pause("Stir mix tube so liquid is homogeneus!")
				check_tip_and_pick(pipette_distributing, variables)
				for position_well_final in positions_final[:reactions_tube]:
					pipette_distributing.transfer(volume_distribute, positions_tube_source[index_tube], position_well_final, new_tip = "never")
					pipette_distributing.blow_out(position_well_final.bottom())
				del positions_final[:reactions_tube]
		
		if pipette_distributing.has_tip == True:
			pipette_distributing.drop_tip()
		if pipette_mixing.has_tip == True:
			pipette_mixing.drop_tip()
		
		return
	
	
	# ------------------------ Transfering plasmids to wells for the construct ------------------------
	"""
	Esta es la parte tricky del script, o la mas tricky y la mas diferente de las anteriormente programadas.
	
	Para ello vamos a usar parte de sus funciones y su script con un poco del nuestro.
	
	Esta va a ser la ultima parte de este protocolo. Vamos a intentar hacer esto primero bien que es la creacion de las placas.
	Por si solo esto puede ser un protocolo, por lo que se puede probar, pero estaria bien hacerlo para muchos y tambien para 1
	para el termociclador. Si hacemos estos 2, que uno seria simplemente lo mismo pero con restriccion de que solo se puede hacer
	1 placa y que habria que introducirle los protocolos del termociclador (lo de los ciclos, las tempreaturas y los tiempos)
	"""
	
	
	def find_dna(name, dna_plate_map_dict, dna_plate_dict, number_rows_plate):
		"""Return a well containing the named DNA
		El dan_plate_map dict es el de los nombres y el dna_plate_dict es el labware, creo"""
		for plate_name, plate_map in dna_plate_map_dict.items():
			for i, row in enumerate(plate_map):
				for j, dna_name in enumerate(row):
					if dna_name == name:
						well_num = number_rows_plate * j + i # Ese 8 habria que cambiarlo por el numero de rows que tenga el labware, pero es minimo el cambio
						return(dna_plate_dict[well_num])
		raise ValueError("Could not find dna piece named \"{0}\"".format(name))
	
	#This function checks if the DNA parts exist in the DNA plates and returns for well locaion of output DNA combinations
	def find_combination(name, combinations_to_make, reaction_plate):
		"""Return a well containing the named combination."""
		for i, combination in enumerate(combinations_to_make):
			if combination["name"] == name:
				return reaction_plate[i] # Este reaction_plate.wells(i) es el i de los que hemos antes distribuido, por lo que la lista la tenemos porque lo hemos hecho el mapa
										 # por columns y opentrons va de esa manera, haciendo una n inversa (una i cirilica) por lo que va bien eso
		raise ValueError("Could not find combination \"{0}\".".format(name))
	
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
	
	def setting_reactives_in_places(number_tubes_reactives, name_reactive, possible_positions):
		"""
		Function that will set tubes of reactives within a type of labware, in the case of PCR sample preparation is setting eppendorfs
		within one or more coldblocks
		
		Because we have already established how many labwares we needed, we can go through this function without having a possible error
		
		This function will be mainly or exclusivelly used in the setting_labware function
		
		We are going to obtain a dictionary with the name of the reactive and its positions
		"""
		# print("+++++++")
		# print(number_tubes_reactives)
		# print(name_reactive)
		# print(len(possible_positions))
		# print("+++++++")
		dict_position_reactive = {}
		for tube in range(int(number_tubes_reactives)):
			try:
				# It is not the first time the reactive is being assigned to a well so we are adding it to the list and not assigning it
				dict_position_reactive[name_reactive] = dict_position_reactive[name_reactive] + [possible_positions[tube]]
			except:
				dict_position_reactive[name_reactive] = [possible_positions[tube]]
		
		return (dict_position_reactive)

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
	
	def mixing_eppendorf_15(location_mixing, volume_mixing, aSuitablePippet_mixing, variables):
		"""
		This is a function to perform an extensive mixing of every eppendorf which should be done before distributing the reactives along the final plates
		
		Mixing is one of the most crucial parts of this workflow and that is why theer is a function only for it
		
		This function will perform the mixing and will return warnings (in case that is needed) and the final pipette that has been used
		"""
	
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
		
		return

	def print_information_user(variables):
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
		print("\t- Name of file where map of DNA parts are: "+str(variables.name_file_source_labware))
		print("\t- Name of file where name|acceptor|DNA parts combinations are: "+str(variables.name_file_combinations))
		print("\t- Name of final plasmid map file: "+str(variables.name_final_maps))
		print("\t- Number of total combinations: "+str(variables.number_reactions))
		print("\t- Number of output labware(s): "+str(len(variables.positions_labware["final labware"])))
		print("\t- Final volume in output labware (uL): "+str(variables.volume_final_reaction))
		if variables.replace_tiprack == True and variables.need_change_tiprack == True:
			print("\t- You should stay near the OT to change the tiprack(s)\n")
		elif variables.replace_tiprack == True and variables.need_change_tiprack == False:
			print("\t- There will be no need to change tipracks\n")
		elif variables.replace_tiprack == False:
			print("\t- You have the replace tiprack option as False so there is no need of replacing none of the tipracks")
		
		# Deck information of the labware
		print("\n--------------------------------------------------------------\nDECK LABWARE POSITIONS\nThis information will be provided by the OT App as well")
		for labware_position in protocol.loaded_labwares.values():
			print("\t"+str(labware_position))
		
		# Coldblock positions for the reactives
		coldblock_positions = variables.positions_labware["cold block"]
		# Printing of the reactives and their volume as it will bein the coldblock(s)
		print("\n--------------------------------------------------------------\nCOLDBLOCK REACTIVES POSITIONS\nReactives in each position of the coldblock(s) and their respective volume (minimum)")
		for coldblock in coldblock_positions:
			coldblock_table = pd.DataFrame(index=["A","B","C"],columns=["1","2","3","4"])
			print("\n--> Coldblock in Slot "+str(coldblock).split(" ")[-1])			
			for name_reactive, info_reactive in variables.positions_reactives.items():
				for index_tube, tube_position in enumerate(info_reactive["positions"]):
					reactive_well = str(tube_position).split(" ")[0]
					if str(tube_position).split(" ")[-1] == str(coldblock).split(" ")[-1]:
						coldblock_table[reactive_well[1:]][reactive_well[0]] = name_reactive.upper()+" - "+str(info_reactive["volume"][index_tube])+str("uL")
			print(coldblock_table.fillna("-")) # for esthetic purposes
			print("\n--------------------------------------------------------------")
		return

	def set_labware(labware_wells_name):
		"""
		Generator of the positions to transfer, it will give you the next element of a list every time it is called
		"""
		for well in labware_wells_name:
			yield well
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
# BODY DEL SCRIPT

	#------------------------- Upload the variables and check them ------------------------
	labware_context = opentrons.protocol_api.labware
	
	current_step = "Reading csv and transforming them to parameters/variables"
	variables_df = pd.read_csv("/data/user_storage/230110_Variables-MoClo-testBiologicalOT.csv", index_col = 0)

	# We are going to convert these parameters into arguments of the class variables and we ar egoing to process some of them so they can be usable (they are going to be dictionaries)
	variables = setted_parameters(variables_df.to_dict(orient = "index"))
	
	# Due to the fact that even in the parameters checking we use the labwares definition, we are going to check that they exist outside that function
	try:
		labware_context.get_labware_definition(variables.name_source_plate)
		labware_context.get_labware_definition(variables.name_final_plate)
		labware_context.get_labware_definition(variables.name_coldblock)
	except:
		raise Exception("One or more of the introduced labwares are not in the custom labware directory of the opentrons. Check for any typo of the API labware name.")
	
	# Check if the provided files exist
	try:
		open(variables.name_file_source_labware, "r")
		open(variables.name_file_combinations, "r")
	except:
		raise Exception("Some of the provided files, either map of modules or combinations, cannot be found")
	
	# Check and set the pipettes, we do it at the same time
	if variables.name_right_pip == "None" or variables.name_right_pip == "-":
		variables.right_pipette = None
	else:
		try:
			variables.right_pipette = protocol.load_instrument(variables.name_right_pip, "right", tip_racks=[])
		except:
			raise Exception("Right pipette cannot be loaded, check for typos in the name, it has to be the API name")
		
	if variables.name_left_pip == "None" or variables.name_left_pip == "-":
		variables.left_pipette = None
	else:
		try:
			variables.left_pipette = protocol.load_instrument(variables.name_left_pip, "left", tip_racks=[])
		except:
			raise Exception("Left pipette cannot be loaded, check for typos in the name, it has to be the API name")
	
	# ------------------------ Creation of Maps and the maps we are going to use ------------------------
	"""
	En esta parte vamos a usar las funciones de DAMP Lab, o al menos una modificacion de ellas
	para generar los diccionarios que usaremos mas tarde para poder transferir los fragmentos correctos a los
	wells finales correctos
	Ademas, en esta seccion tambien se exportan los mapas finales de los labwares, asi no tengo qu epensar como hacerlo
	yo, uso sus fucniones, aunque esa justamente ha tenido que ser modificada para generalizarla un poco a cualquier labware
	"""

	
	# First we create the map of input to play with it and the dictionary of combinations
	variables.map_dna_input = generate_plate_maps(variables.name_file_source_labware)
	variables.maps_reactives = generate_combinations(variables.name_file_combinations)
	
	# After setting the pipettes, the maps and the labware check out we do the checking of the rest of the parameters
	errors_variables = check_setted_parameters(variables)
	if len(errors_variables) > 0:
		print("------------------------------------------------\nERROR(S):")
		for error in errors_variables:
			print("\t- "+error)
		print("------------------------------------------------\nBefore proceeding this error(s) should be fixed")
		quit()
	# Now that we have checked the parameters and we know that it will, in theory, work, we do the final map
	# Now we create the map of final combinations (with the name of combinations in the cells)
	generate_and_save_output_plate_maps(variables.maps_reactives, variables.name_final_maps, variables.name_final_plate)
	
	# ------------------------ Creation of MoClo Reactives Mix------------------------
	"""
	En esta parte lo que vamos a crear es lo necesario del mix para crear todos las combinaciones del MoClo
	Vamos, que lo que vamos a crear es un mix en eppendorfs, tipo la PCR que hacemos,
	un mix teorico (qu e alo mejor luego habra que hacer a mano)
	
	Esta vamos a intentar usar el de eppendorfs de 24 (el labware) con nuestras funciones de la PCR,
	en el caso de que algo de esto cambie, por ejemplo, que tenga que ser en el cold block, eso sera facil de cambiar
	ya que lo que cambiara sera el numero de labwares que se necesitaran, algo que ya tenemos programado,
	y tambien cambiar lo del mix (el mix teorico ese) por un pause o un delay con un comentario.
	Al final no vamos a hacer aqui ni el mix teorico aunque obviamente en el real si haria falta un mix.
	
	Vamos a hacerlo para el Short protocol in restriction buffer primero, para enfocar un poco el problema,
	luego generalizaremos ya porque hay 4 casos que se pueden hacer: short restriction Level 0, short restriction Level 1,
	el long ligase Level 0 y el long ligase Level 1
	"""
	
	
	variables.number_reactions = len(variables.maps_reactives)
	# =============== ESTA PARTE SERA MUY VARIABLE PORQUE NO SE EXACTAMENTE EL VOLUMEN DE CADA COSA ===============
	# reactives for the short protocol in restriction buffer for Level 0
	variables.mix_reactives = variables.volume_restriction_enzyme + variables.volume_buffer + variables.volume_ligase + variables.volume_serum
	
	vol_max_tube = labware_context.get_labware_definition(variables.name_coldblock)["wells"]["A1"]["totalLiquidVolume"]
	# ===============
	
	# tubes_restriction_enzyme = number_tubes_needed (variables.volume_restriction_enzyme, number_reactions, vol_max_tube)
	# tubes_buffer = number_tubes_needed (variables.volume_buffer, number_reactions, vol_max_tube)
	# tubes_ligase = number_tubes_needed (variables.volume_ligase, number_reactions, vol_max_tube)
	# tubes_atp = number_tubes_needed (variables.volume_serum, number_reactions, vol_max_tube)
	# tubes_mix = number_tubes_needed (variables.mix_reactives, number_reactions, vol_max_tube)
	number_tubes_enzyme, variables.positions_reactives["restriction enzyme"]["reactions"], variables.positions_reactives["restriction enzyme"]["volume"] = number_tubes_needed (variables.volume_restriction_enzyme*variables.extra_pip_factor, variables.number_reactions, vol_max_tube)
	number_tubes_buffer, variables.positions_reactives["buffer"]["reactions"], variables.positions_reactives["buffer"]["volume"] = number_tubes_needed (variables.volume_buffer*variables.extra_pip_factor, variables.number_reactions, vol_max_tube)
	number_tubes_ligase, variables.positions_reactives["ligase"]["reactions"], variables.positions_reactives["ligase"]["volume"] = number_tubes_needed (variables.volume_ligase*variables.extra_pip_factor, variables.number_reactions, vol_max_tube)
	number_tubes_serum, variables.positions_reactives["serum"]["reactions"], variables.positions_reactives["serum"]["volume"] = number_tubes_needed (variables.volume_serum*variables.extra_pip_factor, variables.number_reactions, vol_max_tube)
	# print(variables.mix_reactives*variables.extra_pip_factor)
	number_tubes_mix, variables.positions_reactives["mix"]["reactions"], variables.positions_reactives["mix"]["volume"] = number_tubes_needed (variables.mix_reactives*variables.extra_pip_factor, variables.number_reactions, vol_max_tube)
	# print("\n++++++++++++++++++++++++++++++")
	# print(number_tubes_mix)
	# print(variables.positions_reactives["mix"]["reactions"])
	# print(variables.positions_reactives["mix"]["volume"])
	
	
	##++++++++++++++++++++++++++++++++++++++++++++++++
	##++++++++++++++++++++++++++++++++++++++++++++++++
	# Aqui calculamos cunato hay que pasar a cada well y los volumenes del tubo o tubos
	
	volume_every_well = []
	for combination in variables.maps_reactives:
		volume_every_well.append(variables.mix_reactives+variables.volume_acceptor+variables.volume_module_plasmid*(len(combination["parts"])-1))
	# volume_every_well = [5, 10, 15, 20, 0, 18, 3] # Esto tambien habra que calcularlo en algun momento, pero, por ahora pongamoslo mas simple
	
	# lets calculate the volume of water in each well
	volume_water_every_well = [variables.volume_final_reaction-i for i in volume_every_well]
	# lets calculate the water needed
	all_water_needed = sum(volume_water_every_well)
	
	# print(volume_water_every_well)
	
	if any(volume+10 > vol_max_tube for volume in volume_water_every_well) == True:
		raise Exception("One of the volumes of water does not fit in tubes, check combinations, maxime volume of tubes or final reaction volume")
	
	if sum(volume_water_every_well) <= vol_max_tube:
		all_tubes = [volume_water_every_well]
		
	else:
		current_tube = []
		all_tubes = []
		for i, element in enumerate(volume_water_every_well):
			if sum(current_tube)+element+10 <= vol_max_tube:
				current_tube.append(element)
			else:
				all_tubes.append(current_tube)
				current_tube = [element]
			
			if i == len(volume_water_every_well)-1:
				all_tubes.append(current_tube)
	
	# Now we will come with volume and number of tubes
	number_tubes_water = len(all_tubes)
	variables.water_tubes["volumes tubes"] = [sum(volumes_tube) for volumes_tube in all_tubes]
	variables.water_tubes["volumes transfer from tube"] = all_tubes
	
	##++++++++++++++++++++++++++++++++++++++++++++++++
	##++++++++++++++++++++++++++++++++++++++++++++++++
	
	
	total_tubes = number_tubes_enzyme + number_tubes_buffer + number_tubes_ligase + number_tubes_serum + number_tubes_mix + number_tubes_water
	
	# print("Number reactions: "+str(number_reactions))
	# print("Volume needed of reactives:")
	# print("Enzyme: "+str(variables.volume_restriction_enzyme*number_reactions))
	# print("Buffer: "+str(variables.volume_buffer*number_reactions))
	# print("Ligasa: "+str(variables.volume_ligase*number_reactions))
	# print("Serum/ATP: "+str(variables.volume_serum*number_reactions))
	# print("Mix: "+str(mix_reactives*number_reactions))
	
	#==============================================================
	"""
	Aqui vamos a setear los labwares en su posicion despues de haber calculado ya todos los tubos que necesitamos de reactivos
	Primero vamos a calcular cuantos final labwares necesitamos y asi simplemente seteamos todos los labwares necesarios
	Esta parte de setear los labwares necesitaran una modificacion cuando se generalice porque solo tenemos por ahora un source labware, pero eso cambiara
	"""
	# Set the position of source labwares
	variables.positions_labware["source labware"] = setting_number_plates(1, variables.name_source_plate)
	
	# Calculate how many final labwares do we need
	number_wells_final_plates = len(labware_context.get_labware_definition(variables.name_final_plate)["wells"])
	variables.positions_labware["final labware"] = setting_number_plates(math.ceil(variables.number_reactions/number_wells_final_plates), variables.name_final_plate)
	
	# First we need to know how many recatives labwares do we need
	number_wells_reactives_labware = len(labware_context.get_labware_definition(variables.name_coldblock)["wells"])
	plates_reactives_needed = math.ceil(total_tubes/number_wells_reactives_labware)
	variables.positions_labware["cold block"] = setting_number_plates(plates_reactives_needed, variables.name_coldblock)
	
	# Now we are going to generate the dictionary in which we are going to have the positions of the reactives
	positions_coldblock_free = []
	for coldblock in variables.positions_labware["cold block"]:
		positions_coldblock_free += coldblock.wells()
	
	# Vamos a hacerlo de esta manera por intentar no cambiar tanto la fuincion que ya teniamos del otro script
	variables.positions_reactives["restriction enzyme"]["positions"] = list(setting_reactives_in_places(number_tubes_enzyme, "enzyme", positions_coldblock_free).values())[0]
	positions_coldblock_free = positions_coldblock_free[number_tubes_enzyme:]
	
	variables.positions_reactives["buffer"]["positions"] = list(setting_reactives_in_places(number_tubes_buffer, "buffer", positions_coldblock_free).values())[0]
	positions_coldblock_free = positions_coldblock_free[number_tubes_buffer:]
	
	variables.positions_reactives["ligase"]["positions"] = list(setting_reactives_in_places(number_tubes_ligase, "ligase", positions_coldblock_free).values())[0]
	positions_coldblock_free = positions_coldblock_free[number_tubes_ligase:]
	
	variables.positions_reactives["serum"]["positions"] = list(setting_reactives_in_places(number_tubes_serum, "serum", positions_coldblock_free).values())[0]
	positions_coldblock_free = positions_coldblock_free[number_tubes_serum:]
	
	variables.water_tubes["positions"] = list(setting_reactives_in_places(number_tubes_water, "water", positions_coldblock_free).values())[0]
	positions_coldblock_free = positions_coldblock_free[number_tubes_water:]
	
	variables.positions_reactives["mix"]["positions"] = list(setting_reactives_in_places(number_tubes_mix, "mix", positions_coldblock_free).values())[0]
	positions_coldblock_free = positions_coldblock_free[number_tubes_mix:]
	
	# print("++++++++++++++++++++++")
	# print(variables.positions_reactives)
	# print("++++++++++++++++++++++")
	
	##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# Ahora lo que vamos a hacer es distribuir el agua a todos los pocillos necesarios
	#Esta parte a continuacion la hemos cogido de mas abajo, es que la necesitamos aqui tambien
	positions_moclo_transform = []
	for labware_final in variables.positions_labware["final labware"]:
		positions_moclo_transform += labware_final.wells()
	positions_moclo_transform = positions_moclo_transform[:variables.number_reactions]
	generator_positions_wells_contructs = set_labware(positions_moclo_transform)
	
	last_pipette_used = variables.right_pipette # only to initialize
	for index_tube, volumes_transfering in enumerate(variables.water_tubes["volumes transfer from tube"]):
		# print(index_tube)
		# print(volumes_transfering)
		position_tube = variables.water_tubes["positions"][index_tube]
		# print(position_tube)
		# volumes_transfering = tube["volumes"]
		for transfer_well_volume in volumes_transfering:
			# print(transfer_well_volume)
			well_transfer = next(generator_positions_wells_contructs)
			# print(well_transfer)
			# print("=====")
			# print(variables.right_pipette.has_tip)
			# print(variables.left_pipette.has_tip)
			# print("=====")
			if transfer_well_volume != 0:
				pipette_transfer = give_me_optimal_pipette_to_use(transfer_well_volume, variables.left_pipette, variables.right_pipette)
				if pipette_transfer == last_pipette_used and pipette_transfer.has_tip == True:
					pipette_transfer.transfer(transfer_well_volume, position_tube, well_transfer, new_tip="never")
					pipette_transfer.blow_out(well_transfer.bottom())
				elif pipette_transfer == last_pipette_used and pipette_transfer.has_tip == False:
					check_tip_and_pick(pipette_transfer, variables)
					pipette_transfer.transfer(transfer_well_volume, position_tube, well_transfer, new_tip="never")
					pipette_transfer.blow_out(well_transfer.bottom())
				elif pipette_transfer != last_pipette_used and last_pipette_used.has_tip == False:
					check_tip_and_pick(pipette_transfer, variables)
					pipette_transfer.transfer(transfer_well_volume, position_tube, well_transfer, new_tip="never")
					pipette_transfer.blow_out(well_transfer.bottom())
					last_pipette_used = pipette_transfer
				else:
					last_pipette_used.drop_tip()
					check_tip_and_pick(pipette_transfer, variables)
					pipette_transfer.transfer(transfer_well_volume, position_tube, well_transfer, new_tip="never")
					pipette_transfer.blow_out(well_transfer.bottom())
					last_pipette_used = pipette_transfer
			else:
				pass
	last_pipette_used.drop_tip()
	##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	# We already have all the tubes, all the labware setted, now we need to do/create the big mix to distribute to all the final construction positions
	"""
	Lo vamos a hacer sin tener en cuenta las reacciones por tubo, porque por ahora solo tenemos 1
	ya luego lo iremos modificando para hacerlo cada vez mas general
	"""
	# We change the flow rate of th epipettes to transfer the proteings (ligase and enzyme restriction) because they are in a viscous element
	
	if variables.right_pipette != None:
		variables.right_pipette.flow_rate.aspirate = 10
		variables.right_pipette.flow_rate.dispense = 10
	if variables.left_pipette != None:
		variables.left_pipette.flow_rate.aspirate = 10
		variables.left_pipette.flow_rate.dispense = 10
	
	# Transfer enzyme
	# print("Enzyme")
	last_pipette_used = transfer_volume_tracking([variables.volume_restriction_enzyme*number_reactions_tube*variables.extra_pip_factor for number_reactions_tube in variables.positions_reactives["mix"]["reactions"]], variables.positions_reactives["restriction enzyme"]["positions"], variables.positions_reactives["mix"]["positions"], variables.positions_reactives["restriction enzyme"]["volume"], variables)
	last_pipette_used.drop_tip()
	
	# Transfer ligase
	# print("Ligase")
	last_pipette_used = transfer_volume_tracking([variables.volume_ligase*number_reactions_tube*variables.extra_pip_factor for number_reactions_tube in variables.positions_reactives["mix"]["reactions"]], variables.positions_reactives["ligase"]["positions"], variables.positions_reactives["mix"]["positions"], variables.positions_reactives["ligase"]["volume"], variables)
	last_pipette_used.drop_tip()
	
	# We RESTORE aspiration/dispensing speed to defaults
	if variables.right_pipette != None:
		variables.right_pipette.flow_rate.aspirate = 92.86
		variables.right_pipette.flow_rate.dispense = 92.86
	if variables.left_pipette != None:
		variables.left_pipette.flow_rate.aspirate = 92.86
		variables.left_pipette.flow_rate.dispense = 92.86
		
	# Transfer buffer
	# print("Buffer")
	last_pipette_used = transfer_volume_tracking([variables.volume_buffer*number_reactions_tube*variables.extra_pip_factor for number_reactions_tube in variables.positions_reactives["mix"]["reactions"]], variables.positions_reactives["buffer"]["positions"], variables.positions_reactives["mix"]["positions"], variables.positions_reactives["buffer"]["volume"], variables)
	last_pipette_used.drop_tip()
	
	# Transfer atp/serum, depending of the protocol
	# print("ATP")
	last_pipette_used = transfer_volume_tracking([variables.volume_serum*number_reactions_tube*variables.extra_pip_factor for number_reactions_tube in variables.positions_reactives["mix"]["reactions"]], variables.positions_reactives["serum"]["positions"], variables.positions_reactives["mix"]["positions"], variables.positions_reactives["serum"]["volume"], variables)
	last_pipette_used.drop_tip()
	
	# ------------------------ Distribution of the mix ------------------------
	"""
	En esta parte lo unico que vamos a hacer es ver a que posiciones hay que distribuir simplemente a todos los pocillos que tengan
	
	Como en este caso solo tenemos un tipo de mix, simplemente hay que poner a tantos wells como combinaciones haya.
	
	Lo primero que vamos a hacer es sacar la lista de posiciones de la placa y luego ya distribuimos ahi sin ningun problema con la funcion de
	distribucion vol tracking
	"""
	
	"""
	Vamos a hacer esta parte con el combinations_by_part pero he de decir que para lo unico que se usa variables.maps_reactives es para generar el mapa final y ese mapa final
	se puede generar de la misma manera con el variables.maps_reactives, por lo que se podria hacer directamente este sin tener que pasar por otro doble for loop
	
	variables.maps_reactives = [{'name': 'E0040m-1', 'parts': ['E0040m_CD-1', 'DVK_CD-1']}, {'name': 'C0012m-1', 'parts': ['C0012m_CD-1', 'DVK_CD-1']}]
	combinations_by_part = {'E0040m-1':['E0040m_CD-1', 'DVK_CD-1'], 'C0012m-1':['C0012m_CD-1', 'DVK_CD-1']}
	
	lo dejo aqui para que no tengas que printearlos de nuevo y se sepa la estructura que tienen ambas variables
	"""
	combinations_by_part = {}
	for combination in variables.maps_reactives:
		name_combination = combination["name"]
		for part_DNA in combination["parts"]:
			if name_combination in combinations_by_part.keys():
				combinations_by_part[name_combination].append(part_DNA)
			else:
				combinations_by_part[name_combination] = [part_DNA]
	
	# pipette_distributing = give_me_optimal_pipette_to_use(variables.mix_reactives, variables.right_pipette, variables.left_pipette)
	distribute_vol_tracking(variables.mix_reactives, variables.positions_reactives["mix"]["positions"], variables.positions_reactives["mix"]["reactions"], copy.deepcopy(positions_moclo_transform), variables)
	
	"""
	A partir de aqui es donde yo disiento de como lo hacen estas personas y yo lo haria de la siguiente manera (esto es lo esquematico, luego habria que hacerlo bien):
	
	for i in combinations_by_part:
		find_combination --> dara el well donde esta la combinacion a partir depositions_moclo_transform, lo que es equivalente al well final en un transfer
		for part in i:
			find_DNA --> nos dara el well source en un transfer porque correspondera al well del plasmido con la estructura
			transfer(find_DNA, find_combination)
	
	Es parecido a lo que hacen ellos pero le ponen mas pasos y se asgueran de que si tiene mas de un numeor de partes, sea diferente, algo que entiendo que sera
	diferente para nosotros
	"""
	
	dna_plate_dict = []
	for source_plate in variables.positions_labware["source labware"]:
		dna_plate_dict += source_plate.wells()
	
	number_rows_plate = len(labware_context.get_labware_definition(variables.name_source_plate)["ordering"][0])
	for name_combination, parts_combination in combinations_by_part.items():
		# print("+++++++++++++++")
		# print(name_combination)
		well_final_combination = find_combination(name_combination, variables.maps_reactives, positions_moclo_transform)
		# print(well_final_combination)
		# print(parts_combinations)
		number_module = 1
		for part_plasmid in parts_combination:
			# print("==========")
			well_source_DNApart = find_dna(part_plasmid, variables.map_dna_input, dna_plate_dict, number_rows_plate)
			# print(well_source_DNApart)
			if number_module == 1: # The first column, the first part to transfer is going always to be the acceptor vector, which has to be in the map as well
				volume_transfer = variables.volume_acceptor
				number_module += 1
			else:
				volume_transfer = variables.volume_module_plasmid
			pipette_transfer_parts = give_me_optimal_pipette_to_use(volume_transfer, variables.right_pipette, variables.left_pipette)
			check_tip_and_pick(pipette_transfer_parts, variables)
			pipette_transfer_parts.transfer(volume_transfer, well_source_DNApart, well_final_combination, new_tip = "never")
			pipette_transfer_parts.blow_out(well_final_combination.bottom())
			pipette_transfer_parts.drop_tip()
		# 	print("==========")
		# print("+++++++++++++++")
	
	# We print all that we need in the coldblocks and teh general information that will used as an user guide 
	print_information_user(variables)
	
