# Cosas que hay que arreglar aqui:
# - Que sea random o los primeros x muestras que coger, a lo mejor se quiere hacer por random o las ultimas x
# - Variables y organizacion mejor
# - Hay que hacer control de errores

import opentrons.execute
from opentrons import protocol_api
import pandas as pd
import logging
import math
import copy
import numpy as np
import random

metadata = {
'apiLevel':'2.12'
}

def run(protocol: protocol_api.ProtocolContext):
	#Before everything because we will need it for the labware setting
	labware_context = opentrons.protocol_api.labware

	#We are going to create a logger to control errors that can occure in the script that are not taked in account and some of ours that will quit the program
	logging.basicConfig(format="""-------------------------------------------------------------------------------
--> SUM ERROR:
%(message)s
-------------------------------------------------------------------------------
--> INFO ERROR:""", level=logging.ERROR)
	logger = logging.getLogger()

# ----------------------------------
# ----------------------------------
	class setted_parameters():
		#def __init__(self, modular_parameters):
		vol_max_tube = 15000 # The same always because the transfer is setted with the heights of this falcon (15mL)
		# def __init__(self, variables_csv):
		def __init__(self,varibles_csv):
			self.number_source_plates = int(varibles_csv._get_value("Number of Source Plates","Value"))
			self.samples_per_plate = varibles_csv._get_value("Samples Per Plate","Value")
			self.name_right_pipette = varibles_csv._get_value("Name Right Pipette","Value")
			self.name_left_pipette = varibles_csv._get_value("Name Left Pipette","Value")
			# First we are going to check the paraameters and then we are going to settle the pipettes
			self.right_pipette = None
			self.left_pipette = None
			self.starting_tip_right_pip = varibles_csv._get_value("Initial Tip Right Pipette","Value")
			self.starting_tip_left_pip = varibles_csv._get_value("Initial Tip Left Pipette","Value")
			if varibles_csv._get_value("Replace Tipracks","Value").lower() == "true":
				self.replace_tiprack = True
			elif varibles_csv._get_value("Replace Tipracks","Value").lower() == "false":
				self.replace_tiprack = False
			self.name_map = variables_csv._get_value("Name maps","Value")
			self.volume_transfer_sample = float(variables_csv._get_value("Volume Colony transfer (uL)","Value"))
			self.volume_transfer_water = float(variables_csv._get_value("Volume Water transfer (uL)","Value"))
			self.water = {"position_tubes":[], "reactions_per_tube":[]}
			self.name_rack_falcon = variables_csv._get_value("Name 15mL Falcon Rack","Value")
			self.name_source_plate = variables_csv._get_value("Name Source Plates","Value")
			self.name_final_plate = variables_csv._get_value("Name Final Plates","Value")
			self.maps_source_plates_names = variables_csv._get_value("Name Source Plates Maps","Value")[1:-1].replace(" ","").split(",") # We have the list of the maps
			self.maps_source_plate = {}
			for index_map, map_source in enumerate(self.maps_source_plates_names):
				#self.maps_source_plate[index_map] = pd.read_csv("/data/user_storage/"+map_source+".csv", index_col = 0)
				self.maps_source_plate[index_map] = pd.read_csv(map_source+".csv", index_col = 0)
			self.type_selection = variables_csv._get_value("Type of Sample Selection", "Value")
			
			self.number_samples_source_plates = variables_csv._get_value("Number of samples in every source plate", "Value")
			self.index_start_source_plate = variables_csv._get_value("Index of start cell in source plate", "Value")
			#Aqui se puede mejorar lo de que tengamos el numero de muetsras que hay en cada plate y coger de manera random de ellas
			return
		
		def proccess_variables(self, name_parameter, value_default_empty, must_entries):
		# Some of the variables need processing to be usable in the script, usually is from a "(,),(,)" format to a dictionary
			dict_final = {}
			
			print(name_parameter)
			print(getattr(self, name_parameter))
			
			
			if name_parameter in ["number_samples_source_plates", "index_start_source_plate"]:
				value_parameter = getattr(self, name_parameter).replace(" ","")[1:-1].split('),(')
				print(value_parameter)
				for i, element in enumerate(value_parameter):
					value_parameter[i] = element.split(",")
				print(value_parameter)
				if len(value_parameter[0][0]) != 0:
					for well, plate in value_parameter:
						print(well)
						dict_final[int(plate)-1] = int(well)
				else:
						pass
			else:
				value_parameter = getattr(self, name_parameter).replace(" ","")[1:-1].split(',')
				if len(value_parameter) != 0:
					for plate, well in enumerate(value_parameter):
							dict_final[int(plate)] = int(well)
				else:
						pass
			
			# 
			# if len(value_parameter) != 0:
			# 		for plate, well in enumerate(value_parameter):
			# 				dict_final[int(plate)] = int(well)
			# else:
			# 		pass
				
			for i in range(must_entries):
				try:
					dict_final[i]
				except:
					dict_final[i] = value_default_empty
			
			setattr(self, name_parameter, dict_final)
			return
		
	def set_labware(labware_wells_name):
			# Generator of the positions to transfer, it will give you the next element of a list every time it is called
			for well in labware_wells_name:
				yield well

	def give_me_optimal_pipette_to_use(aVolume, pipette_r, pipette_l):
		# This is a function that will give us the pipette that we should use for the volume
		
		#For this function to work at least 1 pipette should be in the OT
		if pipette_r == None and pipette_l == None: #no pipettes attached
			raise Exception("There is not a pippette attached")
			#Error control
		
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
				# Error control
		return
	
	def check_tip_and_pick(pipette_used, variables):
		# This functions is used to pick tips and in case of need, replace the tiprack or add one to the labware
		# This way we add labwares in the process of the simulation and we do not need to calculate it a priori
		# In the OT-App it will appear directly in the deck but it has been added with this function
		try:
			pipette_used.pick_up_tip()
		except:
			if len(pipette_used.tip_racks) == 0: # There are no tipracks associated to this pipette yet (it is the beginning of the protocol)
				try:
					first_position_free = [position for position, labware in protocol.deck.items() if labware == None][0] # Check which is the first position free in the deck
				except:
					raise Exception("There is not enough space in the deck for this protocol, try less samples")
				
				pipette_used.tip_racks += [protocol.load_labware(define_tiprack(pipette_used.name), first_position_free)]
				
				#We establish now the starting tip, it will only be with the first addition, the rest will be establish that the first tip is in A1 directly
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
	
	def define_tiprack (pipette):
		# Define the correspondant tiprack to the pipette
		# This function should be updated as more pipettes are included in the possible ones
		
		if pipette == "p20_single_gen2" or pipette == "p20_multi_gen2":
			return "opentrons_96_tiprack_20ul"
		elif pipette == "p300_single_gen2" or pipette == "p300_multi_gen2":
			return "opentrons_96_tiprack_300ul"
		elif pipette == "p1000_single_gen2":
			return "opentrons_96_tiprack_1000ul"
		else:
			raise Exception("The pippette that you are using is not established.")
		return
	
	def number_tubes_needed (vol_reactive_per_reaction_factor, number_reactions, vol_max_tube):
		# With this function we are going to establish how many tubes we need for a volume of reactive
		# This function does not garante the lower number of tubes (for that check versions lower than 10) but it assures that everything can be picked with the pipettes that we have
		
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
		# We are going to set how many labwares we need (in this case how much coldblocks and final plates we need)
		# We are not going to calculate the number (number_plates), we are only setting them
		try:
			position_plates = [position for position, labware in protocol.deck.items() if labware == None] # We obtain the positions in which there are not labwares
			all_plates = []
			for i in range (number_plates):
				plate = protocol.load_labware(labware_name, position_plates[i])
				all_plates.append(plate)
			return all_plates
		
		except: # There were not enough items in the position_plate to fit the number of labwares that we needed
			# This exception can only work properly if the labwares can be loaded with load_labware() but we have made sure of it previously in this script
			raise Exception("There is not enough space in the deck for this protocol, try less samples")
	
	def setting_labware(variables):
		# In this function we are going to calculate and set the labwares in the deck according to the csv provided, with the exception of the tipracks
		# Also, the number of tubes of each reactive will be calculated here with their respectives number of reactions
		
		# We start settign the source labware which number has been provided
		labware_source = setting_number_plates(variables.number_source_plates, variables.name_source_plate)
		
		# Set the final plates
		# We need to calculat ehow many plates do we need by sum (cells by plate) / wells per final plate
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
	
	def distribute_z_tracking(pipette_used, vol_source, vol_distribute_well, pos_source, pos_final): #modified from the IPTG script
		while len(pos_final) != 0:
			if position_dispense_aspirate(vol_source, pos_source).point == position_dispense_aspirate((vol_source-(len(pos_final)*vol_distribute_well)), pos_source).point:
				pipette_used.distribute(vol_distribute_well, position_dispense_aspirate(vol_source, pos_source), pos_final, new_tip = "never", disposal_volume = 0)
				vol_source -= vol_distribute_well*len(pos_final)
				pos_final = []
			else:
				pos_final_original = pos_final
				for number_positions in range(1, len(pos_final_original)+1): #vamos mirando hasta cuanto puede distribuir sin tener que cambair la altura
					vol_needed = vol_distribute_well*len(pos_final[:number_positions])
					if position_dispense_aspirate(vol_source, pos_source) == position_dispense_aspirate((vol_source-vol_needed), pos_source):
						next
					else:
						pipette_used.distribute(vol_distribute_well, position_dispense_aspirate(vol_source, pos_source), pos_final[:number_positions], new_tip = "never", disposal_volume = 0)
						vol_source -= vol_distribute_well*len(pos_final[:number_positions])
						pos_final = pos_final[number_positions:]
						break
		return
	
	def position_dispense_aspirate (vol_falcon, theory_position):
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
	
	def distribute_z_tracking(pipette_used, vol_source, vol_distribute_well, pos_source, pos_final): #modified from the IPTG script
		while len(pos_final) != 0:
			if position_dispense_aspirate(vol_source, pos_source).point == position_dispense_aspirate((vol_source-(len(pos_final)*vol_distribute_well)), pos_source).point:
				pipette_used.distribute(vol_distribute_well, position_dispense_aspirate(vol_source, pos_source), pos_final, new_tip = "never", disposal_volume = 0)
				vol_source -= vol_distribute_well*len(pos_final)
				pos_final = []
			else:
				pos_final_original = pos_final
				for number_positions in range(1, len(pos_final_original)+1): #vamos mirando hasta cuanto puede distribuir sin tener que cambair la altura
					vol_needed = vol_distribute_well*len(pos_final[:number_positions])
					if position_dispense_aspirate(vol_source, pos_source) == position_dispense_aspirate((vol_source-vol_needed), pos_source):
						next
					else:
						pipette_used.distribute(vol_distribute_well, position_dispense_aspirate(vol_source, pos_source), pos_final[:number_positions], new_tip = "never", disposal_volume = 0)
						vol_source -= vol_distribute_well*len(pos_final[:number_positions])
						pos_final = pos_final[number_positions:]
						break
		return
	
	def pds_labware_sources(labware_source):
		# This function is needed for the creation of the map, we are creating here the pandas dataframe with the dimensions of the labware
		# This function is setted to labware which has rows with alphabetical name (letters) and columns with numbers
		
		labware_wells = labware_context.get_labware_definition(labware_source[0].name)["groups"][0]["wells"]
		name_rows = sorted(list(set([well[0] for well in labware_wells])))
		number_rows = len(name_rows)
		number_columns = len(set([well[1:] for well in labware_wells]))
		
		dataframes_source_with_index = {}
		
		for labware in labware_source:
			dataframes_source_with_index[int(str(labware).split(" ")[-1])] = pd.DataFrame(np.full((number_rows,number_columns),None),columns=list(range(1,number_columns+1)),index=name_rows)
			dataframes_source_with_index[int(str(labware).split(" ")[-1])].index.name = "Row/Column"
	
		return dataframes_source_with_index

	def map_tracking(pos_source, pos_final, dataframes_with_index, map_source_plate):
		# This function process the name of the final well and inserts it in the source plate
		# The results is a map of the source plate filled with the well names in which that colony is sittuated
		
		pos_source_deck = str(pos_source).split(" ")[-1]
		pos_source_labware = str(pos_source).split(" ")[0]
		
		column_pos_cource = str(pos_source_labware[1:])
		row_pos_source = str(pos_source_labware[0])
		
		pos_final_deck = str(pos_final).split(" ")[-1]
		pos_final_labware = str(pos_final).split(" ")[0]
		
		
		identifier_sample = map_source_plate[int(pos_source_deck)-1][column_pos_cource][row_pos_source]
	
		
		dataframes_with_index[int(pos_final_deck)][int(pos_final_labware[1:])][pos_final_labware[0]] = str(identifier_sample+" on "+pos_source_deck)
		#We are going to still put the position on the deck because 2 samples can have the same identifier in different plates
		
		return
	
#++++++++++++++++++++++++++++++++
# Por ahora estamos asumiendo que estan llenos los plates, luego lo haremos mas complejo
	def wells_transfer(source_plates, final_plates, variables):
		# We state which are going to be the final wells, the ones that we are going to transfer to
		# final_wells = []
		# for index_labware, labware in final_plates:
		# 	final_wells += labware.wells()[variables.index_start_source_plate[index_labware]:variables.number_samples_source_plate[index_labware]]
	
		final_wells = sum([labware.wells() for labware in final_plates],[])[:sum(variables.samples_per_plate.values())]
		
		# Now lets take the source wells that we are going to aspirate from
		source_well_distribute = []
		for index_labware, samples_labware in variables.samples_per_plate.items():
			possible_wells_plate = source_plates[index_labware].wells()[variables.index_start_source_plate[index_labware]:variables.number_samples_source_plates[index_labware]]
			if variables.type_selection == "first":
				source_well_distribute += possible_wells_plate[:samples_labware]
			elif variables.type_selection == "last":
				source_well_distribute += list(reversed(possible_wells_plate))[:samples_labware]
			elif variables.type_selection == "random":
				source_well_distribute += random.choices(source_plates[index_labware].wells(), k = samples_labware)
		
		return (source_well_distribute, final_wells)
#++++++++++++++++++++++++++++++++
		
		
		


	def print_information_user(variables, tuberack_positions):
		# The position and the reactions per tube are in the variables (variables.antibiotic)
		# This function purpose is to print all the information needed for the user to perform the protocol
		# This will be printed if performed with opentrons_simulate, it will not be an output of the OT-App
		
		print("""--------------------------------------------------------------\n--------------------------------------------------------------\nINFORMATION FILE
	Data the users need to set the OT and perform the protocol
	It is the output of the python file for the protocol when executed with the opentrons_simulate option\n""")
		#Some general information that the user has setted
		print("--------------------------------------------------------------\nGENERAL INFORMATION\nThis details are set by the user, this is only a remainder of some variables")
		print("\t- Number total of samples (all source plates): "+str(sum(variables.samples_per_plate.values())))
		print("\t- Final volume of output labware (uL): "+str(variables.volume_transfer_water+variables.volume_transfer_sample))
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
				print(tuberack_table.fillna("-")) # for esthetic purposes
				print("\n--------------------------------------------------------------")
		else:
			print("\n--------------------------------------------------------------")

		return
	
# ----------------------------------
# ----------------------------------
	try:
		current_step = "Reading csv and transforming them to parameters/variables"
		# Setting variables and calculating others
		#variables_csv = pd.read_csv("/data/user_storage/Variables_Generation_MasterPlatePCR-version.csv", index_col = 0)
		variables_csv = pd.read_csv("Variables_Generation_MasterPlatePCR-version.csv", index_col = 0)
		print(variables_csv)
		# We are going to convert these parameters into arguments of the class variables and we are going to process some of them so they can be usable (they are going to be dictionaries)
		variables = setted_parameters(variables_csv)
		setted_parameters.proccess_variables(variables, "samples_per_plate", 0, variables.number_source_plates) # We take 0 as default
		
		number_wells_source_plate = len(labware_context.get_labware_definition(variables.name_source_plate)["wells"])
		setted_parameters.proccess_variables(variables, "number_samples_source_plates", number_wells_source_plate, variables.number_source_plates)
		setted_parameters.proccess_variables(variables, "index_start_source_plate", 0, variables.number_source_plates)
		
		print(vars(variables))
		
		
		variables.right_pipette = protocol.load_instrument(variables.name_right_pipette, mount = "right")
		variables.left_pipette = protocol.load_instrument(variables.name_left_pipette, mount = "left")
		
		labware_source, labware_final, labware_falcons = setting_labware(variables)
		
		# We take the wells that we are going to need and the ones taht we are going to dispense to
		# wells_total = []
		# for labware in labware_final:
		# 	wells_total += labware.wells()
		# wells_dispense = wells_total[:sum(variables.samples_per_plate.values())]
		
		source_well_take, wells_dispense = wells_transfer(labware_source, labware_final, variables)

		# Create the dataframes for the mapping
		# First we create the maps for filling, we are creating it even if we dont want it
		maps_final_plates = pds_labware_sources(labware_final)
		
		# Distribute and track in the dataframes corresponding to the final labwares
		# First we create the list of wells from which we are going to take a sample from ++++++++++++++++++++++ ESTA ES LA PARTE EN LA CUAL PUEDO CAMBIAR POR RANDOM, LAST, FIRST
		# source_well_distribute = []
		# for index_labware, samples_labware in variables.samples_per_plate.items():
		# 	source_well_distribute += labware_source[index_labware].wells()[:samples_labware]
		source_well_distribute = set_labware(source_well_take)
		final_well_distribute = set_labware(wells_dispense)
		used_pipette = give_me_optimal_pipette_to_use(variables.volume_transfer_sample, variables.right_pipette, variables.left_pipette)
		for i in range(len(wells_dispense)):
			# transfer
			well_source = next(source_well_distribute)
			well_final = next(final_well_distribute)
			check_tip_and_pick(used_pipette, variables)
			used_pipette.transfer(variables.volume_transfer_sample, well_source, well_final, new_tip="never")
			used_pipette.drop_tip()
			#map
			map_tracking(well_source, well_final, maps_final_plates, variables.maps_source_plate)
		
		# Export the maps
		current_step = "Importing the map of cells and primers for each final plate"
		for map_source in maps_final_plates:
			maps_final_plates[map_source] = maps_final_plates[map_source].fillna(value="-")
			maps_final_plates[map_source].to_csv(str(variables.name_map)+"_"+str(map_source)+".csv")
		
		# Printing user's info
		print_information_user(variables, labware_falcons)
		
		# Homing
		protocol.home()
		
	except Exception as e:
		print("-------------------------------------------------------------------------------")
		print("--> ERROR STEP:\n" + current_step)
		logger.critical(e, exc_info = True)
		print("-------------------------------------------------------------------------------")

# ----------------------------------
# ----------------------------------

#Hay que anadir una funcionq ue sea decirle de donde empieza en las placas iniciales, lo de random y otras cosillas, este hay que mejorarlo
