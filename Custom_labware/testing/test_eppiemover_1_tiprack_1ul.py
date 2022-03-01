import json
from opentrons import protocol_api, types


TEST_TIPRACK_SLOT = '5'

RATE = 0.25  # % of default speeds
SLOWER_RATE = 0.1

PIPETTE_MOUNT = 'right'
PIPETTE_NAME = 'p1000_single_gen2'


TIPRACK_DEF_JSON = """{"ordering":[["A1"]],"brand":{"brand":"eppie_mover","brandId":["3d_printed"]},"metadata":{"displayName":"Eppie_mover 1 Tip Rack 1 µL","displayCategory":"tipRack","displayVolumeUnits":"µL","tags":[]},"dimensions":{"xDimension":127.76,"yDimension":85.47,"zDimension":106.8},"wells":{"A1":{"depth":20,"totalLiquidVolume":1,"shape":"circular","diameter":7.4,"x":13,"y":69.67,"z":86.8}},"groups":[{"metadata":{},"wells":["A1"]}],"parameters":{"format":"irregular","quirks":[],"isTiprack":true,"tipLength":20,"isMagneticModuleCompatible":false,"loadName":"eppiemover_1_tiprack_1ul"},"namespace":"custom_beta","version":1,"schemaVersion":2,"cornerOffsetFromSlot":{"x":0,"y":0,"z":0}}"""
TIPRACK_DEF = json.loads(TIPRACK_DEF_JSON)
TIPRACK_LABEL = TIPRACK_DEF.get('metadata', {}).get(
    'displayName', 'test labware')

metadata = {'apiLevel': '2.0'}


def run(protocol: protocol_api.ProtocolContext):
    tiprack = protocol.load_labware_from_definition(TIPRACK_DEF, TEST_TIPRACK_SLOT, TIPRACK_LABEL)
    pipette = protocol.load_instrument(
        PIPETTE_NAME, PIPETTE_MOUNT, tip_racks=[tiprack])

    num_cols = len(TIPRACK_DEF.get('ordering', [[]]))
    num_rows = len(TIPRACK_DEF.get('ordering', [[]])[0])


    def set_speeds(rate):
        protocol.max_speeds.update({
            'X': (600 * rate),
            'Y': (400 * rate),
            'Z': (125 * rate),
            'A': (125 * rate),
        })

        speed_max = max(protocol.max_speeds.values())

        for instr in protocol.loaded_instruments.values():
            instr.default_speed = speed_max

    set_speeds(RATE)
    firstwell = tiprack.well('A1')
    pipette.move_to(firstwell.top())
    protocol.pause("If the pipette is accurate click 'resume'")
    pipette.pick_up_tip()
    protocol.pause("If the pipette went into the center of the tip, click 'resume'")
    pipette.return_tip()
    protocol.pause("If the pipette successfully picked up the tip(s) but does not eject succesfully, pull the tip(s) off by hand and click 'resume'. Do not worry about tip ejection yet")

    last_col = (num_cols * num_rows) - num_rows
    if (PIPETTE_NAME == 'p20_multi_gen2' or PIPETTE_NAME == 'p300_multi_gen2'):
        well = tiprack.well(last_col)
        pipette.move_to(well.top())
        protocol.pause("If the position is accurate click 'resume'")
        pipette.pick_up_tip(well)
    else:
        last_well = (num_cols) * (num_rows)
        well = tiprack.well(last_well-1)
        pipette.move_to(well.top())
        protocol.pause("If the position is accurate click 'resume'")
        pipette.pick_up_tip(well)

    protocol.pause("If the pipette went to the center of the tip, click 'resume'")
    pipette.return_tip()
    protocol.comment("If the pipette successfully picked up the tip(s) but does not eject succesfully, pull the tip(s) off by hand and click 'resume'. Do not worry about tip ejection yet")

