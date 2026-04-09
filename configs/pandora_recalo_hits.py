from makedf.makedf import *

DFS = [make_pandora_df_calo_update, make_hdrdf, make_potdf_bnb, make_mcnuwgtdf_slim, make_opflashdf, make_trkhitdf_plane0, make_trkhitdf_plane1, make_trkhitdf_plane2]
NAMES = ["evt", "hdr", "pot", "mcnu", "opflash", "hit0", "hit1", "hit2"]
#DFS = [make_pandora_df_calo_update, make_hdrdf, make_potdf_bnb, make_mcnudf, make_opflashdf]
#NAMES = ["evt", "hdr", "pot", "mcnu", "opflash"]
