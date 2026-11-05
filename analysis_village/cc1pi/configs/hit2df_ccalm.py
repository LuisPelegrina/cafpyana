from analysis_village.cc1pi.makedf.make_cc1pidf import *

#DFS = [make_cc1pidf, make_cc1pinudf, make_hdrdf, make_potdf_bnb, make_trkhitdf_plane0, make_trkhitdf_plane1, make_trkhitdf_plane2]
#NAMES = ["cc1pi", "nudf", "hdr", "pot", "hit0", "hit1", "hit2"]
DFS = [make_pandora_df, make_cc1pinudf_no_syst, make_hdrdf, make_trkhitdf_plane2_ccal_m, make_potdf_bnb]
NAMES = ["pandora", "nudf", "hdr", "hit2", "pot"]
