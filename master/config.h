// config file for master.C
//

// consts
const bool debug             = false;
const bool validate_eid      = false;
const bool debug_phi_cc      = false;
const bool debug_pp_mass     = false;
const bool genListOfGoodRuns = false;
const bool run_pp_time_corr  = false;
const bool run_pp_beta       = false;
const bool run_proton_beta   = false;

// number of strictness cases
Int_t n_strict        = 5;
Float_t strict[5]     = {1.00, 1.75, 2.50, 3.25, 5.00};

// ec_edep_inner, ec_geometric, ec_sampling, cc_theta, vertex
Int_t eid_stricts[5] = {2, 2, 2, 4, 2};

// parameters of binning for pion beta
Int_t n_p_bins_pion_beta     = 50;
Float_t p_max_pion_beta      = 3.5;
Float_t p_min_pion_beta      = 0.2;

// parameters of binning for proton beta
Int_t n_p_bins_proton_beta     = 50;
Float_t p_max_proton_beta      = 3.5;
Float_t p_min_proton_beta      = 0.2;
