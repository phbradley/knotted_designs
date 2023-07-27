// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief

#ifndef INCLUDED_apps_pilot_phil_phil_options_HH
#define INCLUDED_apps_pilot_phil_phil_options_HH

#include <apps/pilot/phil/types.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
// #include <core/scoring/methods/EnvElecInfo.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/conformation/symmetry/command_line_hack.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <apps/pilot/phil/phil.hh>
#include <cstdio>
#include <utility/file/file_sys_util.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/phil.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

using namespace core;
using  basic::options::StringOptionKey;
using  basic::options::IntegerOptionKey;
using  basic::options::RealOptionKey;
using  basic::options::BooleanOptionKey;
using  basic::options::StringVectorOptionKey;
using  basic::options::IntegerVectorOptionKey;
using  basic::options::RealVectorOptionKey;
using  basic::options::BooleanVectorOptionKey;

using  basic::options::option;
namespace OK = basic::options::OptionKeys;

basic::Tracer TR_OPTIONS( "apps.pilot.phil.phil_options" );

namespace simtime {
//clock_t start_clock, last_clock;
time_t start_time, last_time, decoy_timer_start_time;
Real max_simtime;
}

void
init_simtime()
{
	simtime::start_time = time( NULL );
	simtime::last_time = simtime::start_time;
	simtime::max_simtime = 0.0;
	//
	simtime::decoy_timer_start_time = 0; // flag failure to start timer
}

/// options common to multiple apps:
namespace my_options {
RealOptionKey grid_spacing("my:grid_spacing");
RealOptionKey cheating_frags_ratio("my:cheating_frags_ratio");
StringOptionKey resfile_dir("my:resfile_dir");
StringOptionKey tmpdir("my:tmpdir");
StringOptionKey star_fragment_rebuild_mintype("my:star_fragment_rebuild_mintype");
StringOptionKey resfile_tag("my:resfile_tag");
RealOptionKey target_twist("my:target_twist");
RealOptionKey target_rise("my:target_rise");
RealOptionKey ftmc_iterations_scaling("my:ftmc_iterations_scaling");
StringOptionKey design_score_function("my:design_score_function");
StringOptionKey final_round_design_score_function("my:final_round_design_score_function");
StringOptionKey pdb_coords_file("my:pdb_coords_file");
StringOptionKey beta_pairings_file("my:beta_pairings_file");
StringOptionKey perfect_helix_file("my:perfect_helix_file");
StringOptionKey r_thetas_file("my:r_thetas_file");
StringVectorOptionKey ab_turns("my:ab_turns" );
StringVectorOptionKey ba_turns("my:ba_turns" );
StringVectorOptionKey ab_turns_up("my:ab_turns_up");
StringVectorOptionKey ab_turns_down("my:ab_turns_down");
StringVectorOptionKey ba_turns_up("my:ba_turns_up");
StringVectorOptionKey ba_turns_down("my:ba_turns_down");
StringOptionKey beta_transforms_file("my:beta_transforms_file");
StringOptionKey layers_file("my:layers_file");
StringOptionKey psipred_directory("my:psipred_directory");
RealOptionKey cart_bonded_weight("my:cart_bonded_weight");
RealOptionKey core_threshold_sse("my:core_threshold_sse");
RealOptionKey core_threshold_loop("my:core_threshold_loop");
RealOptionKey surface_threshold_sse("my:surface_threshold_sse");
RealOptionKey surface_threshold_loop("my:surface_threshold_loop");
BooleanOptionKey use_layers_sequence_bonus("my:use_layers_sequence_bonus");
BooleanOptionKey cleanup_tmp("my:cleanup_tmp");
RealOptionKey favor_native_residue("my:favor_native_residue");
BooleanOptionKey layer_design("my:layer_design");
BooleanOptionKey verbose("my:verbose");
BooleanOptionKey new_unsats("my:new_unsats");
BooleanOptionKey show_pymol("my:show_pymol");
BooleanOptionKey single_sim("my:single_sim");
BooleanOptionKey make_trajectory("my:make_trajectory");
BooleanOptionKey only_minimize("my:only_minimize");
BooleanOptionKey fastrelax_before_design("my:fastrelax_before_design");
BooleanOptionKey force_helix_capping("my:force_helix_capping");
StringOptionKey mr_data_directory("my:mr_data_directory");
IntegerOptionKey nrepeat_for_refolding_rmsd("my:nrepeat_for_refolding_rmsd" );
StringOptionKey mode("my:mode");
StringOptionKey helix_params_file("my:helix_params_file");
StringOptionKey strand_params_file("my:strand_params_file");
BooleanOptionKey disallow_large_exposed_nonpolar("my:disallow_large_exposed_nonpolar");
RealOptionKey exposed_nonpolar_sasa_fraction("my:exposed_nonpolar_sasa_fraction");
BooleanOptionKey force_dna_symmetry( "my:force_dna_symmetry" );
BooleanOptionKey dump_hbonds( "my:dump_hbonds" );
BooleanOptionKey use_tp1_waters( "my:use_tp1_waters" );
BooleanOptionKey bhlh_mode( "my:bhlh_mode" );
BooleanOptionKey test_resampling("my:test_resampling");
BooleanOptionKey flexy_pose( "my:flexy_pose" );
BooleanOptionKey cpdb_input( "my:cpdb_input" );
BooleanOptionKey use_bhlh_constraints( "my:use_bhlh_constraints" );
BooleanOptionKey ramp_fa_dun( "my:ramp_fa_dun" );
BooleanOptionKey dna_root_strand_flip_hack( "my:dna_root_strand_flip_hack" );
IntegerVectorOptionKey ssdna_score_positions( "my:ssdna_score_positions" );
IntegerOptionKey dna_root_shift_hack( "my:dna_root_shift_hack" );
IntegerOptionKey dna_root_shift_window( "my:dna_root_shift_window" );
IntegerOptionKey outer_nstruct_iterations( "my:outer_nstruct_iterations" );
IntegerOptionKey n_cycles( "my:n_cycles" );
IntegerOptionKey n_cycles_rottrials( "my:n_cycles_rottrials" );
StringOptionKey culled_chains( "culled_chains" );
IntegerOptionKey extra_cycles( "my:extra_cycles" );
RealOptionKey buried_nonpolar_bonus( "my:buried_nonpolar_bonus" );
RealOptionKey exposed_charged_bonus( "my:exposed_charged_bonus" );
RealOptionKey exposed_nonpolar_penalty( "my:exposed_nonpolar_penalty" );
RealOptionKey exposed_thr_penalty( "my:exposed_thr_penalty" );
RealOptionKey at_bias_correction( "my:at_bias_correction" );
RealOptionKey extra_temperature( "my:extra_temperature" );
RealOptionKey target_runtime( "my:target_runtime" ); // in minutes
StringVectorOptionKey protocols("my:protocols");
StringVectorOptionKey dna_target_defs("my:dna_target_defs");
StringVectorOptionKey dna_flex_defs("my:dna_flex_defs");
StringVectorOptionKey protein_interface_defs("my:protein_interface_defs");
RealVectorOptionKey tp3_params("my:tp3_params");
RealOptionKey tp3_radius("my:tp3_radius");
IntegerOptionKey zf_fast_relax_inner_cycles_factor( "my:zf_fast_relax_inner_cycles_factor" );
IntegerOptionKey zf_fast_relax_outer_cycles( "my:zf_fast_relax_outer_cycles" );
IntegerVectorOptionKey flatten_pwm_columns( "my:flatten_pwm_columns" );
IntegerVectorOptionKey dna_motif( "my:dna_motif" );
StringOptionKey pack_mode( "my:pack_mode" );
BooleanOptionKey rescale_tp3_charges_as_hbonds( "my:rescale_tp3_charges_as_hbonds" );
BooleanOptionKey remove_water_at_start( "my:remove_water_at_start" );
BooleanOptionKey start_with_native_sidechains( "my:start_with_native_sidechains" );
BooleanOptionKey idealize_dna( "my:idealize_dna" );
BooleanOptionKey hydrate_dna( "my:hydrate_dna" );
BooleanOptionKey hydrate_minor_groove( "my:hydrate_minor_groove" );
BooleanOptionKey fast_relax( "my:fast_relax" );
BooleanOptionKey hydrate_phosphate( "my:hydrate_phosphate" );
BooleanOptionKey deriv_check( "my:deriv_check" );
BooleanOptionKey output_pdb_files( "output_pdb_files" );
BooleanOptionKey use_intra_dna_cp_crossover_4( "my:use_intra_dna_cp_crossover_4" );
RealOptionKey output_pdb_fraction( "output_pdb_fraction" );
RealOptionKey score_filter_acceptance_rate("my:score_filter_acceptance_rate" );
BooleanOptionKey score_filter_pass_early("my:score_filter_pass_early" );
BooleanOptionKey no_mutation_moves_at_edges("my:no_mutation_moves_at_edges" );
RealOptionKey random_filter_acceptance_rate("my:random_filter_acceptance_rate" );
//StringOptionKey fa_scorefxn("my:fa_scorefxn");
StringOptionKey centroid_score_function("my:centroid_score_function");
StringOptionKey prefix("my:prefix");
StringOptionKey shared_output_tag("my:shared_output_tag" );
StringVectorOptionKey score_functions("my:score_functions");
StringVectorOptionKey set_weight("my:set_weight");
StringVectorOptionKey set_atom_property("my:set_atom_property");
StringVectorOptionKey set_extra_atom_property("my:set_extra_atom_property");
StringVectorOptionKey params("my:params");
StringOptionKey params_file("my:params_file");
RealOptionKey cluster_fraction("my:cluster_fraction");
RealVectorOptionKey cluster_fractions( "my:cluster_fractions" );
RealVectorOptionKey clustering_params( "my:clustering_params" );
RealOptionKey profile_output_fraction("my:profile_output_fraction");
RealOptionKey rescale_phosphate_charges("my:rescale_phosphate_charges");
IntegerOptionKey max_decoys_per_cluster_pdb( "my_options:max_decoys_per_cluster_pdb" );
IntegerOptionKey ramping_cycles( "my:ramping_cycles" );
StringOptionKey refseq_file( "my:refseq_file" );
StringOptionKey min_type( "my:min_type" );
//RealOptionKey fullatom_dna_dihedral_weight("my:fullatom_dna_dihedral_weight");
//RealOptionKey fullatom_dna_dihedral_sugar_weight("my:fullatom_dna_dihedral_sugar_weight");
//RealOptionKey fullatom_dna_chainbreak_weight("my:fullatom_dna_chainbreak_weight");
BooleanOptionKey flex_dna( "my:flex_dna" );
BooleanOptionKey flex_dna_sugar( "my:flex_dna_sugar" );
BooleanOptionKey flex_dna_backbone( "my:flex_dna_backbone" );
BooleanOptionKey flex_protein_backbone( "my:flex_protein_backbone" );
BooleanOptionKey flex_protein_dna_jump( "my:flex_protein_dna_jump" );
BooleanOptionKey superimpose_decoys( "my:superimpose_decoys" );
BooleanOptionKey use_hb_env_dep( "my:use_hb_env_dep" );
BooleanOptionKey use_hb_env_dep_DNA( "my:use_hb_env_dep_DNA" );
BooleanOptionKey allow_file_stealing( "my:allow_file_stealing" );
BooleanOptionKey fa_phosphate_constraints( "my:fa_phosphate_constraints" );
BooleanOptionKey use_dna_chainbreak_constraints( "my:use_dna_chainbreak_constraints" );
BooleanOptionKey hacking( "my:hacking" );
RealOptionKey h2o_ref_wt( "my:h2o_ref_wt" );
IntegerOptionKey mutations_per_mutable_base_pair( "my:mutations_per_mutable_base_pair" );
StringVectorOptionKey adjust_ref_weights("my:adjust_ref_weights" );
BooleanOptionKey use_softrep_for_early_design("my:use_softrep_for_early_design" );
BooleanOptionKey use_softrep_for_design("my:use_softrep_for_design" );
IntegerOptionKey limit_design_positions_by_chain_contacts("my:limit_design_positions_by_chain_contacts" );
StringOptionKey unfolded_sasas("my:unfolded_sasas");
RealOptionKey distance_tolerance("my:distance_tolerance");
RealOptionKey cheating_frags_min_torsion_dev("my:cheating_frags_min_torsion_dev");
RealOptionKey cheating_frags_max_torsion_dev("my:cheating_frags_max_torsion_dev");
RealOptionKey autobuild_rwork_threshold("my:autobuild_rwork_threshold");
BooleanOptionKey generate_r_free_flags("my:generate_r_free_flags");
BooleanOptionKey sgalt_all("my:sgalt_all");
BooleanOptionKey sgalt_none("my:sgalt_none");
BooleanOptionKey skip_tncs("my:skip_tncs");
BooleanOptionKey full_search("my:full_search");
StringOptionKey mtz_file("my:mtz_file");
IntegerOptionKey composition("my:composition");
IntegerOptionKey num_models("my:num_models");
RealOptionKey llg_threshold_for_refinement("my:llg_threshold_for_refinement");
StringVectorOptionKey space_groups("my:space_groups");
StringOptionKey space_group("my:space_group");
BooleanOptionKey local_simfile("my:local_simfile");
}

void
add_phil_options()
{
	option.add( my_options::mode, "mode" );
	option.add( my_options::cheating_frags_ratio, "cheating_frags_ratio" ).def(3.0); // eg 100 to 33
	option.add( my_options::verbose, "verbose" );
	option.add( my_options::tmpdir, "tmpdir" ).def("./");
	option.add( my_options::psipred_directory, "psipred_directory" ).def("/home/pbradley/download/psipred/psipred/");
	option.add( my_options::cleanup_tmp, "cleanup_tmp" );
	option.add( my_options::grid_spacing, "grid_spacing").def(2.75); // 3.48 gives same density as water; 2.75 is double
	option.add( my_options::star_fragment_rebuild_mintype, "star_fragment_rebuild_mintype" );
	option.add( my_options::final_round_design_score_function, "final_round_design_score_function" );
	option.add( my_options::design_score_function, "design_score_function" );
	option.add( my_options::ftmc_iterations_scaling, "ftmc_iterations_scaling" ).def(250);
	option.add( my_options::resfile_dir, "resfile_dir" );
	option.add( my_options::resfile_tag, "resfile_tag" );
	option.add( my_options::perfect_helix_file, "perfect_helix_file" ).def("/home/pbradley/csdat/input/perfect_helix_1elwA_A108_A115.pdb");
	option.add( my_options::distance_tolerance, "distance_tolerance" );
	option.add( my_options::cart_bonded_weight, "cart_bonded_weight" ).def( 0.5 );
	option.add( my_options::favor_native_residue, "favor_native_residue" );
	option.add( my_options::new_unsats, "new_unsats" ).def( true ); // DANGER
	option.add( my_options::pdb_coords_file, "pdb_coords_file" );
	option.add( my_options::show_pymol, "show_pymol" );
	option.add( my_options::limit_design_positions_by_chain_contacts, "limit_design_positions_by_chain_contacts" );
	option.add( my_options::r_thetas_file, "r_thetas_file" );
	option.add( my_options::target_twist, "target_twist" );
	option.add( my_options::target_rise, "target_rise" );
	option.add( my_options::beta_pairings_file, "beta_pairings_file" );
	option.add( my_options::ab_turns_up, "ab_turns_up" );
	option.add( my_options::ab_turns_down, "ab_turns_down" );
	option.add( my_options::ba_turns_up, "ba_turns_up" );
	option.add( my_options::ba_turns_down, "ba_turns_down" );
	option.add( my_options::ab_turns, "ab_turns" );
	option.add( my_options::ba_turns, "ba_turns" );
	option.add( my_options::beta_transforms_file, "beta_transforms_file" );
	option.add( my_options::force_helix_capping, "force_helix_capping" );
	option.add( my_options::core_threshold_sse, "core_threshold_sse" ).def( 15 );
	option.add( my_options::core_threshold_loop, "core_threshold_loop" ).def( 25 );
	option.add( my_options::surface_threshold_sse, "surface_threshold_sse" ).def( 60 );
	option.add( my_options::surface_threshold_loop, "surface_threshold_loop" ).def( 40 );
	option.add( my_options::layer_design, "layer_design" );
	option.add( my_options::use_layers_sequence_bonus, "use_layers_sequence_bonus" );
	option.add( my_options::layers_file, "layers_file" );
	option.add( my_options::disallow_large_exposed_nonpolar, "disallow_large_exposed_nonpolar" );
	option.add( my_options::exposed_nonpolar_sasa_fraction, "exposed_nonpolar_sasa_fraction" );
	option.add( my_options::make_trajectory, "make_trajectory" );
	option.add( my_options::test_resampling, "test_resampling" );
	option.add( my_options::buried_nonpolar_bonus   ,    "buried_nonpolar_bonus" ).def( -0.5  );
	option.add( my_options::exposed_nonpolar_penalty, "exposed_nonpolar_penalty" ).def(  0.5  );
	option.add( my_options::exposed_charged_bonus   ,    "exposed_charged_bonus" ).def( -0.25 );
	option.add( my_options::exposed_thr_penalty     ,      "exposed_thr_penalty" ).def(  0.25 );
	option.add( my_options::nrepeat_for_refolding_rmsd, "nrepeat_for_refolding_rmsd" );
	option.add( my_options::helix_params_file, "helix_params_file" );
	option.add( my_options::strand_params_file, "strand_params_file" );
	option.add( my_options::clustering_params, "clustering_params" );
	option.add( my_options::only_minimize, "only_minimize" );
	option.add( my_options::dump_hbonds, "dump_hbonds" ).def( true );
	option.add( my_options::target_runtime, "target_runtime" );
	option.add( my_options::single_sim, "single_sim" );
	option.add( my_options::unfolded_sasas, "unfolded_sasas" ).def("/home/pbradley/gitrepos/symdes/input/reference_fully_exposed_sasa_values.txt"); // cheesy hardcoded default
	option.add( my_options::use_softrep_for_design, "use_softrep_for_design" );
	option.add( my_options::use_softrep_for_early_design, "use_softrep_for_early_design" );
	option.add( my_options::fastrelax_before_design, "fastrelax_before_design" );
	option.add( my_options::adjust_ref_weights, "adjust_ref_weights" );
	option.add( my_options::force_dna_symmetry, "force_dna_symmetry" );
	option.add( my_options::ssdna_score_positions, "ssdna_score_positions" );
	option.add( my_options::params_file, "params_file" );
	option.add( my_options::cpdb_input, "cpdb_input" );
	option.add( my_options::flexy_pose, "flexy_pose" );
	option.add( my_options::use_tp1_waters, "use_tp1_waters" );
	option.add( my_options::mr_data_directory, "mr_data_directory" );
	option.add( my_options::dna_root_strand_flip_hack, "dna_root_strand_flip_hack" );
	option.add( my_options::dna_root_shift_hack, "dna_root_shift_hack" ).def(0);
	option.add( my_options::dna_root_shift_window, "dna_root_shift_window" );
	option.add( my_options::no_mutation_moves_at_edges, "no_mutation_moves_at_edges" );
	option.add( my_options::bhlh_mode, "bhlh_mode" );
	option.add( my_options::idealize_dna, "idealize_dna" );
	option.add( my_options::use_bhlh_constraints, "use_bhlh_constraints" );
	option.add( my_options::dna_motif, "dna_motif" );
	option.add( my_options::flatten_pwm_columns, "flatten_pwm_columns" );
	option.add( my_options::protein_interface_defs, "protein_interface_defs" );
	option.add( my_options::dna_target_defs, "dna_target_defs" );
	option.add( my_options::dna_flex_defs, "dna_flex_defs" );
	option.add( my_options::at_bias_correction, "at_bias_correction" );
	option.add( my_options::outer_nstruct_iterations, "outer_nstruct_iterations" ).def( 7 );
	option.add( my_options::mutations_per_mutable_base_pair, "mutations_per_mutable_base_pair" ).def( 25 );
	option.add( my_options::ramp_fa_dun, "ramp_fa_dun" );
	option.add( my_options::hacking, "hacking" );
	option.add( my_options::rescale_tp3_charges_as_hbonds, "rescale_tp3_charges_as_hbonds" );
	option.add( my_options::tp3_params, "tp3_params" );
	option.add( my_options::tp3_radius, "tp3_radius" );
	option.add( my_options::use_dna_chainbreak_constraints, "use_dna_chainbreak_constraints" );
	option.add( my_options::remove_water_at_start, "remove_water_at_start" );
	option.add( my_options::start_with_native_sidechains, "start_with_native_sidechains" );
	option.add( my_options::n_cycles, "n_cycles" ).def( 0 );
	option.add( my_options::n_cycles_rottrials, "n_cycles_rottrials" ).def( 0 );
	option.add( my_options::culled_chains, "culled_chains" );
	option.add( my_options::extra_cycles, "extra_cycles" ).def( 0 );
	option.add( my_options::extra_temperature, "extra_temperature" ).def( 0.5 );
	option.add( my_options::protocols, "protocols" );
	option.add( my_options::params, "params" );
	option.add( my_options::fast_relax, "fast_relax" );
	option.add( my_options::zf_fast_relax_inner_cycles_factor, "zf_fast_relax_inner_cycles_factor" ).def(30);
	option.add( my_options::zf_fast_relax_outer_cycles, "zf_fast_relax_outer_cycles" ).def(3);
	option.add( my_options::ramping_cycles, "ramping_cycles" ).def(3);
	option.add( my_options::pack_mode, "pack_mode" ).def("pack_rotamers");
	option.add( my_options::min_type, "min_type" ).def("dfpmin_armijo_atol");
	option.add( my_options::fa_phosphate_constraints, "fa_phosphate_constraints" );
	option.add( my_options::h2o_ref_wt, "h2o_ref_wt" ).def( 0.0 );
	option.add( my_options::hydrate_phosphate, "hydrate_phosphate" ).def( false ); // was def true until 4/22/11
	option.add( my_options::hydrate_dna, "hydrate_dna" );
	option.add( my_options::hydrate_minor_groove, "hydrate_minor_groove" );
	option.add( my_options::allow_file_stealing, "allow_file_stealing" ).def( true );
	option.add( my_options::use_hb_env_dep, "use_hb_env_dep" ).def( true );
	option.add( my_options::use_hb_env_dep_DNA, "use_hb_env_dep_DNA" ); // default is FALSE !!!
	option.add( my_options::superimpose_decoys, "superimpose_decoys" );
	option.add( my_options::flex_protein_backbone, "flex_protein_backbone" );
	option.add( my_options::flex_protein_dna_jump, "flex_protein_dna_jump" );
	option.add( my_options::flex_dna, "flex_dna" );
	option.add( my_options::flex_dna_sugar, "flex_dna_sugar" );
	option.add( my_options::flex_dna_backbone, "flex_dna_backbone" );
	option.add( my_options::use_intra_dna_cp_crossover_4, "use_intra_dna_cp_crossover_4" ).def( true ); // make this deflt
	option.add( my_options::set_weight, "set_weight" );
	option.add( my_options::set_atom_property, "set_atom_property" );
	option.add( my_options::set_extra_atom_property, "set_extra_atom_property" );
	//  option.add( my_options::fullatom_dna_dihedral_weight, "fullatom_dna_dihedral_weight" ).def( 0.25 );
	//  option.add( my_options::fullatom_dna_dihedral_sugar_weight, "fullatom_dna_dihedral_sugar_weight" );
	//  option.add( my_options::fullatom_dna_chainbreak_weight, "fullatom_dna_chainbreak_weight" ).def( 1.0 );
	option.add( my_options::output_pdb_fraction, "output_pdb_fraction" ).def( 1.0 );
	option.add( my_options::refseq_file, "refseq_file" ).def( std::string() ); // default is empty string
	option.add( my_options::deriv_check, "deriv_check" );
	option.add( my_options::output_pdb_files, "output_pdb_files" );
	option.add( my_options::score_filter_acceptance_rate, "score_filter_acceptance_rate" );
	option.add( my_options::score_filter_pass_early, "score_filter_pass_early" );
	option.add( my_options::random_filter_acceptance_rate, "random_filter_acceptance_rate" ).def( 0.0 );
	//option.add( my_options::fa_scorefxn, "fa_scorefxn" );
	option.add( my_options::centroid_score_function, "centroid_score_function" );
	option.add( my_options::prefix, "prefix" );
	option.add( my_options::shared_output_tag, "shared_output_tag" );
	option.add( my_options::score_functions, "score_functions");
	option.add( my_options::cluster_fraction, "cluster_fraction" );
	option.add( my_options::cluster_fractions, "cluster_fractions" );
	option.add( my_options::profile_output_fraction, "profile_output_fraction" ).def( 0.05 );
	option.add( my_options::rescale_phosphate_charges, "rescale_phosphate_charges" );
	option.add( my_options::max_decoys_per_cluster_pdb, "max_decoys_per_cluster_pdb" ).def( 10 );
	option.add( my_options::cheating_frags_min_torsion_dev, "cheating_frags_min_torsion_dev" ).def(15); // was 30
	option.add( my_options::cheating_frags_max_torsion_dev, "cheating_frags_max_torsion_dev" ).def(90);
	option.add( my_options::full_search, "full_search" );
	option.add( my_options::autobuild_rwork_threshold, "autobuild_rwork_threshold" );
	option.add( my_options::skip_tncs, "skip_tncs" );
	option.add( my_options::sgalt_all, "sgalt_all" );
	option.add( my_options::sgalt_none, "sgalt_none" );
	option.add( my_options::generate_r_free_flags, "generate_r_free_flags" );
	option.add( my_options::mtz_file, "mtz_file" );
	option.add( my_options::num_models, "num_models" );
	option.add( my_options::composition, "composition" );
	option.add( my_options::llg_threshold_for_refinement, "llg_threshold_for_refinement" ).def( 160.0 );
	option.add( my_options::space_groups, "space_groups" );
	option.add( my_options::space_group, "space_group" );
	option.add( my_options::local_simfile, "local_simfile" );

	init_simtime();
}

void
start_decoy_timer()
{
	simtime::decoy_timer_start_time = time( NULL );
}

// in minutes
Real
check_decoy_timer()
{
	return ( time(NULL) - simtime::decoy_timer_start_time )/60.0;
}

void
check_simtime()
{
	using namespace simtime;

	// check for a STOP order
	if ( option[ my_options::shared_output_tag ].user() ) {
		string const stopfile( option[ my_options::shared_output_tag ]() + ".stop" );
		if ( utility::file::file_exists( stopfile ) ) exit(0);
	}

	if ( !option[ my_options::target_runtime ].user() ) return; // nothing to do

	//++nstruct_completed;
	Real const desired_runtime_in_minutes( option[ my_options::target_runtime ] );

	//clock_t const this_clock( clock() );
	time_t const this_time( time( NULL ) ); // in seconds

	Real const simtime      ( ( Real( this_time -  last_time )/60.0 )); // in minutes
	Real const total_runtime( ( Real( this_time - start_time )/60.0 )); // in minutes

	max_simtime = max( simtime, max_simtime );

	// exit if we'll be out of time
	std::cout << "check_simtime:" <<
		" simtime: " << simtime <<
		" total_runtime: " << total_runtime <<
		" max_simtime: " << max_simtime <<
		" desired_runtime_in_minutes: " << desired_runtime_in_minutes << std::endl;
	fflush( stdout );

	if ( total_runtime + max_simtime > desired_runtime_in_minutes ) exit(0);

	last_time = this_time;
}


bool
hacking()
{
	if ( option[ my_options::hacking ] ) {
		std::cout << "HACKING!" << std::endl;
		std::cerr << "HACKING!" << std::endl;
	}
	return option[ my_options::hacking ];
}

bool
dry_run()
{
	return option[ OK::out::dry_run ] || option[ OK::run::dry_run ];
}

bool
cmdline_verbose()
{
	return option[ my_options::verbose ];
}

bool
cpdb_input()
{
	return option[ my_options::cpdb_input ];
}

string
output_tag()
{
	return option[ OK::out::output_tag ];
}

string
shared_output_tag()
{
	return option[ my_options::shared_output_tag ];
}

/// note has a trailing /
string
shared_output_dir()
{
	return shared_output_tag()+"_output/";
}

/// note has a trailing /
string
shared_output_dir_create_if_necessary()
{
	string const outdir( shared_output_dir() );
	if ( !utility::file::file_exists( outdir ) ) utility::file::create_directory( outdir );
	return outdir;
}

Size
nstruct()
{
	return option[ OK::out::nstruct ];
}

Size
n_inner()
{
	return option[ OK::dna::specificity::n_inner ];
}

Size
n_outer()
{
	return option[ OK::dna::specificity::n_outer ];
}

std::map< std::string, utility::vector1< Real > >
get_params_from_command_line()
{
	using namespace utility;
	using namespace std;
	//using namespace options;
	using devel::blab::split_to_vector1;

	vector1< string > const allparams( basic::options::option[ my_options::params ] );
	map< string, vector1< Real > > params;

	for ( Size i=1; i<= allparams.size(); ++i ) {
		string const p( allparams[i] );
		Size const pos( p.find(':' ) );
		pbassert( pos != string::npos );
		string const tag( p.substr(0,pos) );
		vector1< string > const l( split_to_vector1( p.substr(pos+1), "," ) );
		for ( Size j=1; j<= l.size(); ++j ) {
			pbassert( is_float( l[j] ) );
			params[ tag ].push_back( float_of( l[j] ) );
		}
	}
	return params;

}
///////////////////////////////////////////////////////////////////////////////

void
set_scorefxn_use_hb_env_dep_DNA_from_commandline( ScoreFunction & scorefxn )
{
	scoring::methods::EnergyMethodOptions emoptions( scorefxn.energy_method_options() ); // make a local copy
	emoptions.hbond_options().use_hb_env_dep_DNA( option[ my_options::use_hb_env_dep_DNA ] );
	scorefxn.set_energy_method_options( emoptions );
}


///////////////////////////////////////////////////////////////////////////////
void
adjust_ref_weights_from_command_line(
	ScoreFunction & scorefxn
)
{
	if ( !option[ my_options::adjust_ref_weights ].user() ) return;
	scoring::methods::EnergyMethodOptions options( scorefxn.energy_method_options() );
	Reals ref_weights( options.method_weights( scoring::ref ) );
	runtime_assert( ref_weights.size() == 20 );
	strings const l( option[ my_options::adjust_ref_weights ]() );
	runtime_assert( l.size()%2 == 0 );
	Size const nwts( l.size()/2 );
	for ( Size ii=0; ii< nwts; ++ii ) {
		AA const aa( aa_from_oneletter_code( l[2*ii+1 ][0] ) );
		Real const adjustment( float_of( l[2*ii+2] ) );
		TR_OPTIONS.Trace << "Updating reference energy for " << aa << " oldval: " << F(9,3,ref_weights[aa]) <<
			" newval: " << F(9,3,ref_weights[aa] + adjustment) << endl;
		ref_weights[ aa ] += adjustment;
	}
	options.set_method_weights( scoring::ref, ref_weights );
	scorefxn.set_energy_method_options( options );
}


scoring::ScoreFunctionOP
setup_score_function( std::string const & wtsfile )
{
	bool const exclude_DNA_DNA( false );

	scoring::ScoreFunctionOP scorefxn(0);
	if ( conformation::symmetry::symmetry_is_on() ) scorefxn = ScoreFunctionOP( new scoring::symmetry::SymmetricScoreFunction() );
	else scorefxn = ScoreFunctionOP( new scoring::ScoreFunction() );

	scoring::methods::EnergyMethodOptions emoptions;
	emoptions.exclude_DNA_DNA( exclude_DNA_DNA );
	emoptions.use_intra_DNA_cp_crossover_4( option[ my_options::use_intra_dna_cp_crossover_4 ] );
	scorefxn->set_energy_method_options( emoptions );

	set_scorefxn_use_hb_env_dep_DNA_from_commandline( *scorefxn ); // NOTE NOTE NOTE
	//  scorefxn->set_energy_method_options
	//   ( scoring::methods::EnergyMethodOptions()
	//    .exclude_DNA_DNA( exclude_DNA_DNA )
	//    .use_hb_env_dep( option[ my_options::use_hb_env_dep ] )
	//    .use_intra_dna_cp_crossover_4( option[ my_options::use_intra_dna_cp_crossover_4 ] ) );

	scorefxn->add_weights_from_file( wtsfile );

	// this is the same horrible hack as in ScoreFunctionFactory.cc
	// if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() ) {
	//  scorefxn = new scoring::symmetry::SymmetricScoreFunction( scorefxn );
	// }

	return scorefxn;
}


scoring::ScoreFunctionOP
get_option_key_score_function( StringOptionKey const & key )
{
	if ( !( option[ key ].user() ) ) utility_exit_with_message("scorefxn key is not present on commandline!");
	return setup_score_function( option[ key ] );
}


// /// this will only do it once, even if called multiple times
// void
// rescale_phosphate_charges(
//  Real const scale_factor
// )
// {
//  using namespace chemical;

//  static bool init( false );

//  if ( init ) return; // already done it once
//  init = true;

//  utility::vector1< std::string > atoms;
//  atoms.push_back( "P" );
//  atoms.push_back( "O1P" );
//  atoms.push_back( "O2P" );
//  atoms.push_back( "O5'" );
//  atoms.push_back( "O3'" );

//  ResidueTypeSet & rsd_set( ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD ) );

//  for ( int i= first_DNA_aa; i<= last_DNA_aa; ++i ) {
//   ResidueTypeCAPs const rsd_types( rsd_set.aa_map( AA(i) ) );
//   for ( Size k=1; k<= rsd_types.size(); ++k ) {
//    ResidueType & rsd_type( rsd_set.nonconst_name_map( rsd_types[k]->name() ) );
//    assert( rsd_type.is_DNA() );
//    for ( Size j=1; j<= atoms.size(); ++j ) {
//     Real const old_charge( rsd_type.atomic_charge( rsd_type.atom_index( atoms[j] ) ) );
//     Real const new_charge( old_charge * scale_factor );
//     rsd_type.set_atomic_charge( atoms[j], new_charge );
//     if ( k==1 ) {
//      std::cout << "rescaling phosphate charge on atom " << rsd_type.name() << ' ' << atoms[j] << " from " <<
//       old_charge << " to " << new_charge << std::endl;
//     }
//    }
//   }
//  }
// }


// /// USES THE FA_STANDARD RESIDUE SET !!!!!!!!!!!
// void
// rescale_tp3_charges( Real const scale_factor )
// {
//  using namespace chemical;

//  static bool init( false );
//  static Real default_charge_O(0.0);
//  static Real default_charge_H(0.0);

//  ResidueTypeSet & rsd_set( ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD ) );
//  ResidueType & rsd_type( rsd_set.nonconst_name_map( "TP3" ) );

//  if ( !init ) {
//   init = true;
//   default_charge_O = rsd_type.atomic_charge( rsd_type.atom_index( "O" ) );
//   default_charge_H = rsd_type.atomic_charge( rsd_type.atom_index( "H1" ) );
//  }

//  rsd_type.set_atomic_charge( "O" , default_charge_O * scale_factor );
//  rsd_type.set_atomic_charge( "H1", default_charge_H * scale_factor );
//  rsd_type.set_atomic_charge( "H2", default_charge_H * scale_factor );

//  std::cout << "rescaling phosphate charge on TP3 atoms from " << default_charge_O << ' ' << default_charge_H <<
//   " to " << default_charge_O * scale_factor << ' ' << default_charge_H * scale_factor << std::endl;
// }



/// this will only do it once, even if called multiple times
void
set_atom_property(
	string const & atom_set_name,
	string const & atom_name,
	string const & param,
	Real const & setting
)
{
	using namespace chemical;

	TR_OPTIONS.Trace << "set_atom_property: " << atom_set_name << ' ' << atom_name << ' ' << param << ' ' << setting
		<< std::endl;

	AtomTypeSet & atom_set( ChemicalManager::get_instance()->nonconst_atom_type_set( atom_set_name ) );
	Size const atom_index( atom_set.atom_type_index( atom_name ) );

	/// I would like to uncomment the following if-check, but right now there is an extra parameter file
	/// that defines a parameter with the name LK_DGFREE (memb_fa_params.txt). That's kind of confusing...
	//  if ( atom_set.has_extra_parameter( param ) ) {
	//   Size const param_index( atom_set.extra_parameter_index( param ) );
	//   atom_set[ atom_index ].set_extra_parameter( param_index, setting );
	//  } else {
	{
		atom_set[ atom_index ].set_parameter( param, setting );
	}

}
/// this will only do it once, even if called multiple times
void
set_extra_atom_property(
	string const & atom_set_name,
	string const & atom_name,
	string const & param,
	Real const & setting
)
{
	using namespace chemical;

	TR_OPTIONS.Trace << "set_extra_atom_property: " << atom_set_name << ' ' << atom_name << ' ' << param << ' ' << setting
		<< std::endl;

	AtomTypeSet & atom_set( ChemicalManager::get_instance()->nonconst_atom_type_set( atom_set_name ) );
	Size const atom_index( atom_set.atom_type_index( atom_name ) );

	runtime_assert( atom_set.has_extra_parameter( param ) );
	Size const param_index( atom_set.extra_parameter_index( param ) );
	atom_set[ atom_index ].set_extra_parameter( param_index, setting );
}


// void
// set_atom_property_used_by_scorefxn(
//                   scoring::ScoreFunction const & scorefxn,
//                   string const & atom_set_name,
//                   string const & atom_name,
//                   string param,
//                   Real const & setting
//                   )
// {
//  TR_OPTIONS.Trace << "set_atom_property_used_by_scorefxn: " << atom_set_name << ' ' << atom_name << ' ' <<
//   param << ' ' << setting << std::endl;

//  using namespace scoring;

//  string const etable_type( scorefxn.energy_method_options().etable_type() );

//  if ( etable_type != FA_STANDARD_DEFAULT ) {
//   AtomTypeSet & atom_set( ChemicalManager::get_instance()->nonconst_atom_type_set( atom_set_name ) );
//   TR_OPTIONS.Trace << "Nonstandard etable: " << etable_type << std::endl;
//   pbassert( etable_type.substr( 0, FA_STANDARD_DEFAULT.size() ) == FA_STANDARD_DEFAULT );
//   string const suffix( etable_type.substr( FA_STANDARD_DEFAULT.size() ) );
//   if ( atom_set.has_extra_parameter( param + suffix ) ) {
//    TR_OPTIONS.Trace << "atom_set " << atom_set_name << " HAS   extra parameter: " << param+suffix << std::endl;
//    param += suffix;
//   } else {
//    TR_OPTIONS.Trace << "atom_set " << atom_set_name << " LACKS extra parameter: " << param+suffix << std::endl;
//   }
//  }
//  set_atom_property( atom_set_name, atom_name, param, setting );

// }

// void
// print_etable_parameters( scoring::ScoreFunction const & scorefxn, std::ostream & os )
// {
//  scoring::etable::Etable const & etable
//   ( *scoring::ScoringManager::get_instance()->etable( scorefxn.energy_method_options().etable_type() ) );
//  etable.print_atomtype_parameters( os );
// }


// void
// set_fa_standard_residue_type_charges_from_dipoles()
// {
//  static bool init( false );
//  if ( init ) return;
//  init = true;

//  ResidueTypeSet & rsd_set( ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD ) );

//  for ( chemical::ResidueTypeCOPs::const_iterator it = rsd_set.residue_types().begin();
//     it != rsd_set.residue_types().end(); ++it ) {
//   ResidueType & rsd_type( rsd_set.nonconst_name_map( (*it)->name() ) );

//   if ( rsd_type.is_protein() || rsd_type.is_DNA() ) {

//    scoring::methods::EnvResidueTypeInfo const rsd_info( rsd_type );

//    for ( Size i=1; i<= rsd_type.natoms(); ++i ) {
//     Real const old_charge( rsd_type.atom(i).charge() ), new_charge( rsd_info.default_charge(i));
//     if ( fabs( old_charge - new_charge )>=0.2 ) {
//      TR_OPTIONS.Trace << "set_fa_standard_residue_type_charges_from_dipoles: " <<
//       " old: " << F(9,3,rsd_type.atom(i).charge() ) <<
//       " new: " << F(9,3,rsd_info.default_charge(i)) << ' ' <<
//       rsd_type.name() << ' ' << rsd_type.atom_name(i) << ' ' << rsd_type.atom_type(i).name() << std::endl;
//     }
//     rsd_type.atom(i).charge( rsd_info.default_charge(i) );
//    }
//   }
//  }
// }


void
modify_atom_properties_from_command_line()
{
	using namespace std;
	using namespace basic::options;
	using namespace utility;

	static bool init( false );

	if ( init ) return; // already done it once
	init = true;

	if ( option[ my_options::set_atom_property ].user() ) {
		vector1< string > const & mods( option[ my_options::set_atom_property ] );
		for ( Size i=1; i<= mods.size(); ++i ) {
			///
			/// mod should look like:  "OOC:LK_RADIUS:4.5"
			///
			string const & mod( mods[i] );
			Size const pos1( mod.find(":") );
			if ( pos1 == string::npos ) utility_exit_with_message("bad format "+mod);
			string const atom_name( mod.substr(0,pos1) );
			Size const pos2( mod.substr(pos1+1).find(":") );
			if ( pos2 == string::npos ) utility_exit_with_message("bad format "+mod);
			string const param( mod.substr(pos1+1).substr(0,pos2) );
			Real const setting( float_of( mod.substr(pos1+1).substr(pos2+1) ) );
			set_atom_property( chemical::FA_STANDARD, atom_name, param, setting );
		}
	}

	if ( option[ my_options::set_extra_atom_property ].user() ) {
		vector1< string > const & mods( option[ my_options::set_extra_atom_property ] );
		for ( Size i=1; i<= mods.size(); ++i ) {
			///
			/// mod should look like:  "OOC:LK_RADIUS:4.5"
			///
			string const & mod( mods[i] );
			Size const pos1( mod.find(":") );
			if ( pos1 == string::npos ) utility_exit_with_message("bad format "+mod);
			string const atom_name( mod.substr(0,pos1) );
			Size const pos2( mod.substr(pos1+1).find(":") );
			if ( pos2 == string::npos ) utility_exit_with_message("bad format "+mod);
			string const param( mod.substr(pos1+1).substr(0,pos2) );
			Real const setting( float_of( mod.substr(pos1+1).substr(pos2+1) ) );
			set_extra_atom_property( chemical::FA_STANDARD, atom_name, param, setting );
		}
	}


	if ( basic::options::option[ basic::options::OptionKeys::dna::specificity::dipole_charges ] ) {
		// update the atomic charges using the dipoles array
		runtime_assert( false );
		//set_fa_standard_residue_type_charges_from_dipoles();
	}

}

void
set_h2o_ref_wt( Real const setting, ScoreFunction & scorefxn )
{

	vector1< core::Real > ref_wts( scorefxn.energy_method_options().method_weights( scoring::ref ) );
	ref_wts.resize( aa_h2o, 0.0 );
	ref_wts[ aa_h2o ] = setting;
	runtime_assert( Size( aa_vrt ) > ref_wts.size() || ref_wts[ aa_vrt ] == 0.0 );
	TR_OPTIONS.Trace << "set_h2o_ref_wt: " << setting << endl;
	scorefxn.set_method_weights( scoring::ref, ref_wts );

}

// void
// set_water_params(
//          Real const h2o_ref_wt,
//          Real const h2o_wdepth,
//          Real const h2o_dgfree,
//          Real const h2o_entropy_factor,
//          Real const h2o_hbond_scale_factor,
//          ScoreFunction & scorefxn
//          )
// {
//  std::cout << "set_water_params: " <<
//   " h2o_ref_wt: " << F(9,3,h2o_ref_wt) <<
//   " h2o_wdepth: " << F(9,3,h2o_wdepth) <<
//   " h2o_dgfree: " << F(9,3,h2o_dgfree) <<
//   " h2o_entropy_factor: " << F(9,3,h2o_entropy_factor) <<
//   " h2o_hbond_scale_factor: " << F(9,3,h2o_hbond_scale_factor) << std::endl;

//  set_h2o_ref_wt( h2o_ref_wt, scorefxn );
//  set_atom_property_used_by_scorefxn( scorefxn, FA_STANDARD, "OW", "LJ_WDEPTH", h2o_wdepth );
//  set_atom_property_used_by_scorefxn( scorefxn, FA_STANDARD, "O1", "LJ_WDEPTH", h2o_wdepth );
//  set_atom_property_used_by_scorefxn( scorefxn, FA_STANDARD, "OW", "LK_DGFREE", h2o_dgfree );
//  set_atom_property_used_by_scorefxn( scorefxn, FA_STANDARD, "O1", "LK_DGFREE", h2o_dgfree );

//  if ( option[ my_options::tp3_radius ].user() ) {
//   Real const rad( option[ my_options::tp3_radius ] );
//   set_atom_property_used_by_scorefxn( scorefxn, FA_STANDARD, "OW", "LJ_RADIUS", rad );
//   set_atom_property_used_by_scorefxn( scorefxn, FA_STANDARD, "O1", "LJ_RADIUS", rad );
//  }

//  option[ OK::dna::specificity::h2o_entropy_factor ].value( h2o_entropy_factor );
//  option[ OK::dna::specificity::h2o_hbond_scale_factor ].value( h2o_hbond_scale_factor );
//  if ( option[ my_options::rescale_tp3_charges_as_hbonds ] ) {
//   rescale_tp3_charges( h2o_hbond_scale_factor );
//  }
//  ScoringManager::get_instance()->nonconst_etable_ref
//   ( scorefxn.energy_method_options().etable_type() ).reinitialize();
//  print_etable_parameters( scorefxn, std::cout );
// }


void
modify_weights_from_command_line( ScoreFunction & scorefxn )
{
	using namespace scoring;
	if ( option[ my_options::set_weight ].user() ) {
		strings const opt_settings( option[ my_options::set_weight ] );
		strings settings; // do this so we can handle negative values
		foreach_ ( string tag, opt_settings ) {
			strings const l( split_to_vector1( tag ) );
			foreach_ ( string ltag, l ) {
				settings.push_back( ltag );
			}
		}
		pbassert( settings.size()%2==0 );
		for ( Size i=0; i< settings.size()/2; ++i ) {
			ScoreType const t( score_type_from_name( settings[ 2*i + 1] ) );
			Real const value( float_of( settings[ 2*i + 2 ] ) );
			TR_OPTIONS.Trace << "setting scorefxn weight from command line: " << t << ' ' << value << std::endl;
			scorefxn.set_weight( t, value );
		}
	}

	/// water reference weight:
	if ( option[ my_options::h2o_ref_wt ].user() ) {
		Real const h2o_ref_wt( option[ my_options::h2o_ref_wt ] );
		std::cout << "setting h2o_ref_wt from command line: " << h2o_ref_wt << endl;
		set_h2o_ref_wt( h2o_ref_wt, scorefxn );
	}

	//  if ( option[ my_options::tp3_params ].user() ) {
	//   utility::vector1< Real > const params( option[ my_options::tp3_params ]() );
	//   if ( params.size() != 5 ) {
	//    utility_exit_with_message("-tp3_params h2o_ref_wt h2o_wdepth h2o_dgfree h2o_entropy_factor h2o_hbscalefactor");
	//   }
	//   set_water_params( params[1], params[2], params[3], params[4], params[5], scorefxn );
	//  }

	// fiddle with the phosphate charges -- not sure if this should be here
	// note that inside this function is a check to see if its been done
	//
	//  if ( option[ my_options::rescale_phosphate_charges ].user() ) {
	//   rescale_phosphate_charges( option[ my_options::rescale_phosphate_charges ] );
	//  }

	//  if ( option[ my_options::at_bias_correction ].user() ) {
	//   Real const deltaE( option[ my_options::at_bias_correction ] );
	//   vector1< core::Real > ref_wts( scorefxn.energy_method_options().method_weights( dna_ref ) );
	//   if ( ref_wts.size() != 16 ) utility_exit_with_message("bad dna ref wts?");

	//   for ( Size i=1; i<= 4; ++i ) {
	//    for ( Size j=1; j<= 4; ++j ) {
	//     Size const n_ats( ( i==1 || i == 4 ) + ( j==1 || j == 4 ) );
	//     Size const k( (i-1)*4+j );
	//     std::cout << "Updating dna ref wt: " << i << ' ' << j << ' ' << n_ats << ' ' << deltaE*n_ats << ' ' << k << std::endl;
	//     ref_wts[ k ] = ref_wts[k] + deltaE * n_ats;
	//    }
	//   }

	//   scorefxn.set_method_weights( dna_ref, ref_wts );

	//  }
}

scoring::ScoreFunctionOP
get_score_function_from_command_line()
{
	ScoreFunctionOP scorefxn;
	//  if ( options::option[ my_options::fa_scorefxn ].user() ) {
	//   scorefxn = get_option_key_score_function( my_options::fa_scorefxn );
	//  } else {
	scorefxn = get_option_key_score_function( OK::dna::specificity::score_function );

	modify_weights_from_command_line( *scorefxn );

	modify_atom_properties_from_command_line();

	return scorefxn;
}

scoring::ScoreFunctionOP
get_score_function_from_file( string const & filename )
{
	ScoreFunctionOP scorefxn;
	scorefxn = setup_score_function( filename );

	modify_weights_from_command_line( *scorefxn );

	modify_atom_properties_from_command_line();

	return scorefxn;
}

scoring::ScoreFunctionOP
get_centroid_score_function_from_command_line()
{
	return get_option_key_score_function( my_options::centroid_score_function );
}

void
parse_clustering_params(
	Real & min_threshold,
	Real & max_threshold,
	Size & min_cluster_size,
	Size & min_top_cluster_size,
	Size & try_top_cluster_size,
	Size & max_top_cluster_size,
	Size & max_clusters,
	Size & max_decoys_per_cluster_pdbfile
)
{
	if ( !option[ my_options::clustering_params ].user() || option[ my_options::clustering_params].size() != 8 ) {
		cout << "-clustering_params min_threshold max_threshold " <<
			" min_cluster_size min_top_cluster_size try_top_cluster_size max_top_cluster_size " <<
			" max_clusters max_decoys_per_cluster_pdbfile" << endl;
		utility_exit();
	}

	Reals const clustering_params( option[ my_options::clustering_params ] );

	runtime_assert( clustering_params.size() == 8 );

	min_threshold = clustering_params[1];
	max_threshold = clustering_params[2];
	/// now casting to ints
	min_cluster_size = clustering_params[3];
	min_top_cluster_size = clustering_params[4];
	try_top_cluster_size = clustering_params[5];
	max_top_cluster_size = clustering_params[6];
	max_clusters = clustering_params[7];
	max_decoys_per_cluster_pdbfile = clustering_params[8];
}


#endif
