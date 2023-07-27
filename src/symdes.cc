// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/pilot/phil/test1.cc
/// @brief  Some simple examples of how to use basic functionality + some DNA functionality
/// @author Phil Bradley (pbradley@fhcrc.org)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// libRosetta headers
#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/phil_options.hh>
#include <apps/pilot/phil/phil_io.hh>

#include <apps/pilot/phil/pdb.hh>
#include <apps/pilot/phil/symdes.hh>
#include <apps/pilot/phil/chains.hh>
#include <apps/pilot/phil/fragments.hh>
#include <apps/pilot/phil/star_rebuild.hh>
#include <apps/pilot/phil/sym_centroid.hh>
#include <apps/pilot/phil/stub_frag.hh>
#include <apps/pilot/phil/interface.hh>
#include <apps/pilot/phil/topology.hh>
#include <apps/pilot/phil/cluster.hh>
#include <apps/pilot/phil/helix_frags.hh>
#include <apps/pilot/phil/helix_turns.hh>
#include <apps/pilot/phil/phenix.hh>

#include <core/scoring/constraints/DihedralPairConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/CircularPowerFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <devel/blab/loops/util.hh>
#include <devel/blab/cluster/simple_cluster.tmpl.hh>
#include <core/io/silent/SilentFileData.hh>
//#include <core/io/pdb/pose_io.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/powell.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <core/pose/datacache/CacheableDataType.hh>

//#include <core/id/AtomID.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/pack/rtmin.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/rigid_body_moves.hh>
#include <protocols/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>

#include <protocols/sasa_scores/sasapack.hh>

#include <devel/blab/classic_frags/TorsionFragment.hh>
#include <devel/dna/util.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <core/pose/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/VariantType.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <devel/init.hh>

#include <numeric/util.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <utility/io/izstream.hh>

//////// option key includes
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/phil.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

using numeric::random::random_permutation;
//using numeric::random::RG;

//static numeric::random::RandomGenerator RG(12321); // <- Magic number, do not change it!!!

static basic::Tracer TR( "apps.pilot.phil.symdes" );

// string const relax_scorefxn_weights_tag( "talaris2013" ); // should put this in places where we have score12prime now
// //string const design_scorefxn_weights_tag( "talaris2013" );
// string const soft_design_scorefxn_weights_tag( "talaris2013_soft" );

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
using devel::blab::classic_frags::FragLib;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using numeric::constants::d::pi;
using numeric::square;
using basic::subtract_degree_angles;
using basic::periodic_range;


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


namespace my_options {
StringOptionKey target_hand("my:target_hand");
StringOptionKey phenixseq("my:phenixseq");
//StringOptionKey phenix_trunc("my:phenix_trunc");
StringOptionKey centroid_repeatseq_aas("my:centroid_repeatseq_aas");
StringOptionKey xray_data_tag("my:xray_data_tag");
StringOptionKey refoldseq("my:refoldseq");
StringOptionKey sasapack_datafile("my:sasapack_datafile");
StringOptionKey symm_type("my:symm_type");
StringOptionKey pdb_energies_datafile("my:pdb_energies_datafile");
StringOptionKey fragfile3("my:fragfile3");
StringOptionKey fragfile9("my:fragfile9");
StringOptionKey reconstruct_bb_weights_file("my:reconstruct_bb_weights_file" );
StringOptionKey template_pdb("my:template_pdb");
StringVectorOptionKey template_pdbs("my:template_pdbs");
StringVectorOptionKey turn_pairs("my:turn_pairs");
StringVectorOptionKey tree_colors("my:tree_colors");
StringVectorOptionKey pose1_pdbs("my:pose1_pdbs");
StringVectorOptionKey pose2_pdbs("my:pose2_pdbs");
StringVectorOptionKey symm_types("my:symm_types");
StringOptionKey interface_transforms_file("my:interface_transforms_file");
BooleanOptionKey trefoil_mode("my:trefoil_mode");
BooleanOptionKey pentafoil_mode("my:pentafoil_mode");
BooleanOptionKey ok_filter("my:ok_filter");
BooleanOptionKey centroid_only("my:centroid_only");
BooleanOptionKey debug_derivs("my:debug_derivs");
BooleanOptionKey ignore_n2cdist("my:ignore_n2cdist");
BooleanOptionKey fold_first_chain("my:fold_first_chain");
// BooleanOptionKey force_beta_braid("my:force_beta_braid");
BooleanOptionKey centroid_jumping("my:centroid_jumping");
BooleanOptionKey force_TIM("my:force_TIM");
BooleanOptionKey use_fullatom_twistrise_energy("my:use_fullatom_twistrise_energy");
BooleanOptionKey free_mode("my:free_mode");
BooleanOptionKey tal_mode("my:tal_mode");
BooleanOptionKey force_triangles("my:force_triangles");
BooleanOptionKey force_biangles("my:force_biangles");
BooleanOptionKey pdb_list_io("my:pdb_list_io");
BooleanOptionKey output_all_ca_coords("my:output_all_ca_coords");
BooleanOptionKey restrict_to_bab("my:restrict_to_bab");
BooleanOptionKey pick_decoys("my:pick_decoys");
BooleanOptionKey helix1_is_a_strand("my:helix1_is_a_strand");
BooleanOptionKey helix2_is_a_strand("my:helix2_is_a_strand");
BooleanOptionKey interpolate_helixlens("my:interpolate_helixlens");
IntegerVectorOptionKey polar_pore("my:polar_pore");
IntegerVectorOptionKey src_helixlens("my:src_helixlens");
IntegerVectorOptionKey orientations("my:orientations");
IntegerOptionKey nrepeat_refold("my:nrepeat_refold");
IntegerOptionKey max_turn_length("my:max_turn_length");
IntegerOptionKey phase_fragment("my:phase_fragment");
BooleanOptionKey force_symmetry("my:force_symmetry");
BooleanOptionKey planar("my:planar");
BooleanOptionKey self_dock_2c("my:self_dock_2c");
BooleanOptionKey break_helix_using_proline_frags("my:break_helix_using_proline_frags");
BooleanOptionKey dont_force_first_helix_in("my:dont_force_first_helix_in");
BooleanOptionKey use_asymmetric_centroid_scoring("my:use_asymmetric_centroid_scoring");
BooleanOptionKey insert_frag_sequence("my:insert_frag_sequence");
BooleanOptionKey skip_bump_check("my:skip_bump_check");
//   BooleanOptionKey use_softrep_for_early_design("my:use_softrep_for_early_design" );
//   BooleanOptionKey use_softrep_for_design("my:use_softrep_for_design" );
BooleanOptionKey break_longer_helix("my:break_longer_helix");
BooleanOptionKey make_interface_moves("my:make_interface_moves");
BooleanOptionKey randomize_aana("my:randomize_aana");
RealVectorOptionKey anchor_position_range("my:anchor_position_range");
RealVectorOptionKey helix1_len_range("my:helix1_len_range");
RealOptionKey relax_fraction("my:relax_fraction");
RealOptionKey max_ala_fraction("my:max_ala_fraction");
RealOptionKey max_sasapack_score("my:max_sasapack_score");
RealOptionKey min_depth("my:min_depth");
RealOptionKey max_helix_angle("my:max_helix_angle");
RealOptionKey star_fraction("my:star_fraction");
RealOptionKey fragq_threshold_for_refolding("my:fragq_threshold_for_refolding");
RealOptionKey centroid_bsasa14_threshold("my:centroid_bsasa14_threshold");
RealOptionKey sometimes_break_longer_helix("my:sometimes_break_longer_helix");
RealOptionKey donut_energy_weight("my:donut_energy_weight");
RealOptionKey target_twist_tpr("my:target_twist_tpr");
RealOptionKey target_rise_tpr("my:target_rise_tpr");
RealOptionKey resample_centroid_seqs_fragweight("my:resample_centroid_seqs_fragweight");
BooleanOptionKey no_ss_wiggle("my:no_ss_wiggle");
BooleanOptionKey reversed("my:reversed");
BooleanOptionKey use_tal_fragments("my:use_tal_fragments");
BooleanOptionKey tpr_mode("my:tpr_mode");
BooleanOptionKey barrel_mode("my:barrel_mode");
BooleanOptionKey freeze_anchorseq("my:freeze_anchorseq");
BooleanOptionKey freeze_helixseq("my:freeze_helixseq");
BooleanOptionKey freeze_centroidseq("my:freeze_centroidseq");
BooleanOptionKey nodesign("my:nodesign");
//BooleanOptionKey fastrelax_before_design("my:fastrelax_before_design");
BooleanOptionKey fastrelax_after_refold("my:fastrelax_after_refold");
BooleanOptionKey randomize_centroid_sequence("my:randomize_centroid_sequence");
// IntegerVectorOptionKey nbraids("my:nbraids" );
IntegerOptionKey helixlen_delta("my:helixlen_delta" );
IntegerVectorOptionKey helixlen_deltas("my:helixlen_deltas" );
IntegerOptionKey ntrim("my:ntrim" );
IntegerOptionKey ntrim_for_rmsd("my:ntrim_for_rmsd" );
IntegerOptionKey ctrim_for_rmsd("my:ctrim_for_rmsd" );
IntegerOptionKey num_chains_to_delete("my:num_chains_to_delete" );
IntegerVectorOptionKey ntrims("my:ntrims" );
IntegerVectorOptionKey ctrims("my:ctrims" );
IntegerOptionKey num_simfiles("my:num_simfiles" );
IntegerOptionKey final_num_models("my:final_num_models");
// IntegerOptionKey nrepeat_per_braid_refold("my:nrepeat_per_braid_refold" );
IntegerOptionKey max_repeatlen("my:max_repeatlen" );
IntegerOptionKey ntrim_sses("my:ntrim_sses" );
IntegerOptionKey ctrim_sses("my:ctrim_sses" );
IntegerOptionKey max_helixlen_delta("my:max_helixlen_delta" );
IntegerOptionKey n_refold("my:n_refold" );
IntegerOptionKey min_nsims("my:min_nsims" );
IntegerOptionKey min_min_nsim("my:min_min_nsim" );
IntegerOptionKey n_refold_sspred("my:n_refold_sspred" );
IntegerOptionKey design_cycles("my:design_cycles" );
IntegerOptionKey nrepeat("my:nrepeat" );
IntegerOptionKey nrepeats_to_fold("my:nrepeats_to_fold" );
IntegerOptionKey nres_to_fold("my:nres_to_fold" );
IntegerOptionKey base_repeat("my:base_repeat" );
IntegerOptionKey nrepeat_sim("my:nrepeat_sim" );
//StringVectorOptionKey adjust_ref_weights("my:adjust_ref_weights" );
StringVectorOptionKey centroid_helixseq("my:centroid_helixseq" );
StringVectorOptionKey turns("my:turns" );
StringVectorOptionKey turn1s("my:turn1s" );
StringVectorOptionKey turn2s("my:turn2s" );
StringVectorOptionKey turn3s("my:turn3s" );
StringVectorOptionKey anchorpos_loop_bbs("my:anchorpos_loop_bbs" );
StringVectorOptionKey resample_bbs("my:resample_bbs" );
StringOptionKey resample_params("my:resample_params" );
StringVectorOptionKey resample_centroid_seqs("my:resample_centroid_seqs" );
StringVectorOptionKey force_sequence_positions("my:force_sequence_positions" );
IntegerVectorOptionKey repeatlens("my:repeatlens" );
IntegerVectorOptionKey helix_lens("my:helix_lens" );
IntegerVectorOptionKey helix1_lens("my:helix1_lens" );
IntegerVectorOptionKey helix2_lens("my:helix2_lens" );
IntegerVectorOptionKey helix3_lens("my:helix3_lens" );
IntegerVectorOptionKey nrepeats("my:nrepeats" );
IntegerVectorOptionKey anchor_positions("my:anchor_positions" );
IntegerVectorOptionKey anchorpos_loop_lengths("my:anchorpos_loop_lengths" );
IntegerVectorOptionKey anchorpos_loop_begins("my:anchorpos_loop_begins" );
IntegerVectorOptionKey anchorpos_loop_ends("my:anchorpos_loop_ends" );
IntegerVectorOptionKey cutpoint_loop_lengths("my:cutpoint_loop_lengths" );
IntegerVectorOptionKey cutpoint_loop_begins("my:cutpoint_loop_begins" );
IntegerVectorOptionKey cutpoint_loop_ends("my:cutpoint_loop_ends" );
IntegerVectorOptionKey design_positions("my:design_positions" );

}

void
add_my_options()
{
	add_phil_options();
 	option.add( my_options::ntrim, "ntrim" );
 	option.add( my_options::phenixseq, "phenixseq" );
 	//option.add( my_options::phenix_trunc, "phenix_trunc" );
 	option.add( my_options::ntrim_for_rmsd, "ntrim_for_rmsd" ).def(2);
 	option.add( my_options::ctrim_for_rmsd, "ctrim_for_rmsd" ).def(5);
 	option.add( my_options::max_ala_fraction, "max_ala_fraction" );
 	option.add( my_options::max_sasapack_score, "max_sasapack_score" );
 	option.add( my_options::centroid_only, "centroid_only" );
 	option.add( my_options::design_positions, "design_positions" );
 	option.add( my_options::turn_pairs, "turn_pairs" );
	option.add( my_options::src_helixlens, "src_helixlens" );
	option.add( my_options::min_depth, "min_depth" );
	option.add( my_options::max_helix_angle, "max_helix_angle" );
	option.add( my_options::trefoil_mode, "trefoil_mode" );
	option.add( my_options::pentafoil_mode, "pentafoil_mode" );
	option.add( my_options::centroid_repeatseq_aas, "centroid_repeatseq_aas" ).def("GGGDESTNVLLAMK");
	option.add( my_options::ok_filter, "ok_filter" );
	option.add( my_options::debug_derivs, "debug_derivs" );
	option.add( my_options::max_turn_length, "max_turn_length" );
	option.add( my_options::nrepeat_refold, "nrepeat_refold" );
	option.add( my_options::ignore_n2cdist, "ignore_n2cdist" );
	option.add( my_options::phase_fragment, "phase_fragment" );
	option.add( my_options::refoldseq, "refoldseq" );
	option.add( my_options::xray_data_tag, "xray_data_tag" );
	option.add( my_options::num_chains_to_delete, "num_chains_to_delete" );
	option.add( my_options::ntrims, "ntrims" );
	option.add( my_options::ctrims, "ctrims" );
	option.add( my_options::final_num_models, "final_num_models" );
	option.add( my_options::num_simfiles, "num_simfiles" ).def(25);
	option.add( my_options::relax_fraction, "relax_fraction" ).def(0);
	option.add( my_options::star_fraction, "star_fraction" ).def(0);
	option.add( my_options::fold_first_chain, "fold_first_chain" );
	// option.add( my_options::force_beta_braid, "force_beta_braid" );
	option.add( my_options::centroid_jumping, "centroid_jumping" );
	option.add( my_options::force_TIM, "force_TIM" );
	// option.add( my_options::nrepeat_per_braid_refold, "nrepeat_per_braid_refold" );
	option.add( my_options::tree_colors, "tree_colors" );
	option.add( my_options::orientations, "orientations" );
	option.add( my_options::fragq_threshold_for_refolding, "fragq_threshold_for_refolding" ).def(1e-6);
	// option.add( my_options::nbraids, "nbraids" );
	option.add( my_options::use_fullatom_twistrise_energy, "use_fullatom_twistrise_energy" );
	option.add( my_options::max_repeatlen, "max_repeatlen" );
	option.add( my_options::tal_mode, "tal_mode" );
	option.add( my_options::free_mode, "free_mode" );
	option.add( my_options::force_triangles, "force_triangles" );
	option.add( my_options::force_biangles, "force_biangles" );
	option.add( my_options::pdb_list_io, "pdb_list_io" );
	option.add( my_options::ntrim_sses, "ntrim_sses" );
	option.add( my_options::ctrim_sses, "trim_sses" );
	option.add( my_options::output_all_ca_coords, "output_all_ca_coords" );
	option.add( my_options::restrict_to_bab, "restrict_to_bab" ).def( true );
	option.add( my_options::pick_decoys, "pick_decoys" );
	option.add( my_options::helix1_is_a_strand, "helix1_is_a_strand" );
	option.add( my_options::helix2_is_a_strand, "helix2_is_a_strand" );
	option.add( my_options::interpolate_helixlens, "interpolate_helixlens" );
	option.add( my_options::helixlen_delta, "helixlen_delta" ).def(0);
	option.add( my_options::helixlen_deltas, "helixlen_deltas" );
	option.add( my_options::max_helixlen_delta, "max_helixlen_delta" ).def(0);
	option.add( my_options::target_hand, "target_hand" );
	option.add( my_options::polar_pore, "polar_pore" );
	option.add( my_options::reconstruct_bb_weights_file, "reconstruct_bb_weights_file" );
	option.add( my_options::force_symmetry, "force_symmetry" );
	option.add( my_options::min_nsims, "min_nsims" ).def(3);
	option.add( my_options::min_min_nsim, "min_min_nsim" ).def(11);
	option.add( my_options::helix_lens, "helix_lens" );
	option.add( my_options::helix1_lens, "helix1_lens" );
	option.add( my_options::helix2_lens, "helix2_lens" );
	option.add( my_options::helix3_lens, "helix3_lens" );
	option.add( my_options::planar, "planar" );
	option.add( my_options::self_dock_2c, "self_dock_2c" );
	option.add( my_options::centroid_bsasa14_threshold, "centroid_bsasa14_threshold" );
	option.add( my_options::pose1_pdbs, "pose1_pdbs" );
	option.add( my_options::pose2_pdbs, "pose2_pdbs" );
	option.add( my_options::symm_types, "symm_types" );
	option.add( my_options::template_pdb, "template_pdb" );
	option.add( my_options::template_pdbs, "template_pdbs" );
	option.add( my_options::symm_type, "symm_type" );
	option.add( my_options::break_helix_using_proline_frags, "break_helix_using_proline_frags" );
	option.add( my_options::anchorpos_loop_bbs, "anchorpos_loop_bbs" );
	option.add( my_options::force_sequence_positions, "force_sequence_positions" );
	option.add( my_options::target_twist_tpr, "target_twist_tpr" ).def( 48.58 ); // mean values for 1n0a (regan tpr)
	option.add( my_options::target_rise_tpr, "target_rise_tpr" ).def( 9.74 ); // mean values for 1n0a (regan tpr)
	option.add( my_options::resample_centroid_seqs_fragweight, "resample_centroid_seqs_fragweight" ).def( 1.0 );
	option.add( my_options::nodesign, "nodesign" );
	option.add( my_options::tpr_mode, "tpr_mode" );
	option.add( my_options::donut_energy_weight, "donut_energy_weight" ).def(0.0);
	option.add( my_options::dont_force_first_helix_in, "dont_force_first_helix_in" );
	option.add( my_options::resample_bbs, "resample_bbs" );
	option.add( my_options::resample_params, "resample_params" );
	option.add( my_options::resample_centroid_seqs, "resample_centroid_seqs" );
	option.add( my_options::use_asymmetric_centroid_scoring, "use_asymmetric_centroid_scoring" );
	option.add( my_options::centroid_helixseq, "centroid_helixseq" );
	option.add( my_options::barrel_mode, "barrel_mode" );
	option.add( my_options::helix1_len_range, "helix1_len_range" );
	option.add( my_options::turns, "turns" );
	option.add( my_options::turn1s, "turn1s" );
	option.add( my_options::turn2s, "turn2s" );
	option.add( my_options::turn3s, "turn3s" );
	option.add( my_options::no_ss_wiggle, "no_ss_wiggle" );
	option.add( my_options::nrepeat, "nrepeat" );
	option.add( my_options::design_cycles, "design_cycles" );
	option.add( my_options::base_repeat, "base_repeat" );
	option.add( my_options::nrepeat_sim, "nrepeat_sim" );
	option.add( my_options::insert_frag_sequence, "insert_frag_sequence" );
	option.add( my_options::skip_bump_check, "skip_bump_check" );
	//  option.add( my_options::use_softrep_for_design, "use_softrep_for_design" );
	//  option.add( my_options::use_softrep_for_early_design, "use_softrep_for_early_design" );
	//option.add( my_options::adjust_ref_weights, "adjust_ref_weights" );
	option.add( my_options::n_refold, "n_refold" ).def(0);
	option.add( my_options::n_refold_sspred, "n_refold_sspred" ).def(0);
	option.add( my_options::interface_transforms_file, "interface_transforms_file" );
	option.add( my_options::randomize_aana, "randomize_aana" );
	option.add( my_options::make_interface_moves, "make_interface_moves" );
	option.add( my_options::anchor_position_range, "anchor_position_range" );
	option.add( my_options::repeatlens, "repeatlens" );
	option.add( my_options::nrepeats, "nrepeats" );
	option.add( my_options::nrepeats_to_fold, "nrepeats_to_fold" );
	option.add( my_options::nres_to_fold, "nres_to_fold" );
	option.add( my_options::sometimes_break_longer_helix, "sometimes_break_longer_helix" ).def(0.0); // prob of doing so
	option.add( my_options::break_longer_helix, "break_longer_helix" );
	option.add( my_options::sasapack_datafile, "sasapack_datafile" );
	option.add( my_options::pdb_energies_datafile, "pdb_energies_datafile" );
	option.add( my_options::fastrelax_after_refold, "fastrelax_after_refold" );
	//  option.add( my_options::fastrelax_before_design, "fastrelax_before_design" );
	option.add( my_options::use_tal_fragments, "use_tal_fragments" );
	option.add( my_options::randomize_centroid_sequence, "randomize_centroid_sequence" );
	option.add( my_options::freeze_anchorseq, "freeze_anchorseq" );
	option.add( my_options::freeze_helixseq, "freeze_helixseq" );
	option.add( my_options::freeze_centroidseq, "freeze_centroidseq" );
	option.add( my_options::anchor_positions, "anchor_positions" );
	option.add( my_options::anchorpos_loop_lengths, "anchorpos_loop_lengths" );
	option.add( my_options::anchorpos_loop_begins, "anchorpos_loop_begins" );
	option.add( my_options::anchorpos_loop_ends, "anchorpos_loop_ends" );
	option.add( my_options::cutpoint_loop_lengths, "cutpoint_loop_lengths" );
	option.add( my_options::cutpoint_loop_begins, "cutpoint_loop_begins" );
	option.add( my_options::cutpoint_loop_ends, "cutpoint_loop_ends" );
	option.add( my_options::reversed, "reversed" );
	option.add( my_options::fragfile3, "fragfile3" );
	option.add( my_options::fragfile9, "fragfile9" );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
Real
get_total_rotation_degrees(
	Vectors const & coords,
	Vector const & axis_center,
	Vector const & axis_vector
)
{
	Real total_rotation(0);
	runtime_assert( fabs( axis_vector.length()-1)<1e-3 );
	// Vectors coords_projected;
	// foreach_( Vector const & v, coords ) {
	//  coords_projected.push_back( v + axis_vector * axis_vector.dot( axis_center - v ) );
	// }
	for ( Size i=1; i< coords.size(); ++i ) {
		total_rotation += numeric::dihedral_degrees( coords[i], axis_center, axis_center+axis_vector, coords[i+1] );
		// runtime_assert( fabs( axis_vector.dot( coords_projected[i] - axis_center ) )<1e-3 );
		// total_rotation += numeric::angle_radians( coords_projected[i], axis_center, coords_projected[i+1] );
	}
	return total_rotation;
	// return numeric::conversions::degrees( total_rotation );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// struct HelixStubFrag {
//  string id;
//  RT rt;
//  Vectors coords;
//  Real raw_strain;
//  Real norm_strain;
//  HelixParams hparams;
// };


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef NOCOMPILE

class GeneralizedHelixBarrelMultifunc : public optimization::Multifunc {
public:
	typedef optimization::Multivec Multivec;

	GeneralizedHelixBarrelMultifunc():
		nrepeat_(0),
		params_scale_(10)
	{}

	virtual
	Real
	operator()( Multivec const & wts ) const;


	virtual
	void
	dfunc( Multivec const & phipsi, Multivec & dE_dphipsi ) const
 	{
 		dE_dphipsi.clear(); dE_dphipsi.resize( phipsi.size(), 0.0 );
 	}

	void
	set_nrepeat_and_frags_and_helix_lengths(
																					Size const nrepeat,
																					StubFrags const & helix_frags,
																					StubFrags const & turn_frags,
																					Sizes const & helixlens // not really necessary?
																					)
	{
		nrepeat_ = nrepeat;
		turn_frags_ = turn_frags;
		helixlens_ = helixlens;
		params_from_helix_frags( helix_frags, target_params_ );
	}

	void
	params_from_helix_frags( StubFrags const & helix_frags, Multivec & params ) const;

	void
	params_from_helix_params( HelixParams const & h1params, HelixParams const & h2params, Multivec & params ) const;

	void
	helix_frags_from_params( Multivec const & params, StubFrag & h1,  StubFrag & h2 ) const;

private:

	void
	unpack_params( Multivec const & params ) const;

private:
	Size nrepeat_;
	StubFrag t1_, t2_;
	Size h1len_, h2len_;
	Multivec target_params_;
	mutable HelixParams h1params_, h2params_;
	Real params_scale_;

};




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
GeneralizedHelixBarrelMultifunc::unpack_params( Multivec const & params ) const
{
	h1params_.rise           = params[1]/params_scale_;
	h1params_.twist          = params[2]/params_scale_;
	h1params_.tilt           = params[3]/params_scale_;
	h1params_.tilt_direction = params[4]/params_scale_;
	h1params_.ca_distance    = params[5]/params_scale_;

	h2params_.rise           = params[6]/params_scale_;
	h2params_.twist          = params[7]/params_scale_;
	h2params_.tilt           = params[8]/params_scale_;
	h2params_.tilt_direction = params[9]/params_scale_;
	h2params_.ca_distance    = params[10]/params_scale_;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
GeneralizedHelixBarrelMultifunc::operator()( Multivec const & params ) const
{
	runtime_assert( nrepeat_ );

	// unpack the multivec params into the helical params arrays h1params_ and h2params_ (they are mutable)
	unpack_params( params );

	///
	Stub const start_stub;

	RT const h1_rt( start_stub, generate_helix_stub( start_stub, h1params_, h1len_ ) );
	RT const h2_rt( start_stub, generate_helix_stub( start_stub, h2params_, h2len_ ) );

	// TR.Trace << "h1_rt_trans: " <<
	// 	F(9,3,h1_rt.get_translation().x()) <<
	// 	F(9,3,h1_rt.get_translation().y()) <<
	// 	F(9,3,h1_rt.get_translation().z()) << endl;

	// TR.Trace << "h2_rt_trans: " <<
	// 	F(9,3,h2_rt.get_translation().x()) <<
	// 	F(9,3,h2_rt.get_translation().y()) <<
	// 	F(9,3,h2_rt.get_translation().z()) << endl;

	Matrix const
		&R1( h1_rt.get_rotation() ), &R2( t1_.rt.get_rotation() ),
		&R3( h2_rt.get_rotation() ), &R4( t2_.rt.get_rotation() );

	Vector const
		&t1( h1_rt.get_translation() ), &t2( t1_.rt.get_translation() ),
		&t3( h2_rt.get_translation() ), &t4( t2_.rt.get_translation() );

	/// what happens to start_stub.v as we go through a single repeat?
	Matrix Rtot( R1 );
	Vector ttot( t1 );

	ttot += Rtot * t2; Rtot.right_multiply_by( R2 );
	ttot += Rtot * t3; Rtot.right_multiply_by( R3 );
	ttot += Rtot * t4; Rtot.right_multiply_by( R4 );


	Real theta, theta_dev, rise;
	{
		Vector center, n, t;
		get_stub_transform_data( start_stub, Stub( Rtot, ttot ), center, n, t, theta );
		theta = numeric::conversions::degrees( theta );
		theta_dev = numeric::square( basic::subtract_degree_angles( theta, 360.0 / nrepeat_ ) );
		//theta_dev = fabs( basic::subtract_degree_angles( theta, 360.0 / nrepeat_ ) );
		rise = n.dot(t);
	}

	Vector vfinal( start_stub.v ); // is 0.0
	for ( Size n=1; n<= nrepeat_; ++n ) {
		vfinal += ttot;
		ttot = Rtot * ttot;
		//vfinal = Rtot * vfinal + ttot;
	}


	if ( false ) { // DEBUGGING
		Stub stub( start_stub );
		for ( Size n=1; n<= nrepeat_; ++n ) {
			stub = t2_.rt.make_jump( h2_rt.make_jump( t1_.rt.make_jump( h1_rt.make_jump( stub ) ) ) );
		}
		Vector recompute_vfinal( stub.v );
		TR.Trace << "recompute_vfinal_distance: " << F(9,3,vfinal.distance( recompute_vfinal ) ) << endl;
		vfinal = recompute_vfinal;
	}

	// constraint term
	Real cstE( 0.0 );
	Real const cstwt( 1000.0 );
	static Reals const cst_weights( make_vector1( 1.0, // rise              Angstroms
																								1.0, // twist             radians
																								1.0, // tilt              radians
																								0.1, // tilt_direction    radians
																								1.0, // ca_distance       Angstroms
																								1.0, // rise
																								1.0, // twist
																								1.0, // tilt
																								0.1, // tilt_direction
																								1.0 // ca_distance
																								));
	for ( Size i=1; i<= 10; ++i ) {
		cstE += cst_weights[i] * numeric::square( ( params[i] - target_params_[i] )/params_scale_ );
	}

	Real const vfinal_dis2( vfinal.distance_squared( start_stub.v ) );

	//Real const func( vfinal_dis2 + cstwt * cstE + theta_dev );
	Real const func( cstwt * cstE + theta_dev + numeric::square( nrepeat_ * rise ) );

	TR.Trace << "barrel_func " << F(9,3,func) <<
		" vfdis: " << F(9,3,sqrt( vfinal_dis2 ) ) <<
		" theta: " << F(9,3,theta ) <<
		" rise: " << F(9,3,rise ) <<
		" wtdcstE: " << F(9,3,cstwt*cstE ) <<
		" theta_dev: " << F(9,3,theta_dev ) << endl;
	TR.Trace << "barrel_params " << F(9,3,func);
	for ( Size i=1; i<= 10; ++i ) TR.Trace << ' ' << F(6,3,params[i]/params_scale_);
	TR.Trace << endl;

	return func;

}

////////////////////////////////////////////
void
GeneralizedHelixBarrelMultifunc::params_from_helix_frags(
																							StubFrag const & h1,
																							StubFrag const & h2,
																							Multivec & params
																							) const
{
	params_from_helix_params( h1.hparams, h2.hparams, params );
}

////////////////////////////////////////////
////////////////////////////////////////////
void
GeneralizedHelixBarrelMultifunc::params_from_helix_params(
																							 HelixParams const & h1params,
																							 HelixParams const & h2params,
																							 Multivec & params
																							 ) const
{
	params.resize( 10 );
	params[1] = params_scale_ * h1params.rise;
	params[2] = params_scale_ * h1params.twist;
	params[3] = params_scale_ * h1params.tilt;
	params[4] = params_scale_ * h1params.tilt_direction;
	params[5] = params_scale_ * h1params.ca_distance;

	params[6] = params_scale_ * h2params.rise;
	params[7] = params_scale_ * h2params.twist;
	params[8] = params_scale_ * h2params.tilt;
	params[9] = params_scale_ * h2params.tilt_direction;
	params[10]= params_scale_ * h2params.ca_distance;

}

////////////////////////////////////////////
void
GeneralizedHelixBarrelMultifunc::helix_frags_from_params(
																							Multivec const & params,
																							StubFrag & h1,
																							StubFrag & h2
																							) const
{
	unpack_params( params );

	Stub const start_stub; // default
	Stub const h1_stop_stub( generate_helix_coords( start_stub, h1params_, h1len_, h1.coords ) );
	h1.rt = RT( start_stub, h1_stop_stub );
	h1.hparams = h1params_; // set by unpack_params

	Stub const h2_stop_stub( generate_helix_coords( start_stub, h2params_, h2len_, h2.coords ) );
	h2.rt = RT( start_stub, h2_stop_stub );
	h2.hparams = h2params_; // set by unpack_params

}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
get_axis_axis_lambdas(
	Vector const & center1,
	Vector const & axis1,
	Vector const & center2,
	Vector const & axis2,
	Real & lambda1,
	Real & lambda2
)
{

	//Real const angle( numeric::conversions::degrees( std::acos( axis1.dot( axis2 ) ) ) );
	Real const v12( axis1.dot( axis2 ) ), v11( axis1.dot(axis1)), v22( axis2.dot(axis2)),
		c11( center1.dot( axis1 ) ), c12( center1.dot( axis2) ), c21( center2.dot( axis1 ) ), c22( center2.dot( axis2) ),
		a( c21 - c11 ), b( c22 - c12 );
	lambda1 = ( a * v22 - b * v12 )/( v11 * v22 - v12 * v12 );
	lambda2 = ( a * v12 - b * v11 )/( -1*( v12 * v12 - v22 * v11 ) );

	// {
	//  Vector const p1( center1 + lambda1 * axis1 ), p2( center2 + l2 * axis2 );
	//  runtime_assert( fabs( (p1-p2).dot( axis1 ) ) < 1e-2 );
	//  runtime_assert( fabs( (p1-p2).dot( axis2 ) ) < 1e-2 );
	// }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Real
get_axis_distance_squared(
	Vector const & center1,
	Vector const & axis1,
	Vector const & center2,
	Vector const & axis2
)
{
	//bool const debug( true );
	Real lambda1, lambda2;

	get_axis_axis_lambdas( center1, axis1, center2, axis2, lambda1, lambda2 );

	return ( center1 + lambda1 * axis1 ).distance_squared( center2 + lambda2 * axis2 );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
get_helix_helix_distance_squared_numeric(
	Vector const & p1,
	Vector const & q1,
	Vector const & p2,
	Vector const & q2
)
{
	Size const nsteps( 100 );
	Vector const step1( (q1-p1)/nsteps ), step2( (q2-p2)/nsteps );
	Real mindis2( 1e6 );
	for ( Size i=0; i<= nsteps; ++i ) {
		for ( Size j=0; j<= nsteps; ++j ) {
			mindis2 = min( mindis2, ( p1 + i * step1 ).distance_squared( p2 + j * step2 ) );
		}
	}
	return mindis2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// note: axis should point from p to q
///

Real
get_helix_helix_distance_squared(
	Vector const & axis1,
	Vector const & p1,
	Vector const & q1,
	Vector const & axis2,
	Vector const & p2,
	Vector const & q2
)
{
	//Real lambda1, lambda2;
	Real const
		p1dot1( p1.dot( axis1 ) ), p1dot2( p1.dot( axis2 ) ),
		q1dot1( q1.dot( axis1 ) ), q1dot2( q1.dot( axis2 ) ),
		p2dot1( p2.dot( axis1 ) ), p2dot2( p2.dot( axis2 ) ),
		q2dot1( q2.dot( axis1 ) ), q2dot2( q2.dot( axis2 ) );

	runtime_assert( p1dot1 < q1dot1 );
	runtime_assert( p2dot2 < q2dot2 );

	Size seg1pos, seg2pos;

	if      ( max( p1dot2, q1dot2 ) < p2dot2 ) seg2pos = 1; // closest point in seg2 will  be p2
	else if ( min( p1dot2, q1dot2 ) > q2dot2 ) seg2pos = 3; // closest point in seg2 will  be q2
	else                                       seg2pos = 2; // closest point in seg2 *may* be in the middle somewhere

	if      ( max( p2dot1, q2dot1 ) < p1dot1 ) seg1pos = 1; // closest point in seg1 will  be p1
	else if ( min( p2dot1, q2dot1 ) > q1dot1 ) seg1pos = 3; // closest point in seg1 will  be q1
	else                                       seg1pos = 2; // closest point in seg1 *may* be in the middle somewhere

	// edge cases
	if      ( seg1pos == 1 && seg2pos == 1 ) return p1.distance_squared( p2 );
	else if ( seg1pos == 1 && seg2pos == 2 ) return p1.distance_squared( p2 + max( Real(0),min(q2dot2,p1dot2)-p2dot2)*axis2 );
	else if ( seg1pos == 1 && seg2pos == 3 ) return p1.distance_squared( q2 );
	else if ( seg1pos == 3 && seg2pos == 1 ) return q1.distance_squared( p2 );
	else if ( seg1pos == 3 && seg2pos == 2 ) return q1.distance_squared( p2 + max( Real(0),min(q2dot2,q1dot2)-p2dot2)*axis2 );
	else if ( seg1pos == 3 && seg2pos == 3 ) return q1.distance_squared( q2 );
	else if ( seg1pos == 2 && seg2pos == 1 ) return p2.distance_squared( p1 + max( Real(0),min(q1dot1,p2dot1)-p1dot1)*axis1 );
	else if ( seg1pos == 2 && seg2pos == 3 ) return q2.distance_squared( p1 + max( Real(0),min(q1dot1,q2dot1)-p1dot1)*axis1 );


	/// call p1 "center1" and p2 "center2"
	Real const v12( axis1.dot( axis2 ) ),// v11( axis1.dot(axis1)), v22( axis2.dot(axis2)),
		c11( p1dot1 ), c12( p1dot2 ), c21( p2dot1 ), c22( p2dot2 ),
		a( c21 - c11 ), b( c22 - c12 ), d( 1 - v12 * v12 );

	if ( d<1e-3 ) { // almost parallel
		// just choose one of the edges
		Real const
			dp1( p1.distance_squared( p2 + max( Real(0),min(q2dot2,p1dot2)-p2dot2)*axis2 ) ),
			dq1( q1.distance_squared( p2 + max( Real(0),min(q2dot2,q1dot2)-p2dot2)*axis2 ) ),
			dp2( p2.distance_squared( p1 + max( Real(0),min(q1dot1,p2dot1)-p1dot1)*axis1 ) ),
			dq2( q2.distance_squared( p1 + max( Real(0),min(q1dot1,q2dot1)-p1dot1)*axis1 ) );
		return min( dp1, min( dq1, min( dp2, dq2 ) ) );

	} else {
		Real const lambda1 = ( a - b * v12 )/ d;
		Real const lambda2 = ( a * v12 - b )/ d;
		// Real const lambda1 = ( a * v22 - b * v12 )/( v11 * v22 - v12 * v12 );
		// Real const lambda2 = ( a * v12 - b * v11 )/( -1*( v12 * v12 - v22 * v11 ) );

		if ( false ) { // DEBUGGING
			Vector const v1( p1 + lambda1 * axis1 ), v2( p2 + lambda2 * axis2 );
			runtime_assert( fabs( (v1-v2).dot( axis1 ) ) < 1e-2 );
			runtime_assert( fabs( (v1-v2).dot( axis2 ) ) < 1e-2 );
		}


		if ( lambda2 < 0 || lambda2 > q2dot2 - p2dot2 ||
				lambda1 < 0 || lambda1 > q1dot1 - p1dot1 ) {
			// need to be careful
			Real const
				dp1( p1.distance_squared( p2 + max( Real(0),min(q2dot2,p1dot2)-p2dot2)*axis2 ) ),
				dq1( q1.distance_squared( p2 + max( Real(0),min(q2dot2,q1dot2)-p2dot2)*axis2 ) ),
				dp2( p2.distance_squared( p1 + max( Real(0),min(q1dot1,p2dot1)-p1dot1)*axis1 ) ),
				dq2( q2.distance_squared( p1 + max( Real(0),min(q1dot1,q2dot1)-p1dot1)*axis1 ) );
			return min( dp1, min( dq1, min( dp2, dq2 ) ) );

		} else {
			return ( p1 + lambda1 * axis1 ).distance_squared( p2 + lambda2 * axis2 );
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct HelixDat {
	Vector axis;
	Vector begin;
	Vector end;
};
typedef vector1< HelixDat > HelixDats;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
get_helix_helix_distance_squared(
	HelixDat const & hd1,
	HelixDat const & hd2
)
{
	return get_helix_helix_distance_squared( hd1.axis, hd1.begin, hd1.end, hd2.axis, hd2.begin, hd2.end );
}


///////////////////////////////////////////////////////////////////////////////


void
get_symm_type_vertices(
	char const symm_type,
	Vectors & vertices,
	Vectors & nbr_vertices
)
{
	vertices.clear();
	nbr_vertices.clear();

	if ( symm_type == 'W' ) {
		// read from special pdb file
		string const filename("/home/pbradley/csdat/trap_cage/ring_centers.pdb");
		vertices = read_CA_coords_from_file( filename );
		Vector centroid(0,0,0);
		for ( Vector v : vertices ) centroid += v;
		centroid /= vertices.size();
		for ( Vector & v : vertices ) v -= centroid;
		for ( Vector & v : vertices ) v.normalize();

	} else if ( symm_type == 'O' ) {
		for ( int i=-1; i<= 1; i+= 2 ) {
			for ( int j=-1; j<= 1; j+= 2 ) {
				for ( int k=-1; k<= 1; k+= 2 ) {
					vertices.push_back( Vector( i, j, k ) );
				}
			}
		}

	} else if ( symm_type == 'T' ) {
		for ( int i=-1; i<= 1; i+= 2 ) {
			vertices.push_back( Vector( i, 0, -1.0 / sqrt(2) ) );
			vertices.push_back( Vector( 0, i,  1.0 / sqrt(2) ) );
		}

	} else if ( symm_type == 'I' ) {
		/// vertices of cube:
		for ( int i=-1; i<= 1; i+= 2 ) {
			for ( int j=-1; j<= 1; j+= 2 ) {
				for ( int k=-1; k<= 1; k+= 2 ) {
					vertices.push_back( Vector( i, j, k ) );
				}
			}
		}
		Real const golden( 0.5 * ( 1 + sqrt(5) ) ), inv_golden( 1.0/ golden );
		for ( int i=-1; i<= 1; i+= 2 ) {
			for ( int j=-1; j<= 1; j+= 2 ) {
				vertices.push_back( Vector( 0, i * inv_golden, j * golden ) );
				vertices.push_back( Vector( i * inv_golden, j * golden, 0 ) );
				vertices.push_back( Vector( i * golden, 0, j * inv_golden ) );
			}
		}
	} else {
		utility_exit_with_message("unrecognized symm_type: "+symm_type );
	}


	/// we are only going to simulate two of the monomers
	/// get the angle between the symmetry axes from


	Size n_monomers( vertices.size() );

	/// get nbr relationships among vertices
	nbr_vertices.clear(); nbr_vertices.resize( n_monomers );
	Size nbr_count1(0);
	for ( Size i=1; i<= n_monomers; ++i ) {
		Vector const & vi( vertices[i] );
		Real mindis2( 1e6 ), epsilon( 1e-2 );
		Size nbr_count(0);
		for ( Size j=1; j<= n_monomers; ++j ) {
			if ( j==i ) continue;
			Vector const & vj( vertices[j] );
			Real const dis2( vi.distance_squared(vj) );
			if ( dis2<mindis2-epsilon ) {
				mindis2 = dis2;
				nbr_vertices[i] = vj;
				nbr_count = 1;
			} else if ( dis2<mindis2+epsilon ) {
				++nbr_count;
			}
		}
		TR.Trace << "nbr_count: " << i << ' ' << nbr_count << F(9,3,sqrt(mindis2)) << endl;
		if ( i==1 ) {
			nbr_count1 = nbr_count;
		} else {
			runtime_assert( nbr_count == nbr_count1 );
		}
	}

}




////////////////////////////////////////////////
void
get_symm_type_vertices(
	char const symm_type,
	Vectors & vertices,
	vector1< Sizes > & nbr_vertices
)
{
	vertices.clear();

	if ( symm_type == 'O' ) {
		for ( int i=-1; i<= 1; i+= 2 ) {
			for ( int j=-1; j<= 1; j+= 2 ) {
				for ( int k=-1; k<= 1; k+= 2 ) {
					vertices.push_back( Vector( i, j, k ) );
				}
			}
		}

	} else if ( symm_type == 'T' ) {
		for ( int i=-1; i<= 1; i+= 2 ) {
			vertices.push_back( Vector( i, 0, -1.0 / sqrt(2) ) );
			vertices.push_back( Vector( 0, i,  1.0 / sqrt(2) ) );
		}

	} else if ( symm_type == 'I' ) {
		/// vertices of cube:
		for ( int i=-1; i<= 1; i+= 2 ) {
			for ( int j=-1; j<= 1; j+= 2 ) {
				for ( int k=-1; k<= 1; k+= 2 ) {
					vertices.push_back( Vector( i, j, k ) );
				}
			}
		}
		Real const golden( 0.5 * ( 1 + sqrt(5) ) ), inv_golden( 1.0/ golden );
		for ( int i=-1; i<= 1; i+= 2 ) {
			for ( int j=-1; j<= 1; j+= 2 ) {
				vertices.push_back( Vector( 0, i * inv_golden, j * golden ) );
				vertices.push_back( Vector( i * inv_golden, j * golden, 0 ) );
				vertices.push_back( Vector( i * golden, 0, j * inv_golden ) );
			}
		}
	} else {
		utility_exit_with_message("unrecognized symm_type: "+symm_type );
	}


	/// we are only going to simulate two of the monomers
	/// get the angle between the symmetry axes from


	Size n_monomers( vertices.size() );

	/// get nbr relationships among vertices
	nbr_vertices.clear(); nbr_vertices.resize( n_monomers );
	Size nbr_count(0);
	for ( Size i=1; i<= n_monomers; ++i ) {
		Vector const & vi( vertices[i] );
		Real mindis2( 1e6 ), epsilon( 1e-2 );
		for ( Size j=1; j<= n_monomers; ++j ) {
			if ( j==i ) continue;
			Vector const & vj( vertices[j] );
			Real const dis2( vi.distance_squared(vj) );
			if ( dis2<mindis2-epsilon ) {
				mindis2 = dis2;
				nbr_vertices[i] = make_vector1( j );
			} else if ( dis2<mindis2+epsilon ) {
				nbr_vertices[i].push_back( j );
			}
		}
		TR.Trace << "nbr_count: " << i << ' ' << nbr_vertices[i].size() << F(9,3,sqrt(mindis2)) << endl;
		if ( i==1 ) {
			nbr_count = nbr_vertices[i].size();
		} else {
			runtime_assert( nbr_count == nbr_vertices[i].size() );
		}
	}

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// this is pretty slow...
//

// void
// get_segment_window_stubs(
//              SegmentType const & segtype,
//              Vectors const & coords,
//              Size const segbegin,
//              Stub & start_stub,
//              Stub & stop_stub
//              )
// {
//  using namespace optimization;
//  /// from helical_params_test
//  Real const
//   ideal_rise_helix( 1.494 ),
//   ideal_twist_helix( numeric::conversions::radians( 99.236 ) ),
//   ideal_radius_helix( 2.292 ),
//   ideal_rise_strand( 3.26 ), // these guys from strand_test output
//   ideal_twist_strand( numeric::conversions::radians( 196.0 ) ),
//   ideal_radius_strand( 0.98 );

//  //Size const helix_window( 7 );
//  // fit ideal helix to the coords, return starting stub
//  static HelicalParamsFitMultifunc func;
//  func.set_target_coords( segbegin, segbegin + get_segment_window( segtype ) -1, coords );
//  func.optimize_tilt( true );

//  Multivec params( 5 );
//  if ( segtype == alpha_type ) {
//   params[1] = ideal_rise_helix;
//   params[2] = ideal_twist_helix;
//   params[3] = 0;
//   params[4] = 0;
//   params[5] = ideal_radius_helix;
//  } else {
//   runtime_assert( segtype == beta_type );
//   params[1] = ideal_rise_strand;
//   params[2] = ideal_twist_strand;
//   params[3] = 0;
//   params[4] = 0;
//   params[5] = ideal_radius_strand;
//  }

//  Real const tolerance( 1e-3 );
//  Size iterations;
//  Real final_func_value;
//  optimization::powell( params, func, tolerance, iterations, final_func_value );

//  func.get_helix_stubs_for_params( params, start_stub, stop_stub );

// }






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
setup_helix_dat(
	Size const helix_begin,
	Size const helix_end,
	Vectors const & coords,
	HelixDat & hd
)
{
	static Size const helix_window( 7 );
	runtime_assert( helix_end - helix_begin + 1 >= helix_window );

	Stub const begin_stub( get_helix_stub( coords, helix_begin ) ), end_stub( get_helix_stub(coords,helix_end-helix_window+1));
	Vector const & begin_axis( begin_stub.M.col_x() ), & end_axis( end_stub.M.col_x() ),
		&begin_center( begin_stub.v ), &end_center( end_stub.v );
	hd.begin = begin_center + ( coords[ helix_begin ] - begin_center ).dot( begin_axis ) * begin_axis;
	hd.end   =   end_center + ( coords[ helix_end   ] -   end_center ).dot(   end_axis ) *   end_axis;
	hd.axis = ( hd.end - hd.begin ).normalized();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
build_helix_library(
	Sizes const & helixlens,
	Size const max_frags,
	map< Size, StubFrags > & all_frags
)
{
	static bool init( false );

	Size const helix_window( 7 );

	static vector1< std::pair< string, Vectors > > all_helix_coords;

	if ( !init ) { //read info on bb coords
		init = true;
		string const pdb_coords_file( option[ my_options::pdb_coords_file ] );
		ifstream data( pdb_coords_file.c_str() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l[1] == "helix_coords:" ) {
				Size const helixlen( int_of( l[2] ) );
				Vectors coords;
				for ( Size i=1; i<= helixlen; ++i ) {
					coords.push_back( Vector( float_of( l[ 3*i+1 ] ), float_of( l[ 3*i+2 ] ), float_of( l[3*i+3 ] ) ) );
				}
				string const tag( filebase( l.back() )+"_"+ l[ l.size()-1 ]+"_"+l[ l.size()-3 ] );
				all_helix_coords.push_back( make_pair( tag, coords ) );
			}
		}
		data.close();
		numeric::random::random_permutation( all_helix_coords, numeric::random::rg() );
		TR.Trace << "Read " << all_helix_coords.size() << " helix coords from file " << pdb_coords_file << endl;
	} // initialize the coordinates /////////////////////////////////////////////////////////////////////

	all_frags.clear();
	for ( Sizes::const_iterator hlen= helixlens.begin(); hlen != helixlens.end(); ++hlen ) {
		Size const helixlen( *hlen );
		all_frags[ helixlen ];
		StubFrags & frags( all_frags.find( helixlen )->second );
		Reals all_strains;
		for ( Size ii=1; ii<= all_helix_coords.size(); ++ii ) {
			Vectors const & pdb_helix_coords( all_helix_coords[ii].second );
			Size const pdb_helixlen( pdb_helix_coords.size() );
			if ( pdb_helixlen >= helixlen ) {
				for ( Size i=1; i<= pdb_helixlen-helixlen+1; ++i ) { // pdb helix start pos
					Size const stub1_helix_begin( i ), stub2_helix_begin( i + helixlen - helix_window );
					Stub const
						stub1( get_helix_stub( pdb_helix_coords, stub1_helix_begin ) ),
						stub2( get_helix_stub( pdb_helix_coords, stub2_helix_begin ) );
					if ( false && helixlen > helix_window ) { // look at offset between stubs
						Vector const v( stub1.global2local( stub2.v ) );
						Matrix const R( stub2.M * stub1.M.transposed() );
						Real theta;
						Vector n( numeric::rotation_axis( R, theta ) );
						if ( n.dot( stub1.M.col_x() ) < 0 ) { n *= -1; theta *= -1; }
						Size const nres_shift( helixlen-helix_window );
						// v.x shift is roughly 1.5A /res
						//
						TR.Trace << "helix_stub_delta: nres_shift: " << nres_shift << " v-shift: " <<
							F(9,3,v.x()/nres_shift) <<
							F(9,3,v.y()/nres_shift) <<
							F(9,3,v.z()/nres_shift) <<
							" theta: " << F(9,3,basic::periodic_range( numeric::conversions::degrees( theta ), 360.0 )) <<
							" xtheta: " << F(9,3,basic::periodic_range( 99.0 * nres_shift, 360.0 )) << endl;
					}
					StubFrag f;
					f.id = all_helix_coords[ii].first+"_I"+string_of(i)+"_L"+string_of(helixlen);
					f.rt = RT( stub1, stub2 );
					for ( Size j=1; j<= helixlen; ++j ) {
						f.coords.push_back( stub1.global2local( pdb_helix_coords[i+j-1] ) );
					}
					if ( false ) {
						Real avg_rmsd, med_rmsd;
						compute_helix_strain( 1, helixlen, f.coords, avg_rmsd, med_rmsd );
						f.raw_strain = avg_rmsd;
						all_strains.push_back( f.raw_strain );
					}
					frags.push_back( f );
					if ( frags.size() >= max_frags ) break;
				} // loop over windows in this pdb helix
			}
			if ( frags.size() >= max_frags ) break;
		} // loop over pdb helices


		if ( false ) {
			Size const nfrags( frags.size() );
			for ( Size i=1; i<= nfrags; ++i ) {
				StubFrag & f( frags[i] );
				Size nbetter(0);
				Real const epsilon( 1e-3 );
				for ( Reals::const_iterator s= all_strains.begin(); s!= all_strains.end(); ++s ) {
					if ( *s + epsilon < f.raw_strain ) ++nbetter;
				}

				f.norm_strain = Real( nbetter );
				// nbetter will range from 0 to nfrags-1
				if ( nfrags > 1 ) f.norm_strain /= ( nfrags-1 ); // now range in [0,1]
				//TR.Trace << "norm_strain_helix: " << helixlen << F(9,3,f.raw_strain) << F(9,3,f.norm_strain) << endl;
			}
		}
		TR.Trace << "build_helix_library:: Created " << frags.size() << " helix frags of length " << helixlen << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// struct HelixParamsLine {
//  HelixParams hparams;
//  string tag;
//  Size helixlen;
// };
// typedef vector1< HelixParamsLine > HelixParamsLines;

// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// void
// build_helix_library_v2(
//             Sizes const & helixlens,
//             Size const max_frags,
//             map< Size, HelixStubFrags > & all_frags
//             )
// {
//  using numeric::conversions::radians;
//  static bool init( false );

//  //Size const helix_window( 7 );

//  static HelixParamsLines all_helix_params_lines;

//  if ( !init ) { //read info on helical params
//   init = true;
//   string const pdb_coords_file( option[ my_options::helix_params_file ] );
//   ifstream data( pdb_coords_file.c_str() );
//   string line;
//   HelixParamsLine hpline;
//   while ( getline( data, line ) ) {
//    strings const l( split_to_vector1( line ) );
//    if ( l[1] == "final_params_helix_full_tilt" ) {
//     hpline.hparams.rise           =          float_of( l[ 7] );
//     hpline.hparams.twist          = radians( float_of( l[ 9] ) );
//     hpline.hparams.tilt           = radians( float_of( l[11] ) );
//     hpline.hparams.tilt_direction = radians( float_of( l[13] ) );
//     hpline.hparams.ca_distance    =          float_of( l[15] );
//     hpline.tag = l.back();
//     hpline.helixlen = int_of( l[3] );
//     all_helix_params_lines.push_back( hpline );
//    }
//   }
//   data.close();
//   numeric::random::random_permutation( all_helix_params_lines, numeric::random::RG );
//   TR.Trace << "Read " << all_helix_params_lines.size() << " helix lines from file " << pdb_coords_file << endl;
//  } // initialize the coordinates /////////////////////////////////////////////////////////////////////

//  all_frags.clear();
//  Vectors coords;
//  for ( Sizes::const_iterator hlen= helixlens.begin(); hlen != helixlens.end(); ++hlen ) {
//   Size const helixlen( *hlen );
//   coords.reserve( helixlen );
//   all_frags[ helixlen ];
//   HelixStubFrags & frags( all_frags.find( helixlen )->second );
//   Reals all_strains;
//   for ( Size ii=1; ii<= all_helix_params_lines.size(); ++ii ) {
//    HelixParamsLine const & hpline( all_helix_params_lines[ii] );
//    if ( hpline.helixlen > Real(helixlen)*1.25 ||
//       hpline.helixlen < Real(helixlen)*0.8 ) continue;

//    // reconstruct some coords
//    HelixStubFrag f;
//    Stub start_stub; // default
//    Stub const stop_stub( generate_helix_coords( start_stub, hpline.hparams, helixlen, f.coords ) );
//    f.id = hpline.tag;
//    f.rt = RT( start_stub, stop_stub );
//    f.hparams = hpline.hparams;
//    frags.push_back( f );
//    if ( frags.size() >= max_frags ) break;
//   } // loop over pdb helices

//   TR.Trace << "build_helix_library_v2:: Created " << frags.size() << " helix frags of length " << helixlen << endl;
//  }
// }


void
build_turn_library(
	strings const & turns,
	map< string, StubFrags > & all_frags
)
{
	static bool init( false );

	static map< string, vector1< std::pair< string, Vectors > > > all_turn_coords;

	Size const turn_buffer( 7 ); // should be the same as helix_window

	if ( !init ) { //read info on bb coords
		init = true;
		string const pdb_coords_file( option[ my_options::pdb_coords_file ] );
		ifstream data( pdb_coords_file.c_str() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l[1] == "turn_coords:" ) {
				string const turn( l[2] );
				if ( all_turn_coords.find( turn ) == all_turn_coords.end() ) all_turn_coords[ turn ];
				Vectors coords;
				Size const coordslen( int_of( l[6] ) );
				runtime_assert( coordslen == turn.size() + 2*turn_buffer );
				for ( Size i=1; i<= coordslen; ++i ) {
					coords.push_back( Vector( float_of( l[ 3*i+4 ] ), float_of( l[ 3*i+5 ] ), float_of( l[3*i+6 ] ) ) );
				}
				string const turntag( filebase( l.back() ) + "_" + l[ l.size()-2 ] + "_" + l[ l.size()-1 ] +"_"+turn );
				all_turn_coords.find( turn )->second.push_back( make_pair( turntag, coords ) );
				//TR.Trace << "read turn: " << turn << ' ' << all_turn_coords.find( turn )->second.size() << endl;
			}

		}
		data.close();
	} // initialize the coordinates /////////////////////////////////////////////////////////////////////


	all_frags.clear();
	for ( strings::const_iterator turn= turns.begin(); turn != turns.end(); ++turn ) {
		if ( all_turn_coords.find( *turn ) == all_turn_coords.end() ) {
			utility_exit_with_message("no PDB turn frags for turn "+(*turn));
		}
		all_frags[ *turn ];
		Size const turnlen( turn->size() );
		StubFrags & frags( all_frags.find( *turn )->second );
		vector1< std::pair< string, Vectors > > const & pdb_turn_coords( all_turn_coords.find( *turn )->second );
		Reals all_strains;
		for ( Size ii=1; ii<= pdb_turn_coords.size(); ++ii ) {
			Vectors const & turn_coords( pdb_turn_coords[ii].second );
			runtime_assert( turn_coords.size() == turnlen + 2*turn_buffer );
			Size const stub1_helix_begin( 1 ), stub2_helix_begin( turn_buffer+turnlen+1 );
			Stub const
				stub1( get_helix_stub( turn_coords, stub1_helix_begin ) ),
				stub2( get_helix_stub( turn_coords, stub2_helix_begin ) );
			StubFrag f;
			f.id = pdb_turn_coords[ii].first;
			f.rt = RT( stub1, stub2 );
			for ( Size j=1; j<= turnlen; ++j ) f.coords.push_back( stub1.global2local( turn_coords[turn_buffer+j] ));
			Real avg_rmsd, med_rmsd;
			compute_turn_strain( *turn, turn_buffer+1, turn_coords, avg_rmsd, med_rmsd );
			f.raw_strain = avg_rmsd;
			all_strains.push_back( f.raw_strain );
			frags.push_back( f );
			//if ( frags.size() >= max_frags ) break;
		} // loop over pdb turns of this type

		Size const nfrags( frags.size() );

		for ( Size i=1; i<= nfrags; ++i ) {
			StubFrag & f( frags[i] );
			Size nbetter(0);
			Real const epsilon( 1e-3 );
			for ( Reals::const_iterator s= all_strains.begin(); s!= all_strains.end(); ++s ) {
				if ( *s + epsilon < f.raw_strain ) ++nbetter;
			}

			f.norm_strain = Real( nbetter );
			// nbetter will range from 0 to nfrags-1
			if ( nfrags > 1 ) f.norm_strain /= ( nfrags-1 ); // now range in [0,1]
			//TR.Trace << "norm_strain_turn: " << *turn << F(9,3,f.raw_strain) << F(9,3,f.norm_strain) << endl;
		}


		TR.Trace << "build_turn_library:: Read " << frags.size() << " turns of type " << *turn << endl;
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// struct TurnCoordsLine {
//  string tag;
//  Vectors coords;
//  vector1< Reals > torsions;
//  string sequence;
// };



// void
// build_new_turn_library(
//             string const & turn_type_tag,
//             strings const & turns,
//             map< string, StubFrags > & all_frags
//             )
// {
//  SegmentType lower_segment_type( alpha_type ), upper_segment_type( alpha_type );
//  if ( turn_type_tag == "alpha_beta" ) {
//   lower_segment_type = alpha_type;
//   upper_segment_type = beta_type;
//  } else if ( turn_type_tag == "beta_alpha" ) {
//   lower_segment_type = beta_type;
//   upper_segment_type = alpha_type;
//  } else {
//   utility_exit_with_message("unrecognized turn_type_tag: "+turn_type_tag);
//  }

//  Size const lower_buffer( get_segment_window( lower_segment_type ) ),
//   upper_buffer( get_segment_window( upper_segment_type ) );
//  static bool init( false );

//  static map< string, vector1< TurnCoordsLine > > all_turn_coords;

//  if ( !init ) { //read info on bb coords
//   init = true;
//   string const pdb_coords_file( option[ my_options::new_turns_file ] );
//   ifstream data( pdb_coords_file.c_str() );
//   string line;
//   while ( getline( data, line ) ) {
//    strings const l( split_to_vector1( line ) );
//    if ( l[1] == "turn_coords:" && l[2] == turn_type_tag ) {
//     TurnCoordsLine tcline;
//     string const turn( l[4] ); // might be a "-"
//     Size const turnlen( turn == "-" ? 0 : turn.size() );
//     runtime_assert( int(turnlen) == int_of( l[3] ) );
//     if ( all_turn_coords.find( turn ) == all_turn_coords.end() ) all_turn_coords[ turn ];
//     //Vectors coords;
//     Size const coordslen( int_of( l[8] ) );
//     runtime_assert( coordslen == turnlen + lower_buffer + upper_buffer );
//     for ( Size i=1; i<= coordslen; ++i ) {
//      tcline.coords.push_back( Vector( float_of( l[ 3*i+4 ] ), float_of( l[ 3*i+5 ] ), float_of( l[3*i+6 ] ) ) );
//     }
//     Size const torsions_start( 3*coordslen+10 );
//     runtime_assert( l[ torsions_start-1 ] == "bb_torsions:" );
//     for ( Size i=0; i< turnlen; ++i ) { //read the bb torsions
//      tcline.torsions.push_back( make_vector1( Real( float_of( l[ torsions_start+3*i   ] ) ),
//                           Real( float_of( l[ torsions_start+3*i+1 ] ) ),
//                           Real( float_of( l[ torsions_start+3*i+2 ] ) ) ) );
//     }
//     tcline.tag = filebase( l.back() ) + "_" + l[ l.size()-2 ] + "_" + l[ l.size()-1 ] +"_"+turn;
//     tcline.sequence = l[5]; // might be a "-"
//     runtime_assert( tcline.sequence.size() == turn.size() );
//     all_turn_coords.find( turn )->second.push_back( tcline );
//     //TR.Trace << "read turn: " << turn << ' ' << all_turn_coords.find( turn )->second.size() << endl;
//    }

//   }
//   data.close();
//  } // initialize the coordinates /////////////////////////////////////////////////////////////////////


//  all_frags.clear();
//  for ( strings::const_iterator turn= turns.begin(); turn != turns.end(); ++turn ) {
//   if ( all_turn_coords.find( *turn ) == all_turn_coords.end() ) {
//    utility_exit_with_message("no PDB turn frags for turn "+(*turn));
//   }
//   all_frags[ *turn ];
//   Size const turnlen( *turn == "-" ? 0 : turn->size() );
//   StubFrags & frags( all_frags.find( *turn )->second );
//   //vector1< std::pair< string, Vectors > > const & pdb_turn_coords( all_turn_coords.find( *turn )->second );
//   vector1< TurnCoordsLine > const & pdb_turn_coords( all_turn_coords.find( *turn )->second );
//   Reals all_strains;
//   for ( Size ii=1; ii<= pdb_turn_coords.size(); ++ii ) {
//    TurnCoordsLine const & tcl( pdb_turn_coords[ii] );
//    Vectors const & turn_coords( tcl.coords );
//    runtime_assert( turn_coords.size() == turnlen + lower_buffer + upper_buffer );
//    Size const stub1_segment_begin( 1 ), stub2_segment_begin( lower_buffer+turnlen+1 );
//    Stub stub1, stub2, tmpstub;
//    // stop_stub of segment before turn
//    get_segment_window_stubs( lower_segment_type, turn_coords, stub1_segment_begin, tmpstub, stub1 );
//    // start_stub of segment after turn
//    get_segment_window_stubs( upper_segment_type, turn_coords, stub2_segment_begin, stub2, tmpstub );
//    StubFrag f;
//    f.rt = RT( stub1, stub2 );
//    for ( Size j=1; j<= turnlen; ++j ) f.coords.push_back( stub1.global2local( turn_coords[lower_buffer+j] ));
//    // Real avg_rmsd, med_rmsd;
//    // compute_turn_strain( *turn, turn_buffer+1, turn_coords, avg_rmsd, med_rmsd );
//    // f.raw_strain = avg_rmsd;
//    // all_strains.push_back( f.raw_strain );
//    f.torsions = tcl.torsions;
//    f.sequence = tcl.sequence;
//    f.id = tcl.tag;
//    f.raw_strain = f.norm_strain = 0.0;
//    frags.push_back( f );

//    { //get some geometry params
//     Real const axis_angle( degrees( acos( stub1.M.col_x().dot( stub2.M.col_x() ) ) ) );
//     // dihedral
//     Vector const p1( stub1.v - stub1.M.col_x() ), p2 ( stub1.v ), p3( stub2.v ), p4( stub2.v + stub2.M.col_x() );
//     Real const axis_twist( numeric::dihedral_degrees( p1,p2,p3,p4) );
//     TR.Trace << "turn_geom " << turn_type_tag << ' ' <<
//      *turn << " axis_enddist: " << F(9,3,stub1.v.distance(stub2.v) ) <<
//      " axis_angle: " << F(9,3,axis_angle) << " axis_twist: " << F(9,3,axis_twist) << endl;
//    }
//    //if ( frags.size() >= max_frags ) break;
//   } // loop over pdb turns of this type

//   Size const nfrags( frags.size() );

//   // for ( Size i=1; i<= nfrags; ++i ) {
//   //  StubFrag & f( frags[i] );
//   //  Size nbetter(0);
//   //  Real const epsilon( 1e-3 );
//   //  for ( Reals::const_iterator s= all_strains.begin(); s!= all_strains.end(); ++s ) {
//   //   if ( *s + epsilon < f.raw_strain ) ++nbetter;
//   //  }

//   //  f.norm_strain = Real( nbetter );
//   //  // nbetter will range from 0 to nfrags-1
//   //  if ( nfrags > 1 ) f.norm_strain /= ( nfrags-1 ); // now range in [0,1]
//   //  //TR.Trace << "norm_strain_turn: " << *turn << F(9,3,f.raw_strain) << F(9,3,f.norm_strain) << endl;
//   // }


//   TR.Trace << "build_turn_library:: Read " << frags.size() << " turns of type " << *turn << endl;
//  }
// }




void
compute_bb_strain_rmsd(
	Vectors const & ca_coords,
	Size const nrepeat,
	Size const repeatlen,
	Size const helix1_len,
	Size const helix2_len,
	string const & turn1,
	string const & turn2,
	Real & helix_avg_rmsd,
	Real & turn_avg_rmsd
)
{
	/// not very robust right now
	/// will fail if helices are too short or the turns are the wrong kind...
	/// does not confirm that turns start and stop with non-alpha
	///

	Size const turn_buffer( 7 );
	//Size const helix_window( 7 );

	helix_avg_rmsd = turn_avg_rmsd = 0.0;

	Sizes
		helixbegins( make_vector1( Size(1), helix1_len + turn1.size() + 1 ) ),
		helixends  ( make_vector1( helix1_len  , repeatlen-turn2.size() ) ),
		turnbegins ( make_vector1( helix1_len+1, repeatlen-turn2.size()+1 ) );
	//turnends   ( make_vector1( helix1_len+turn1.size(), repeatlen ) );

	runtime_assert( ca_coords.size() == repeatlen * nrepeat );
	runtime_assert( helix1_len >= turn_buffer && helix2_len >= turn_buffer );


	// HELICES
	for ( Size r=1; r<= 2; ++r ) {
		Real avg_avg_rmsd( 0 ), avg_med_rmsd( 0 );
		compute_helix_strain( helixbegins[r], helixends[r], ca_coords, avg_avg_rmsd, avg_med_rmsd );
		helix_avg_rmsd += avg_avg_rmsd;
	} // r=1,2
	helix_avg_rmsd /= 2;


	// TURNS
	for ( Size r=1; r<= 2; ++r ) {
		string const turn( r==1 ? turn1 : turn2 );
		Real avg_rmsd, med_rmsd;
		compute_turn_strain( turn, turnbegins[r], ca_coords, avg_rmsd, med_rmsd );
		turn_avg_rmsd += avg_rmsd;
	}
	turn_avg_rmsd /= 2;

}







///////////////////////////////////////////////////////////////////////////////
void
get_sasas_by_atomtype(
	Pose const & pose,
	bools const & subset,
	Real & polar_sasa,
	Real & nonpolar_sasa
)
{
	Reals rsd_polar_sasa, rsd_nonpolar_sasa;
	get_rsd_sasas_by_atomtype( pose, subset, polar_sasa, nonpolar_sasa, rsd_polar_sasa, rsd_nonpolar_sasa );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////



Size
repeatlen_from_symminfo( Pose const & pose )
{
	//runtime_assert( pose::symmetry::is_symmetric( pose ) );
	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );
	return repeatlen;
}
Size
nrepeat_from_symminfo( Pose const & pose )
{
	//runtime_assert( pose::symmetry::is_symmetric( pose ) );
	Size nrepeat, repeatlen, base_repeat;
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );
	return repeatlen;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void
// add_silly_vrt_residues( Pose & pose )
// {
//  Size nrepeat, repeatlen, base_repeat;
//  parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

// }




///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// skips non-protein
// void
// find_buried_unsatisfied_polars(
//                 vector1< id::AtomID > & buried_unsatisfied_donors,
//                 vector1< id::AtomID > & buried_unsatisfied_acceptors,
//                 Pose const & pose_in,
//                 string const tag = string("")
//                 )
// {

//  Pose pose( pose_in );

//  pose::symmetry::make_asymmetric_pose( pose );

//  Real const hbond_energy_threshold( -0.01 );
//  Real const burial_threshold( 0.01 );
//  Real const probe_radius( 1.0 ); // less than 1.4 --> fewer things "buried" --> more conservative


//  id::AtomID_Map< Real > atom_sasa, total_hbond_energy;
//  utility::vector1< Real > rsd_sasa;
//  bool const use_big_polar_H( true );
//  scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, use_big_polar_H );

//  // HACK!!!!!!!!!!!!!! TO TRIGGER 10A nbr graph update
//  scoring::hbonds::HBondSet hbond_set;
//  {
//   using namespace scoring;
//   ScoreFunction sf;
//    sf.set_weight( hbond_sc, 1.0 );
//   sf(pose);
//   pose.update_residue_neighbors();
//   scoring::hbonds::fill_hbond_set( pose, false, hbond_set );
//  }

//  core::pose::initialize_atomid_map( total_hbond_energy, pose, 0.0 );

//  for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
//   if ( true || hbond_set.allow_hbond(i) ) { // EVEN UNALLOWED HBONDS...
//    scoring::hbonds::HBond const & hb( hbond_set.hbond(i) );
//    id::AtomID const hatm( hb.don_hatm(), hb.don_res() );
//    id::AtomID const aatm( hb.acc_atm() , hb.acc_res() );
//    total_hbond_energy[ hatm ] += hb.energy(); // unweighted
//    total_hbond_energy[ aatm ] += hb.energy(); // unweighted
//   }
//  }

//  // now find unsat+buried
//  buried_unsatisfied_donors.clear();
//  buried_unsatisfied_acceptors.clear();

//  for ( Size i=1; i<= pose.total_residue(); ++i ) {
//   conformation::Residue const & rsd( pose.residue(i) );
//   if ( !rsd.is_protein() ) continue; /// SKIPPING NON-PROTEIN RIGHT NOW !!!!!!!!!!!!!!!!!!!!!!!!

//   // donors
//   for ( chemical::AtomIndices::const_iterator
//       hnum  = rsd.Hpos_polar().begin(),
//       hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
//    id::AtomID const hatm( *hnum, i );
//    if ( total_hbond_energy[ hatm ] >= hbond_energy_threshold &&
//       atom_sasa[ hatm ] <= burial_threshold ) {
//     TR.Trace << "UNSAT_DONOR " << I(4,i) << ' ' << rsd.name3() << ' ' << rsd.atom_name(*hnum) << ' ' <<
//      tag << std::endl;
//     buried_unsatisfied_donors.push_back( hatm );
//     //++buried_unsatisfied_donors;
//    }
//   }

//   // acceptors
//   for ( chemical::AtomIndices::const_iterator
//       anum  = rsd.accpt_pos().begin(),
//       anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
//    id::AtomID const aatm( *anum, i );
//    if ( total_hbond_energy[ aatm ] >= hbond_energy_threshold &&
//       atom_sasa[ aatm ] <= burial_threshold ) {
//     TR.Trace << "UNSAT_ACCEPTOR " << I(4,i) << ' ' << rsd.name3() << ' ' << rsd.atom_name(*anum) << ' ' <<
//      tag << std::endl;
//     buried_unsatisfied_acceptors.push_back( aatm );
//     //++buried_unsatisfied_acceptors;
//    }
//   }
//  } // i
// }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void
// count_bb_and_sc_atomids(
//             vector1< id::AtomID > const & atomids,
//             Pose const & pose,
//             Size & n_bb,
//             Size & n_sc
//             )
// {
//  n_bb = n_sc = 0;

//  for ( Size i=1; i<= atomids.size(); ++i ) {
//   if ( pose.residue( atomids[i].rsd() ).atom_is_backbone( atomids[i].atomno() ) ) ++n_bb;
//   else ++n_sc;
//  }
// }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// skips non-protein
// string
// get_buried_unsatisfied_string(
//                Pose const & pose
//                )
// {
//  bools const subset( pose.total_residue(), true );
//  return get_buried_unsatisfied_string( subset,pose);
// using id::AtomID;

// vector1< id::AtomID > buried_unsatisfied_acceptors, buried_unsatisfied_donors;

// find_buried_unsatisfied_polars( subset, pose, buried_unsatisfied_donors, buried_unsatisfied_acceptors );

// Size n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
//  n_buried_unsatisfied_acceptors_sc;

// count_bb_and_sc_atomids( buried_unsatisfied_donors, pose,
//              n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc );

// count_bb_and_sc_atomids( buried_unsatisfied_acceptors, pose,
//              n_buried_unsatisfied_acceptors_bb, n_buried_unsatisfied_acceptors_sc );

// ostringstream out;
// out << " n_buried_unsatisfied_donors: " << buried_unsatisfied_donors.size() <<
//  " n_buried_unsatisfied_donors_bb: " << n_buried_unsatisfied_donors_bb <<
//  " n_buried_unsatisfied_donors_sc: " << n_buried_unsatisfied_donors_sc <<
//  " buried_unsatisfied_donors: ";
// if ( buried_unsatisfied_donors.empty() ) out << '-';
// for ( Size i=1; i<= buried_unsatisfied_donors.size(); ++i ) {
//  if ( i>1 ) out << ',';
//  AtomID const & id( buried_unsatisfied_donors[i] );
//  Residue const & rsd( pose.residue(id.rsd() ) );
//  out << id.rsd() << "." << rsd.name1() << "." << stripped( rsd.atom_name(id.atomno()));
// }

// out << " n_buried_unsatisfied_acceptors: " << buried_unsatisfied_acceptors.size() <<
//  " n_buried_unsatisfied_acceptors_bb: " << n_buried_unsatisfied_acceptors_bb <<
//  " n_buried_unsatisfied_acceptors_sc: " << n_buried_unsatisfied_acceptors_sc <<
//  " buried_unsatisfied_acceptors: ";
// if ( buried_unsatisfied_acceptors.empty() ) out << '-';
// for ( Size i=1; i<= buried_unsatisfied_acceptors.size(); ++i ) {
//  if ( i>1 ) out << ',';
//  AtomID const & id( buried_unsatisfied_acceptors[i] );
//  Residue const & rsd( pose.residue(id.rsd() ) );
//  out << id.rsd() << "." << rsd.name1() << "." << stripped( rsd.atom_name(id.atomno()));
// }

// return out.str();

// }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
dna_test()
{
	// get familiar with using the symmetry code
	//
	// read in a dna molecule
	//
	// declare symmetry
	//
	// try minimization
	//
	//
	string const filename( start_file() );

	Pose pose;
	pose_from_pdb( pose, filename );

	set_base_partner( pose );

	Size const nb( pose.total_residue()/2 );

	/// try to figure out the z-axis (normal and location)
	///
	Vector center, n, t;
	Real theta;
	{
		Stub const stub1( scoring::dna::get_base_pair_stub_slow( pose.residue(1), pose.residue(2*nb) ) );
		Stub const stub2( scoring::dna::get_base_pair_stub_slow( pose.residue(2), pose.residue(2*nb-1) ) );

		// two matrices are related by a rotation about the axis
		// can find the normal using numeric code
		//
		Matrix const R( stub2.M * stub1.M.transposed() );
		n = numeric::rotation_axis( R, theta );

		Vector const v1( pose.residue(1).xyz("C1'" ) ), v2( pose.residue(2).xyz("C1'") );
		t = (v2-v1).dot(n) * n;
		Vector const v2p( v2 - t );


		cout << "theta: " << F(9,3,numeric::conversions::degrees(theta)) << " n: " <<
			F(9,3,n.x()) << F(9,3,n.y()) << F(9,3,n.z()) << endl;

		// confirm that these guys are now in the same plane
		runtime_assert( is_small( (v1-v2p).dot( n ) ) );


		// now try to figure out what the x,y coords of the axis are
		Real const y( v1.distance(v2p)/2 ), x( y / tan( theta/2 ) );
		Vector const midpoint( 0.5 * (v1+v2p) ), ihat( ( v2p - v1 ).normalized() ), khat( n ),
			jhat( khat.cross( ihat ) );
		center = midpoint + x * jhat;

		// check
		Real const dev( v2.distance( ( R*( v1-center) + t + center ) ) );
		cout << "dev: " << F(9,3,dev) << endl;
		runtime_assert( dev<1e-3 );

	}

	{ // add some virtual residues with the right geometry
		ResidueOP vrtrsd
			( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRT" ) ) );

		for ( Size i=1; i<= nb; ++i ) {
			// a virtual residue with origin along the axis
			// x pointing toward C1'
			// y perp to x and axis normal
			Vector const v( pose.residue(i).xyz("C1'") ), orig( v - (   ( v - center ) - (v-center).dot(n)*n ) );
			runtime_assert( is_small( orig.dot(n) - v.dot(n) ) );
			runtime_assert( is_small( (orig-center).dot(orig-v) ) );
			Vector const z(n), x( (v-orig).normalized() ), y( z.cross(x) );

			vrtrsd->set_xyz( "ORIG", orig );
			vrtrsd->set_xyz( "X", orig + x );
			vrtrsd->set_xyz( "Y", orig + y );

			pose.append_residue_by_jump( *vrtrsd, 1 );
		}
	}
	runtime_assert( pose.total_residue() == 2*nb + nb );


	{ // setup a simple fold tree
		FoldTree f( pose.total_residue() );
		for ( Size i=1; i<= nb; ++i ) {
			Size const ip( retrieve_base_partner_from_pose( pose )[i] );
			f.new_jump( i, ip, ip-1 );
		}
		for ( Size i=1; i<= nb; ++i ) {
			Size const ivrt( 2*nb+i );
			Size const cut( i==nb ? 2*nb : i );
			f.new_jump( i, ivrt, cut );
		}
		for ( Size i=2*nb+1; i<pose.total_residue(); ++i ) {
			f.new_jump( i, i+1, i );
		}

		f.reorder( 2*nb + 1 ); // root at first of the virtual residues
		pose.fold_tree(f);
	}

	Pose asympose( pose ); // save asymmetric pose

	// these two guys together are the independent monomer
	Size const pos1( 5 ), pos1p( retrieve_base_partner_from_pose( pose )[ pos1 ] ), pos1v( 2*nb+pos1 );

	// try to setup some symmetry info
	//  Size const jump1a( pose.fold_tree().get_jump_that_builds_residue( pos1+1 ) ),

	//
	conformation::symmetry::SymmetryInfo symminfo;

	for ( Size i=1; i<= nb; ++i ) {
		if ( i == pos1 ) continue;
		Size const ip( retrieve_base_partner_from_pose( pose )[i] );

		symminfo.add_bb_clone( pos1, i );
		symminfo.add_chi_clone( pos1, i );

		symminfo.add_bb_clone( pos1p, ip );
		symminfo.add_chi_clone( pos1p, ip );

	}

	/// now the jumps

	/// this is the base "intra-monomer" jump
	Size const jump1( pose.fold_tree().get_jump_that_builds_residue( pos1p ) );
	/// this is the base "vrt-->monomer" jump
	Size const jump2( pose.fold_tree().get_jump_that_builds_residue( pos1 ) );
	/// this is the base "intra-vrt jump"
	Size const jump3( pose.fold_tree().get_jump_that_builds_residue( pos1v+1 ) );

	runtime_assert( pose.fold_tree().upstream_jump_residue(jump1) == pos1 );
	runtime_assert( pose.fold_tree().upstream_jump_residue(jump2) == pos1v );
	runtime_assert( pose.fold_tree().upstream_jump_residue(jump3) == pos1v );

	for ( Size i=1; i<= nb; ++i ) {
		Size const ip( retrieve_base_partner_from_pose( pose )[i] );
		Size const iv( 2*nb+i );
		if ( i != pos1 ) {
			symminfo.add_jump_clone( jump1, pose.fold_tree().get_jump_that_builds_residue( ip  ), 0.0 );
			symminfo.add_jump_clone( jump2, pose.fold_tree().get_jump_that_builds_residue( i   ), 0.0 );
			if ( iv+1<= pose.total_residue() ) {
				symminfo.add_jump_clone( jump3, pose.fold_tree().get_jump_that_builds_residue( iv+1 ), 0.0 );
			}
		}
	}

	/// now need to set the symdofs
	/// jumps that build the monomers have all dofs flexible except z-trans and z-rot
	/// jumps that build vrt rsds have only z-trans and z-rot
	///
	using core::conformation::symmetry::SymDof;

	map< Size, SymDof > symdofs;
	if ( true ) {
		SymDof symdof2, symdof3;
		symdof2.read( "x y angle_x angle_y" );
		symdof3.read( "z angle_z" );
		symdofs[ jump2 ] = symdof2;
		symdofs[ jump3 ] = symdof3;
	}

	symminfo.set_dofs( symdofs );


	symminfo.num_virtuals( nb );
	symminfo.set_use_symmetry( true );
	//symminfo.contiguous_monomers( false );
	symminfo.set_flat_score_multiply( pose.total_residue(), 1 ); // was nb

	for ( Size i=1; i<pos1; ++i ) {
		Size const ip( retrieve_base_partner_from_pose( pose )[i] );
		symminfo.set_score_multiply( i, 0.0 );
		symminfo.set_score_multiply( ip, 0.0 );
	}

	symminfo.update_score_multiply_factor();


	//cout << symminfo << endl;

	pose::symmetry::make_symmetric_pose( pose, symminfo );


	// try scoring with a symmetric scoring function
	scoring::symmetry::SymmetricScoreFunctionOP scorefxn( new scoring::symmetry::SymmetricScoreFunction() );
	scoring::methods::EnergyMethodOptions options;
	options.exclude_DNA_DNA( false );
	scorefxn->set_energy_method_options( options );

	scorefxn->set_weight(  fa_atr, 0.8 );
	scorefxn->set_weight(  fa_rep, 0.44 );
	scorefxn->set_weight(  fa_sol, 0.65 );
	//scorefxn->set_weight(  hbond_sc, 1.1 );
	//scorefxn->set_weight( lk_ball, 0.325 );
	//scorefxn->set_weight( lk_ball_iso, -0.325 );

	Real score( (*scorefxn)( pose ) );

	cout << "final_scores: " << F(9,3,score) << " premin " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	// now minimize
	MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
	movemap->set_bb(false);
	movemap->set_chi(false);
	movemap->set_jump(true);
	movemap->set_jump(jump1,false); // get rid of the intra-monomer jump for the moment


	if ( false ) { // asymmetric minimization

		Pose & pose( asympose );

		scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction() );
		scoring::methods::EnergyMethodOptions options;
		options.exclude_DNA_DNA( false );
		scorefxn->set_energy_method_options( options );

		scorefxn->set_weight(  fa_atr, 0.8 );
		scorefxn->set_weight(  fa_rep, 0.44 );
		scorefxn->set_weight(  fa_sol, 0.65 );
		//scorefxn->set_weight(  hbond_sc, 1.1 );
		//scorefxn->set_weight( lk_ball, 0.325 );
		//scorefxn->set_weight( lk_ball_iso, -0.325 );

		Real score( (*scorefxn)( pose ) );

		cout << "final_scores: " << F(9,3,score) << " premin " <<
			pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

		movemap->set_jump( false );
		movemap->set_jump( jump2, true );
		movemap->set_jump( jump3, true );

		protocols::minimization_packing::MinMoverOP min_mover
			( new protocols::minimization_packing::MinMover(movemap, scorefxn, "dfpmin", 0.00001, true,
			true, true ));

		min_mover->apply( pose );

		pose.dump_pdb("aftermin_asym.pdb");

		score = ( (*scorefxn)( pose ) );

		cout << "final_scores: " << F(9,3,score) << " postmin_asym " <<
			pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

		exit(0);

	}


	protocols::minimization_packing::symmetry::SymMinMoverOP min_mover
		( new protocols::minimization_packing::symmetry::SymMinMover(movemap, scorefxn, "dfpmin_armijo_nonmonotone", 0.00001, true,
		true, true ));

	min_mover->apply( pose );

	pose.dump_pdb("aftermin.pdb");

	score = ( (*scorefxn)( pose ) );

	cout << "final_scores: " << F(9,3,score) << " postmin " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	exit(0);

	// great! now let's try packing...
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	task->or_include_current( true );
	task->restrict_to_repacking();

	Size const nloop( 25 );
	protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( scorefxn, task, nloop );

	packmover.apply( pose );

	score = ( (*scorefxn)( pose ) );

	cout << "final_scores: " << F(9,3,score) << " postpack " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	pose.dump_pdb("afterpack.pdb");

	protocols::minimization_packing::symmetry::SymRotamerTrialsMover rtmover( scorefxn, *task );
	rtmover.apply( pose );
	score = ( (*scorefxn)( pose ) );

	cout << "final_scores: " << F(9,3,score) << " postrt " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	pose.dump_pdb("afterrt.pdb");

}

void
tal_test()
{
	Size const anchorpos( 13 ); // might want to make this configurable

	//
	string const filename( start_file() );

	Pose pose;
	pose_from_pdb( pose, filename );

	set_base_partner( pose );

	runtime_assert( pose.conformation().num_chains() == 3 );

	Size const nrepeat( pose.conformation().chain_end(2) - pose.conformation().chain_begin(2) + 1 );
	Size const repeatlen( 34 );
	Size base_repeat(0);
	if ( nrepeat == 7 ) base_repeat = 3;
	else if ( nrepeat == 8 ) base_repeat = 3;
	else utility_exit_with_message("bad nrepeat "+string_of(nrepeat));

	runtime_assert( pose.conformation().chain_end(1) == nrepeat*repeatlen );
	runtime_assert( nrepeat == ( pose.conformation().chain_end(3) - pose.conformation().chain_begin(3) + 1 ) );

	/// try to figure out the z-axis (normal and location)
	///
	Vector center, n, t;
	Real theta;
	{
		Size const dpos1( pose.conformation().chain_begin(2)+1 ), dpos2( dpos1+1 ),
			dpos1p( retrieve_base_partner_from_pose( pose )[dpos1 ] ),
			dpos2p( retrieve_base_partner_from_pose( pose )[dpos2 ] );
		Stub const stub1( scoring::dna::get_base_pair_stub_slow( pose.residue(dpos1), pose.residue(dpos1p) ) );
		Stub const stub2( scoring::dna::get_base_pair_stub_slow( pose.residue(dpos2), pose.residue(dpos2p) ) );

		// two matrices are related by a rotation about the axis
		// can find the normal using numeric code
		//
		Matrix const R( stub2.M * stub1.M.transposed() );
		n = numeric::rotation_axis( R, theta );

		Vector const v1( pose.residue(dpos1).xyz("C1'" ) ), v2( pose.residue(dpos2).xyz("C1'") );
		t = (v2-v1).dot(n) * n;
		Vector const v2p( v2 - t );


		cout << "theta: " << F(9,3,numeric::conversions::degrees(theta)) << " rise: " << F(9,3,t.dot(n)) <<
			" n: " << F(9,3,n.x()) << F(9,3,n.y()) << F(9,3,n.z()) << endl;

		// confirm that these guys are now in the same plane
		runtime_assert( is_small( (v1-v2p).dot( n ) ) );


		// now try to figure out what the x,y coords of the axis are
		Real const y( v1.distance(v2p)/2 ), x( y / tan( theta/2 ) );
		Vector const midpoint( 0.5 * (v1+v2p) ), ihat( ( v2p - v1 ).normalized() ), khat( n ),
			jhat( khat.cross( ihat ) );
		center = midpoint + x * jhat;

		// check
		Real const dev( v2.distance( ( R*( v1-center) + t + center ) ) );
		cout << "dev: " << F(9,3,dev) << endl;
		runtime_assert( dev<1e-3 );

		exit(0);

	}

	{ // add some virtual residues with the right geometry
		ResidueOP vrtrsd
			( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRT" ) ) );

		for ( Size i=0; i< nrepeat; ++i ) {
			// a virtual residue with origin along the axis
			// x pointing toward C1'
			// y perp to x and axis normal
			Size const dpos( pose.conformation().chain_begin(2)+i );
			Vector const v( pose.residue(dpos).xyz("C1'") ), orig( v - (   ( v - center ) - (v-center).dot(n)*n ) );
			runtime_assert( is_small( orig.dot(n) - v.dot(n) ) );
			runtime_assert( is_small( (orig-center).dot(orig-v) ) );
			Vector const z(n), x( (v-orig).normalized() ), y( z.cross(x) );

			vrtrsd->set_xyz( "ORIG", orig );
			vrtrsd->set_xyz( "X", orig + x );
			vrtrsd->set_xyz( "Y", orig + y );

			pose.append_residue_by_jump( *vrtrsd, 1 );
		}
	}
	Size const nres( pose.total_residue() ), nres_protein( repeatlen*nrepeat ), nres_dna( 2*nrepeat ),
		nres_nonvirtual( nres_protein + nres_dna );
	runtime_assert( nres == repeatlen*nrepeat + 2*nrepeat + nrepeat );


	{ // setup a simple fold tree
		FoldTree f( pose.total_residue() );
		for ( Size i=0; i< nrepeat; ++i ) {
			Size const dpos( pose.conformation().chain_begin(2)+i ), dposp( retrieve_base_partner_from_pose( pose )[dpos ] ),
				ppos( i*repeatlen + anchorpos ), pcut( (i+1)*repeatlen ), vrtpos( nres_nonvirtual+i+1 );
			f.new_jump( dpos, dposp, dposp-1 ); // base-pair jump, cut before dposp
			f.new_jump( ppos, dpos, pcut ); // dna-->prot jump, cut after prot repeat
			if ( i==nrepeat-1 ) f.new_jump( dpos, vrtpos, nres_nonvirtual ); // last dna-vrt jump: cut before 1st virtual
			else f.new_jump( dpos, vrtpos, dpos ); // dna-vrt jump: cut after dpos
		}
		for ( Size i=nres_nonvirtual+1; i<nres; ++i ) {
			f.new_jump( i, i+1, i ); // jumps between the virtuals
		}

		f.reorder( nres_nonvirtual + 1 ); // root at first of the virtual residues
		pose.fold_tree(f);
	}

	Pose asympose( pose ); // save asymmetric pose


	//
	conformation::symmetry::SymmetryInfo symminfo;

	Size const base_protoffset( (base_repeat-1)*repeatlen ), base_dpos( nres_protein + base_repeat ),
		base_dposp( retrieve_base_partner_from_pose( pose )[ base_dpos ] ), base_protanchor( base_protoffset+anchorpos ),
		base_vpos( nres_nonvirtual+base_repeat );

	cout << "base_info " << base_repeat << ' ' << base_dpos << ' ' << base_dposp << ' ' << base_protoffset << ' ' <<
		base_vpos << ' ' << base_protanchor << endl;

	/// add all the bb/chi clone information
	for ( Size i=1; i<= nrepeat; ++i ) {
		if ( i == base_repeat ) continue;

		Size const dpos( nres_protein+i ), dposp( retrieve_base_partner_from_pose( pose )[dpos ] ),
			protoffset( (i-1)*repeatlen );

		/// dna clones
		symminfo.add_bb_clone ( base_dpos, dpos );
		symminfo.add_chi_clone( base_dpos, dpos );

		symminfo.add_bb_clone ( base_dposp, dposp );
		symminfo.add_chi_clone( base_dposp, dposp );

		/// protein clones
		for ( Size j=1; j<= repeatlen; ++j ) {
			symminfo.add_bb_clone ( base_protoffset+j, protoffset+j );
			symminfo.add_chi_clone( base_protoffset+j, protoffset+j );
		}
	}

	/// now the jumps

	/// this is the base "intra-dna" jump
	Size const jump1( pose.fold_tree().get_jump_that_builds_residue( base_dposp ) );
	/// this is the base "prot-dna" jump
	Size const jump2( pose.fold_tree().get_jump_that_builds_residue( base_protanchor ) );
	/// this is the base "vrt-->monomer" jump
	Size const jump3( pose.fold_tree().get_jump_that_builds_residue( base_dpos ) );
	/// this is the base "intra-vrt jump"
	Size const jump4( pose.fold_tree().get_jump_that_builds_residue( base_vpos+1 ) );

	runtime_assert( pose.fold_tree().upstream_jump_residue(jump1) == base_dpos );
	runtime_assert( pose.fold_tree().upstream_jump_residue(jump2) == base_dpos );
	runtime_assert( pose.fold_tree().upstream_jump_residue(jump3) == base_vpos );
	runtime_assert( pose.fold_tree().upstream_jump_residue(jump4) == base_vpos );

	for ( Size i=1; i<= nrepeat; ++i ) {
		if ( i == base_repeat ) continue;

		Size const dpos( nres_protein+i ), dposp( retrieve_base_partner_from_pose( pose )[dpos ] ),
			protoffset( (i-1)*repeatlen ), protanchor( protoffset + anchorpos ), vpos( nres_nonvirtual+i );

		symminfo.add_jump_clone( jump1, pose.fold_tree().get_jump_that_builds_residue( dposp ), 0.0 );
		symminfo.add_jump_clone( jump2, pose.fold_tree().get_jump_that_builds_residue( protanchor ), 0.0 );
		symminfo.add_jump_clone( jump3, pose.fold_tree().get_jump_that_builds_residue( dpos ), 0.0 );
		if ( vpos+1 <= pose.total_residue() ) {
			symminfo.add_jump_clone( jump4, pose.fold_tree().get_jump_that_builds_residue( vpos+1 ), 0.0 );
		}
	}

	/// now need to set the symdofs
	/// jumps that build the monomers have all dofs flexible except z-trans and z-rot
	/// jumps that build vrt rsds have only z-trans and z-rot
	///
	using core::conformation::symmetry::SymDof;

	map< Size, SymDof > symdofs;
	if ( false ) { /// EXCLUDE VRT JUMP MOVEMENT
		SymDof symdof3, symdof4;
		symdof3.read( "x y angle_x angle_y" );
		symdof4.read( "z angle_z" );
		symdofs[ jump3 ] = symdof3;
		symdofs[ jump4 ] = symdof4;
	}

	symminfo.set_dofs( symdofs );


	symminfo.num_virtuals( nrepeat );
	symminfo.set_use_symmetry( true );
	//symminfo.contiguous_monomers( false );
	symminfo.set_flat_score_multiply( pose.total_residue(), 1 ); // was nb

	/// repeats prior to base_repeat have score_multiply ==> 0
	for ( Size i=1; i<base_repeat; ++i ) {
		Size const dpos( nres_protein+i ), dposp( retrieve_base_partner_from_pose( pose )[dpos] ),
			protoffset( (i-1)*repeatlen );
		symminfo.set_score_multiply( dpos, 0 );
		symminfo.set_score_multiply( dposp, 0 );
		for ( Size j=1; j<= repeatlen; ++j ) symminfo.set_score_multiply( protoffset+j, 0 );
	}

	symminfo.update_score_multiply_factor();


	//cout << symminfo << endl;

	pose::symmetry::make_symmetric_pose( pose, symminfo );



	// try scoring with a symmetric scoring function
	scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	//create_score_function( scoring::STANDARD_WTS, make_vector1( SCORE12_PATCH ) ) );


	//  scoring::symmetry::SymmetricScoreFunctionOP scorefxn
	//   ( dynamic_cast< SymmetricScoreFunctionOP >(

	//   ( new scoring::symmetry::SymmetricScoreFunction() );
	scoring::methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn->set_energy_method_options( options );

	//  scorefxn->set_weight(  fa_atr, 0.8 );
	//  scorefxn->set_weight(  fa_rep, 0.44 );
	//  scorefxn->set_weight(  fa_sol, 0.65 );
	//scorefxn->set_weight(  hbond_sc, 1.1 );
	//scorefxn->set_weight( lk_ball, 0.325 );
	//scorefxn->set_weight( lk_ball_iso, -0.325 );

	scorefxn->show( cout );

	Real score( (*scorefxn)( pose ) );

	cout << "final_scores: " << F(9,3,score) << " start " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;



	/// try "symmetrizing" the dofs

	pose.dump_pdb("start.pdb");

	pose.set_jump( jump4, pose.jump( jump4 ) ); // vrt-vrt
	pose.set_jump( jump3, pose.jump( jump3 ) ); // vrt-dna
	pose.set_jump( jump2, pose.jump( jump2 ) ); // prot-d
	pose.set_jump( jump1, pose.jump( jump1 ) ); // bp

	pose.dump_pdb("aftersym.pdb");

	score = (*scorefxn)( pose );

	cout << "final_scores: " << F(9,3,score) << " postsym " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;


	{ /// try fast relax
		protocols::relax::FastRelax fastrelax( scorefxn, 0 );

		MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
		movemap->set_bb(true);
		movemap->set_chi(true);
		movemap->set_jump(true);
		movemap->set_jump(jump1,false); // turn off base-pair optimization

		fastrelax.set_movemap( movemap );
		fastrelax.apply( pose );

		score = (*scorefxn)( pose );
	}


	// now minimize
	MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
	movemap->set_bb(false);
	movemap->set_chi(false);
	movemap->set_jump(true);
	movemap->set_jump(jump1,false); // turn off base-pair optimization

	//movemap->set_jump(jump1,false); // get rid of the intra-monomer jump for the moment




	protocols::minimization_packing::symmetry::SymMinMoverOP min_mover
		( new protocols::minimization_packing::symmetry::SymMinMover(movemap, scorefxn, "dfpmin_armijo_nonmonotone", 0.00001, true,
		true, true ));
	//false, false ) );

	min_mover->apply( pose );

	pose.dump_pdb("aftermin.pdb");

	score = ( (*scorefxn)( pose ) );

	cout << "final_scores: " << F(9,3,score) << " postmin " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	//exit(0);

	// great! now let's try packing...
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	task->or_include_current( true );

	if ( false ) { // only repacking
		task->restrict_to_repacking();
	} else { // still allow design
		bools is_protein( pose.total_residue(), false );
		for ( Size i=1; i<= nres_protein; ++i ) is_protein[i] = true;
		task->restrict_to_residues( is_protein );
	}

	if ( false ) { // hacking for crazy rottrials bug
		runtime_assert( pose.residue(118).aa() == aa_lys );
		bools posn_118( pose.total_residue(), false ); posn_118[118] = true;
		task->restrict_to_residues( posn_118 );
	}




	Size const nloop( 25 );
	protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( scorefxn, task, nloop );

	packmover.apply( pose );

	score = ( (*scorefxn)( pose ) );

	cout << "final_scores: " << F(9,3,score) << " postpack " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	pose.dump_pdb("afterpack.pdb");

	protocols::minimization_packing::symmetry::SymRotamerTrialsMover rtmover( scorefxn, *task );
	rtmover.apply( pose );
	score = ( (*scorefxn)( pose ) );

	cout << "final_scores: " << F(9,3,score) << " postrt " <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	pose.dump_pdb("afterrt.pdb");

}


////////////////////////////////////////////////////////////////////////////////
void
setup_symmetric_pose(
	Size const nrepeat,
	Size const repeatlen,
	Size const base_repeat,
	Size const anchorpos,
	bool const flex_symdofs,
	Size & bp_jump,
	Size & dnaprot_jump,
	Size & vrtdna_jump,
	Size & vrtvrt_jump,
	Pose & pose,
	bool const reversed = false,
	Size const anchorpos_in = 13
)
{
	runtime_assert( num_chains( pose ) == 3 );
	runtime_assert( pose.conformation().chain_end(1) == nrepeat*repeatlen );
	runtime_assert( nrepeat == ( pose.conformation().chain_end(3) - pose.conformation().chain_begin(3) + 1 ) );
	runtime_assert( nrepeat == ( pose.conformation().chain_end(2) - pose.conformation().chain_begin(2) + 1 ) );

	/// try to figure out the z-axis (normal and location)
	///
	Vector center, n, t;
	Real theta;
	{
		Size const dpos1( pose.conformation().chain_begin(2)+1 ), dpos2( dpos1+1 ),
			dpos1p( retrieve_base_partner_from_pose( pose )[dpos1 ] ),
			dpos2p( retrieve_base_partner_from_pose( pose )[dpos2 ] );
		Stub const stub1( scoring::dna::get_base_pair_stub_slow( pose.residue(dpos1), pose.residue(dpos1p) ) );
		Stub const stub2( scoring::dna::get_base_pair_stub_slow( pose.residue(dpos2), pose.residue(dpos2p) ) );

		// two matrices are related by a rotation about the axis
		// can find the normal using numeric code
		//
		Matrix const R( stub2.M * stub1.M.transposed() );
		n = numeric::rotation_axis( R, theta );

		Vector const v1( pose.residue(dpos1).xyz("C1'" ) ), v2( pose.residue(dpos2).xyz("C1'") );
		t = (v2-v1).dot(n) * n;
		Vector const v2p( v2 - t );


		//    cout << "theta: " << F(9,3,numeric::conversions::degrees(theta)) << " n: " <<
		//     F(9,3,n.x()) << F(9,3,n.y()) << F(9,3,n.z()) << endl;

		// confirm that these guys are now in the same plane
		runtime_assert( is_small( (v1-v2p).dot( n ) ) );


		// now try to figure out what the x,y coords of the axis are
		Real const y( v1.distance(v2p)/2 ), x( y / tan( theta/2 ) );
		Vector const midpoint( 0.5 * (v1+v2p) ), ihat( ( v2p - v1 ).normalized() ), khat( n ),
			jhat( khat.cross( ihat ) );
		center = midpoint + x * jhat;

		// check
		Real const dev( v2.distance( ( R*( v1-center) + t + center ) ) );
		//    cout << "dev: " << F(9,3,dev) << endl;
		runtime_assert( dev<1e-3 );

	}

	{ // add some virtual residues with the right geometry
		ResidueOP vrtrsd
			( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRT" ) ) );

		for ( Size i=0; i< nrepeat; ++i ) {
			// a virtual residue with origin along the axis
			// x pointing toward C1'
			// y perp to x and axis normal
			Size const dpos( pose.conformation().chain_begin(2)+i );
			Vector const v( pose.residue(dpos).xyz("C1'") ), orig( v - (   ( v - center ) - (v-center).dot(n)*n ) );
			runtime_assert( is_small( orig.dot(n) - v.dot(n) ) );
			runtime_assert( is_small( (orig-center).dot(orig-v) ) );
			Vector const z(n), x( (v-orig).normalized() ), y( z.cross(x) );

			vrtrsd->set_xyz( "ORIG", orig );
			vrtrsd->set_xyz( "X", orig + x );
			vrtrsd->set_xyz( "Y", orig + y );

			pose.append_residue_by_jump( *vrtrsd, 1 );
		}
	}
	Size const nres( pose.total_residue() ), nres_protein( repeatlen*nrepeat ), nres_dna( 2*nrepeat ),
		nres_nonvirtual( nres_protein + nres_dna );
	runtime_assert( nres == repeatlen*nrepeat + 2*nrepeat + nrepeat );

	set_base_partner( pose );


	{ // setup a simple fold tree
		FoldTree f( pose.total_residue() );
		for ( Size i=0; i< nrepeat; ++i ) {
			Size const protrepeat( reversed ? nrepeat - i : i+1 );
			Size const dpos( pose.conformation().chain_begin(2)+i ), dposp( retrieve_base_partner_from_pose( pose )[dpos ] ),
				ppos( (protrepeat-1)*repeatlen + anchorpos ), pcut( protrepeat*repeatlen ), vrtpos( nres_nonvirtual+i+1 );
			f.new_jump( dpos, dposp, dposp-1 ); // base-pair jump, cut before dposp
			f.new_jump( ppos, dpos, pcut ); // dna-->prot jump, cut after prot repeat
			if ( i==nrepeat-1 ) f.new_jump( dpos, vrtpos, nres_nonvirtual ); // last dna-vrt jump: cut before 1st virtual
			else f.new_jump( dpos, vrtpos, dpos ); // dna-vrt jump: cut after dpos
		}
		for ( Size i=nres_nonvirtual+1; i<nres; ++i ) {
			f.new_jump( i, i+1, i ); // jumps between the virtuals
		}

		f.reorder( nres_nonvirtual + 1 ); // root at first of the virtual residues

		/// guarantee that dna-->protein jumps will remain within residue
		for ( Size i=0; i<nrepeat; ++i ) {
			Size const protpos( i*repeatlen + anchorpos );
			Size const jumpno( f.get_jump_that_builds_residue( protpos ) );
			f.set_jump_atoms( jumpno, "C1'", "N", false ); // keep in residue
		}


		pose.fold_tree(f);
	}

	// make some cutpoint variants
	for ( Size i=1; i<nrepeat; ++i ) {
		Size const cut( i*repeatlen );
		pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cut );
		pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cut+1 );
	}



	//
	conformation::symmetry::SymmetryInfo symminfo;

	Size const base_protoffset( reversed ? (nrepeat-base_repeat)*repeatlen : (base_repeat-1)*repeatlen ),
		base_dpos( nres_protein + base_repeat ), base_dposp( retrieve_base_partner_from_pose( pose )[ base_dpos ] ),
		base_protanchor( base_protoffset+anchorpos ),
		base_vpos( nres_nonvirtual+base_repeat );

	//   cout << "base_info " << base_repeat << ' ' << base_dpos << ' ' << base_dposp << ' ' << base_protoffset << ' ' <<
	//    base_vpos << ' ' << base_protanchor << endl;

	/// add all the bb/chi clone information
	for ( Size i=1; i<= nrepeat; ++i ) {
		if ( i == base_repeat ) continue;

		Size const dpos( nres_protein+i ), dposp( retrieve_base_partner_from_pose( pose )[dpos ] ),
			protoffset( reversed ? (nrepeat-i)*repeatlen : (i-1)*repeatlen );

		/// dna clones
		symminfo.add_bb_clone ( base_dpos, dpos );
		symminfo.add_chi_clone( base_dpos, dpos );

		symminfo.add_bb_clone ( base_dposp, dposp );
		symminfo.add_chi_clone( base_dposp, dposp );

		/// protein clones
		for ( Size j=1; j<= repeatlen; ++j ) {
			symminfo.add_bb_clone ( base_protoffset+j, protoffset+j );
			symminfo.add_chi_clone( base_protoffset+j, protoffset+j );
		}
	}

	/// now the jumps

	/// this is the base "intra-dna" jump
	Size const jump1( pose.fold_tree().get_jump_that_builds_residue( base_dposp ) );
	/// this is the base "prot-dna" jump
	Size const jump2( pose.fold_tree().get_jump_that_builds_residue( base_protanchor ) );
	/// this is the base "vrt-->monomer" jump
	Size const jump3( pose.fold_tree().get_jump_that_builds_residue( base_dpos ) );
	/// this is the base "intra-vrt jump"
	Size const jump4( pose.fold_tree().get_jump_that_builds_residue( base_vpos+1 ) );

	runtime_assert( pose.fold_tree().upstream_jump_residue(jump1) == base_dpos );
	runtime_assert( pose.fold_tree().upstream_jump_residue(jump2) == base_dpos );
	runtime_assert( pose.fold_tree().upstream_jump_residue(jump3) == base_vpos );
	runtime_assert( pose.fold_tree().upstream_jump_residue(jump4) == base_vpos );


	for ( Size i=1; i<= nrepeat; ++i ) {
		if ( i == base_repeat ) continue;

		Size const dpos( nres_protein+i ), dposp( retrieve_base_partner_from_pose( pose )[dpos ] ),
			protoffset( reversed ? (nrepeat-i)*repeatlen : (i-1)*repeatlen ),
			protanchor( protoffset + anchorpos ), vpos( nres_nonvirtual+i );

		symminfo.add_jump_clone( jump1, pose.fold_tree().get_jump_that_builds_residue( dposp ), 0.0 );
		symminfo.add_jump_clone( jump2, pose.fold_tree().get_jump_that_builds_residue( protanchor ), 0.0 );
		symminfo.add_jump_clone( jump3, pose.fold_tree().get_jump_that_builds_residue( dpos ), 0.0 );
		if ( vpos+1 <= pose.total_residue() ) {
			symminfo.add_jump_clone( jump4, pose.fold_tree().get_jump_that_builds_residue( vpos+1 ), 0.0 );
		}
	}

	/// now need to set the symdofs
	/// jumps that build the monomers have all dofs flexible except z-trans and z-rot
	/// jumps that build vrt rsds have only z-trans and z-rot
	///
	using core::conformation::symmetry::SymDof;

	map< Size, SymDof > symdofs;
	if ( flex_symdofs ) {
		SymDof symdof3, symdof4;
		symdof3.read( "x y angle_x angle_y" );
		symdof4.read( "z angle_z" );
		symdofs[ jump3 ] = symdof3;
		symdofs[ jump4 ] = symdof4;
	}

	symminfo.set_dofs( symdofs );

	symminfo.num_virtuals( nrepeat );
	symminfo.set_use_symmetry( true );
	symminfo.set_flat_score_multiply( pose.total_residue(), 1 );

	/// repeats prior to base_repeat have score_multiply ==> 0
	for ( Size i=1; i<base_repeat; ++i ) {
		Size const dpos( nres_protein+i ), dposp( retrieve_base_partner_from_pose( pose )[dpos] ),
			protoffset( reversed ? (nrepeat-i)*repeatlen : (i-1)*repeatlen );
		symminfo.set_score_multiply( dpos, 0 );
		symminfo.set_score_multiply( dposp, 0 );
		for ( Size j=1; j<= repeatlen; ++j ) symminfo.set_score_multiply( protoffset+j, 0 );
	}

	symminfo.update_score_multiply_factor();

	pose::symmetry::make_symmetric_pose( pose, symminfo );

	bp_jump = jump1;
	dnaprot_jump = jump2;
	vrtdna_jump = jump3;
	vrtvrt_jump = jump4;

	/// now symmetrize the dofs
	//pose.dump_pdb("before_symjumps.pdb");

	if ( anchorpos != anchorpos_in ) {
		using namespace id;
		// tricky: get the desired transform
		Size const srcprotanchor( (base_repeat-1)*repeatlen + anchorpos_in );
		Residue const & srcprotrsd( pose.residue( srcprotanchor ) ), &dnarsd( pose.residue( base_dpos ) ),
			destprotrsd( pose.residue( base_protanchor ) );
		StubID const
			upstub( AtomID( dnarsd.chi_atoms(1)[4], base_dpos ),
			AtomID( dnarsd.chi_atoms(1)[3], base_dpos ),
			AtomID( dnarsd.chi_atoms(1)[2], base_dpos ) ),
			//    upstub( AtomID( dnarsd.atom_index( "P"   ), base_dpos ),
			//        AtomID( dnarsd.atom_index( "O5'" ), base_dpos ),
			//        AtomID( dnarsd.atom_index( "C5'" ), base_dpos ) ),
			src_downstub(  AtomID(  srcprotrsd.atom_index( "N"  ), srcprotanchor ),
			AtomID(  srcprotrsd.atom_index( "CA" ), srcprotanchor ),
			AtomID(  srcprotrsd.atom_index( "C"  ), srcprotanchor ) ),
			dest_downstub( AtomID( destprotrsd.atom_index( "N"  ), base_protanchor ),
			AtomID( destprotrsd.atom_index( "CA" ), base_protanchor ),
			AtomID( destprotrsd.atom_index( "C"  ), base_protanchor ) );
		kinematics::RT const rt( pose.conformation().get_stub_transform( upstub, src_downstub ) );

		//pose.dump_pdb("before.pdb");

		// this routine is not symmetrized:
		pose.conformation().set_stub_transform( upstub, dest_downstub, rt );

		//pose.dump_pdb("after1.pdb");

		/// this will now symmetrize the conformation:
		pose.set_jump( jump2, pose.jump( jump2 ) ); // prot-d

		//pose.dump_pdb("after2.pdb");

	} else { // using pos13 as anchorpos

		if ( reversed ) {
			// tricky...
			// we want to preserve the geometry that currently exists between
			Size const altprotanchor( (base_repeat-1)*repeatlen + anchorpos ),
				altjump2( pose.fold_tree().get_jump_that_builds_residue( altprotanchor ) );
			kinematics::Stub const upstub( pose.conformation().upstream_jump_stub( jump2 ) ),
				downstub( pose.conformation().downstream_jump_stub( altjump2 ) );
			pose.set_jump( jump2, kinematics::Jump( upstub, downstub ) );
		} else {
			pose.set_jump( jump2, pose.jump( jump2 ) ); // prot-d
		}

	}
	pose.set_jump( jump1, pose.jump( jump1 ) ); // bp
	pose.set_jump( jump3, pose.jump( jump3 ) ); // vrt-dna
	pose.set_jump( jump4, pose.jump( jump4 ) ); // vrt-vrt

	//pose.dump_pdb("after_symjumps.pdb");

}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// simple protocol: repeat 5 times: {relax, design}
//
void
design_test()
{
	Size const anchorpos( 13 ); // might want to make this configurable
	Size const cycles(5);

	string const output_tag( basic::options::option[ basic::options::OptionKeys::out::output_tag ] );
	Size const nstruct( basic::options::option[ basic::options::OptionKeys::out::nstruct] );

	//
	strings files( start_files() );

	numeric::random::random_permutation( files, numeric::random::rg() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {

		Pose pose;
		pose_from_pdb( pose, files[fi] );

		set_base_partner( pose );

		runtime_assert( pose.conformation().num_chains() == 3 );

		Size const nrepeat( pose.conformation().chain_end(2) - pose.conformation().chain_begin(2) + 1 );
		Size const repeatlen( 34 );
		Size base_repeat(0);
		if ( nrepeat == 7 ) base_repeat = 3;
		else if ( nrepeat == 8 ) base_repeat = 3;
		else utility_exit_with_message("bad nrepeat "+string_of(nrepeat));

		Size jump1, jump2, jump3, jump4;
		bool const flex_symdofs( false );
		setup_symmetric_pose( nrepeat, repeatlen, base_repeat, anchorpos, flex_symdofs, jump1, jump2, jump3, jump4, pose );


		// try scoring with a symmetric scoring function
		//
		// this is so bad: uses the symmetry_definition command line argument to decide whether to make the scorefxn
		// symmetric...
		//
		scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		runtime_assert( pose::symmetry::is_symmetric( *scorefxn ) );
		scorefxn->set_weight( chainbreak, 3.0 );

		scoring::methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
		options.exclude_DNA_DNA( false );
		scorefxn->set_energy_method_options( options );


		Real score( (*scorefxn)( pose ) );


		/// try "symmetrizing" the dofs

		pose.set_jump( jump4, pose.jump( jump4 ) ); // vrt-vrt
		pose.set_jump( jump3, pose.jump( jump3 ) ); // vrt-dna
		pose.set_jump( jump2, pose.jump( jump2 ) ); // prot-d
		pose.set_jump( jump1, pose.jump( jump1 ) ); // bp

		score = (*scorefxn)( pose );

		Pose const start_pose( pose );

		for ( Size n=1; n<= nstruct; ++n ) {

			pose = start_pose;

			for ( Size m=1; m<= cycles; ++m ) {

				{ /// try fast relax
					protocols::relax::FastRelax fastrelax( scorefxn, 0 );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					movemap->set_bb(true);
					movemap->set_chi(true);
					movemap->set_jump(true);
					movemap->set_jump(jump1,false); // turn off base-pair optimization

					fastrelax.set_movemap( movemap );

					fastrelax.apply( pose );

					score = (*scorefxn)( pose );

					cout << "cycle_scores: " << m << ' ' << F(9,3,score) << " postrelax " <<
						pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;
				}


				{
					// great! now let's try packing...
					pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
					task->initialize_from_command_line();
					task->or_include_current( true );

					bools is_protein( pose.total_residue(), false );
					for ( Size i=1; i<= chain_end(1,pose); ++i ) is_protein[i] = true;
					task->restrict_to_residues( is_protein );


					Size const nloop( 25 );
					protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( scorefxn, task, nloop );

					packmover.apply( pose );

					score = ( (*scorefxn)( pose ) );

					cout << "cycle_scores: " << m << ' ' << F(9,3,score) << " postpack " <<
						pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;
				}
			} // cycles

			score = (*scorefxn)(pose);

			string const outfilename( output_tag + "relaxdesign_" + filebase( files[fi] ) + "_N"+
				lead_zero_string_of( n, 4 )+".pdb" );

			string const repeatseq( pose.sequence().substr(0,repeatlen) );

			cout << "final_scores: " << F(9,3,score) << ' ' << outfilename << ' ' <<
				" repeatseq: " << repeatseq <<
				" weighted_energies: " << pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

			pose.dump_pdb(outfilename);
		} // nstruct
	} // fi
}


///////////////////////////////////////////////////////////////////////////////
void
rescore_test()
{
	strings const files( start_files() );

	scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	scoring::methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn->set_energy_method_options( options );


	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		Real const score( (*scorefxn)( pose ) );

		/// now delete the dna, rescore
		///
		Pose dnapose( pose ), protpose( pose );
		dnapose.conformation().delete_residue_range_slow( pose.conformation().chain_begin(1),
			pose.conformation().chain_end(1) );
		protpose.conformation().delete_residue_range_slow( pose.conformation().chain_begin(2),
			pose.conformation().chain_end(3) );
		Real const dnascore( (*scorefxn)( dnapose ) ), protscore( (*scorefxn)( protpose ) );

		runtime_assert( pose.residue(1).is_protein() );
		runtime_assert( pose.residue( pose.conformation().chain_begin(2)).is_DNA() );
		runtime_assert( pose.residue( pose.conformation().chain_begin(3)).is_DNA() );

		cout << "rescore: " << F(9,3,score) << ' ' << files[fi] <<
			" protscore: " << F(9,3,protscore) <<
			" dnascore: " << F(9,3,dnascore) <<
			" bindscore: " << F(9,3,score - protscore - dnascore ) <<
			" nondnascore: " << F(9,3,score - dnascore ) <<
			std::endl;
	}

}

///////////////////////////////////////////////////////////////////////////////
typedef std::pair<AA, AA > AA_Pair;
typedef map< AA_Pair, RTs > AllRTs;

void
load_rts_library(
	string const & filename,
	AllRTs & rts_map
)
{
	if ( !utility::file::file_exists( filename ) ) utility_exit_with_message("missing "+filename);

	ifstream data( filename.c_str() );

	string line;

	while ( getline( data, line ) ) {
		istringstream l( line );
		string tag;
		l >> tag;
		if ( !l.fail() && tag == "PROT_DNA_RT:" ) {
			RT rt;
			AA aa, na;
			l >> na >> aa >> rt;
			if ( !l.fail() ) {
				AA_Pair aana( aa, na );
				if ( rts_map.find( aana ) == rts_map.end() ) {
					TR.Trace << "New pair: " << aa << ' ' << na << endl;
					rts_map[ aana ];
				}
				rts_map.find( aana )->second.push_back( rt );
			} else {
				utility_exit_with_message("bad line: "+line);
			}
		}
	}

	for ( AllRTs::const_iterator it = rts_map.begin(); it != rts_map.end(); ++it ) {
		TR.Trace << "load_rts_library: " << it->first.first << ':' << it->first.second << I(6,it->second.size() ) << endl;
	}

	data.close();
}


///////////////////////////////////////////////////////////////////////////////////////
void
insert_interface_fragment_symmetrically(
	AllRTs const & all_rts,
	Pose & pose
)
{
	Size nrepeat, repeatlen, base_repeat, nres_protein( chain_end( 1, pose ) );
	parse_symminfo( pose, nrepeat, repeatlen, base_repeat );

	FoldTree const & f( pose.fold_tree() );

	bool found_jump( false );

	RT rt;
	AA_Pair aana;

	for ( Size i=1; i<= f.num_jump(); ++i ) {
		Size const pos1( f.upstream_jump_residue(i) ), pos2( f.downstream_jump_residue(i) );
		Residue const & rsd1( pose.residue( pos1 ) ), &rsd2( pose.residue( pos2 ) );
		if ( rsd1.is_DNA() && rsd2.is_protein() ) {
			Size const prot_repeatno( (pos2-1)/repeatlen+1 ), dna_repeatno( pos1-nres_protein );

			AA_Pair const aana_thisrepeat = make_pair( rsd2.aa(), rsd1.aa() );
			TR.Trace << "insert_interface_fragment_symmetrically: " <<
				aana_thisrepeat.first << ' ' << aana_thisrepeat.second <<
				I(4,dna_repeatno) << I(4,prot_repeatno) << I(4,pos1) << I(4,pos2) << endl;

			if ( !found_jump ) {
				aana = aana_thisrepeat;
				found_jump = true;
				runtime_assert( all_rts.count(aana) && all_rts.find( aana )->second.size() );
				rt = random_element( all_rts.find( aana )->second );
			} else {
				runtime_assert( aana == aana_thisrepeat );
			}

			id::StubID const stubid1( id::AtomID( rsd1.chi_atoms(1)[4], pos1 ),
				id::AtomID( rsd1.chi_atoms(1)[3], pos1 ),
				id::AtomID( rsd1.chi_atoms(1)[2], pos1 ) );
			id::StubID const stubid2( id::AtomID( rsd2.atom_index("N" ), pos2 ),
				id::AtomID( rsd2.atom_index("CA"), pos2 ),
				id::AtomID( rsd2.atom_index("C" ), pos2 ) );
			pose.conformation().set_stub_transform( stubid1, stubid2, rt );
		}
	}
	runtime_assert( found_jump );
}
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
void
simple_frag_sim(
	Size const cycles,
	ScoreFunction const & scorefxn_in,
	FragLib const & fraglib,
	AllRTs const & all_rts,
	bool const make_interface_moves,
	Sizes const & fragseq_poslist,
	//bool const freeze_anchorseq,
	Pose & pose
)
{
	protocols::moves::PyMOLMoverOP pymol(0);
	if ( option[ my_options::show_pymol ] ) {
		static Size counter(0);
		++counter;
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol->pymol_name("simple_frag_sim_"+lead_zero_string_of(counter,2));
		pymol->keep_history(true);
	}

	runtime_assert( pose::symmetry::is_symmetric( scorefxn_in ) );
	runtime_assert( pose::symmetry::is_symmetric( pose ) );

	using namespace protocols::moves;

	Real const jump_move_frequency( 0.25 );
	Real init_temp( 10.0 ), final_temp( 1.0 );
	Real init_protein_chainbreak_weight( 0.0 ), final_protein_chainbreak_weight( 1.0 );
	Real init_coordinate_constraint_weight( 0.0 ), final_coordinate_constraint_weight( 1.0 );


	Real const gamma = std::pow( (final_temp/init_temp), Real(1.0/( cycles-1 ) ) );
	Real temperature = init_temp / gamma;
	ScoreFunctionOP scorefxn( scorefxn_in.clone() );
	MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, temperature ) );

	Size outer_cycles( 5 ), inner_cycles( cycles/ outer_cycles );
	if ( dry_run() ) {
		outer_cycles = 2; // no division by zero...
		inner_cycles = 1;
	}

	mc->reset_counters();
	for ( Size m=1; m<= outer_cycles; ++m ) {
		Real const protein_chainbreak_weight( init_protein_chainbreak_weight +
			( m-1 ) * ( final_protein_chainbreak_weight -
			init_protein_chainbreak_weight ) / ( outer_cycles-1 ) );

		Real const coordinate_constraint_weight( init_coordinate_constraint_weight +
			( m-1 ) * ( final_coordinate_constraint_weight -
			init_coordinate_constraint_weight ) / ( outer_cycles-1 ) );


		scorefxn->set_weight( atom_pair_constraint, protein_chainbreak_weight );
		//scorefxn->set_weight( chainbreak, protein_chainbreak_weight );
		scorefxn->set_weight( coordinate_constraint, coordinate_constraint_weight );

		mc->score_function( *scorefxn );

		//devel::blab::attach_trajectory_observer( "tal_simple_frag_sim", *mc );

		mc->reset( pose ); // OK since we recovered low at the end of the previous loop

		for ( Size n=1; n<= inner_cycles; ++n ) {
			temperature *= gamma;
			mc->set_temperature( temperature );

			/// protein frag insertion
			Size const fragsize( ( m<=(2*outer_cycles)/3 ) ? 9 : 3 );
			insert_protein_fragment( fragsize, fragseq_poslist, fraglib, pose );
			string movename( "PROT"+string_of( fragsize ) );
			bool pass_mc( mc->boltzmann( pose, movename ) );
			if ( pymol && mc->mc_accepted() == protocols::moves::MCA_accepted_score_beat_low ) pymol->apply( pose );

			TR.Trace << "MC " << I(3,m) << I(6,n) << " mc_temp: " << F(9,3,mc->temperature() ) <<
				" cbweight: " << F(9,3,mc->score_function().get_weight( atom_pair_constraint ) ) <<
				" cbscore: " << F(9,3,pose.energies().total_energies()[ atom_pair_constraint ] ) <<
				" pass: " << pass_mc << ' ' <<
				F(9,3,mc->last_accepted_score()) << F(9,3,mc->lowest_score() ) << ' ' << movename << endl;

			if ( make_interface_moves && uniform() < jump_move_frequency ) {
				insert_interface_fragment_symmetrically( all_rts, pose );
				string movename( "INT" );
				bool pass_mc( mc->boltzmann( pose, movename ) );
				if ( pymol && mc->mc_accepted() == protocols::moves::MCA_accepted_score_beat_low ) pymol->apply( pose );

				TR.Trace << "MC " << I(3,m) << I(6,n) << " mc_temp: " << F(9,3,mc->temperature() ) <<
					" cbweight: " << F(9,3,mc->score_function().get_weight( chainbreak ) ) <<
					" cbscore: " << F(9,3,pose.energies().total_energies()[ chainbreak ] ) <<
					" pass: " << pass_mc << ' ' <<
					F(9,3,mc->last_accepted_score()) << F(9,3,mc->lowest_score() ) << ' ' << movename << endl;
			}

			if ( n%100 == 0 ) {
				mc->show_counters();
				//mc->show_counters( TR.Trace, false );
			}
		} // inner cycles
		mc->recover_low( pose );
	} // outer cycles

	mc->show_counters();
	//  mc->show_counters( TR.Trace, false );
	mc->recover_low( pose );

}


// ///////////////////////////////////////////////////////////////////////////////
// void
// count_atom_occupancy(
//            string const & filename,
//            Pose const & pose,
//            id::AtomID_Map< Size > & count
//            )
// {
//  core::pose::initialize_atomid_map( count, pose, Size(0) );
//  ifstream data( filename.c_str() );

//  string line;
//  while ( getline( data, line ) ) {
//   if ( line.substr(0,6) == "ATOM  " ) {
//    string const atomname( line.substr(12,4) );
//    char const chain( line[21] );
//    Size const resnum( int_of( line.substr(22,4 ) ) );
//    char const insertcode( line[26] );
//    Size const pos( pose.pdb_info()->pdb2pose( chain, resnum, insertcode  ) );
//    if ( !pos ) {
//     TR.Trace << "no pose pos: "<< chain << ' ' << resnum << ' '<< insertcode << endl;
//     continue;
//    }
//    Residue const & rsd( pose.residue( pos  ) );
//    if ( !rsd.has( atomname ) ) {
//     TR.Trace << "rsd missing atom " << rsd.name() << ' ' << atomname << endl;
//     continue;
//    }
//    ++count[ id::AtomID( rsd.atom_index( atomname ), pos ) ];
//   }
//  }
//  data.close();
// }



// ///////////////////////////////////////////////////////////////////////////////
// bools
// identify_residues_with_all_single_occupancy(
//                       Pose const & pose,
//                       string const filename
//                       )
// {
//  bools goodrsd( pose.total_residue(), false );

//  id::AtomID_Map< Size > atom_occupancy;
//  count_atom_occupancy( filename, pose, atom_occupancy );

//  for ( Size i=1; i<= pose.total_residue(); ++i ) {
//   Residue const & rsd( pose.residue(i) );
//   bool good( true );
//   for ( Size ii=1; ii<= rsd.nheavyatoms(); ++ii ) {
//    if ( rsd.atom_type(ii).is_virtual() ) continue;
//    if ( atom_occupancy[ id::AtomID( ii, i) ] != 1 ) good = false;
//   }
//   goodrsd[i] = good;
//   TR.Trace << "goodrsd: " << I(4,i) << ' ' << rsd.name1() << ' ' << good << ' ' << filename << endl;
//  }

//  return goodrsd;
// }

///////////////////////////////////////////////////////////////////////////////
void
sasa_fractions_test()
{
	strings const files( start_files() );

	Sizes aa_counts( num_canonical_aas, 0 );
	vector1< Reals > all_rsd_sasas;
	all_rsd_sasas.resize( num_canonical_aas );

	Reals avg_sasas;
	load_avg_residue_sasas( avg_sasas );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );
		Pose pose;
		pose_from_pdb( pose, filename );

		cout << "numchains: " << num_chains( pose ) << ' ' << filename << endl;

		// NOT USING GOODRSD FOR THIS !!!
		//bools const goodrsd( identify_residues_with_all_single_occupancy( pose, filename ) );
		set_ss_from_dssp( filename, pose ); // phil_io.hh

		/// compute rsd sasas
		Real const probe_radius( 1.4 );
		Reals rsd_sasas;
		protocols::sasa_scores::compute_residue_sasas_for_sasa_scores( probe_radius, pose, rsd_sasas );

		{
			bools subset( pose.total_residue(), false );


			for ( Size ch=1; ch<= num_chains( pose ); ++ch ) {
				Size const chbegin( chain_begin( ch, pose ) ), chlen( chain_end( ch, pose ) - chbegin + 1 );
				string const chainss( pose.secstruct().substr(chbegin-1,chlen) );

				if ( chainss.find('H') == string::npos ) continue;

				Size const first_H( chainss.find_first_of("H")+chbegin ), last_H( chainss.find_last_of("H") + chbegin );
				runtime_assert( pose.secstruct( first_H ) == 'H' );
				runtime_assert( pose.secstruct( last_H ) == 'H' );


				for ( Size i=first_H; i<= last_H; ++i ) {
					if ( pose.secstruct(i) == 'H' ) subset[i] = true;
					else if ( pose.secstruct(i) == 'L' ) {
						Size j(i);
						while ( pose.secstruct(j) != 'H' ) --j;
						Size const loopbegin( j+1 );
						j = i;
						while ( pose.secstruct(j) != 'H' ) ++j;
						Size const loopend( j-1 ), looplen( loopend - loopbegin+1 );
						subset[i] = ( looplen < 11 ); // hack
					}
				}
			}


			string substring;
			Size nres_subset( 0 );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const & rsd( pose.residue(i) );
				if ( subset[i] ) {
					runtime_assert( rsd.is_protein() );
					if ( !rsd.is_lower_terminus() && !rsd.is_upper_terminus() ) {
						++nres_subset;
						++aa_counts[ rsd.aa() ];
						all_rsd_sasas[ rsd.aa() ].push_back( rsd_sasas[i] );
					}
					substring.push_back(pose.secstruct(i) );
				} else substring.push_back( '-' );
			}

			if ( !nres_subset ) continue;

			cout << "ss " << pose.secstruct() << ' ' << filename << endl;
			cout << "su " << substring << ' ' << filename << endl;

			//get_sasas_by_atomtype( pose, subset, polar_sasa, nonpolar_sasa );
			Real polar_sasa(0.0), nonpolar_sasa(0.0), buried_polar_sasa( 0.0 ), buried_nonpolar_sasa( 0.0 );
			get_exposed_and_buried_sasas_by_atomtype( pose, subset, polar_sasa, nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa );

			/// normalize to nbr of residues
			polar_sasa /= nres_subset;
			nonpolar_sasa /= nres_subset;
			buried_polar_sasa /= nres_subset;
			buried_nonpolar_sasa /= nres_subset;

			cout << "SASA exposed_polar_sasa: " << F(9,3,polar_sasa) <<
				" exposed_nonpolar_sasa: "<< F(9,3,nonpolar_sasa) <<
				" buried_polar_sasa: "<< F(9,3,buried_polar_sasa) <<
				" buried_nonpolar_sasa: "<< F(9,3,buried_nonpolar_sasa) <<
				" nres_subset: " << I(4,nres_subset) << ' ' << filename << endl;
		}

	}


	/// show amino acid frequencies
	Size total(0);
	for ( Size i=1; i<= num_canonical_aas; ++i ) total += aa_counts[i];
	for ( Size i=1; i<= num_canonical_aas; ++i ) {
		cout << "AA_FREQ: " << AA(i) << F(9,6,Real( aa_counts[i] )/total) << I(12,aa_counts[i]) << I(12,total) << endl;
	}


	for ( Size i=1; i<= num_canonical_aas; ++i ) {
		AA const aa = AA(i);
		Reals const & rsd_sasas( all_rsd_sasas[aa] );

		if ( rsd_sasas.empty() ) continue;

		/// compute avg
		Real avg(0.0);
		for ( Size i=1; i<= rsd_sasas.size(); ++i ) avg += rsd_sasas[i];
		avg /= rsd_sasas.size();

		cout << "delta_sasa: " << F(9,3,avg-avg_sasas[aa]) << ' ' << aa << ' ' << I(10,rsd_sasas.size()) <<
			" avg: " << F(9,3,avg) << " pdbavg: " << F(9,3,avg_sasas[aa]) << endl;
	}
}


///////////////////////////////////////////////////////////////////////////////////////
Vector
helix_axis(
	Size const begin,
	Size const end,
	Pose const & pose
)
{
	Size const len( end - begin + 1 );
	runtime_assert( len >= 7 );

	runtime_assert( pose.secstruct().substr( begin-1, len ) == string( len, 'H' ) );

	Vector bxyz( 0,0,0 ), exyz( 0,0,0 );
	for ( Size i=0; i<4; ++i ) {
		bxyz += 0.25 * pose.residue( begin+i ).xyz("CA");
		exyz += 0.25 * pose.residue( end-i ).xyz("CA");
	}
	return ( exyz - bxyz ).normalized();

}

///////////////////////////////////////////////////////////////////////////////
void
helix_loop_test()
{
	Size const h1_len( 7 ), h2_len( 7 ), loop_maxlen( 7 );

	strings const files( start_files() );

	for ( Size fi =1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );

		Pose pose;
		//pose_from_pdb( pose, filename );
		cenpose_from_pdb( pose, filename );

		set_ss_from_dssp( filename, pose );

		add_termini_at_protein_chainbreaks( pose );

		string const ss( pose.secstruct() );

		string h1l( h1_len, 'H' ); h1l.push_back( 'L' );
		string lh2( h2_len, 'H' ); lh2 = string("L") + lh2;

		Size const nres( pose.total_residue() );


		for ( Size i=1; i<= nres - h1_len - h2_len; ++i ) {
			if ( ss.substr(i-1, h1_len+1) == h1l ) {
				Size looppos( i + h1_len );
				pbassert( ss[looppos-1] == 'L' );
				while ( looppos <= nres-h2_len && ss[looppos] == 'L' ) ++looppos;
				if ( looppos > nres-h2_len ) break;
				pbassert( ss[looppos-1] == 'L' && ss[looppos] != 'L' );
				Size const loopbegin( i + h1_len ), loopend( looppos ), looplen( loopend - loopbegin+1 );
				//h1_begin( i ), h2_begin( loopend+1 );
				if ( looplen <= loop_maxlen &&
						looppos < nres - h2_len && ss.substr(looppos-1,h2_len+1) == lh2 &&
						pose.chain(i) == pose.chain( loopend + h2_len ) ) {
					bool is_term( false );
					for ( Size j=i; j<= loopend+h2_len; ++j ) {
						if ( pose.residue(j).is_terminus() ) is_term = true;
					}
					if ( is_term ) continue;

					// found one
					Vector const h1_axis( helix_axis( i, i+h1_len-1, pose ) ),
						h2_axis( helix_axis( loopend+1, loopend+h2_len, pose ) );
					Real const dotprod( h1_axis.dot( h2_axis ) );
					cout << "helix_loop: " << F(9,6,dotprod ) << ' ' << I(3,looplen) << ' ' <<
						ss.substr(i-1,h1_len + h2_len + looplen ) << ' ' <<
						torsion2big_bin_string( i, loopend+h2_len, pose ) << ' ' <<
						pose.pdb_info()->chain(i) << ' ' << pose.pdb_info()->number(i) << ' ' <<
						filename << endl;
				}
			}
		}
	}
}


void
rescore_logfile_test()
{
	//Size const min_overlap( 30 );
	bool const show_sasa_scores( true );
	Size const nrepeat( 8 ), default_repeatlen( 34 ), base_repeat( 3 );

	Reals avg_sasas;
	load_avg_residue_sasas( avg_sasas );

	strings const files( start_files() );

	Size counter(0), filecounter(0);

	/// read the logfile
	string const logfile( start_file() );

	ifstream data( logfile.c_str() );

	// rhino1 symdes$ showfields tmp9 | more
	// 1 ./output/job9_1_000_hyraxB50_br2norm/job9_1_000_0.log:final_scores
	// 2 -70.994
	// 3 job9_1_000_0symfragrebuild_symmpose_02.pdb_Srevcendesign_N0001.pdb
	// 4 relax_rmsd:
	// 5 13.652
	// 6 passed_score_filter:
	// 7 0

	vector1< Reals > all_rsd_sasas;
	all_rsd_sasas.resize( num_canonical_aas );

	string line;
	while ( getline( data, line ) ) {
		++counter;
		if ( counter%10000 == 0 ) cerr << counter << ' ' << filecounter << endl;

		//istringstream l( line );
		strings const l( split_to_vector1( line ) );
		string const dirtag( l[1] ), filename( dirtag.substr( 0, dirtag.find_last_of( "/" ) ) +"/" + l[3] );
		bool const passed_score_filter( l[7] == "1" );

		if ( !passed_score_filter ) continue;

		if ( !utility::file::file_exists( filename ) ) {
			cout << "MISSING file: " << filename << endl;
			continue;
		}

		Size repeatlen( default_repeatlen );
		if ( l[46] == "repeatlen:" ) {
			repeatlen = int_of( l[47] );
		} else runtime_assert( line.find( "repeatlen:" ) == string::npos );

		bool const reversed( filename.find( "rev" ) != string::npos );
		Size const base_repeat_offset( !reversed ? (base_repeat-1)*repeatlen : (nrepeat-base_repeat)*repeatlen );

		cerr << "read " << filename << endl;
		++filecounter;

		Pose pose;
		pose_from_pdb( pose, filename );

		runtime_assert( chain_end( 1, pose ) == nrepeat * repeatlen );

		{
			Real const probe_radius( 1.4 );
			Reals rsd_sasas;
			protocols::sasa_scores::compute_residue_sasas_for_sasa_scores( probe_radius, pose, rsd_sasas );

			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) {
				AA const aa( pose.residue(i).aa() );
				all_rsd_sasas[ aa ].push_back( rsd_sasas[i] );
			}
		}

		Real sasapack_score_single_repeat(0.0), norme_score_single_repeat( 0.0 ), sasa14_normalized_single_repeat(0.0);
		if ( show_sasa_scores ) {
			Real average_sasapack, average_normsasa, average_avge;
			Reals rsd_sasapack_scores, rsd_norme_scores, rsd_sasa14_normalized;
			protocols::sasa_scores::compute_sasapack_scores( pose, rsd_sasapack_scores, rsd_sasa14_normalized,
				average_sasapack, average_normsasa );
			protocols::sasa_scores::compute_avge_scores( pose, rsd_norme_scores, rsd_sasa14_normalized, average_avge,
				average_normsasa );

			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) {
				sasapack_score_single_repeat += rsd_sasapack_scores[i];
				norme_score_single_repeat += rsd_norme_scores[i];
				sasa14_normalized_single_repeat += rsd_sasa14_normalized[i];
			}
			sasapack_score_single_repeat /= repeatlen;
			norme_score_single_repeat /= repeatlen;
			sasa14_normalized_single_repeat /= repeatlen;
		}

		Real polar_sasa(0.0), nonpolar_sasa(0.0), buried_polar_sasa(0.0), buried_nonpolar_sasa(0.0),
			unbound_polar_sasa(0.0), unbound_nonpolar_sasa(0.0),
			unbound_buried_polar_sasa(0.0), unbound_buried_nonpolar_sasa(0.0);
		{
			bools subset( pose.total_residue(), false );
			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;
			get_exposed_and_buried_sasas_by_atomtype( pose, subset, polar_sasa, nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa );
			polar_sasa /= repeatlen;
			nonpolar_sasa /= repeatlen;
			buried_polar_sasa /= repeatlen;
			buried_nonpolar_sasa /= repeatlen;

			Pose unbound_pose( pose );
			unbound_pose.conformation().delete_residue_range_slow( chain_end( 1, pose )+1, pose.total_residue() );

			get_exposed_and_buried_sasas_by_atomtype( unbound_pose, subset, unbound_polar_sasa, unbound_nonpolar_sasa,
				unbound_buried_polar_sasa, unbound_buried_nonpolar_sasa );
			unbound_polar_sasa /= repeatlen;
			unbound_nonpolar_sasa /= repeatlen;
			unbound_buried_polar_sasa /= repeatlen;
			unbound_buried_nonpolar_sasa /= repeatlen;
		}

		// abego string, use internal repeat to get well-defined torsions at cutpoints
		string const repeatbb( torsion2big_bin_string( base_repeat_offset+1, base_repeat_offset+repeatlen, pose ) );
		runtime_assert( repeatbb.size() == repeatlen );

		cout << line;

		if ( show_sasa_scores ) {
			cout << " sasapack_score_single_repeat: " << F(9,3,sasapack_score_single_repeat) <<
				" norme_score_single_repeat: " << F(9,3,norme_score_single_repeat) <<
				" sasa14_normalized_single_repeat: " << F(9,3,sasa14_normalized_single_repeat) <<
				" repeatbb: " << repeatbb; // silly that this is in here...
		}
		cout << " polar_sasa: " << F(9,3,polar_sasa) <<
			" nonpolar_sasa: "<< F(9,3,nonpolar_sasa) <<
			" buried_polar_sasa: " << F(9,3,buried_polar_sasa) <<
			" buried_nonpolar_sasa: "<< F(9,3,buried_nonpolar_sasa) <<
			" unbound_polar_sasa: " << F(9,3,unbound_polar_sasa) <<
			" unbound_nonpolar_sasa: "<< F(9,3,unbound_nonpolar_sasa) <<
			" unbound_buried_polar_sasa: " << F(9,3,unbound_buried_polar_sasa) <<
			" unbound_buried_nonpolar_sasa: "<< F(9,3,unbound_buried_nonpolar_sasa);

		/// newline
		cout << endl;
	}

	/// loop over all the aa types
	for ( Size i=1; i<= num_canonical_aas; ++i ) {
		AA const aa = AA(i);

		Reals const & rsd_sasas( all_rsd_sasas[aa] );

		if ( rsd_sasas.empty() ) continue;

		/// compute avg
		Real avg(0.0);
		for ( Size i=1; i<= rsd_sasas.size(); ++i ) avg += rsd_sasas[i];
		avg /= rsd_sasas.size();

		cout << "delta_sasa: " << F(9,3,avg-avg_sasas[aa]) << ' ' << aa << ' ' << I(10,rsd_sasas.size()) <<
			" avg: " << F(9,3,avg) << " pdbavg: " << F(9,3,avg_sasas[aa]) << endl;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
rescore_unbound_logfile_test()
{
	Size const base_repeat( 3 );
	Size const turn_buffer( 7 );

	bool const recompute_psipred( false ), show_sasa_scores( false ), show_helix_twists( false ),
		show_unsatisfied( true );

	Reals avg_sasas;
	load_avg_residue_sasas( avg_sasas );

	strings const files( start_files() );

	Size counter(0), filecounter(0);

	/// read the logfile
	string const logfile( start_file() );

	ifstream data( logfile.c_str() );


	vector1< Reals > all_rsd_sasas;
	all_rsd_sasas.resize( num_canonical_aas );


	vector1< Vectors > all_helix_coords;
	map< string, vector1< Vectors > > all_turn_coords;

	{ //read info on bb coords
		string const pdb_coords_file( option[ my_options::pdb_coords_file ] );
		ifstream data( pdb_coords_file.c_str() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			if ( l[1] == "helix_coords:" ) {
				Size const helixlen( int_of( l[2] ) );
				Vectors coords;
				for ( Size i=1; i<= helixlen; ++i ) {
					coords.push_back( Vector( float_of( l[ 3*i+1 ] ), float_of( l[ 3*i+2 ] ), float_of( l[3*i+3 ] ) ) );
				}
				all_helix_coords.push_back( coords );
			} else if ( l[1] == "turn_coords:" ) {
				string const turn( l[2] );
				if ( all_turn_coords.find( turn ) == all_turn_coords.end() ) all_turn_coords[ turn ];
				Vectors coords;
				Size const coordslen( int_of( l[6] ) );
				runtime_assert( coordslen == turn.size() + 2*turn_buffer );
				for ( Size i=1; i<= coordslen; ++i ) {
					coords.push_back( Vector( float_of( l[ 3*i+4 ] ), float_of( l[ 3*i+5 ] ), float_of( l[3*i+6 ] ) ) );
				}
				all_turn_coords.find( turn )->second.push_back( coords );
				//TR.Trace << "read turn: " << turn << ' ' << all_turn_coords.find( turn )->second.size() << endl;

			}

		}
		data.close();
		TR.Trace << "Read " << all_helix_coords.size() << " helix coords from file " << pdb_coords_file << endl;

	}



	string line;
	while ( getline( data, line ) ) {
		++counter;
		if ( counter%10000 == 0 ) cerr << counter << ' ' << filecounter << endl;

		//istringstream l( line );
		strings const l( split_to_vector1( line ) );
		string const dirtag( l[1] ), filename( dirtag.substr( 0, dirtag.find_last_of( "/" ) ) +"/" + l[3] );

		bool passed_score_filter( false ), relaxed( true );
		Size repeatlen(0), nrepeat(0), helix1_len(0), helix2_len(0);
		string repeatss, repeatbb, turn1, turn2;
		Real handedness(0), helix1_dist(0), helix2_dist(0);


		for ( Size i=1; i< l.size(); ++i ) {
			if      ( l[i] == "passed_score_filter:" ) passed_score_filter = ( l[i+1] == "1" );
			else if ( l[i] == "passed_fullatom_score_filter:" ) passed_score_filter = ( l[i+1] == "1" );
			else if ( l[i] == "handedness:" ) handedness = float_of( l[i+1] );
			else if ( l[i] == "repeatlen:" ) repeatlen = int_of( l[i+1] );
			else if ( l[i] == "nrepeat:" ) nrepeat = int_of( l[i+1] );
			else if ( l[i] == "relaxed:" ) relaxed = ( l[i+1] == "1" );
			else if ( l[i] == "repeatss:" ) repeatss = l[i+1];
			else if ( l[i] == "repeatbb:" ) repeatbb = l[i+1];
			else if ( l[i] == "helix1_len:" ) helix1_len = int_of( l[i+1] );
			else if ( l[i] == "helix2_len:" ) helix2_len = int_of( l[i+1] );
			else if ( l[i] == "helix1_dist:" ) helix1_dist = float_of( l[i+1] );
			else if ( l[i] == "helix2_dist:" ) helix2_dist = float_of( l[i+1] );
			else if ( l[i] == "turn1:" ) turn1 = l[i+1];
			else if ( l[i] == "turn2:" ) turn2 = l[i+1];
		}

		if ( !relaxed ) passed_score_filter = false;

		if ( !passed_score_filter ) continue;


		runtime_assert( repeatlen * nrepeat > 0 );

		if ( !utility::file::file_exists( filename ) ) {
			cout << "MISSING file: " << filename << endl;
			continue;
		}

		{ //hacking-- test the pdb strain code

			while ( true ) {
				// shrink turn1 if necessary
				turn1 = repeatbb.substr( helix1_len, turn1.size() );
				if ( turn1[0] == 'A' ) {
					TR.Trace << "trimbegin: " << turn1 << endl;
					helix1_len += 1;
					turn1.erase( 0, 1 );
					continue;
				}
				if ( turn1[turn1.size()-1 ] == 'A' ) {
					TR.Trace << "trimend: " << turn1 << endl;
					helix2_len += 1;
					turn1.erase( turn1.end()-1 );
					continue;
				}
				break;
			}

			//bool fail( false );
			while ( true ) {
				// shrink turn2 if necessary
				turn2 = repeatbb.substr( helix1_len+turn1.size()+helix2_len, turn2.size() );
				if ( turn2[0] == 'A' ) {
					TR.Trace << "trimbegin: " << turn2 << endl;
					helix2_len += 1;
					turn2.erase( 0, 1 );
					continue;
				}
				break;
			}

			//
			Sizes
				helixbegins( make_vector1( Size(1), helix1_len + turn1.size() + 1 ) ),
				helixends  ( make_vector1( helix1_len  , repeatlen-turn2.size() ) ),
				turnbegins ( make_vector1( helix1_len+1, repeatlen-turn2.size()+1 ) ),
				turnends   ( make_vector1( helix1_len+turn1.size(), repeatlen ) );

			{
				string tmpbb( repeatbb.substr( repeatlen-1 ) + repeatbb + repeatbb.substr(0,1) ); // 1-indexed now
				runtime_assert( tmpbb.size() == repeatlen+2 );

				bool fail( false );
				for ( Size r=1; r<= 2; ++r ) {
					Size const helixlen( helixends[r] - helixbegins[r] + 1 );
					if ( tmpbb.substr( helixbegins[r], helixlen ) != string( helixlen, 'A' ) ||
							tmpbb[ helixbegins[r]-1] == 'A' ||
							tmpbb[ helixends[r]+1 ] == 'A' ) {
						fail = true;
						TR.Trace << "badbb " << r << ' ' << tmpbb.substr( helixbegins[r], helixlen ) << ' ' <<
							tmpbb[ helixbegins[r]-1 ] << ' ' << tmpbb[ helixends[r]+1 ] << ' ' << repeatbb << ' ' <<
							turn1 << ' ' << turn2 << endl;
					}
				}
				if ( fail ) {
					continue;
				}

				string oldturn1( turn1 ), oldturn2( turn2 );
				turn1 = tmpbb.substr( turnbegins[1], turn1.size() );
				turn2 = tmpbb.substr( turnbegins[2], turn2.size() );
				if ( oldturn1 != turn1 || oldturn2 != turn2  ) {
					TR.Trace << "turnchange " << turn1 << ' '<< oldturn1 << ' ' << turn2 << ' ' <<  oldturn2 << endl;
				}
			}

			if ( turn1 != "GBB" || turn2 != "GBB" ) continue; // HACKING

			TR.Trace << "read " << filename << endl;
			vector1< Vectors > chain_coords( read_CA_coords_from_complex_file( filename ) );
			Vectors const & ca_coords( chain_coords[1] );
			runtime_assert( ca_coords.size() == repeatlen * nrepeat );

			Size const helix_window( 7 );
			// compute rmsds for helices and turns
			if ( !( helix1_len >= turn_buffer && helix2_len >= turn_buffer &&
					helix1_len >= helix_window && helix2_len >= helix_window ) ) {
				TR.Trace << "bad helix lens " << line << endl;
				continue;
			}
			runtime_assert( helix1_len >= turn_buffer && helix2_len >= turn_buffer );


			// look at helix rmsd
			Real helix_avg_rmsd( 0 );
			for ( Size r=1; r<= 2; ++r ) {
				Size const helixbegin( helixbegins[r] ), helixend( helixends[r] ), helixlen( helixend - helixbegin + 1 );

				Real avg_avg_rmsd( 0 ), avg_med_rmsd( 0 );

				for ( Size i=1; i<= helixlen - helix_window+1; ++i ) {
					Vectors coords;
					for ( Size pos = helixbegin+i-1; pos<= helixbegin+i-1+helix_window-1; ++pos ) {
						coords.push_back( ca_coords[pos] );
					}
					Size const max_rmsds( 250 );
					Reals rmsds;
					for ( Size ii=1; ii<= all_helix_coords.size(); ++ii ) {
						if ( all_helix_coords[ii].size() >= i+helix_window-1 ) {
							Vectors pdb_coords;
							for ( Size j=i; j<= i+helix_window-1; ++j ) pdb_coords.push_back( all_helix_coords[ii][j] );
							runtime_assert( coords.size() == pdb_coords.size() );
							rmsds.push_back( numeric::model_quality::calc_rms( coords, pdb_coords ) );
							if ( rmsds.size() >= max_rmsds ) break;
						}
					}

					std::sort( rmsds.begin(), rmsds.end() );

					Real const median_rmsd( rmsds[ rmsds.size()/2 ] );

					Size const n_outliers( 15 );
					runtime_assert( n_outliers < rmsds.size() );
					for ( Size i=1; i<= n_outliers; ++i ) rmsds.pop_back();

					Real avg_rmsd( 0 );
					for ( Size i=1; i<= rmsds.size(); ++i ) avg_rmsd += rmsds[i];
					avg_rmsd/= rmsds.size();

					avg_avg_rmsd += avg_rmsd;
					avg_med_rmsd += median_rmsd;
				} // loop over helix windows
				Size const nwindows( helixlen - helix_window + 1 );
				avg_avg_rmsd /= nwindows;
				avg_med_rmsd /= nwindows;

				helix_avg_rmsd += avg_avg_rmsd;

				cout << "helix_rmsd: " << r << ' ' << helixlen << F(9,3,avg_avg_rmsd) << F(9,3,avg_med_rmsd) << ' ' <<
					filename << endl;
			} // r=1,2



			Real turn_avg_rmsd( 0. );
			for ( Size r=1; r<= 2; ++r ) {
				string const turn( r==1 ? turn1 : turn2 );
				Size const turnbegin( r == 1 ? helix1_len+1 : helix1_len+turn1.size()+helix2_len+1 ),
					coordsbegin( turnbegin-turn_buffer ), coordsend( turnbegin+turn.size()+turn_buffer-1);
				Vectors coords;
				for ( Size i= coordsbegin; i<= coordsend; ++i ) coords.push_back( ca_coords[i] );

				// compute rmsds to all turn coords
				vector1< Vectors > const & pdb_coords( all_turn_coords.find( turn )->second );
				TR.Trace << "pdb_coords: " << pdb_coords.size() << endl;
				Reals rmsds;
				for ( Size i=1; i<= pdb_coords.size(); ++i ) {
					runtime_assert( coords.size() == pdb_coords[i].size() );
					rmsds.push_back( numeric::model_quality::calc_rms( coords, pdb_coords[i] ) );
				}
				std::sort( rmsds.begin(), rmsds.end() );

				Real const median_rmsd( rmsds[ rmsds.size()/2 ] );

				Size const n_outliers( 5 );
				runtime_assert( n_outliers < rmsds.size() );
				for ( Size i=1; i<= n_outliers; ++i ) rmsds.pop_back();

				Real avg_rmsd( 0 );
				for ( Size i=1; i<= rmsds.size(); ++i ) avg_rmsd += rmsds[i];
				avg_rmsd/= rmsds.size();

				turn_avg_rmsd += avg_rmsd;

				cout << "turn_rmsd: " << r << ' ' << turn << F(9,3,avg_rmsd) << F(9,3,median_rmsd) << ' ' <<
					filename << endl;
			}

			char const hand( handedness<0 ? 'L' : 'R' );

			cout << "all_rmsd " << F(9,3,helix_avg_rmsd+turn_avg_rmsd) <<
				" helix: " << F(9,3,helix_avg_rmsd) <<
				" turn: " << F(9,3,turn_avg_rmsd) << ' ';
			if ( helix1_dist < helix2_dist ) cout << helix1_len << turn1 << helix2_len << turn2 << hand;
			else cout << helix2_len << turn2 << helix1_len << turn1 << hand;
			cout << ' ' << filename << endl;


			continue;
		} // scope

		Size const base_repeat_offset( (base_repeat-1)*repeatlen );

		cerr << "read " << filename << endl;
		++filecounter;

		Pose pose;
		pose_from_pdb( pose, filename );

		runtime_assert( num_chains( pose ) == 1 );
		runtime_assert( chain_end( 1, pose ) == nrepeat * repeatlen );

		{
			Real const probe_radius( 1.4 );
			Reals rsd_sasas;
			protocols::sasa_scores::compute_residue_sasas_for_sasa_scores( probe_radius, pose, rsd_sasas );

			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) {
				AA const aa( pose.residue(i).aa() );
				all_rsd_sasas[ aa ].push_back( rsd_sasas[i] );
			}
		}

		Real sasapack_score_single_repeat(0.0), norme_score_single_repeat( 0.0 ), sasa14_normalized_single_repeat(0.0);

		Real polar_sasa(0.0), nonpolar_sasa(0.0), buried_polar_sasa(0.0), buried_nonpolar_sasa(0.0),
			buried_nonpolar_sasa_sc(0.0), buried_nonpolar_rsd_sasa_sc(0.0), packstat_score(0.0);
		bools subset( pose.total_residue(), false );
		for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;
		compute_sasa_scores_for_subset( subset, pose, sasapack_score_single_repeat, norme_score_single_repeat,
			sasa14_normalized_single_repeat, polar_sasa, nonpolar_sasa,
			buried_polar_sasa, buried_nonpolar_sasa,
			buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
			packstat_score );

		//Real const n2cdist( pose.residue(1).xyz("N").distance( pose.residue(chain_end(1,pose)).xyz("C") ) );

		/// compute psipred ss using single sequence
		string repeatsspred( "-" );
		Real sspred_match(0.0);
		if ( recompute_psipred ) {
			{ // run psipred to re-predict the secondary structure
				vector1< Reals > pred_HEL;
				string const sequence( pose.sequence().substr(0,chain_end(1,pose) ) );
				//string const secstruct( pose.secstruct().substr(0,chain_end(1,pose) ) );
				string sspred;

				run_psipred( sequence, sspred, pred_HEL );

				if ( sspred.size() == sequence.size() ) { // success
					repeatsspred = sspred.substr( (base_repeat-1)*repeatlen, repeatlen );

					if ( repeatss.size() == repeatlen ) { // logfile had repeatss info
						runtime_assert( sspred.size() == nrepeat * repeatlen );
						for ( Size i=1; i<= sspred.size(); ++i ) {
							char const ss( repeatss[ (i-1)%repeatlen ] );
							Real const prediction( ss == 'H' ? pred_HEL[i][1] : ( ss == 'E' ? pred_HEL[i][2] : pred_HEL[i][3] ) );
							sspred_match += prediction;
						}
						sspred_match /= sspred.size();
					}
				}
			}
		}

		Real helix1_twist(0), helix2_twist(0);
		//Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);
		if ( helix1_len ) {
			runtime_assert( helix1_len + turn1.size() + helix2_len + turn2.size() == repeatlen );

			compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
				helix1_dist, helix1_twist, helix2_dist, helix2_twist );
		}


		string const unsatstring( show_unsatisfied ? get_buried_unsatisfied_string( pose ) : string("unused") );


		cout << line;
		if ( show_sasa_scores ) {
			cout <<
				" sasapack_score_single_repeat: " << F(9,3,sasapack_score_single_repeat) <<
				" norme_score_single_repeat: " << F(9,3,norme_score_single_repeat) <<
				" sasa14_normalized_single_repeat: " << F(9,3,sasa14_normalized_single_repeat) <<
				" polar_sasa: " << F(9,3,polar_sasa) <<
				" nonpolar_sasa: "<< F(9,3,nonpolar_sasa) <<
				" buried_polar_sasa: " << F(9,3,buried_polar_sasa) <<
				" buried_nonpolar_sasa: "<< F(9,3,buried_nonpolar_sasa) <<
				" buried_nonpolar_sasa_sc: "<< F(9,3,buried_nonpolar_sasa_sc) <<
				" buried_nonpolar_rsd_sasa_sc: "<< F(9,3,buried_nonpolar_rsd_sasa_sc) <<
				" packstat_score: " << F(9,3,packstat_score);
		}
		//" n2cdist: " << F(9,3,n2cdist) <<
		if ( recompute_psipred ) {
			cout <<
				" repeatsspred: " << repeatsspred <<
				" sspred_match: " << F(9,3,sspred_match);
		}

		if ( show_helix_twists ) {
			cout <<
				" helix1_dist: " << F(9,3,helix1_dist) <<
				" helix1_twist: " << F(9,3,helix1_twist) <<
				" helix2_dist: " << F(9,3,helix2_dist) <<
				" helix2_twist: " << F(9,3,helix2_twist);
		}

		if ( show_unsatisfied ) cout << ' ' << unsatstring;

		cout << endl;
	}

	return; ///////// HACKING


	/// loop over all the aa types
	for ( Size i=1; i<= num_canonical_aas; ++i ) {
		AA const aa = AA(i);

		Reals const & rsd_sasas( all_rsd_sasas[aa] );

		if ( rsd_sasas.empty() ) continue;

		/// compute avg
		Real avg(0.0);
		for ( Size i=1; i<= rsd_sasas.size(); ++i ) avg += rsd_sasas[i];
		avg /= rsd_sasas.size();

		cout << "delta_sasa: " << F(9,3,avg-avg_sasas[aa]) << ' ' << aa << ' ' << I(10,rsd_sasas.size()) <<
			" avg: " << F(9,3,avg) << " pdbavg: " << F(9,3,avg_sasas[aa]) << endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
redesign_unbound_logfile_test()
{
	Size const design_cycles( option[ my_options::design_cycles ] );
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( get_score_function_from_command_line() );
	//  ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );

	adjust_ref_weights_from_command_line( *fa_scorefxn );
	//Size const base_repeat( option[ my_options::base_repeat ] );

	/// read the logfile
	strings all_lines;
	Sizes goodlines;

	{ // read the logfile, store all the goodlines
		string const logfile( start_file() );

		ifstream data( logfile.c_str() );

		string line;
		while ( getline( data, line ) ) {
			all_lines.push_back( line ); // every single line, not just final_scores lines

			//istringstream l( line );
			bool passed_score_filter( false );
			strings const l( split_to_vector1( line ) );
			for ( Size i=1; i< l.size(); ++i ) {
				if ( l[i] == "passed_score_filter:" ) passed_score_filter = ( l[i+1] == "1" );
			}
			if ( passed_score_filter ) goodlines.push_back( all_lines.size() );
		}
	}

	numeric::random::random_permutation( goodlines, numeric::random::rg() );


	for ( Size gi=1; gi<= goodlines.size(); ++gi ) {
		Size const line_index( goodlines[gi] );

		string const simfile( shared_output_tag() + "_line"+string_of( line_index)+".work" );
		string const worktag( "redesign" );

		Size const first_n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );

		if ( first_n > nstruct() ) continue;

		//istringstream l( line );
		strings const l( split_to_vector1( all_lines[ line_index ] ) );
		string const dirtag( l[1] ), filename( dirtag.substr( 0, dirtag.find_last_of( "/" ) ) +"/" + l[3] );

		bool passed_score_filter( false );
		Size repeatlen(0), nrepeat(0), nrepeat_asym(0);
		Sizes design_poslist, repack_poslist, polar_poslist;

		Size base_repeat( 0 );
		for ( Size i=1; i< l.size(); ++i ) {
			if      ( l[i] == "passed_score_filter:" ) passed_score_filter = ( l[i+1] == "1" );
			else if ( l[i] == "repeatlen:" ) repeatlen = int_of( l[i+1] );
			else if ( l[i] == "nrepeat:" ) nrepeat = int_of( l[i+1] );
			else if ( l[i] == "nrepeat_asym:" ) nrepeat_asym = int_of( l[i+1] );
			else if ( l[i] == "base_repeat:" ) base_repeat = int_of( l[i+1] );
			else if ( l[i] == "design_poslist:" ) {
				strings const ll( split_to_vector1( l[i+1], "," ) );
				for ( Size ii=1; ii<= ll.size(); ++ii ) {
					runtime_assert( is_int( ll[ii] ) );
					design_poslist.push_back( int_of( ll[ii] ) );
				}
			} else if ( l[i] == "repack_poslist:" ) {
				strings const ll( split_to_vector1( l[i+1], "," ) );
				for ( Size ii=1; ii<= ll.size(); ++ii ) {
					runtime_assert( is_int( ll[ii] ) );
					repack_poslist.push_back( int_of( ll[ii] ) );
				}
			} else if ( l[i] == "polar_poslist:" ) {
				strings const ll( split_to_vector1( l[i+1], "," ) );
				for ( Size ii=1; ii<= ll.size(); ++ii ) {
					runtime_assert( is_int( ll[ii] ) );
					polar_poslist.push_back( int_of( ll[ii] ) );
				}
			}
		}

		if ( !base_repeat ) base_repeat = option[ my_options::base_repeat ];

		runtime_assert( passed_score_filter );
		runtime_assert( repeatlen * nrepeat * nrepeat_asym > 0 );
		runtime_assert( !design_poslist.empty() );
		runtime_assert( nrepeat%nrepeat_asym == 0 ); // e.g. 6 and 2

		if ( !utility::file::file_exists( filename ) ) {
			cout << "MISSING file: " << filename << endl;
			continue;
		}

		//Size const base_repeat_offset( (base_repeat-1)*repeatlen );

		Pose pose;
		pose_from_pdb( pose, filename );

		runtime_assert( num_chains( pose ) == 1 );
		runtime_assert( chain_end( 1, pose ) == nrepeat * repeatlen );

		{
			//Residue const & lastrsd( pose.residue( nres_protein ) );
			remove_upper_terminus_type_from_pose_residue( pose, pose.total_residue() );
			ResidueOP vrtrsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
			pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...
			pose.conformation().insert_chain_ending( pose.total_residue()-1 );

			kinematics::FoldTree f( pose.total_residue() );
			f.reorder( pose.total_residue() );
			pose.fold_tree( f );
		}


		Pose const startpose( pose );

		bool first_time_through( true );

		while ( true ) {
			Size const n( first_time_through ?
				first_n : get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( n > nstruct() ) break;
			first_time_through = false;

			pose = startpose;

			{ /// set up the symm info object
				///
				conformation::symmetry::SymmetryInfo symminfo;

				//Size const nclones( nrepeat/nrepeat_asym - 1 ); // e.g. 2, if nrepeat = 6 and nrepeat_asym = 2
				for ( Size i=1; i<= nrepeat; ++i ) {
					// is this a clone or base?
					if ( i >= base_repeat && i <base_repeat+nrepeat_asym ) {
						// base
					} else {
						Size const ibase( (i-1)%nrepeat_asym + base_repeat );
						Size const offset( (i-1)*repeatlen ), base_offset( (ibase-1)*repeatlen );

						for ( Size j=1; j<= repeatlen; ++j ) {
							symminfo.add_bb_clone ( base_offset+j, offset+j );
							symminfo.add_chi_clone( base_offset+j, offset+j );
						}
					}
				}


				symminfo.num_virtuals( 1 );
				symminfo.set_use_symmetry( true );
				symminfo.set_flat_score_multiply( pose.total_residue(), 1 );

				/// repeats prior to base_repeat have score_multiply ==> 0
				for ( Size i=1; i<base_repeat; ++i ) {
					Size const offset( (i-1)*repeatlen );
					for ( Size j=1; j<= repeatlen; ++j ) symminfo.set_score_multiply( offset+j, 0 );
				}

				symminfo.update_score_multiply_factor();

				pose::symmetry::make_symmetric_pose( pose, symminfo );
			} // setup symmetry info object

			/// now design/relax
			clock_t starttime( clock() );
			for ( Size m=1; m<= design_cycles; ++m ) {

				{ // design
					pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
					task->initialize_from_command_line();
					task->or_include_current( true );

					runtime_assert( pose.total_residue() == nrepeat * repeatlen+1 ); // sanity
					task->nonconst_residue_task( pose.total_residue() ).prevent_repacking();
					for ( Size i=1; i<= nrepeat*repeatlen; ++i ) {
						Size const rpos( (i-1)%repeatlen+1 );
						bool const is_designable( has_element( design_poslist, rpos ) );
						bool const is_repackable( has_element( repack_poslist, rpos ) );
						bool const is_polar( has_element( polar_poslist, rpos ) );

						if ( !is_designable && !is_repackable ) {
							task->nonconst_residue_task( i ).prevent_repacking();
						} else if ( !is_designable ) {
							task->nonconst_residue_task( i ).restrict_to_repacking();
						} else if ( is_polar ) {
							///
							bools allowed_aas( num_canonical_aas, false );
							string const polaraas( "STDENQKRH" );
							for ( string::const_iterator c=polaraas.begin(); c!= polaraas.end(); ++c ) {
								allowed_aas[ aa_from_oneletter_code( *c ) ] = true;
							}
							task->nonconst_residue_task( i ).restrict_absent_canonical_aas( allowed_aas );
						}
					} // i=1,nres

					if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
						protocols::task_operations::LimitAromaChi2Operation lp_op;
						lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negatives)
						lp_op.apply( pose, *task );
					}
					Size const nloop( 25 );
					ScoreFunctionOP design_scorefxn(0);
					if ( ( option[ my_options::use_softrep_for_early_design ] && m<= design_cycles/2 ) ||
							option[ my_options::use_softrep_for_design ] ) {
						design_scorefxn = ScoreFunctionFactory::create_score_function( "soft_rep_design.wts" );
						adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
					} else {
						design_scorefxn = fa_scorefxn; // already adjusted refwts
					}
					protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
					if ( !dry_run() ) packmover.apply( pose );
				}

				{ // now sidechain relax
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					movemap->set_bb (false);
					movemap->set_chi(false);

					for ( Size i=1; i<= nrepeat*repeatlen; ++i ) {
						Size const rpos( (i-1)%repeatlen+1 );
						bool const is_designable( has_element( design_poslist, rpos ) );
						bool const is_repackable( has_element( repack_poslist, rpos ) );
						if ( is_designable || is_repackable ) movemap->set_chi( i, true );
					}
					fastrelax.set_movemap( movemap );
					if ( !dry_run() ) fastrelax.apply( pose );
				}
			} // cycles

			Real const final_score( (*fa_scorefxn)( pose ) );

			clock_t stoptime = clock();
			Real const simtime = ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes


			Real sasapack_score(0.0), norme_score( 0.0 ), normsasa_score( 0.0 ), exposed_polar_sasa(0.0),
				exposed_nonpolar_sasa(0.0), buried_polar_sasa(0.0), buried_nonpolar_sasa(0.0),
				buried_nonpolar_sasa_sc(0.0), buried_nonpolar_rsd_sasa_sc(0.0), packstat_score(0.0);
			bools subset( pose.total_residue(), false );
			Size const base_repeat_offset( (base_repeat-1)*repeatlen );
			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+nrepeat_asym*repeatlen; ++i ) subset[i] = true;
			compute_sasa_scores_for_subset( subset, pose, sasapack_score, norme_score, normsasa_score,
				exposed_polar_sasa, exposed_nonpolar_sasa, buried_polar_sasa, buried_nonpolar_sasa,
				buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score );


			string const unsatstring( get_buried_unsatisfied_string( pose ) );

			bool const passed_score_filter
				( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );
			string const outfilename( output_tag() + "redesign_L"+ string_of( line_index )+
				"_" + string_of( nrepeat ) +"_"+ string_of( repeatlen ) +
				"_N"+ lead_zero_string_of( n, 4 )+".pdb" );

			cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
				" nrepeat: " << nrepeat <<
				" repeatlen: " << repeatlen <<
				" fullseq: " << pose.sequence() <<
				" start_file: " << filename <<
				" passed_score_filter: " << passed_score_filter <<
				" simtime: " << F(9,3,simtime) <<
				" sasapack_score: " << F(9,3,sasapack_score) <<
				" norme_score: " << F(9,3,norme_score) <<
				" normsasa_score: " << F(9,3,normsasa_score) <<
				" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
				" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
				" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
				" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
				" packstat_score: " << F(9,3,packstat_score) <<
				' ' << unsatstring <<
				endl;

			if ( passed_score_filter ) {
				pose.dump_pdb( outfilename );
			}

			fflush( stdout );
			check_simtime();
		} // nstruct
	} // goodlines
	signal_that_job_is_done();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void
// pick_nnmake_fragments_from_single_sequence(
//                       string const & seq,
//                       FragLib & fraglib,
//                       bool const use_nohoms_flag = false
//                       )
// {
//  string const prefix( "/tmp/pbradley_"+filebase( output_tag())+"_"+string_of(uniform())),
//   small_frags_filename( prefix+"_3mers.txt" ),
//   large_frags_filename( prefix+"_9mers.txt" ),
//   logfile( prefix+"_log.txt" ),
//   errfile( prefix+"_err.txt" );

//  /// NOTE: we are using -nohoms flag:
//  ///
//  string extra_flags;
//  if ( use_nohoms_flag ) extra_flags += " -nohoms ";
//  if ( option[ my_options::ssblast ] ) extra_flags += " -ssblast ";
//  string const cmd( "python /home/pbradley/nnmake/create_fragment_files_from_sequence.py "+
//           small_frags_filename+" "+large_frags_filename +" "+seq+" "+extra_flags+
//           " > " + logfile + " 2> " + errfile );
//  run_command( cmd );

//  if ( !utility::file::file_exists( small_frags_filename ) || !utility::file::file_exists( small_frags_filename ) ) {
//   cout << "nnmake failed " << get_hostname() << endl;
//   run_command( "date" );
//   run_command( "cat "+logfile );
//   run_command( "cat "+errfile );
//  } else {
//   fraglib.library( 3 ).read_file( small_frags_filename, 3, 3 );
//   fraglib.library( 9 ).read_file( large_frags_filename, 9, 3 );
//  }

//  run_command( "rm "+prefix+"*" );
// }



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
refold_logfile_test()
{
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	bool const use_asymmetric_centroid_scoring( option[ my_options::use_asymmetric_centroid_scoring ] );
	Size const extra_repeats_for_centroid_folding( use_asymmetric_centroid_scoring ? 0 : 2 );
	Size const base_repeat( 3 );

	/// read all the lines in the file, noting which are the "good" ones
	strings all_lines;
	Sizes good_lines;
	{
		string const logfile( start_file() );
		string line;
		ifstream data( logfile.c_str() );
		runtime_assert( data.good() );
		while ( getline( data, line ) ) {
			all_lines.push_back( line );
			strings const l( split_to_vector1( line ) );
			bool passed_score_filter( false );
			for ( Size i=1; i< l.size(); ++i ) if ( l[i] == "passed_score_filter:" ) passed_score_filter = ( l[i+1] == "1" );
			if ( passed_score_filter ) good_lines.push_back( all_lines.size () );
		}
		data.close();
	}


	numeric::random::random_permutation( good_lines, numeric::random::rg() );

	for ( Size gi=1; gi<= good_lines.size(); ++gi ) {
		Size const line_index( good_lines[gi] );

		string const simfile( shared_output_tag() + "_line"+string_of( line_index)+".work" );
		string const worktag( "refold" );

		Size const first_n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );

		if ( first_n > nstruct() ) continue;

		//string const & line( all_lines[ line_index ] );
		strings const l( split_to_vector1( all_lines[ line_index ] ) );

		string const dirtag( l[1] ), filename( dirtag.substr( 0, dirtag.find_last_of( "/" ) ) +"/" + l[3] );

		Size repeatlen(0), nrepeat(0);
		string repeatseq;

		for ( Size i=1; i< l.size(); ++i ) {
			if      ( l[i] == "repeatlen:" ) repeatlen = int_of( l[i+1] );
			else if ( l[i] == "nrepeat:" ) nrepeat = int_of( l[i+1] );
			else if ( l[i] == "repeatseq:" ) repeatseq = l[i+1];
		}


		if ( !utility::file::file_exists( filename ) ) {
			cout <<"missing " << filename << endl;
			cerr <<"missing " << filename << endl;
			continue;
		}

		Pose pdb_pose;
		pose_from_pdb( pdb_pose, filename );

		runtime_assert( num_chains( pdb_pose ) == 1 );
		runtime_assert( pdb_pose.total_residue() == nrepeat * repeatlen );
		runtime_assert( pdb_pose.sequence().substr(0,repeatlen) == repeatseq );


		/// compute psipred ss using single sequence, for frag picking
		string sspred_for_frags;
		{ // run psipred to re-predict the secondary structure
			vector1< Reals > pred_HEL;
			string sequence;
			for ( Size i=1; i<= nrepeat + extra_repeats_for_centroid_folding; ++i ) sequence += repeatseq;
			runtime_assert( sequence.size() == repeatlen * ( nrepeat + extra_repeats_for_centroid_folding ) );
			run_psipred( sequence, sspred_for_frags, pred_HEL );
			runtime_assert( sspred_for_frags.size() == sequence.size() );
		}


		bool first_time_through( true );

		while ( true ) {
			Size const n( first_time_through ?
				first_n : get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( n > nstruct() ) break;
			first_time_through = false;

			Pose pose;
			bool const abrelax( true );

			refold_repeat_pose( repeatseq, sspred_for_frags, nrepeat+extra_repeats_for_centroid_folding, base_repeat, pose,
				abrelax, fa_scorefxn, extra_repeats_for_centroid_folding,
				use_asymmetric_centroid_scoring );
			runtime_assert( chain_end( 1,pose ) == pdb_pose.total_residue() ); // back to the same size
			runtime_assert( pose.total_residue() == pdb_pose.total_residue()+1 ); // the virtual residue

			Real const final_score( ( *fa_scorefxn )( pose ) );

			/// calc rmsd to pdb_pose
			Real rmsd( 0.0 ), rmsd_single_repeat( 0.0 );
			{
				using namespace core::id;
				AtomID_Map< AtomID > atom_map;
				initialize_atomid_map( atom_map, pose, id::GLOBAL_BOGUS_ATOM_ID );
				for ( Size i=1; i<= repeatlen; ++i ) {
					atom_map[ AtomID( pose.residue(i).atom_index("CA"), i ) ] =
						AtomID( pdb_pose.residue(i).atom_index("CA"),i);
				}
				rmsd_single_repeat = rmsd_by_mapping( pose, pdb_pose, atom_map );
				for ( Size i=1; i<= pdb_pose.total_residue(); ++i ) {
					atom_map[ AtomID( pose.residue(i).atom_index("CA"), i ) ] =
						AtomID( pdb_pose.residue(i).atom_index("CA"),i);
				}
				rmsd = rmsd_by_mapping( pose, pdb_pose, atom_map );
				conformation::symmetry::SymmetryInfo const symminfo( *pose::symmetry::symmetry_info( pose ) );
				pose::symmetry::make_asymmetric_pose( pose );
				superimpose_pose( pose, pdb_pose, atom_map ); // map pose into the same reference frame as pdb_pose
				pose::symmetry::make_symmetric_pose( pose, symminfo );
			}

			/// more stats
			///
			bools subset( pose.total_residue(), false );
			Size const base_repeat_offset( repeatlen*(base_repeat-1)), nres_protein( nrepeat*repeatlen );
			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;

			Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa, packstat_score;
			compute_sasa_scores_for_subset( subset, pose, sasapack_score, norme_score, normsasa_score,
				exposed_polar_sasa, exposed_nonpolar_sasa, buried_polar_sasa, buried_nonpolar_sasa,
				packstat_score );



			Real const n2cdist( pose.residue(1).xyz("N").distance( pose.residue(chain_end(1,pose)).xyz("C") ) );

			/// RE-compute transform: twist, axis, after relax
			Real rise, twist, min_radius, max_radius, com_radius, handedness;
			compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius );
			{
				Vectors coords;
				for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
				handedness = get_chirality( repeatlen, coords );
			}

			string const repeatbb( torsion2big_bin_string( base_repeat_offset+1, base_repeat_offset+repeatlen, pose ) );
			runtime_assert( repeatbb.size() == repeatlen );

			string const outfilename( output_tag() + "refold_L"+ string_of( line_index )+
				"_" + string_of( nrepeat ) +"_"+ string_of( repeatlen ) +
				"_N"+ lead_zero_string_of( n, 4 )+".pdb" );

			bool const passed_score_filter
				( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );

			cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
				" nrepeat: " << nrepeat <<
				" repeatlen: " << repeatlen <<
				" repeatseq: " << repeatseq <<
				" repeatss: " << pose.secstruct().substr(0,repeatlen) <<
				" repeatbb: " << repeatbb <<
				" rmsd: " << F(9,3,rmsd) <<
				" rmsd_single_repeat: " << F(9,3,rmsd_single_repeat) <<
				" start_file: " << filename <<
				" handedness: " << F(9,3,handedness) <<
				" max_radius: " << F(9,3,max_radius )<<
				" com_radius: " << F(9,3,com_radius )<<
				" min_radius: " << F(9,3,min_radius ) <<
				" rise: " << F(9,3,rise ) <<
				" twist: " << F(9,3,twist) <<
				" n2cdist: " << F(9,3,n2cdist ) <<
				" passed_score_filter: " << passed_score_filter <<
				" sasapack_score: " << F(9,3,sasapack_score) <<
				" norme_score: " << F(9,3,norme_score) <<
				" normsasa_score: " << F(9,3,normsasa_score) <<
				" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
				" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
				" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
				" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
				" packstat_score: " << F(9,3,packstat_score) <<
				endl;

			fflush( stdout );

			if ( passed_score_filter ) {
				pose.dump_pdb( outfilename );
			}
		}

	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
unbound_frag_extend_test()
{
	strings const files( start_files() );

	Size const nrepeats_to_append( option[ my_options::nrepeat ] );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		if ( !pose.residue( pose.total_residue() ).is_protein() ) {
			pose.conformation().delete_residue_slow( pose.total_residue());
		}

		remove_upper_terminus_type_from_pose_residue( pose, pose.total_residue() );
		remove_lower_terminus_type_from_pose_residue( pose, 1 );

		Size const repeatlen( deduce_repeatlen( pose.sequence() ) ), nrepeat( pose.total_residue()/repeatlen );
		runtime_assert( pose.total_residue() == repeatlen * nrepeat );

		// what is the transform that maps between repeats?
		Stub const istub( pose.residue(1).xyz("CA"), pose.residue(2).xyz("CA"), pose.residue(3).xyz("CA") ),
			jstub( pose.residue(repeatlen+1).xyz("CA"),
			pose.residue(repeatlen+2).xyz("CA"),
			pose.residue(repeatlen+3).xyz("CA") );

		// what's the transform (R*x + v) that transforms istub onto jstub
		// R*istub.M = jstub.M
		Matrix const R( jstub.M * istub.M.transposed() );
		// R*istub.v + v = jstub.v
		Vector const v( jstub.v - R * istub.v );

		Pose movpose( pose );

		for ( Size n=1; n<= nrepeats_to_append; ++n ) {

			movpose.apply_transform_Rx_plus_v( R, v );

			for ( Size i=1; i<=repeatlen; ++i ) {
				pose.append_residue_by_bond( movpose.residue( (nrepeat-1)*repeatlen+i) );
			}
		}

		pose.dump_pdb(output_tag() + "ext_"+filebase( files[fi] )+".pdb" );

	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
regan_tpr_test()
{
	Size const repeatlen( 34 ), nrepeat_fix( 3 );

	string const filename( "1na0A_3rep.pdb" );

	Pose fixpose;
	pose_from_pdb( fixpose, filename );

	set_ss_from_dssp( filename, fixpose ); // phil_io.hh

	Pose movpose( fixpose );

	{ /// look at repeat params
		string const bb( torsion2big_bin_string( 1, fixpose.total_residue(), fixpose ) );
		for ( Size i=5; i<= fixpose.total_residue()-repeatlen-2; ++i ) {
			Residue const & rsd1( fixpose.residue( i ) ), &rsd2( fixpose.residue(i+repeatlen) );
			Stub const
				stub1( rsd1.xyz("N"), rsd1.xyz("CA"), rsd1.xyz("C") ),
				stub2( rsd2.xyz("N"), rsd2.xyz("CA"), rsd2.xyz("C") );

			Real twist0, rise0;
			{
				Real theta;
				Vector n,t,cen;
				get_stub_transform_data( stub1, stub2, cen, n, t, theta );
				twist0 = numeric::conversions::degrees( theta );
				rise0 = t.dot( n );
			}
			cout << "twist_rise: " << F(9,3,twist0) << F(9,3,rise0) << I(5,i) << ' ' << bb[i-1] << ' ' <<
				fixpose.secstruct(i) << endl;
		}


		Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);
		Size const helix1_len( 14 ), turn1_len( 2 ), helix2_len( 14 ), turn2_len( 4 ), base_repeat(2);

		runtime_assert( helix1_len + turn1_len + helix2_len + turn2_len == repeatlen );

		compute_helix_axis_angles( fixpose, base_repeat, repeatlen, helix1_len, turn1_len, helix2_len,
			helix1_dist, helix1_twist, helix2_dist, helix2_twist );
		cout << "helix_angles: " <<
			" helix1_dist: " << F(9,3,helix1_dist) <<
			" helix1_twist: " << F(9,3,helix1_twist) <<
			" helix2_dist: " << F(9,3,helix2_dist) <<
			" helix2_twist: " << F(9,3,helix2_twist) << endl;
		exit(0);
	}



	for ( Size n_extra=1; n_extra<= 7; ++n_extra ) {
		Size nrepeat_mov( nrepeat_fix+n_extra-1 );
		runtime_assert( fixpose.total_residue() == repeatlen * nrepeat_fix );
		runtime_assert( movpose.total_residue() == repeatlen * nrepeat_mov );

		// superimpose the last two repeats of movpose onto the first 2 repeats of fixpose
		// then append the last repeat from fixpose onto movpose
		Size const fixbegin( 1 ), movbegin( movpose.total_residue() - 2*repeatlen+1 );

		using namespace core::id;
		AtomID_Map< AtomID > atom_map;
		initialize_atomid_map( atom_map, movpose, id::GLOBAL_BOGUS_ATOM_ID );
		for ( Size i=0; i< 2*repeatlen; ++i ) {
			atom_map[ AtomID( movpose.residue( movbegin+i ).atom_index("CA"), movbegin+i ) ] =
				AtomID( fixpose.residue( fixbegin+i ).atom_index("CA"), fixbegin+i );
		}
		superimpose_pose( movpose, fixpose, atom_map );

		remove_upper_terminus_type_from_pose_residue( movpose, movpose.total_residue() );
		for ( Size i=1; i<= repeatlen; ++i ) {
			movpose.append_residue_by_bond( fixpose.residue( fixpose.total_residue() - repeatlen + i ) );
		}

		cout << "chainbreak? " << n_extra << F(9,3,movpose.residue( movpose.total_residue()-34 ).xyz("C").distance
			( movpose.residue( movpose.total_residue()-33 ).xyz("N") )) << endl;

		++nrepeat_mov;
		movpose.dump_pdb("1na0A_extended_nr"+string_of(nrepeat_mov)+".pdb");

	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
w3b_test()
{
	Size const repeatlen( 34 ), nrepeat( 3 ); // how many repeats to write out

	Pose fullpose;
	pose_from_pdb( fullpose, "1w3b.pdb" );

	// where does the 1st repeat start?
	Size const r1begin( fullpose.pdb_info()->pdb2pose( 'A', 45 ) );
	Size const rNend( fullpose.pdb_info()->pdb2pose( 'A', 384 ) );
	runtime_assert( fullpose.residue( r1begin ).name1() == 'T' );

	fullpose.conformation().delete_residue_range_slow( rNend+1, fullpose.total_residue() );
	fullpose.conformation().delete_residue_range_slow( 1, r1begin-1 );

	runtime_assert( fullpose.total_residue()%repeatlen == 0 );

	Size const nrepeat_full( fullpose.total_residue() / repeatlen );

	for ( Size r=2; r< nrepeat_full; ++r ) {
		bools subset( fullpose.total_residue(), false );
		Size const rbegin( (r-1)*repeatlen+1 ), rend( r*repeatlen );
		for ( Size i=rbegin; i<= rend; ++i ) subset[i] = true;

		Real sasapack_score, norme_score, normsasa_score,
			exposed_polar_sasa, exposed_nonpolar_sasa,
			buried_polar_sasa, buried_nonpolar_sasa,
			buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
			packstat_score;

		compute_sasa_scores_for_subset( subset, fullpose,
			sasapack_score, norme_score, normsasa_score,
			exposed_polar_sasa, exposed_nonpolar_sasa,
			buried_polar_sasa, buried_nonpolar_sasa,
			buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
			packstat_score );

		cout << "repeat_sasa_scores: " << I(2,r) <<
			" sasapack_score: " << F(9,3,sasapack_score) <<
			" norme_score: " << F(9,3,norme_score) <<
			" normsasa_score: " << F(9,3,normsasa_score) <<
			" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa) <<
			" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa) <<
			" buried_polar_sasa: " << F(9,3,buried_polar_sasa) <<
			" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa) <<
			" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc) <<
			" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc) <<
			" packstat_score: " << F(9,3,packstat_score) << endl;
	}



	Size const numwindows( nrepeat_full - nrepeat + 1 );
	for ( Size i=0; i< numwindows; ++i ) {
		Size const begin( i*repeatlen+1 ), end( (i+nrepeat)*repeatlen );
		Pose pose( fullpose );
		if ( end < pose.total_residue() ) pose.conformation().delete_residue_range_slow( end+1, pose.total_residue() );
		if ( begin > 1 ) pose.conformation().delete_residue_range_slow( 1, begin-1 );
		runtime_assert( pose.total_residue() == nrepeat * repeatlen );
		pose.dump_pdb("1w3b_3rep"+lead_zero_string_of( i,2)+".pdb" );
	}


}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
run_phaser_and_phenix_refinement(
	string const modelfile,
	Real const llg_threshold_for_refinement,
	Real & llg,
	Real & tfz,
	Real & pak,
	Real & rphaser,
	Real & rwork,
	Real & rfree,
	Real & coord_error,
	Real & xray_target,
	string & space_group,
	string & refined_pdbfile
)
{

	/// some possible useful commands for refining at low-resolution, from frank
	// > optimize_wxc=True
	// > ramachandran_restraints=True
	// > main.secondary_structure_restraints=True


	// if ( simulated_annealing ) {
	//  refinement_flags += " simulated_annealing=True strategy=rigid_body+individual_sites+individual_adp main.number_of_macro_cycles=10 ";
	// }

	run_command( "python /home/pbradley/python/phenix/run_toroid_phaser_and_phenix_refinement.py "+modelfile+" "+
		string_of( llg_threshold_for_refinement ) );

	string const scoresfile( modelfile +"_MR.scores" );
	if ( utility::file::file_exists( scoresfile ) ) {
		ifstream data( scoresfile.c_str() );
		string line;
		getline( data, line );
		istringstream l(line );
		l >> llg >> tfz >> pak >> rphaser >> rwork >> rfree >> space_group >> coord_error >> xray_target;
		if ( l.fail() ) {
			rphaser = rwork = rfree = 1.0;
			llg = tfz = pak = 0.0; coord_error = xray_target = 1000.0;
		}
		data.close();
		refined_pdbfile = modelfile +"_MR_refine_001.pdb";
		//new_solfile = modelfile +"_MR.sol";
	} else {
		rphaser = rwork = rfree = 1.0;
		llg = tfz = pak = 0.0;
		space_group.clear();
		refined_pdbfile.clear();
		//new_solfile.clear();
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void
run_phaser_and_phenix_refinement_9x2(
	string const modelfile,
	Real const llg_threshold_for_refinement,
	Real & llg,
	Real & tfz,
	Real & rwork,
	Real & rfree,
	string & space_group,
	string & refined_pdbfile
)
{

	run_command( "python /home/pbradley/python/phenix/run_9x2_phaser_and_phenix_refinement.py "+modelfile+" "+
		string_of( llg_threshold_for_refinement ) );

	string const scoresfile( modelfile +"_MR.scores" );
	if ( utility::file::file_exists( scoresfile ) ) {
		ifstream data( scoresfile.c_str() );
		string line, tmp;
		getline( data, line );
		istringstream l(line );
		l >> tmp >> llg >> tmp >> tfz >> tmp >> rwork >> tmp >> rfree >> tmp >> space_group;
		if ( l.fail() ) {
			rwork = rfree = 1.0;
			llg = tfz = 0.0;
		}
		data.close();
		refined_pdbfile = modelfile +"_MR_refine_001.pdb";
	} else {
		rwork = rfree = 1.0;
		llg = tfz  = 0.0;
		space_group.clear();
		refined_pdbfile.clear();
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void
run_phaser_and_phenix_refinement_3x2(
	string const modelfile,
	Real const llg_threshold_for_refinement,
	Real & llg,
	Real & tfz,
	Real & rwork,
	Real & rfree,
	string & space_group,
	string & refined_pdbfile
)
{
	run_command( "python /home/pbradley/python/phenix/run_3x2_phaser_and_phenix_refinement.py "+modelfile+" "+
		string_of( llg_threshold_for_refinement ) );

	string const scoresfile( modelfile +"_MR.scores" );
	if ( utility::file::file_exists( scoresfile ) ) {
		ifstream data( scoresfile.c_str() );
		string line, tmp;
		getline( data, line );
		istringstream l(line );
		l >> tmp >> llg >> tmp >> tfz >> tmp >> rwork >> tmp >> rfree >> tmp >> space_group;
		if ( l.fail() ) {
			rwork = rfree = 1.0;
			llg = tfz = 0.0;
		}
		data.close();
		refined_pdbfile = modelfile +"_MR_refine_001.pdb";
	} else {
		rwork = rfree = 1.0;
		llg = tfz  = 0.0;
		space_group.clear();
		refined_pdbfile.clear();
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void
run_phaser_and_phenix_refinement_tordim1(
	string const modelfile,
	Real const llg_threshold_for_refinement,
	Real & llg,
	Real & tfz,
	Real & rwork,
	Real & rfree,
	string & space_group,
	string & refined_pdbfile
)
{
	run_command( "python /home/pbradley/python/phenix/run_tordim1_phaser_and_phenix_refinement.py "+modelfile+" "+
		string_of( llg_threshold_for_refinement ) );

	string const scoresfile( modelfile +"_MR.scores" );
	if ( utility::file::file_exists( scoresfile ) ) {
		ifstream data( scoresfile.c_str() );
		string line, tmp;
		getline( data, line );
		istringstream l(line );
		l >> tmp >> llg >> tmp >> tfz >> tmp >> rwork >> tmp >> rfree >> tmp >> space_group;
		if ( l.fail() ) {
			rwork = rfree = 1.0;
			llg = tfz = 0.0;
		}
		data.close();
		refined_pdbfile = modelfile +"_MR_refine_001.pdb";
	} else {
		rwork = rfree = 1.0;
		llg = tfz  = 0.0;
		space_group.clear();
		refined_pdbfile.clear();
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void
run_phaser_and_phenix_refinement_tri(
	string const target_space_group, // maybe ''
	string const modelfile,
	Real const llg_threshold_for_refinement,
	Size & actual_num_models,
	Size & num_models, // number searched for
	Real & llg,
	Real & tfz,
	Real & rwork,
	Real & rfree,
	string & space_group,
	string & refined_pdbfile,
	Size const target_num_models = 0,
	Size const composition = 0,
	string const fixed_file = ""
)
{
	string cmd( "python /home/pbradley/python/phenix/run_tri_phaser_and_phenix_refinement.py "+modelfile+" "+
		string_of( llg_threshold_for_refinement ) );
	if ( !target_space_group.empty() ) cmd += " --space_group "+target_space_group;
	if ( target_num_models ) cmd += " --num_models "+string_of(target_num_models);
	if ( composition ) cmd += " --composition "+string_of(composition);
	if ( dry_run() ) cmd += " --dry_run ";
	if ( option[ my_options::skip_tncs ] ) cmd += " --skip_tncs ";
	if ( option[ my_options::sgalt_all ] ) cmd += " --sgalt_all ";
	if ( option[ my_options::sgalt_none ] ) cmd += " --sgalt_none ";
	if ( option[ my_options::full_search ] ) cmd += " --full_search ";
	if ( option[ my_options::xray_data_tag ].user() ) cmd += " --xray_data_tag "+option[ my_options::xray_data_tag ]();
	if ( !fixed_file.empty() ) {
		cmd += " --fixed_file "+fixed_file;
	}
	if ( option[ my_options::autobuild_rwork_threshold ].user() ) {
		cmd += " --autobuild --autobuild_rwork_threshold "+string_of( option[ my_options::autobuild_rwork_threshold ] );
	}
	run_command( cmd );

	string const scoresfile( modelfile +"_MR.scores" );
	if ( utility::file::file_exists( scoresfile ) ) {
		ifstream data( scoresfile.c_str() );
		string line, tmp;
		getline( data, line );
		istringstream l(line );
		l >> tmp >> llg >> tmp >> tfz >> tmp >> rwork >> tmp >> rfree >> tmp >> space_group >> tmp >> actual_num_models >> tmp >> num_models;
		if ( l.fail() ) {
			rwork = rfree = 1.0;
			llg = tfz = 0.0;
			space_group = "PARSE_ERROR";
		}
		data.close();
		refined_pdbfile = modelfile +"_MR_refine_001.pdb";
	} else {
		rwork = rfree = 1.0;
		llg = tfz  = 0.0;
		space_group = "NO_SCORESFILE";
		refined_pdbfile  = "NONE";
	}
}


/// some possible useful commands for refining at low-resolution, from frank
// > optimize_wxc=True
// > ramachandran_restraints=True
// > main.secondary_structure_restraints=True


// if ( simulated_annealing ) {
//  refinement_flags += " simulated_annealing=True strategy=rigid_body+individual_sites+individual_adp main.number_of_macro_cycles=10 ";
// }



void
abrelax_sequence(
	string const & fullseq,
	ScoreFunction const & fa_scorefxn_in,
	Pose & pose
)
{
	ScoreFunctionOP fa_scorefxn( fa_scorefxn_in.clone() );

	// create a starting pose from the sequence
	//Pose pose;
	pose.clear();
	for ( Size i=0; i< fullseq.size(); ++i ) {
		pose.append_residue_by_bond( *get_vanilla_protein_residue( fullseq[i] ), true );
	}
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

	core::fragment::FragSetOP small_frags(0), large_frags(0);
	pick_nnmake_fragments_from_single_sequence( pose.sequence(), small_frags, large_frags );

	devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		pose.set_phi  ( i, init_phi   );
		pose.set_psi  ( i, init_psi   );
		pose.set_omega( i, init_omega );
		pose.set_secstruct( i, 'L' );
	}
	for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

	{
		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		mm->set_bb( true );
		protocols::abinitio::ClassicAbinitio abinitio( small_frags, large_frags, mm );
		abinitio.init( pose );
		if ( !dry_run() ) abinitio.apply( pose );
	}


	/// now fastrelax
	devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
	{
		protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

		MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
		movemap->set_bb (true);
		movemap->set_chi(true);
		fastrelax.set_movemap( movemap );
		if ( !dry_run() ) fastrelax.apply( pose );
	}
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
strand_test()
{
	/// read library of pdb files
	strings const files( start_files() );

	Size const strand_window( 4 );

	/// store downstream helix-atom coords in defined reference frames
	//strings const turns( option[ my_options::turns ]() );

	HelicalParamsFitMultifunc func;




	for ( Size fi=1; fi<= files.size(); ++fi ) {

		string const filename( files[fi] );

		Pose pose;
		pose_from_pdb( pose, filename );

		set_ss_from_dssp( filename, pose );

		add_termini_at_protein_chainbreaks( pose );

		string const allbb( torsion2big_bin_string( 1, pose.total_residue(), pose, true ) );
		runtime_assert( allbb.size() == pose.total_residue() );

		bools is_beta( pose.total_residue(), false );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			is_beta[i] = ( pose.residue(i).is_protein() && pose.secstruct(i) == 'E' && allbb[i-1] == 'B' );
		}

		/// first identify bb NH--OOC hbonds
		Sizes
			bbn_hbond_partner ( pose.total_residue(), 0 ),
			bbn_hbond_partner2( pose.total_residue(), 0 ),
			bbo_hbond_partner ( pose.total_residue(), 0 ),
			bbo_hbond_partner2( pose.total_residue(), 0 );

		Real const bb_hbond_maxdis2( 3.4 * 3.4 );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const irsd( pose.residue(i));
			if ( !is_beta[i] ) continue;
			if ( irsd.aa() == aa_pro ) continue; // skip prolines here !!!!!!!!!!!!!!!!!!!!
			for ( Size j=1; j<= pose.total_residue(); ++j ) {
				Residue const jrsd( pose.residue(j));
				if ( !is_beta[j] ) continue;
				if ( j>=i-2 && i>= j-2 ) continue;

				if ( irsd.xyz("N").distance_squared( jrsd.xyz("O") ) <= bb_hbond_maxdis2 ) {
					if ( !bbn_hbond_partner[i] ) bbn_hbond_partner[i] = j;
					else bbn_hbond_partner2[i] = j; // bifurcation
					if ( !bbo_hbond_partner[j] ) bbo_hbond_partner[j] = i;
					else bbo_hbond_partner2[j] = i; // bifurcation
				}
			}
		}


		Sizes beta_partner( pose.total_residue(), 0 );

		/// look for i NH-->OC(j-1)
		/// look for i OC<--NH(j+1)
		bools is_beta_paired( pose.total_residue(), false );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const irsd( pose.residue(i));
			if ( !is_beta[i] ) continue;
			if ( bbn_hbond_partner[i] && bbo_hbond_partner[i] && bbn_hbond_partner[i] + 2== bbo_hbond_partner[i] &&
					is_beta[ bbn_hbond_partner[i]+1 ] ) {
				Size const j( bbn_hbond_partner[i] + 1 ); // rsd that is paired with i
				runtime_assert( is_beta[j-1] );
				runtime_assert( is_beta[j+1] );
				// parallel beta bridge
				is_beta_paired[i] = is_beta_paired[ j-1] = is_beta_paired[j] = is_beta_paired[j+1] = true;

				if ( i>1                    && is_beta[i-1] ) is_beta_paired[i-1] = true;
				if ( i<pose.total_residue() && is_beta[i+1] ) is_beta_paired[i+1] = true;
			}
		}

		// now look for continuous strand segments
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( !is_beta_paired[i] ) continue;
			Size j(i);
			while ( j<pose.total_residue() && is_beta_paired[j+1] ) ++j;
			// found a strand
			Size const strandlen( j-i + 1);
			if ( strandlen >= strand_window ) {
				// these are guesses
				Reals params( 3 );
				params[1] = 3.3;
				params[2] = radians( 180.0 );
				params[3] = 0.9; // ca-distance

				static bool init( false );
				if ( !init ) {
					init = true;
					HelixParams hparams;
					hparams.rise = params[1];
					hparams.twist = params[2];
					hparams.tilt = 0;
					hparams.tilt_direction = 0;
					hparams.ca_distance = params[3]; // guess

					Vectors coords;
					Stub start_stub;
					generate_helix_coords( start_stub, hparams, 10, coords );

					write_ca_pdbfile( coords, "strand10.pdb" );
				}

				Vectors coords;
				for ( Size k=i; k<= j; ++k ) coords.push_back( pose.residue(k).xyz("CA") );

				func.set_target_coords( 1, strandlen, coords );
				func.optimize_tilt( false );

				Size iterations;
				Real const tolerance( option[ basic::options::OptionKeys::run::min_tolerance ] );
				Real final_func_value;
				Real start_func_value( func( params ) );
				optimization::powell( params, func, tolerance, iterations, final_func_value );

				cout << "final_params_notilt " << I(3,strandlen) << ' ' << I(4,i) << ' ' << I(4,j) <<
					F(9,3,start_func_value) << F(9,3,final_func_value) <<
					F(9,3,params[1]) << F(9,3,degrees(params[2])) <<F(9,3,params[3]) << ' ' << files[fi] << endl;

				{// try optimizing twist
					func.optimize_tilt( true );
					Reals params2( 5, 0.0 );
					params2[1] = params[1];
					params2[2] = params[2];
					params2[3] = 0.0;
					params2[4] = 0.0;
					params2[5] = params[3];

					start_func_value = func( params2 );
					optimization::powell( params2, func, tolerance, iterations, final_func_value );

					cout << "final_params_tilt " << I(3,strandlen) << ' ' << I(4,i) << ' ' << I(4,j) <<
						F(9,3,start_func_value) << F(9,3,final_func_value) <<
						F(9,3,params2[1]) <<
						F(9,3,degrees(params2[2])) <<
						F(9,3,degrees(params2[3])) <<
						F(9,3,degrees(params2[4])) <<
						F(9,3,params2[5]) << ' ' << files[fi] << endl;
				}

			}
			i = j+1;
		}
	} // files
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
strand_turn_test()
{
	/// read library of pdb files
	strings const files( start_files() );

	Size const strand_window( 4 ), helix_window( 7 );

	bool const restrict_to_bab( option[ my_options::restrict_to_bab ] );
	Size const max_bab_seqsep( 40 );
	Size const max_turnlen( 7 );

	for ( Size fi=1; fi<= files.size(); ++fi ) {

		string const filename( files[fi] );

		Pose pose;
		pose_from_pdb( pose, filename );

		set_ss_from_dssp( filename, pose );

		add_termini_at_protein_chainbreaks( pose );

		string allbb( string("X")+torsion2big_bin_string( 1, pose.total_residue(), pose, true ) );

		runtime_assert( allbb.size() == pose.total_residue()+1 ); // so we can 1-index

		bools is_beta( pose.total_residue(), false ), is_alpha( pose.total_residue(), false );
		Sizes strand_num( pose.total_residue(), 0 );
		{
			Size strand(0);
			bool in_strand( false );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				is_beta [i] = ( pose.residue(i).is_protein() && pose.secstruct(i) == 'E' && allbb[i] == 'B' );
				is_alpha[i] = ( pose.residue(i).is_protein() && pose.secstruct(i) == 'H' && allbb[i] == 'A' );
				if ( pose.secstruct(i) == 'E' ) {
					if ( !in_strand ) {
						in_strand = true;
						++strand;
					}
					strand_num[i] = strand;
					TR.Trace << "strand_num: " << i << ' ' << strand << endl;
				} else in_strand = false;
			}
		}

		bools is_bab_helix( pose.total_residue(), false );
		{ // try identifying beta-alpha-beta units

			/// first identify bb NH--OOC hbonds
			Sizes
				bbn_hbond_partner ( pose.total_residue(), 0 ),
				bbn_hbond_partner2( pose.total_residue(), 0 ),
				bbo_hbond_partner ( pose.total_residue(), 0 ),
				bbo_hbond_partner2( pose.total_residue(), 0 );

			Real const bb_hbond_maxdis2( 3.4 * 3.4 );

			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const irsd( pose.residue(i));
				if ( !is_beta[i] ) continue;
				if ( irsd.aa() == aa_pro ) continue; // skip prolines here !!!!!!!!!!!!!!!!!!!!
				for ( Size j=1; j<= pose.total_residue(); ++j ) {
					Residue const jrsd( pose.residue(j));
					if ( !is_beta[j] ) continue;
					if ( j>=i-2 && i>= j-2 ) continue;

					if ( irsd.xyz("N").distance_squared( jrsd.xyz("O") ) <= bb_hbond_maxdis2 ) {
						if ( !bbn_hbond_partner[i] ) bbn_hbond_partner[i] = j;
						else bbn_hbond_partner2[i] = j; // bifurcation
						if ( !bbo_hbond_partner[j] ) bbo_hbond_partner[j] = i;
						else bbo_hbond_partner2[j] = i; // bifurcation
					}
				}
			}


			//Sizes beta_partner( pose.total_residue(), 0 );

			/// look for i NH-->OC(j-1)
			/// look for i OC<--NH(j+1)
			//bools is_beta_paired( pose.total_residue(), false );
			Vectors coords;
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( pose.residue(i).is_protein() ) coords.push_back( pose.residue(i).xyz("CA"));
				else coords.push_back( Vector(0,0,0));
			}
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const irsd( pose.residue(i));
				if ( !is_beta[i] ) continue;
				if ( bbn_hbond_partner[i] && bbo_hbond_partner[i] && bbn_hbond_partner[i] + 2== bbo_hbond_partner[i] &&
						is_beta[ bbn_hbond_partner[i]+1 ] ) {
					Size const j( bbn_hbond_partner[i] + 1 ); // rsd that is paired with i
					runtime_assert( is_beta[j-1] );
					runtime_assert( is_beta[j+1] );

					// found a parallel beta bridge
					// is this  a bab unit?
					Size const pos1( min(i,j) ), pos2( max(i,j) );
					TR.Trace << "beta_bridge: " << pos1 << ' ' << pos2 << ' '<<
						strand_num[pos1] << ' ' << strand_num[pos2] << ' ' <<
						pose.pdb_info()->chain(pos1) << ' ' << pose.pdb_info()->number(pos1) << ' ' <<
						filename << endl;

					if ( is_beta[i-1] && is_beta[i+1] ) { // show the geometry
						// what is the chirality
						Reals distances;
						Real orientation;
						get_beta_pairing_geometry( coords, pos1-1, pos2-1, orientation, distances );
						cout << "parallel_beta_pairing " << F(9,3,orientation);
						for ( Size ii=1; ii<= distances.size(); ++ii ) cout << F(9,3,distances[ii]);
						cout << endl;
						get_beta_pairing_geometry( coords, pos2-1, pos1-1, orientation, distances );
						cout << "parallel_beta_pairing " << F(9,3,orientation);
						for ( Size ii=1; ii<= distances.size(); ++ii ) cout << F(9,3,distances[ii]);
						cout << endl;
						{
							Residue const & jrsd( pose.residue(j) );
							Stub const
								istub( irsd.xyz("N"), irsd.xyz("CA"), irsd.xyz("C") ),
								jstub( jrsd.xyz("N"), jrsd.xyz("CA"), jrsd.xyz("C") );
							cout << "parallel_beta_transform " << RT( istub, jstub ) << endl;
						}
					}

					if ( strand_num[ pos2 ] == strand_num[ pos1 ]+1 && pose.chain(pos1) == pose.chain(pos2) &&
							pos2 - pos1 <= max_bab_seqsep ) {
						// sequential strands, parallel paired
						for ( Size k=pos1; k<= pos2; ++k ) {
							if ( is_alpha[k] ) {
								if ( !is_bab_helix[k] ) {
									TR.Trace << "is_bab: " << pos1 << ' ' << pos2 << ' '<<
										strand_num[pos1] << ' ' << strand_num[pos2] << ' ' << k <<
										' ' << pose.pdb_info()->chain(k) << ' ' << pose.pdb_info()->number(k) << ' ' <<
										filename << endl;
								}
								is_bab_helix[k] = true;
							}
						}
					}
				}
			}
		} // scope: finding bab units


		// find beta-alpha motifs
		for ( Size i=2; i< pose.total_residue(); ++i ) {
			if ( is_beta[i-1] && !is_beta[i] ) {
				Size strandbegin(i-1);
				Size const cb( chain_begin( pose.chain(i), pose ) );
				while ( strandbegin>cb && is_beta[strandbegin-1] ) --strandbegin;
				// end of a strand
				// advance i until we're not B backbone anymore
				Size const ce( chain_end( pose.chain(i), pose ) );
				while ( i<ce && allbb[i] == 'B' ) ++i;
				if ( i==ce ) { ++i; continue; }
				runtime_assert( allbb[i] != 'B' );
				runtime_assert( allbb[i-1] == 'B' );
				// now go forward until we hit an alpha helix
				Size j(i);
				bool ok( true );
				while ( j<ce && !is_alpha[j] ) {
					if ( is_beta[j] ) { ok = false; break; }
					++j;
				}
				if ( !ok || j==ce ) { ++i; continue; }
				runtime_assert( is_alpha[j] );
				Size helixend( j );
				if ( restrict_to_bab && !is_bab_helix[j] ) { ++i; continue; } // check for bab unit
				while ( helixend<ce && is_alpha[ helixend+1] ) ++helixend;
				// got to an alpha helix
				// back up until start of A segment
				while ( j>i && allbb[j-1]=='A' ) --j;
				runtime_assert( allbb[j-1] != 'A' );
				runtime_assert( allbb[j  ] == 'A' );
				// show the turn
				Size const turnbegin( i ), turnend( j-1 ), turnlen( j-i ), strandlen( i-strandbegin ),
					helixlen( helixend-turnend );
				TR.Trace << "strandlen: " << strandlen << " helixlen: " << helixlen <<
					I(4,turnbegin) << ' ' << I(4,turnend) << ' ' <<
					pose.pdb_info()->chain(turnbegin) << ' ' << pose.pdb_info()->number( turnbegin ) << ' ' <<
					filename << endl;

				if ( strandlen >= strand_window && helixlen >= helix_window && turnlen <= max_turnlen ) {
					// write out the turn
					string const turnbb( turnlen>0 ? allbb.substr(turnbegin,turnlen) : string("-") );
					string const turnseq( turnlen>0 ? pose.sequence().substr(turnbegin-1,turnlen) : string("-") );
					string const beforeseq( pose.sequence().substr(turnbegin-1-strand_window,strand_window) );
					string const afterseq( pose.sequence().substr(turnend,helix_window) );
					Size const coordslen( turnlen + strand_window + helix_window );
					cout << "turn_coords: beta_alpha " << turnlen << ' ' << turnbb << ' ' <<
						beforeseq << ' ' << turnseq << ' ' << afterseq << ' ' << coordslen;
					// show ca-coords for flanking sses and turn
					for ( Size i=turnbegin-strand_window; i<= turnend + helix_window; ++i ) {
						cout <<
							F(9,3,pose.residue(i).xyz("CA").x()) <<
							F(9,3,pose.residue(i).xyz("CA").y()) <<
							F(9,3,pose.residue(i).xyz("CA").z());
					}
					// show the torsion angles:
					cout << " bb_torsions: ";
					for ( Size i=turnbegin; i<= turnend; ++i ) {
						cout << F(9,3,pose.phi(i)) << F(9,3,pose.psi(i)) << F(9,3,pose.omega(i) );
					}
					cout << I(4,turnbegin) << ' ' << I(4,turnend) << ' ' <<
						pose.pdb_info()->chain(turnbegin) << ' ' << pose.pdb_info()->number( turnbegin ) << ' ' <<
						filename << endl;
					{ /// write out the turn transform
						Residue const & irsd( pose.residue( turnbegin-1 ) ), &jrsd( pose.residue( turnend+1 ) ); // no reln to i,j
						Stub const
							istub( irsd.xyz("N"), irsd.xyz("CA"), irsd.xyz("C") ),
							jstub( jrsd.xyz("N"), jrsd.xyz("CA"), jrsd.xyz("C") );
						cout << "turn_transform: beta_alpha " << turnlen << ' ' << turnbb << ' ' << turnseq << ' ' <<
							RT( istub, jstub ) << ' ' << pose.pdb_info()->chain(turnbegin) << ' ' <<
							pose.pdb_info()->number( turnbegin ) << ' ' << filename << endl;
					} // scope
				}
			}
		}
		// find alpha-beta motifs
		for ( Size i=2; i< pose.total_residue(); ++i ) {
			if ( is_alpha[i-1] && !is_alpha[i] ) {
				if ( restrict_to_bab && !is_bab_helix[i-1] ) { ++i; continue; } // check for bab unit
				Size helixbegin(i-1);
				Size const cb( chain_begin( pose.chain(i), pose ) );
				while ( helixbegin>cb && is_alpha[helixbegin-1] ) --helixbegin;
				// end of a helix
				// advance i until we're not A backbone anymore
				Size const ce( chain_end( pose.chain(i), pose ) );
				while ( i<ce && allbb[i] == 'A' ) ++i;
				if ( i==ce ) { ++i; continue; }
				runtime_assert( allbb[i] != 'A' );
				runtime_assert( allbb[i-1] == 'A' );
				// now go forward until we hit a beta strand
				Size j(i);
				bool ok( true );
				while ( j<ce && !is_beta[j] ) {
					if ( is_alpha[j] ) { ok = false; break; }
					++j;
				}
				if ( !ok || j==ce ) { ++i; continue; }
				runtime_assert( is_beta[j] );
				Size strandend( j );
				while ( strandend<ce && is_beta[ strandend+1] ) ++strandend;
				// got to a beta strand
				// back up until start of B segment
				while ( j>i && allbb[j-1]=='B' ) --j;
				runtime_assert( allbb[j-1] != 'B' );
				runtime_assert( allbb[j  ] == 'B' );
				// show the turn
				Size const turnbegin( i ), turnend( j-1 ), turnlen( j-i ), helixlen( i-helixbegin ),
					strandlen( strandend-turnend );
				TR.Trace << "strandlen: " << strandlen << " helixlen: " << helixlen <<
					I(4,turnbegin) << ' ' << I(4,turnend) << ' ' <<
					pose.pdb_info()->chain(turnbegin) << ' ' << pose.pdb_info()->number( turnbegin ) << ' ' <<
					filename << endl;

				if ( strandlen >= strand_window && helixlen >= helix_window && turnlen <= max_turnlen ) {
					// write out the turn
					string const turnbb( turnlen>0 ? allbb.substr(turnbegin,turnlen) : string("-") );
					string const turnseq( turnlen>0 ? pose.sequence().substr(turnbegin-1,turnlen) : string("-") );
					string const beforeseq( pose.sequence().substr(turnbegin-1-helix_window,helix_window) );
					string const afterseq( pose.sequence().substr(turnend,strand_window) );
					Size const coordslen( turnlen + strand_window + helix_window );
					cout << "turn_coords: alpha_beta " << turnlen << ' ' << turnbb << ' ' <<
						beforeseq << ' ' << turnseq << ' ' << afterseq << ' ' << coordslen;
					// show ca-coords for flanking sses and turn
					for ( Size i=turnbegin-helix_window; i<= turnend + strand_window; ++i ) {
						cout <<
							F(9,3,pose.residue(i).xyz("CA").x()) <<
							F(9,3,pose.residue(i).xyz("CA").y()) <<
							F(9,3,pose.residue(i).xyz("CA").z());
					}
					// show the torsion angles:
					cout << " bb_torsions: ";
					for ( Size i=turnbegin; i<= turnend; ++i ) {
						cout << F(9,3,pose.phi(i)) << F(9,3,pose.psi(i)) << F(9,3,pose.omega(i) );
					}
					cout << ' ' << I(4,turnbegin) << ' ' << I(4,turnend) << ' ' <<
						pose.pdb_info()->chain(turnbegin) << ' ' << pose.pdb_info()->number( turnbegin ) << ' ' <<
						filename << endl;
					{ /// write out the turn transform
						Residue const & irsd( pose.residue( turnbegin-1 ) ), &jrsd( pose.residue( turnend+1 ) ); // no reln to i,j
						Stub const
							istub( irsd.xyz("N"), irsd.xyz("CA"), irsd.xyz("C") ),
							jstub( jrsd.xyz("N"), jrsd.xyz("CA"), jrsd.xyz("C") );
						cout << "turn_transform: alpha_beta " << turnlen << ' ' << turnbb << ' ' << turnseq << ' ' <<
							RT( istub, jstub ) << ' ' << pose.pdb_info()->chain(turnbegin) << ' ' <<
							pose.pdb_info()->number( turnbegin ) << ' ' << filename << endl;
					} // scope
				}
			}
		}

	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
hairpin_test()
{
	/// read library of pdb files
	strings const files( start_files() );

	Size const strand_window( 4 );//, helix_window( 7 );

	//bool const restrict_to_bab( true );
	//Size const max_bab_seqsep( 40 );
	Size const max_turnlen( 4 );

	for ( Size fi=1; fi<= files.size(); ++fi ) {

		string const filename( files[fi] );

		Pose pose;
		pose_from_pdb( pose, filename );

		set_ss_from_dssp( filename, pose );

		add_termini_at_protein_chainbreaks( pose );

		string allbb( string("X")+torsion2big_bin_string( 1, pose.total_residue(), pose, true ) );

		runtime_assert( allbb.size() == pose.total_residue()+1 ); // so we can 1-index

		bools is_beta( pose.total_residue(), false );//, is_alpha( pose.total_residue(), false );
		Sizes strand_num( pose.total_residue(), 0 );
		{
			Size strand(0);
			bool in_strand( false );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				is_beta [i] = ( pose.residue(i).is_protein() && pose.secstruct(i) == 'E' && allbb[i] == 'B' );
				//is_alpha[i] = ( pose.residue(i).is_protein() && pose.secstruct(i) == 'H' && allbb[i] == 'A' );
				if ( pose.secstruct(i) == 'E' ) {
					if ( !in_strand ) {
						in_strand = true;
						++strand;
					}
					strand_num[i] = strand;
					TR.Trace << "strand_num: " << i << ' ' << strand << endl;
				} else in_strand = false;
			}
		}

		//bools is_hairpin_turn( pose.total_residue(), false );

		/// first identify bb NH--OOC hbonds
		Sizes
			bbn_hbond_partner ( pose.total_residue(), 0 ),
			bbn_hbond_partner2( pose.total_residue(), 0 ),
			bbo_hbond_partner ( pose.total_residue(), 0 ),
			bbo_hbond_partner2( pose.total_residue(), 0 );

		Real const bb_hbond_maxdis2( 3.4 * 3.4 );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const irsd( pose.residue(i));
			if ( !is_beta[i] ) continue;
			if ( irsd.aa() == aa_pro ) continue; // skip prolines here !!!!!!!!!!!!!!!!!!!!
			for ( Size j=1; j<= pose.total_residue(); ++j ) {
				Residue const jrsd( pose.residue(j));
				if ( !is_beta[j] ) continue;
				if ( j>=i-2 && i>= j-2 ) continue;

				if ( irsd.xyz("N").distance_squared( jrsd.xyz("O") ) <= bb_hbond_maxdis2 ) {
					if ( !bbn_hbond_partner[i] ) bbn_hbond_partner[i] = j;
					else bbn_hbond_partner2[i] = j; // bifurcation
					if ( !bbo_hbond_partner[j] ) bbo_hbond_partner[j] = i;
					else bbo_hbond_partner2[j] = i; // bifurcation
				}
			}
		}



		/// antiparallel pairing
		/// look for i NH-->OC(j)
		/// look for i OC<--NH(j)
		//bools is_beta_paired( pose.total_residue(), false );
		Vectors coords;
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_protein() ) coords.push_back( pose.residue(i).xyz("CA"));
			else coords.push_back( Vector(0,0,0));
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const irsd( pose.residue(i));
			if ( !is_beta[i] ) continue;
			if ( bbo_hbond_partner[i] > i && bbo_hbond_partner[i] == bbn_hbond_partner[i] ) {
				Size j( bbn_hbond_partner[i] ); // rsd that is paired with i, j>i

				Size const cb( chain_begin( pose.chain(i), pose ) ), ce( chain_end( pose.chain(i), pose ) );
				if ( i>cb && i<ce && j>cb && j < ce &&
						is_beta[i-1] && is_beta[i+1] && is_beta[j-1] && is_beta[j+1] ) {

					Reals distances;
					Real orientation;
					get_beta_pairing_geometry( coords, i-1, j-1, orientation, distances );
					cout << "antiparallel_beta_pairing " << F(9,3,orientation);
					for ( Size ii=1; ii<= distances.size(); ++ii ) cout << F(9,3,distances[ii]);
					cout << endl;
					get_beta_pairing_geometry( coords, j-1, i-1, orientation, distances );
					cout << "antiparallel_beta_pairing " << F(9,3,orientation);
					for ( Size ii=1; ii<= distances.size(); ++ii ) cout << F(9,3,distances[ii]);
					cout << endl;

					/// show the transform
					{
						Residue const & jrsd( pose.residue(j) );
						Stub const
							istub( irsd.xyz("N"), irsd.xyz("CA"), irsd.xyz("C") ),
							jstub( jrsd.xyz("N"), jrsd.xyz("CA"), jrsd.xyz("C") );
						cout << "antiparallel_beta_transform " << RT( istub, jstub ) << endl;
						cout << "antiparallel_beta_transform " << RT( jstub, istub ) << endl;
					}

					{ // show the other antiparallel transform?
						if ( i+2<=ce && j-2<=ce && is_beta[i+2] && is_beta[j-2] &&
								bbo_hbond_partner[i+2] == bbn_hbond_partner[i+2] &&
								bbo_hbond_partner[i+2] == j-2 ) {
							Residue const & irsd( pose.residue(i+1) ); // same name!
							Residue const & jrsd( pose.residue(j-1) );
							Stub const
								istub( irsd.xyz("N"), irsd.xyz("CA"), irsd.xyz("C") ),
								jstub( jrsd.xyz("N"), jrsd.xyz("CA"), jrsd.xyz("C") );
							cout << "antiparallel2_beta_transform " << RT( istub, jstub ) << endl;
							cout << "antiparallel2_beta_transform " << RT( jstub, istub ) << endl;
						}
					}
				}

				if ( strand_num[ j ] == strand_num[ i ]+1 && pose.chain(i) == pose.chain(j) ) {
					// sequential strands, anti-parallel paired
					Size strand1begin(i);
					while ( strand1begin>cb && is_beta[strand1begin-1] ) --strand1begin;
					// end of a strand
					// advance i until we're not B backbone anymore
					while ( i<ce && allbb[i] == 'B' ) ++i;
					if ( i==ce ) { ++i; continue; }
					runtime_assert( allbb[i] != 'B' );
					runtime_assert( allbb[i-1] == 'B' ); // i is the beginning of the turn, big-bin defn

					Size strand2end( j );
					while ( strand2end<ce && is_beta[strand2end+1] ) ++strand2end;

					while ( j>i && allbb[j-1]=='B' ) --j;
					runtime_assert( allbb[j-1] != 'B' );
					runtime_assert( allbb[j  ] == 'B' );

					Size const turnbegin( i ), turnend( j-1 ), turnlen( j-i ), strand1len( i-strand1begin ),
						strand2len( strand2end-turnend );

					TR.Trace << "strand1len: " << strand1len << " strand2len: " << strand2len <<
						I(4,turnbegin) << ' ' << I(4,turnend) << ' ' <<
						pose.pdb_info()->chain(turnbegin) << ' ' << pose.pdb_info()->number( turnbegin ) << ' ' <<
						filename << endl;

					if ( strand1len >= strand_window && strand2len >= strand_window && turnlen <= max_turnlen ) {
						// write out the turn
						string const turnbb( turnlen>0 ? allbb.substr(turnbegin,turnlen) : string("-") );
						string const turnseq( turnlen>0 ? pose.sequence().substr(turnbegin-1,turnlen) : string("-") );
						string const beforeseq( pose.sequence().substr(turnbegin-1-strand_window,strand_window) );
						string const afterseq( pose.sequence().substr(turnend,strand_window) );
						Size const coordslen( turnlen + strand_window + strand_window );
						cout << "turn_coords: beta_beta " << turnlen << ' ' << turnbb << ' ' <<
							beforeseq << ' ' << turnseq << ' ' << afterseq << ' ' << coordslen;
						// show ca-coords for flanking sses and turn
						for ( Size i=turnbegin-strand_window; i<= turnend + strand_window; ++i ) {
							cout <<
								F(9,3,pose.residue(i).xyz("CA").x()) <<
								F(9,3,pose.residue(i).xyz("CA").y()) <<
								F(9,3,pose.residue(i).xyz("CA").z());
						}
						// show the torsion angles:
						cout << " bb_torsions: ";
						for ( Size i=turnbegin; i<= turnend; ++i ) {
							cout << F(9,3,pose.phi(i)) << F(9,3,pose.psi(i)) << F(9,3,pose.omega(i) );
						}
						cout << I(4,turnbegin) << ' ' << I(4,turnend) << ' ' <<
							pose.pdb_info()->chain(turnbegin) << ' ' << pose.pdb_info()->number( turnbegin ) << ' ' <<
							filename << endl;
					}
					i = j; // skip past the turn
				}
			} // antiparallel beta-pairing
		} // i=1,nres
	} // files
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
bb_strain_test()
{
	/// read library of pdb files
	strings const files( start_files() );

	Size const min_turn_buffer( 7 );
	Size const helix_window( 7 );

	Size const max_turn_length( option[ my_options::max_turn_length ].user() ?
		option[ my_options::max_turn_length ] : 7 );

	/// store downstream helix-atom coords in defined reference frames
	strings turns;
	if ( option[ my_options::turns ].user() ) turns = option[ my_options::turns ]();

	for ( Size fi=1; fi<= files.size(); ++fi ) {

		string const filename( files[fi] );

		Pose pose;

		cout << "try to read " << filename << endl;
		fflush(stdout);

		try {
			pose_from_pdb( pose, filename );

			set_ss_from_dssp( filename, pose );

			add_termini_at_protein_chainbreaks( pose );

		} catch ( utility::excn::Exception const & e ) {
			std::cout << "caught exception " << e.msg() << ' ' << filename << std::endl;
			std::cerr << "caught exception " << e.msg() << ' ' << filename << std::endl;
			continue;
		}


		if ( pose.total_residue()<min_turn_buffer ) continue; // mostly filter out empty poses


		for ( Size i=1; i< pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_protein() && ( !pose.residue(i+1).is_protein() ||
					pose.residue(i).is_upper_terminus() ) ) {
				pose.conformation().insert_chain_ending(i);
			}
			if ( pose.residue(i+1).is_protein() && !pose.residue(i).is_protein() ) {
				pose.conformation().insert_chain_ending(i);
			}
		}

		{ // add some coords using the same SS definitions used during some other design code
			SizePairs helices, strands, antiparallel_bridges, parallel_bridges;
			bools subset( pose.total_residue(), false );
			for ( Size i=1; i<= pose.total_residue(); ++i ) subset[i] = ( pose.residue(i).is_protein() );
			get_backbone_hbond_interactions( subset, pose, helices, strands, antiparallel_bridges, parallel_bridges );

			foreach_ ( SizePair const & h, helices ) {
				Size const i( h.first ), j( h.second ), helixlen( j-i+1 );
				cout << "my_helix_coords: " << I(4,helixlen) << ' ' << pose.sequence().substr(i-1,j-i+1);
				for ( Size ii=i; ii<= j; ++ii ) {
					cout <<
						F(9,3,pose.residue(ii).xyz("CA").x()) <<
						F(9,3,pose.residue(ii).xyz("CA").y()) <<
						F(9,3,pose.residue(ii).xyz("CA").z());
				}
				cout << ' ' << pose.pdb_info()->number(i) << ' ' << pose.pdb_info()->chain(i) <<
					pose.pdb_info()->number(j) << ' ' << pose.pdb_info()->chain(j) << ' ' << filename << endl;
			}

			foreach_ ( SizePair const & h, strands ) {
				Size const i( h.first ), j( h.second ), strandlen( j-i+1 );
				cout << "my_strand_coords: " << I(4,strandlen) << ' ' << pose.sequence().substr(i-1,j-i+1);
				for ( Size ii=i; ii<= j; ++ii ) {
					cout <<
						F(9,3,pose.residue(ii).xyz("CA").x()) <<
						F(9,3,pose.residue(ii).xyz("CA").y()) <<
						F(9,3,pose.residue(ii).xyz("CA").z());
				}
				cout << ' ' << pose.pdb_info()->number(i) << ' ' << pose.pdb_info()->chain(i) <<
					pose.pdb_info()->number(j) << ' ' << pose.pdb_info()->chain(j) << ' ' << filename << endl;
			}
		}



		// get alpha helices
		bools bb_hbond_ip4( pose.total_residue(), false );
		for ( Size i=1; i<= pose.total_residue()-4; ++i ) {
			Size const j(i+4);
			if ( pose.residue(i).is_protein() && pose.residue(j).is_protein() && pose.chain(i) == pose.chain(j) &&
					pose.residue(i).xyz("O").distance_squared( pose.residue(j).xyz("N") ) < 3.4*3.4 ) {
				bb_hbond_ip4[i] = true;
			}
		}

		// these are going to be 1-indexed now
		string const all_bb( string("X")+torsion2big_bin_string( 1, pose.total_residue(), pose, true ) );
		string const all_ss( string("X")+pose.secstruct() );

		vector1< Sizes > all_helices;

		Vectors ca_coords;
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_protein() ) ca_coords.push_back( pose.residue(i).xyz("CA") );
			else ca_coords.push_back( Vector(0,0,0) );
		}

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_protein() && all_ss[i] == 'H' && all_bb[i] == 'A' ) {
				// start of an alpha helix
				// look for the end
				Size j(i);
				while ( j<= pose.total_residue() && pose.residue(j).is_protein() && pose.chain(i) == pose.chain(j) &&
						all_ss[j] == 'H' && all_bb[j] == 'A' ) ++j;
				--j;
				runtime_assert( all_bb[j] == 'A' && all_ss[j] == 'H' && pose.chain(i) == pose.chain(j) );

				Size const helixlen( j-i+1 );
				{ // write out the helix coords
					cout << "helix_coords: " << I(4,helixlen) << ' ' << pose.sequence().substr(i-1,j-i+1);
					for ( Size ii=i; ii<= j; ++ii ) {
						cout <<
							F(9,3,pose.residue(ii).xyz("CA").x()) <<
							F(9,3,pose.residue(ii).xyz("CA").y()) <<
							F(9,3,pose.residue(ii).xyz("CA").z());
					}
					cout << ' ' << pose.pdb_info()->number(i) << ' ' << pose.pdb_info()->chain(i) <<
						pose.pdb_info()->number(j) << ' ' << pose.pdb_info()->chain(j) << ' ' << filename << endl;
				}

				if ( option[ my_options::pdb_coords_file ].user() && helixlen >= helix_window ) {
					Real avg_rmsd, med_rmsd;
					compute_helix_strain( i, j, ca_coords, avg_rmsd, med_rmsd );
					cout << "helix_strain: " << I(4,helixlen) << F(9,3,avg_rmsd)<<F(9,3,med_rmsd)<< ' ' <<
						pose.pdb_info()->number(i) << ' ' << pose.pdb_info()->chain(i) <<
						pose.pdb_info()->number(j) << ' ' << pose.pdb_info()->chain(j) << ' ' << filename << endl;
				}

				all_helices.push_back( make_vector1( i,j ) );
				i=j; // gets incremented
			}
		}

		/// look for beta strands, too
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_protein() && all_ss[i] == 'E' && all_bb[i] == 'B' ) {
				// start of an alpha helix
				// look for the end
				Size j(i);
				while ( j<= pose.total_residue() && pose.residue(j).is_protein() && pose.chain(i) == pose.chain(j) &&
						all_ss[j] == 'E' && all_bb[j] == 'B' ) ++j;
				--j;
				runtime_assert( all_bb[j] == 'B' && all_ss[j] == 'E' && pose.chain(i) == pose.chain(j) );

				Size const strandlen( j-i+1 );
				{ // write out the strand coords
					cout << "strand_coords: " << I(4,strandlen) << ' ' << pose.sequence().substr(i-1,j-i+1);
					for ( Size ii=i; ii<= j; ++ii ) {
						cout <<
							F(9,3,pose.residue(ii).xyz("CA").x()) <<
							F(9,3,pose.residue(ii).xyz("CA").y()) <<
							F(9,3,pose.residue(ii).xyz("CA").z());
					}
					cout << ' ' << pose.pdb_info()->number(i) << ' ' << pose.pdb_info()->chain(i) <<
						pose.pdb_info()->number(j) << ' ' << pose.pdb_info()->chain(j) << ' ' << filename << endl;
				}
				i=j; // gets incremented
			}
		}

		//// find "clean helices" with all bb-hbonds
		for ( Size i=1; i<= pose.total_residue()-4; ++i ) {
			if ( pose.residue(i).is_protein() &&
					all_ss.substr(i,5) == "HHHHH" &&
					all_bb.substr(i,5) == "AAAAA" &&
					pose.residue(i+4).is_protein() &&
					pose.chain(i+4) == pose.chain(i) &&
					bb_hbond_ip4[i] ) {
				// start of an alpha helix
				// look for the end
				Size j(i);
				while ( j<= pose.total_residue()-4 &&
						pose.residue(j+4).is_protein() &&
						pose.chain(j+4) == pose.chain(i) &&
						all_ss.substr(j,5) == "HHHHH" &&
						all_bb.substr(j,5) == "AAAAA" &&
						bb_hbond_ip4[j] ) ++j;
				--j;
				Size const helixbegin(i), helixend(j+4), helixlen( helixend-helixbegin+1 );
				runtime_assert( all_bb.substr(helixbegin,helixlen) == string(helixlen,'A') );
				runtime_assert( all_ss.substr(helixbegin,helixlen) == string(helixlen,'H') );
				runtime_assert( pose.chain(helixbegin) == pose.chain(helixend) );

				{ // write out the helix coords
					cout << "clean_helix_coords: " << I(4,helixlen) << ' ' << pose.sequence().substr(helixbegin-1,helixlen);
					for ( Size ii=helixbegin; ii<= helixend; ++ii ) {
						cout <<
							F(9,3,pose.residue(ii).xyz("CA").x()) <<
							F(9,3,pose.residue(ii).xyz("CA").y()) <<
							F(9,3,pose.residue(ii).xyz("CA").z());
					}
					// show the torsion angles:
					cout << " bb_torsions: ";
					for ( Size ii=helixbegin; ii<= helixend; ++ii ) {
						cout << F(9,3,pose.phi(ii)) << F(9,3,pose.psi(ii)) << F(9,3,pose.omega(ii) );
					}
					cout << ' ' << pose.pdb_info()->number(helixbegin) << ' ' << pose.pdb_info()->chain(helixbegin) <<
						pose.pdb_info()->number(helixend) << ' ' << pose.pdb_info()->chain(helixend) << ' ' << filename << endl;
				}

				if ( option[ my_options::pdb_coords_file ].user() && helixlen >= helix_window ) {
					Real avg_rmsd, med_rmsd;
					compute_helix_strain( helixbegin, helixend, ca_coords, avg_rmsd, med_rmsd );
					cout << "clean_helix_strain: " << I(4,helixlen) << F(9,3,avg_rmsd)<<F(9,3,med_rmsd)<< ' ' <<
						pose.pdb_info()->number(i) << ' ' << pose.pdb_info()->chain(i) <<
						pose.pdb_info()->number(j) << ' ' << pose.pdb_info()->chain(j) << ' ' << filename << endl;
				}
				i=j; // gets incremented
			}
		}

		for ( Size ii=1; ii< all_helices.size(); ++ii ) {
			Size h1_end( all_helices[ii][2] ), h2_begin( all_helices[ii+1][1] );
			if ( pose.chain(h1_end) != pose.chain(h2_begin) ) continue;
			while ( all_bb[   h1_end+1 ] == 'A' && h1_end<h2_begin ) ++h1_end;
			while ( all_bb[ h2_begin-1 ] == 'A' && h1_end<h2_begin ) --h2_begin;
			if ( h1_end == h2_begin ) {
				cout << "whoah " << all_helices[ii][2] << ' ' << h1_end << ' ' << all_helices[ii+1][1] << ' ' << h2_begin <<
					endl;
				continue;
			}
			string const turn( all_bb.substr(h1_end+1,h2_begin-h1_end-1));
			string const turnseq( pose.sequence().substr(h1_end,h2_begin-h1_end-1));
			Size const h1_len( h1_end - all_helices[ii][1]+1 ), h2_len( all_helices[ii+1][2] - h2_begin+1 );
			runtime_assert( turn[0] !='A' && turn[ turn.size()-1 ] != 'A' );
			TR.Trace << "found turn " << turn << ' ' << filename << endl;
			if ( h1_len >= min_turn_buffer && h2_len >= min_turn_buffer &&
					turn.size() <= max_turn_length &&
					( turns.empty() || has_element( turns, turn ) ) ) {
				string const beforeseq( pose.sequence().substr( h1_end-min_turn_buffer,min_turn_buffer) ),
					afterseq( pose.sequence().substr( h2_begin-1, min_turn_buffer ) );

				cout << "turn_coords: " << turn << ' ' << beforeseq << ' ' << turnseq << ' '<< afterseq << ' ' <<
					turn.size() + 2*min_turn_buffer;
				for ( Size i=h1_end-min_turn_buffer+1; i <= h2_begin+min_turn_buffer-1; ++i ) {
					cout <<
						F(9,3,pose.residue(i).xyz("CA").x()) <<
						F(9,3,pose.residue(i).xyz("CA").y()) <<
						F(9,3,pose.residue(i).xyz("CA").z());
				}
				// show the torsion angles:
				cout << " bb_torsions: ";
				for ( Size i=h1_end+1; i< h2_begin; ++i ) {
					cout << F(9,3,pose.phi(i)) << F(9,3,pose.psi(i)) << F(9,3,pose.omega(i) );
				}
				cout << ' ' << pose.pdb_info()->chain(h1_end) << ' ' << pose.pdb_info()->number( h1_end ) <<' ' <<
					filename << endl;

				if ( option[ my_options::pdb_coords_file ].user() ) {
					Real avg_rmsd, med_rmsd;
					compute_turn_strain( turn, h1_end+1, ca_coords, avg_rmsd, med_rmsd );
					cout << "turn_strain: " << turn << F(9,3,avg_rmsd)<<F(9,3,med_rmsd)<< ' ' <<
						pose.pdb_info()->number(h1_end+1) << ' ' << pose.pdb_info()->chain(h1_end+1) << filename << endl;
				}

			}
		}
	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
refold_toroid_test()
{
	// the new correct starting sequence
	string const fullseq("GVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLG");


	//string const fullseq("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLG");
	//string const fullseq("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGGW");
	//string const fullseq("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLG");
	// string const fullseq("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGGWLEH");

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );

	string const simfile( shared_output_tag()+"_refold_toroid1" ), worktag("tmp1");

	// create a starting pose from the sequence
	Pose pose;
	for ( Size i=0; i< fullseq.size(); ++i ) {
		pose.append_residue_by_bond( *get_vanilla_protein_residue( fullseq[i] ), true );
	}
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

	core::fragment::FragSetOP small_frags(0), large_frags(0);
	pick_nnmake_fragments_from_single_sequence( pose.sequence(), small_frags, large_frags );

	Pose const start_pose( pose );

	while ( true ) {
		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		pose = start_pose;
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.set_phi  ( i, init_phi   );
			pose.set_psi  ( i, init_psi   );
			pose.set_omega( i, init_omega );
			pose.set_secstruct( i, 'L' );
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

		{
			kinematics::MoveMapOP mm( new kinematics::MoveMap() );
			mm->set_bb( true );
			protocols::abinitio::ClassicAbinitio abinitio( small_frags, large_frags, mm );
			abinitio.init( pose );
			if ( !dry_run() ) abinitio.apply( pose );
		}


		/// now fastrelax
		devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
		{
			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
			movemap->set_bb (true);
			movemap->set_chi(true);
			fastrelax.set_movemap( movemap );
			if ( !dry_run() ) fastrelax.apply( pose );
		}


		Real const final_score( ( *fa_scorefxn )( pose ) );


		string const outfilename( output_tag() + "refold_toroid1_N"+ lead_zero_string_of( n, 4 )+".pdb" );


		bool const passed_score_filter
			( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		// run phenix/phaser
		Real llg(0.0), tfz(0.0), pak(0.0), rphaser( 1.0 ), rwork(1.0), rfree(1.0), coord_error(0.0), xray_target(0.0);
		string space_group("X1"), refined_pdbfile("none");
		if ( passed_score_filter ) {
			pose.dump_pdb( outfilename );
			run_phaser_and_phenix_refinement( outfilename, option[ my_options::llg_threshold_for_refinement ],
				llg, tfz, pak, rphaser,
				rwork, rfree, coord_error, xray_target, space_group, refined_pdbfile );
		}


		cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" passed_score_filter: " << passed_score_filter <<
			" llg: " << F(9,3,llg) <<
			" tfz: " << F(9,3,tfz) <<
			" rwork: " << F(9,3,rwork) <<
			" rfree: " << F(9,3,rfree) <<
			" coord_error: " << F(9,3,coord_error) <<
			" xray_target: " << F(9,3,xray_target) <<
			" refined_pdbfile: " << refined_pdbfile <<
			endl;

		fflush( stdout );

	}

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
refold_9x2_test()
{
	runtime_assert( !option[ basic::options::OptionKeys::out::file::output_virtual ] ); // no NV in pdb files
	confirm_phenix_version();
	runtime_assert( option[ my_options::llg_threshold_for_refinement ].user() );

	// the new correct starting sequence
	string const fullseq("GISVEELLKLAKAAYYSGTTVEEAYKLALKLGISVEELLKLAEAAYYSGTTVEEAYKLALKLGISVEELLKLAKAAYYSGTTVEEAYKLALKLG");
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ) );

	string const simfile( shared_output_tag()+"_refold_9x2.work" ), worktag("tmp1");

	// create a starting pose from the sequence
	Pose pose;
	for ( Size i=0; i< fullseq.size(); ++i ) {
		pose.append_residue_by_bond( *get_vanilla_protein_residue( fullseq[i] ), true );
	}
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

	core::fragment::FragSetOP small_frags(0), large_frags(0);
	pick_nnmake_fragments_from_single_sequence( pose.sequence(), small_frags, large_frags );

	Pose const start_pose( pose );

	while ( true ) {
		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		pose = start_pose;
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.set_phi  ( i, init_phi   );
			pose.set_psi  ( i, init_psi   );
			pose.set_omega( i, init_omega );
			pose.set_secstruct( i, 'L' );
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

		{
			kinematics::MoveMapOP mm( new kinematics::MoveMap() );
			mm->set_bb( true );
			protocols::abinitio::ClassicAbinitio abinitio( small_frags, large_frags, mm );
			abinitio.init( pose );
			if ( !dry_run() ) abinitio.apply( pose );
		}


		/// now fastrelax
		devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
		{
			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
			movemap->set_bb (true);
			movemap->set_chi(true);
			fastrelax.set_movemap( movemap );
			if ( !dry_run() ) fastrelax.apply( pose );
		}


		Real const final_score( ( *fa_scorefxn )( pose ) );


		string const outfilename( output_tag() + "refold_toroid1_N"+ lead_zero_string_of( n, 4 )+".pdb" );


		bool const passed_score_filter
			( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		// run phenix/phaser
		Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
		string space_group("X1"), refined_pdbfile("none");
		if ( passed_score_filter ) {
			pose.dump_pdb( outfilename );
			run_phaser_and_phenix_refinement_9x2( outfilename, option[ my_options::llg_threshold_for_refinement ],
				llg, tfz, rwork, rfree,
				space_group, refined_pdbfile );
		}


		cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" passed_score_filter: " << passed_score_filter <<
			" llg: " << F(9,3,llg) <<
			" tfz: " << F(9,3,tfz) <<
			" rwork: " << F(9,3,rwork) <<
			" rfree: " << F(9,3,rfree) <<
			" refined_pdbfile: " << refined_pdbfile <<
			endl;

		fflush( stdout );
		check_simtime();

	}

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
refold_3x2_test()
{
	runtime_assert( !option[ basic::options::OptionKeys::out::file::output_virtual ] ); // no NV in pdb files
	confirm_phenix_version();
	runtime_assert( option[ my_options::llg_threshold_for_refinement ].user() );

	// the new correct starting sequence
	string const fullseq("GKSPTEVLLELIAEASGTTREEVKEKFLKELRKGKSPTEVLLELIAEASGTTKEEVKEKFLKELSFGKSPTEVLLELIAEASGTTKEEVKKKFWKELSL");
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ) );

	string const simfile( shared_output_tag()+"_refold_3x2.work" ), worktag("tmp1");

	// create a starting pose from the sequence
	Pose pose;
	for ( Size i=0; i< fullseq.size(); ++i ) {
		pose.append_residue_by_bond( *get_vanilla_protein_residue( fullseq[i] ), true );
	}
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

	core::fragment::FragSetOP small_frags(0), large_frags(0);
	pick_nnmake_fragments_from_single_sequence( pose.sequence(), small_frags, large_frags );

	Pose const start_pose( pose );

	while ( true ) {
		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		pose = start_pose;
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.set_phi  ( i, init_phi   );
			pose.set_psi  ( i, init_psi   );
			pose.set_omega( i, init_omega );
			pose.set_secstruct( i, 'L' );
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

		{
			kinematics::MoveMapOP mm( new kinematics::MoveMap() );
			mm->set_bb( true );
			protocols::abinitio::ClassicAbinitio abinitio( small_frags, large_frags, mm );
			abinitio.init( pose );
			if ( !dry_run() ) abinitio.apply( pose );
		}


		/// now fastrelax
		devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
		{
			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
			movemap->set_bb (true);
			movemap->set_chi(true);
			fastrelax.set_movemap( movemap );
			if ( !dry_run() ) fastrelax.apply( pose );
		}


		Real const final_score( ( *fa_scorefxn )( pose ) );


		string const outfilename( output_tag() + "refold_toroid1_N"+ lead_zero_string_of( n, 4 )+".pdb" );


		bool const passed_score_filter
			( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		// run phenix/phaser
		Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
		string space_group("X1"), refined_pdbfile("none");
		if ( passed_score_filter ) {
			pose.dump_pdb( outfilename );
			run_phaser_and_phenix_refinement_3x2( outfilename, option[ my_options::llg_threshold_for_refinement ],
				llg, tfz, rwork, rfree,
				space_group, refined_pdbfile );
		}


		cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" passed_score_filter: " << passed_score_filter <<
			" llg: " << F(9,3,llg) <<
			" tfz: " << F(9,3,tfz) <<
			" rwork: " << F(9,3,rwork) <<
			" rfree: " << F(9,3,rfree) <<
			" refined_pdbfile: " << refined_pdbfile <<
			endl;

		fflush( stdout );
		check_simtime();

	}

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
refold_tordim1_test()
{
	runtime_assert( !option[ basic::options::OptionKeys::out::file::output_virtual ] ); // no NV in pdb files
	confirm_phenix_version();
	runtime_assert( option[ my_options::llg_threshold_for_refinement ].user() );

	// the new correct starting sequence
	string const fullseq("GVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLG");
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ) );

	string simfile, worktag("tmp1");

	if ( option[ my_options::local_simfile ] ) {
		runtime_assert( shared_output_tag().substr(0,3) == "../" );
		simfile = shared_output_tag().substr(1)+".work";
	} else {
		simfile = shared_output_tag()+"_refold_3x2.work";
	}

	// create a starting pose from the sequence
	Pose pose;
	for ( Size i=0; i< fullseq.size(); ++i ) {
		pose.append_residue_by_bond( *get_vanilla_protein_residue( fullseq[i] ), true );
	}
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

	core::fragment::FragSetOP small_frags(0), large_frags(0);
	pick_nnmake_fragments_from_single_sequence( pose.sequence(), small_frags, large_frags );

	Pose const start_pose( pose );

	while ( true ) {
		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		pose = start_pose;
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.set_phi  ( i, init_phi   );
			pose.set_psi  ( i, init_psi   );
			pose.set_omega( i, init_omega );
			pose.set_secstruct( i, 'L' );
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

		{
			kinematics::MoveMapOP mm( new kinematics::MoveMap() );
			mm->set_bb( true );
			protocols::abinitio::ClassicAbinitio abinitio( small_frags, large_frags, mm );
			abinitio.init( pose );
			if ( !dry_run() ) abinitio.apply( pose );
		}


		/// now fastrelax
		devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
		{
			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
			movemap->set_bb (true);
			movemap->set_chi(true);
			fastrelax.set_movemap( movemap );
			if ( !dry_run() ) fastrelax.apply( pose );
		}


		Real const final_score( ( *fa_scorefxn )( pose ) );


		string const outfilename( output_tag() + "refold_toroid1_N"+ lead_zero_string_of( n, 4 )+".pdb" );


		bool const passed_score_filter
			( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		// run phenix/phaser
		Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
		string space_group("X1"), refined_pdbfile("none");
		if ( passed_score_filter ) {
			pose.dump_pdb( outfilename );
			run_phaser_and_phenix_refinement_tordim1( outfilename, option[ my_options::llg_threshold_for_refinement ],
				llg, tfz, rwork, rfree,
				space_group, refined_pdbfile );
			run_command("rm "+outfilename);
			run_command("gzip "+outfilename+"*");
		}


		cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" passed_score_filter: " << passed_score_filter <<
			" llg: " << F(9,3,llg) <<
			" tfz: " << F(9,3,tfz) <<
			" rwork: " << F(9,3,rwork) <<
			" rfree: " << F(9,3,rfree) <<
			" refined_pdbfile: " << refined_pdbfile <<
			endl;

		fflush( stdout );
		check_simtime();

	}

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
refold_tri_test()
{
	runtime_assert( !option[ basic::options::OptionKeys::out::file::output_virtual ] ); // no NV in pdb files
	confirm_phenix_version();
	runtime_assert( option[ my_options::llg_threshold_for_refinement ].user() );

	// the new correct starting sequence
	string const fullseq("PELDKILARNPELKKILERNPELAKILERNPELAKILERNPELAKILERNPELAKILERNPELAKILSVNPELAKILERNPDLAAVLERNPEAALELEKN"); // removed "GN" at Nterminus
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ) );

	string const simfile( shared_output_tag()+"_refold_tri_"+
		lead_zero_string_of( numeric::random::random_range( 1, option[my_options::num_simfiles]() ), 3 )+".work" ), worktag("tmp1");

	// create a starting pose from the sequence
	Pose pose;
	for ( Size i=0; i< fullseq.size(); ++i ) {
		pose.append_residue_by_bond( *get_vanilla_protein_residue( fullseq[i] ), true );
	}
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

	core::fragment::FragSetOP small_frags(0), large_frags(0);
	pick_nnmake_fragments_from_single_sequence( pose.sequence(), small_frags, large_frags );

	Pose const start_pose( pose );

	while ( true ) {
		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		pose = start_pose;
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.set_phi  ( i, init_phi   );
			pose.set_psi  ( i, init_psi   );
			pose.set_omega( i, init_omega );
			pose.set_secstruct( i, 'L' );
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

		{
			kinematics::MoveMapOP mm( new kinematics::MoveMap() );
			mm->set_bb( true );
			protocols::abinitio::ClassicAbinitio abinitio( small_frags, large_frags, mm );
			abinitio.init( pose );
			if ( !dry_run() ) abinitio.apply( pose );
		}


		/// now fastrelax
		devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
		{
			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
			movemap->set_bb (true);
			movemap->set_chi(true);
			fastrelax.set_movemap( movemap );
			if ( !dry_run() ) fastrelax.apply( pose );
		}


		Real const final_score( ( *fa_scorefxn )( pose ) );


		string const outfilename( output_tag() + "refold_toroid1_N"+ lead_zero_string_of( n, 4 )+".pdb" );


		bool const passed_score_filter
			( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		// run phenix/phaser
		Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
		string space_group("X1"), refined_pdbfile("none");
		Size actual_num_models(0), num_models(0);
		Size const target_num_models( option[ my_options::num_models ].user() ? option[ my_options::num_models ]() : 0 );
		Size const composition( option[ my_options::composition ].user() ? option[ my_options::composition ]() : 0 );

		if ( passed_score_filter ) {
			pose.dump_pdb( outfilename );
			run_phaser_and_phenix_refinement_tri( string(""), outfilename, option[ my_options::llg_threshold_for_refinement ],
				actual_num_models, num_models, llg, tfz, rwork, rfree,
				space_group, refined_pdbfile, target_num_models, composition );
			run_command("gzip "+outfilename+"*");
		}


		cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" passed_score_filter: " << passed_score_filter <<
			" llg: " << F(9,3,llg) <<
			" tfz: " << F(9,3,tfz) <<
			" rwork: " << F(9,3,rwork) <<
			" rfree: " << F(9,3,rfree) <<
			" actual_num_models: " << actual_num_models <<
			" num_models: " << num_models <<
			" space_group: " << space_group <<
			" refined_pdbfile: " << refined_pdbfile <<
			endl;

		fflush( stdout );
		check_simtime();

	}

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
refold_generic_test()
{
	runtime_assert( !option[ basic::options::OptionKeys::out::file::output_virtual ] ); // no NV in pdb files
	runtime_assert( option[ my_options::local_simfile ] );

	confirm_phenix_version();

	runtime_assert( utility::file::file_exists( option[ my_options::mtz_file ] ) );
	runtime_assert( option[ my_options::num_models ].user() );
	runtime_assert( option[ my_options::composition ].user() );
	runtime_assert( option[ my_options::llg_threshold_for_refinement ].user());
	runtime_assert( option[ my_options::space_group ].user() ||
		option[ my_options::space_groups ].user());

	// the new correct starting sequence
	string const fullseq( option[ my_options::refoldseq ] );

	string phenixseq( fullseq );
	if ( option[my_options::phenixseq].user() ) {
		phenixseq = option[my_options::phenixseq];
		// runtime_assert( phenixseq.size() == fullseq.size() );
		runtime_assert( phenixseq != fullseq ); // no point otherwise!
	}

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ) );

	string simfile, worktag("tmp1");
	if ( option[ my_options::local_simfile ] ) {
		runtime_assert( shared_output_tag().substr(0,3) == "../" );
		simfile = shared_output_tag().substr(1)+".work";
	} else {
		simfile = shared_output_tag()+"_refold_generic_"+
			lead_zero_string_of( numeric::random::random_range( 1, option[my_options::num_simfiles]() ), 3 )+".work";
	}

	// create a starting pose from the sequence
	Pose pose;
	for ( Size i=0; i< fullseq.size(); ++i ) {
		pose.append_residue_by_bond( *get_vanilla_protein_residue( fullseq[i] ), true );
	}
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

	core::fragment::FragSetOP small_frags(0), large_frags(0);
	pick_nnmake_fragments_from_single_sequence( pose.sequence(), small_frags, large_frags );

	Pose const start_pose( pose );

	while ( true ) {
		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		pose = start_pose;
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.set_phi  ( i, init_phi   );
			pose.set_psi  ( i, init_psi   );
			pose.set_omega( i, init_omega );
			pose.set_secstruct( i, 'L' );
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

		{
			kinematics::MoveMapOP mm( new kinematics::MoveMap() );
			mm->set_bb( true );
			protocols::abinitio::ClassicAbinitio abinitio( small_frags, large_frags, mm );
			abinitio.init( pose );
			if ( !dry_run() ) abinitio.apply( pose );
		}


		/// now fastrelax
		devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
		{
			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
			movemap->set_bb (true);
			movemap->set_chi(true);
			fastrelax.set_movemap( movemap );
			if ( !dry_run() ) fastrelax.apply( pose );
		}


		Real const final_score( ( *fa_scorefxn )( pose ) );


		string const outfilename( output_tag() + "refold_toroid1_N"+ lead_zero_string_of( n, 4 )+".pdb" );


		bool const passed_score_filter
			( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		// run phenix/phaser
		Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
		string actual_space_group("X1"), refined_pdbfile("none");
		Size actual_num_models(0);

		string target_space_group;
		if ( option[ my_options::space_group ].user() ) {
			target_space_group = option[ my_options::space_group ]();
		} else {
			runtime_assert( option[ my_options::space_groups ].user() );
			target_space_group = random_element( option[ my_options::space_groups ]() );
		}

		if ( passed_score_filter ) {
			if ( phenixseq != fullseq ) {
				runtime_assert( pose.sequence() == fullseq );

				if ( phenixseq.size() < fullseq.size() ) {
					pose.conformation().delete_residue_range_slow( phenixseq.size()+1, pose.total_residue() );
				}

				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					if ( phenixseq[i-1] != fullseq[i-1] ) { // need to make a mutation
						make_sequence_change( i, aa_from_oneletter_code(phenixseq[i-1]), pose );
					}
				}
				runtime_assert( pose.sequence() == phenixseq );
			}

			pose.dump_pdb( outfilename );
			string const fasta_file( outfilename+"_tmp.fasta" );
			{
				ofstream out( fasta_file.c_str() );
				out << ">tmp\n" << phenixseq << endl;
				out.close();
			}
			run_phaser_and_phenix_refinement_generic( outfilename,
				option[ my_options::mtz_file ],
				fasta_file,
				target_space_group,
				option[ my_options::num_models ],
				option[ my_options::composition ],
				option[ my_options::llg_threshold_for_refinement ],
				actual_num_models, llg, tfz, rwork, rfree,
				actual_space_group, refined_pdbfile );
			run_command("gzip "+outfilename+"*");
		}


		cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" passed_score_filter: " << passed_score_filter <<
			" llg: " << F(9,3,llg) <<
			" tfz: " << F(9,3,tfz) <<
			" rwork: " << F(9,3,rwork) <<
			" rfree: " << F(9,3,rfree) <<
			" actual_num_models: " << actual_num_models <<
			" target_num_models: " << option[ my_options::num_models ] <<
			" space_group: " << actual_space_group <<
			" target_space_group: " << target_space_group <<
			" refined_pdbfile: " << refined_pdbfile <<
			endl;

		fflush( stdout );
		check_simtime();

	}

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// read in decoys, compute rmsds to the tri design
//
void
tri_rmsd_test()
{

	Vectors const tri_coords( read_CA_coords_from_file( "/home/pbradley/csdat/phenix/tri/chainA.pdb" ) );
	runtime_assert( tri_coords.size() == 100 );

	strings const logfiles( start_files() );

	foreach_ ( string const logfile, logfiles ) {

		ifstream data( logfile.c_str() );

		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			string const dirtag( l[1] );
			string filename = dirtag.substr( 0, dirtag.find_last_of( "/" ) ) +"/" + l[3];

			if ( !utility::file::file_exists(filename) &&
					utility::file::file_exists(filename+".gz") ) filename += ".gz";

			if ( !utility::file::file_exists(filename) ) {
				cerr << "missing: " << filename << ' ' << logfile << endl;
				continue;
			}

			Vectors const coords( read_CA_coords_from_file( filename ) );

			Size const nres( coords.size() );
			if ( nres%10 != 0 || nres>100 ) {
				cerr << "bad nres: " << nres << ' '<< filename << endl;
				continue;
			}

			Vectors ref_coords( tri_coords );
			ref_coords.resize( nres );
			Real const rmsd( numeric::model_quality::calc_rms( coords, ref_coords ) );

			cout << line << " nres: " << nres << " rmsd: " << F(9,3,rmsd) << endl;
		}
		data.close();
	}


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
symmetric_refold_tri_test()
{
	runtime_assert( !option[ basic::options::OptionKeys::out::file::output_virtual ] ); // no NV in pdb files
	confirm_phenix_version();
	runtime_assert( option[ my_options::llg_threshold_for_refinement ].user() );

	strings space_groups;
	if ( option[ my_options::space_groups ].user() ) space_groups = option[ my_options::space_groups ]();

	// the new correct starting sequence
	string const fullseq("PELDKILARNPELKKILERNPELAKILERNPELAKILERNPELAKILERNPELAKILERNPELAKILSVNPELAKILERNPDLAAVLERNPEAALELEKN"); // removed "GN" at Nterminus

	runtime_assert( fullseq.size() == 100 );
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ) );

	string simfile, worktag("tmp1");
	if ( option[ my_options::local_simfile ] ) {
		runtime_assert( shared_output_tag().substr(0,3) == "../" );
		simfile = shared_output_tag().substr(1)+".work";
	} else {
		simfile = shared_output_tag()+"_refold_tri_"+
			lead_zero_string_of( numeric::random::random_range( 1, option[my_options::num_simfiles]() ), 3 )+".work";
	}

	// create a starting pose from the sequence
	Pose pose;
	for ( Size i=0; i< fullseq.size(); ++i ) {
		pose.append_residue_by_bond( *get_vanilla_protein_residue( fullseq[i] ), true );
	}
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

	// core::fragment::FragSetOP small_frags(0), large_frags(0);
	devel::blab::classic_frags::FragLibOP fraglib( new devel::blab::classic_frags::FragLib() );
	pick_nnmake_fragments_from_single_sequence( pose.sequence(), *fraglib );

	Pose const start_pose( pose );

	while ( true ) {
		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		pose = start_pose;
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.set_phi  ( i, init_phi   );
			pose.set_psi  ( i, init_psi   );
			pose.set_omega( i, init_omega );
			pose.set_secstruct( i, 'L' );
		}
		for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

		{
			runtime_assert( pose.total_residue() == 100 );
			symminfo_hack::nrepeat_ = 10;
			symminfo_hack::repeatlen_ = 10;
			symminfo_hack::base_repeat_ = 5;
			Sizes const fragseq_poslist;
			if ( !dry_run() ) simple_fold_abinitio( *fraglib, fragseq_poslist, pose );
			symminfo_hack::nrepeat_ = 0;
			symminfo_hack::repeatlen_ = 0;
			symminfo_hack::base_repeat_ = 0;
		}


		/// now fastrelax
		devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
		{
			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
			movemap->set_bb (true);
			movemap->set_chi(true);
			fastrelax.set_movemap( movemap );
			if ( !dry_run() ) fastrelax.apply( pose );
		}


		Real const final_score( ( *fa_scorefxn )( pose ) );


		string const outfilename( output_tag() + "refold_toroid1_N"+ lead_zero_string_of( n, 4 )+".pdb" );


		bool const passed_score_filter
			( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		// run phenix/phaser
		Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
		string space_group("X1"), refined_pdbfile("none");
		Size actual_num_models(0), num_models(0);
		Size const target_num_models( option[ my_options::num_models ].user() ? option[ my_options::num_models ]() : 0 );
		Size const composition( option[ my_options::composition ].user() ? option[ my_options::composition ]() : 0 );
		Size const ntrim( option[ my_options::ntrims ].user() ? random_element( option[ my_options::ntrims ]() ) : 0 );
		Size const ctrim( option[ my_options::ctrims ].user() ? random_element( option[ my_options::ctrims ]() ) : 0 );
		start_decoy_timer();
		if ( passed_score_filter ) {
			if ( ctrim ) pose.conformation().delete_residue_range_slow( pose.total_residue()-ctrim+1, pose.total_residue() );
			if ( ntrim ) pose.conformation().delete_residue_range_slow( 1, ntrim );
			pose.dump_pdb( outfilename );
			string desired_space_group("");
			if ( !space_groups.empty() ) desired_space_group = random_element( space_groups );
			run_phaser_and_phenix_refinement_tri( desired_space_group, outfilename,
				option[ my_options::llg_threshold_for_refinement ],
				actual_num_models, num_models, llg, tfz, rwork, rfree,
				space_group, refined_pdbfile, target_num_models, composition );
			run_command("rm "+outfilename+"*mtz");
			run_command("gzip "+outfilename+"*");
		}
		Real const phaser_simtime( check_decoy_timer() );


		cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" passed_score_filter: " << passed_score_filter <<
			" llg: " << F(9,3,llg) <<
			" tfz: " << F(9,3,tfz) <<
			" rwork: " << F(9,3,rwork) <<
			" rfree: " << F(9,3,rfree) <<
			" actual_num_models: " << actual_num_models <<
			" num_models: " << num_models <<
			" ntrim: " << ntrim <<
			" ctrim: " << ctrim <<
			" space_group: " << space_group <<
			" phaser_simtime: " << F(9,3,phaser_simtime) <<
			" refined_pdbfile: " << refined_pdbfile <<
			endl;

		fflush( stdout );
		check_simtime();

	}

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
rephase_tri_test()
{
	runtime_assert( !option[ basic::options::OptionKeys::out::file::output_virtual ] ); // no NV in pdb files
	confirm_phenix_version();
	runtime_assert( option[ my_options::llg_threshold_for_refinement ].user() );

	Size const ntrim( option[ my_options::ntrim ].user() ? option[ my_options::ntrim]() : 0 );

	strings space_groups( option[ my_options::space_groups ] );
	strings files( start_files() );
	numeric::random::random_permutation( space_groups, numeric::random::rg() );

	foreach_ ( string const space_group, space_groups ) {

		string const simfile( shared_output_tag()+"_rephase_tri_"+space_group+".work" );
		if ( simfile_is_done(simfile) ) {
			check_if_job_is_done();
			continue;
		}

		numeric::random::random_permutation( files, numeric::random::rg() );

		foreach_ ( string filename, files  ) {
			Size const nstruct(1);
			Size const n( get_next_decoy_number_and_reserve_if_not_done( filebase(filename), nstruct, simfile ) );
			if ( n > nstruct ) {
				if ( simfile_is_done(simfile) ) break;
				continue;
			}

			if ( !utility::file::file_exists( filename ) && utility::file::file_exists( filename+".gz") ) filename += ".gz";

			if ( !utility::file::file_exists( filename ) ) {
				cerr << "missing " << filename << endl;
				continue;
			}
			string fb( filebase( filename ) );
			run_command( "cp "+filename+" ./");
			if ( fb.substr(fb.size()-3)==".gz" ) {
				run_command( "gunzip "+fb );
				fb = fb.substr(0,fb.size()-3);
			}
			if ( ntrim ) {
				string const fbnew( fb+"_trim"+string_of(ntrim)+".pdb" );
				run_command( "python /home/pbradley/python/trim_pdb.py "+fb+" "+string_of(ntrim)+" > "+fbnew );
				run_command( "rm "+fb);
				fb = fbnew;
			}
			Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
			string space_group_from_phaser("X1"), refined_pdbfile("none");
			Size actual_num_models(0), num_models(0);
			Size const target_num_models( option[ my_options::num_models ].user() ? option[ my_options::num_models ]() : 0 );
			Size const composition( option[ my_options::composition ].user() ? option[ my_options::composition ]() : 0 );

			run_phaser_and_phenix_refinement_tri( space_group, fb, option[ my_options::llg_threshold_for_refinement ],
				actual_num_models, num_models, llg, tfz, rwork, rfree,
				space_group_from_phaser, refined_pdbfile, target_num_models, composition );


			run_command( "rm "+fb );
			run_command( "rm "+fb+"*mtz" );
			run_command( "gzip "+fb+"*" );
			cout << "final_scores " << F(9,3,-1*llg) << ' ' << fb << "_MR.1.pdb" <<
				" llg: " << F(9,3,llg) <<
				" tfz: " << F(9,3,tfz) <<
				" rwork: " << F(9,3,rwork) <<
				" rfree: " << F(9,3,rfree) <<
				" actual_num_models: " << actual_num_models <<
				" num_models: " << num_models <<
				" ntrim: " << ntrim <<
				" target_space_group: " << space_group <<
				" space_group: " << space_group_from_phaser <<
				" refined_pdbfile: " << refined_pdbfile <<
				endl;

			fflush( stdout );
			check_simtime();
		} // files
		signal_that_simfile_is_done( simfile );
		check_if_job_is_done();
	} // space_groups

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
extend_tri_test()
{
	Size const redundancy( 2 );

	runtime_assert( !option[ basic::options::OptionKeys::out::file::output_virtual ] ); // no NV in pdb files
	confirm_phenix_version();
	runtime_assert( option[ my_options::llg_threshold_for_refinement ].user() );

	// Size const target_num_models( option[ my_options::num_models ]() ); // for the single chain guy
	Size const composition( option[ my_options::composition ] );
	// Size const ntrim( option[ my_options::ntrim ].user() ? option[ my_options::ntrim]() : 0 );

	// strings space_groups( option[ my_options::space_groups ] );
	strings files( start_files() );
	// numeric::random::random_permutation( space_groups, numeric::random::rg() );

	Size const ndel( option[ my_options::num_chains_to_delete ] );

	{

		numeric::random::random_permutation( files, numeric::random::rg() );

		foreach_ ( string filename, files  ) {

			string const simfile( shared_output_tag()+"_rephase_tri_"+filebase(filename)+".work" );
			if ( simfile_is_done(simfile) ) {
				check_if_job_is_done();
				continue;
			}

			while ( true ) {
				Size const n( get_next_decoy_number_and_reserve_if_not_done( "tmp", nstruct(), simfile ) );
				if ( n > nstruct() ) {
					check_if_job_is_done();
					break;
				}

				if ( !utility::file::file_exists( filename ) && utility::file::file_exists( filename+".gz") ) filename += ".gz";

				if ( !utility::file::file_exists( filename ) ) {
					cerr << "missing " << filename << endl;
					break;
				}

				Pose pose;
				pose_from_pdb( pose, filename );

				// figure out the space group
				string space_group_from_pdbfile;
				{
					utility::io::izstream data( filename );
					string line;
					while ( getline( data,line ) ) {
						strings const l( split_to_vector1( line ) );
						if ( l.size()>6 && l[1] == "CRYST1" ) {
							for ( Size i=8; i<l.size(); ++i ) {
								space_group_from_pdbfile += l[i];
							}
							break;
						}
					}
					TR.Trace << "space_group_from_pdbfile: " << space_group_from_pdbfile << endl;
				}
				add_termini_at_protein_chainbreaks( pose );

				Size const nchains( num_chains( pose ) );

				vector1< Sizes > all_chaindels;
				if ( ndel == 0 ) all_chaindels.push_back( Sizes() );
				else if ( ndel == 1 ) {
					for ( Size i=1; i<= nchains; ++i ) all_chaindels.push_back( make_vector1(i) );

				} else if ( ndel == 2 ) {
					for ( Size i=1; i< nchains; ++i ) {
						for ( Size j=i+1; j<=nchains; ++j ) {
							all_chaindels.push_back( make_vector1(i,j) );
						}
					}

				} else if ( ndel == 3 ) {
					for ( Size i=1; i<= nchains; ++i ) {
						for ( Size j=i+1; j<=nchains; ++j ) {
							for ( Size k=j+1; k<=nchains; ++k ) {
								all_chaindels.push_back( make_vector1(i,j,k) );
							}
						}
					}
				} else {
					utility_exit_with_message("too many ndel");
				}


				if ( n> redundancy * all_chaindels.size() ) break;

				if ( !all_chaindels.empty() ) std::reverse( all_chaindels.begin(), all_chaindels.end() ); // prefer to delete later chains

				Sizes chaindels( all_chaindels[1 + (n-1)%(all_chaindels.size()) ] );
				string chaindels_tag;
				Pose small_pose( pose ), big_pose( pose );
				small_pose.conformation().delete_residue_range_slow( chain_begin(2,pose), pose.total_residue() );

				if ( chaindels.empty() ) {
					chaindels_tag = "-";
				} else {
					std::sort( chaindels.begin(), chaindels.end() );
					std::reverse( chaindels.begin(), chaindels.end() ); // now in decreasing order

					foreach_ ( Size ch, chaindels ) chaindels_tag += string_of(ch)+".";
					chaindels_tag.erase( chaindels_tag.size()-1 );
					foreach_ ( Size ch, chaindels ) {
						big_pose.conformation().delete_residue_range_slow( chain_begin( ch, big_pose ), chain_end( ch, big_pose ) );
					}
				}

				randomly_shift_and_tilt_pose( small_pose );
				randomly_shift_and_tilt_pose( big_pose );

				string const prefix( filebase( filename )+"_dc_"+chaindels_tag+"_N"+lead_zero_string_of(n,4) ),
					bigfile( prefix+"_big.pdb" ), smallfile( prefix+"_small.pdb" );

				big_pose.dump_pdb( bigfile );
				small_pose.dump_pdb( smallfile );

				Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
				string space_group_from_phaser("X1"), refined_pdbfile("none");
				Size actual_num_models(0), num_models(0);

				Size target_num_models( nchains < composition ? ndel+1 : ndel );
				if ( option[my_options::num_models].user() ) target_num_models = option[ my_options::num_models]();
				if ( option[my_options::final_num_models].user() ) target_num_models = option[ my_options::final_num_models] - ( nchains-ndel );

				TR.Trace << "chaindels_tag: " << chaindels_tag << ' ' << n << ' ' << nchains << ' ' << chain_end(1,pose) <<
					' ' << pose.total_residue() << ' ' << target_num_models << ' ' << composition << ' ' << ndel << endl;

				run_phaser_and_phenix_refinement_tri( space_group_from_pdbfile, smallfile, option[ my_options::llg_threshold_for_refinement ],
					actual_num_models, num_models, llg, tfz, rwork, rfree,
					space_group_from_phaser, refined_pdbfile, target_num_models, composition,
					bigfile );


				run_command( "rm "+bigfile );
				run_command( "rm "+smallfile );
				run_command( "rm "+prefix+"*mtz" );
				run_command( "gzip "+prefix+"*" );

				cout << "final_scores " << F(9,3,-1*llg) << ' ' << smallfile << "_MR.1.pdb" <<
					" chaindels: " << chaindels_tag <<
					" llg: " << F(9,3,llg) <<
					" tfz: " << F(9,3,tfz) <<
					" rwork: " << F(9,3,rwork) <<
					" rfree: " << F(9,3,rfree) <<
					" actual_num_models: " << actual_num_models <<
					" num_models: " << num_models <<
					" start_num_models: " << nchains  <<
					" final_num_models: " << ( nchains - ndel + actual_num_models - 1 ) << // actual_num_models includes bigfile and one or more smallfile's
					" target_space_group: " << space_group_from_pdbfile <<
					" space_group: " << space_group_from_phaser <<
					" refined_pdbfile: " << refined_pdbfile <<
					endl;

				fflush( stdout );
				check_simtime();
			} // nstruct loop
			signal_that_simfile_is_done( simfile );
			check_if_job_is_done();
		} // files
		check_if_job_is_done();
	} // was space_groups

	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
refold_logfile_asym_test()
{
	Real const relax_fraction( option[ my_options::relax_fraction ] ),
		star_fraction( option[ my_options::star_fraction ] ), centroid_fraction( 1.0 - relax_fraction - star_fraction );

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	//bool const use_asymmetric_centroid_scoring( option[ my_options::use_asymmetric_centroid_scoring ] );
	//Size const extra_repeats_for_centroid_folding( use_asymmetric_centroid_scoring ? 0 : 2 );
	//Size const base_repeat( 3 );

	/// read all the lines in the file, noting which are the "good" ones
	strings all_lines;
	Sizes good_lines;
	bool pdb_list_io( option[ my_options::pdb_list_io ] );
	if ( pdb_list_io || start_files().size()> 1 ) { // we're reading a list file
		pdb_list_io = true;
		all_lines = start_files();
		for ( Size i=1; i<= all_lines.size(); ++i ) good_lines.push_back( i );

	} else {
		string const logfile( start_file() );
		string line;
		ifstream data( logfile.c_str() );
		runtime_assert( data.good() );
		while ( getline( data, line ) ) {
			all_lines.push_back( line );
			strings const l( split_to_vector1( line ) );
			bool passed_score_filter( false );
			for ( Size i=1; i< l.size(); ++i ) if ( l[i] == "passed_score_filter:" ) passed_score_filter = ( l[i+1] == "1" );
			if ( passed_score_filter ) good_lines.push_back( all_lines.size () );
		}
		data.close();
	}


	numeric::random::random_permutation( good_lines, numeric::random::rg() );

	string const outdir( shared_output_dir_create_if_necessary() );

	for ( Size gi=1; gi<= good_lines.size(); ++gi ) {
		Size const line_index( good_lines[gi] );

		string const simfile( outdir + "refold_line"+string_of( line_index)+".work" );
		if ( simfile_is_done( simfile ) ) continue;
		string const worktag( "refold" );

		// Size const first_n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );

		// if ( first_n > nstruct() ) continue;


		Size repeatlen(0), nrepeat(0);
		string repeatseq("none");

		string filename;
		if ( pdb_list_io ) {
			filename = all_lines[ line_index ];

		} else {

			//string const & line( all_lines[ line_index ] );
			strings const l( split_to_vector1( all_lines[ line_index ] ) );

			string const dirtag( l[1] );
			filename = dirtag.substr( 0, dirtag.find_last_of( "/" ) ) +"/" + l[3];

			// Size repeatlen(0), nrepeat(0);
			// string repeatseq("none");

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "repeatlen:" ) repeatlen = int_of( l[i+1] );
				else if ( l[i] == "nrepeat:" ) nrepeat = int_of( l[i+1] );
				else if ( l[i] == "repeatseq:" ) repeatseq = l[i+1];
			}
		}

		if ( !utility::file::file_exists( filename ) &&
				utility::file::file_exists( filename+".gz") ) filename += ".gz";

		if ( !utility::file::file_exists( filename ) &&
				utility::file::file_exists( string("../.") + filename ) ) filename = string("../.")+filename;

		if ( !utility::file::file_exists( filename ) &&
				utility::file::file_exists( string("../.") + filename+".gz" ) ) filename = string("../.")+filename+".gz";

		if ( !utility::file::file_exists( filename ) ) {
			cout <<"missing " << filename << endl;
			cerr <<"missing " << filename << endl;
			continue;
		}

		Pose pdb_pose;
		pose_from_pdb( pdb_pose, filename );

		if ( star_fraction > 1e-6 ) {
			set_ss_from_dssp( filename, pdb_pose );
		}

		if ( option[ my_options::fold_first_chain ] ) {
			if ( num_chains( pdb_pose ) >1  ) {
				pdb_pose.conformation().delete_residue_range_slow( chain_begin( 2, pdb_pose ), pdb_pose.total_residue() );
			}
		}

		for ( Size i=1; i<= pdb_pose.total_residue(); ++i ) {
			TR.Trace << "pdb_pose: " << i << ' ' << pdb_pose.residue(i).name() << ' ' << pdb_pose.residue(i).chain() << endl;
		}

		while ( !pdb_pose.residue( pdb_pose.total_residue() ).is_protein() ) {
			pdb_pose.conformation().delete_residue_slow( pdb_pose.total_residue() );
		}
		// if ( pdb_pose.residue( pdb_pose.total_residue() ).name()=="VRT" ) {
		// 	pdb_pose.conformation().delete_residue_slow( pdb_pose.total_residue() );
		// }

		for ( Size i=1; i<= pdb_pose.total_residue(); ++i ) { runtime_assert( pdb_pose.residue(i).is_protein() ); }

		if ( option[ my_options::nrepeats_to_fold ].user() ) {

			if ( !nrepeat ) { // try to figure it out
				string const seq( pdb_pose.sequence() );
				repeatlen = deduce_repeatlen( seq );
				if ( !repeatlen ) utility_exit_with_message("repeatlen problem "+filename );
				runtime_assert( seq.size()%repeatlen == 0 );
				nrepeat = ( seq.size() / repeatlen );
				repeatseq = ( seq.substr(0,repeatlen) );
			}

			Size const nrepeats_to_fold( option[ my_options::nrepeats_to_fold ] );
			runtime_assert( nrepeat >= nrepeats_to_fold );
			runtime_assert( pdb_pose.total_residue() == nrepeat * repeatlen );
			Size const nrepeats_to_trim( nrepeat - nrepeats_to_fold );
			if ( nrepeats_to_trim ) {
				pdb_pose.conformation().delete_residue_range_slow( pdb_pose.total_residue() - nrepeats_to_trim*repeatlen + 1,
					pdb_pose.total_residue() );
			}
			nrepeat = nrepeats_to_fold;
			runtime_assert( pdb_pose.total_residue() == nrepeat * repeatlen );
		}

		if ( option[ my_options::nres_to_fold ].user() ) {
			Size const nres( option[ my_options::nres_to_fold ]  );
			runtime_assert( pdb_pose.total_residue() >= nres );
			if ( pdb_pose.total_residue() > nres ) {
				pdb_pose.conformation().delete_residue_range_slow( nres+1, pdb_pose.total_residue() );
				add_upper_terminus_type_to_pose_residue( pdb_pose, nres );
			}
		}


		if ( nrepeat ) {
			runtime_assert( num_chains( pdb_pose ) == 1 );
			runtime_assert( pdb_pose.total_residue() == nrepeat * repeatlen );
			runtime_assert( pdb_pose.sequence().substr(0,repeatlen) == repeatseq );
		} else {
			// 12/13/13: now that we are using jumping to build designs, lets not do this:
			// add_termini_at_protein_chainbreaks( pdb_pose, 3.5 ); // bigger threshold, just in case
			// runtime_assert( num_chains( pdb_pose ) == 1 );
			nrepeat = 1;
			repeatlen = pdb_pose.sequence().size();
			repeatseq = pdb_pose.sequence();
		}

		Size const ntrim_sses( option[ my_options::ntrim_sses ].user() ? option[ my_options::ntrim_sses ] : 0 );
		Size const ctrim_sses( option[ my_options::ctrim_sses ].user() ? option[ my_options::ctrim_sses ] : 0 );
		if ( ntrim_sses || ctrim_sses ) {
			Pose pose( pdb_pose );

			set_simple_fold_tree_that_respects_cutpoint_variants( pose );

			//Size repeatlen, nrepeat;
			string repeatbb, repeatss, topology; //, repeatseq;
			Sizes sse_lens;
			strings turns;
			vector1< SizePairs > potential_jump_points;
			vector1< Sizes > potential_cutpoints;
			SegmentTypes sse_types;
			BetaPairingTypes jump_types_single_repeat;
			bool success
				( parse_topology_parameters_from_pose( pose, topology, nrepeat, repeatlen,
				repeatseq, repeatbb, repeatss,
				sse_lens, sse_types, turns,
				potential_jump_points, potential_cutpoints, jump_types_single_repeat,
				filename ) );
			if ( !success ) {
				cout << "refold_logfile_asym_test:: parse_topology_parameters_from_pose FAILED 1 " << filename << endl;
				signal_that_simfile_is_done( simfile );
				continue;
			}

			string const bbtag( get_bbtag_from_topology_sses_and_nrepeat( topology, sse_lens, turns, nrepeat ) );
			TR.Trace << "refold_logfile_asym_test:: bbtag: " << bbtag << ' ' << filename << endl;

			if ( ctrim_sses ) {
				// trim some sses off the c-terminus
				// also remove trailing turn
				Size const num_sses( sse_lens.size() );
				if ( ctrim_sses < num_sses ) {
					Size ctrim(0);
					for ( Size i= num_sses-ctrim_sses+1; i<= num_sses; ++i ) ctrim += sse_lens[i] + get_turnlen( turns[i] );
					TR.Trace << "refold_logfile_asym_test: ctrim " << ctrim_sses << ' ' << ctrim << endl;
					pose.conformation().delete_residue_range_slow( pose.total_residue() - ctrim+1, pose.total_residue() );
				} else success = false; // PROBLEMO
			}

			if ( ntrim_sses ) {
				// trim some sses off the n-terminus
				// also remove trailing turn
				if ( ntrim_sses < sse_lens.size() ) {
					Size ntrim(0);
					for ( Size i= 1; i<= ntrim_sses; ++i ) ntrim += sse_lens[i] + get_turnlen( turns[i] );
					TR.Trace << "refold_logfile_asym_test: ntrim " << ntrim_sses << ' ' << ntrim << endl;
					pose.conformation().delete_residue_range_slow( 1, ntrim );
				} else success = false; // PROBLEMO
			}

			if ( !success ) {
				cout << "refold_logfile_asym_test:: parse_topology_parameters_from_pose FAILED 2 " << filename << endl;
				signal_that_simfile_is_done( simfile );
				continue;
			}
			/// update params
			pdb_pose = pose; // now pdb_pose will have wonky fold tree
			nrepeat = 1;
			repeatlen = pdb_pose.sequence().size();
			repeatseq = pdb_pose.sequence();
		}

		if ( option[ my_options::test_resampling ] ) {
			TR.Trace << "test_resampling: success " << filename << endl;
			continue;
		}

		core::fragment::FragSetOP small_frags(0), large_frags(0);
		if ( !dry_run() && centroid_fraction > 1e-3 ) {
			pick_nnmake_fragments_from_single_sequence( pdb_pose.sequence(), small_frags, large_frags );
		}

		// cheating frags if we're doing star-rebuilding
		devel::blab::classic_frags::FragLibOP cheating_fraglib(0);
		if ( star_fraction > 1e-6 ) {
			set_simple_fold_tree_that_respects_cutpoint_variants( pdb_pose ); // allows torsion calculation properly
			if ( !dry_run() ) cheating_fraglib = setup_cheating_fragments( pdb_pose, Sizes() );
		}

		while ( true ) {
			ScoreFunctionOP fa_scorefxn(0);
			if ( option[ OptionKeys::dna::specificity::score_function ].user() ) {
				fa_scorefxn = get_score_function_from_command_line();
			} else {
				fa_scorefxn = ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag );
			}

			Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( n > nstruct() ) break;

			Pose pose; // new approach (1/3/14)

			pose.constraint_set(0);

			// possibly try something else

			string simtag;
			Real const prob( uniform() );
			if ( prob < relax_fraction ) {
				simtag = "relax";
				pose = pdb_pose; // with foldtree, variants, etc...
				set_simple_fold_tree_that_respects_cutpoint_variants( pose );

			} else if ( prob < relax_fraction + star_fraction ) {
				simtag = "star";
				pose = pdb_pose;

				// pick segments
				vector1< Sizes > protein_segments;
				setup_protein_segments( pose, protein_segments );

				append_virtual_residue( pose );

				devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID ); // was CENTROID_DNA

				// fold tree
				setup_star_fold_tree_cutpoint_variants_and_virtual_residue( false/*adjib*/, protein_segments, pose );

				// constraints
				Real const distol(0);
				add_simple_cartesian_constraints( pdb_pose, pose, distol );

				// frag rebuild
				ScoreFunctionOP cen_scorefxn( get_centroid_score_function_from_command_line() );
				cen_scorefxn->set_weight( coordinate_constraint, 1.0 );
				//cen_scorefxn->set_weight( atom_pair_constraint, 1.0 ); // used for disulfides
				fa_scorefxn->set_weight( coordinate_constraint, 1.0 );
				fa_scorefxn->set_weight( chainbreak, 1.0 ); // centroid chainbreaks handled in star_fragment_rebuild

				if ( !dry_run() ) {
					// do some centroid rebuilding
					bool const use_superimpose_segments( true );
					bools const is_flexible_dna( pose.total_residue(), false );
					bools is_flexible_protein( pose.total_residue(), true ); is_flexible_protein[ pose.total_residue() ] = false;
					Size const big_frag_size(9);
					star_fragment_rebuild( *cen_scorefxn, is_flexible_dna, is_flexible_protein, cheating_fraglib, pose,
						10, 3, 8, use_superimpose_segments, big_frag_size );
				}

				// to fullatom
				devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
				pose.constraint_set(0);

				// fold tree
				setup_star_fold_tree_cutpoint_variants_and_virtual_residue( false/*adjib*/, protein_segments, pose );

				// new constraints
				add_simple_cartesian_constraints( pdb_pose, pose, distol );


			} else { // centroid simulation
				simtag = "cen";
				for ( Size i=1; i<= pdb_pose.total_residue(); ++i ) {
					pose.append_residue_by_bond( *get_vanilla_protein_residue( pdb_pose.residue(i).name1() ), true );
				}

				devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID ); // was CENTROID_DNA

				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					pose.set_phi  ( i, init_phi   );
					pose.set_psi  ( i, init_psi   );
					pose.set_omega( i, init_omega );
					pose.set_secstruct( i, 'L' );
				}
				// not necessary now:
				for ( Size i=1; i<= pose.total_residue(); ++i ) conformation::idealize_position( i, pose.conformation() );

				if ( !dry_run() ) {
					kinematics::MoveMapOP mm( new kinematics::MoveMap() );
					mm->set_bb( true );
					protocols::abinitio::ClassicAbinitio abinitio( small_frags, large_frags, mm );
					abinitio.init( pose );
					abinitio.apply( pose );
				}
				devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
			}

			/// now fastrelax
			{
				protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

				MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
				movemap->set_bb (true);
				movemap->set_chi(true);
				movemap->set_jump(true);
				if ( !pose.residue( pose.total_residue() ).is_protein() ) {
					movemap->set_bb ( pose.total_residue(), false );
					movemap->set_chi( pose.total_residue(), false );
				}
				fastrelax.set_movemap( movemap );
				if ( !dry_run() ) fastrelax.apply( pose );
			}

			// dont include constraints in final score
			pose.constraint_set(0);


			Real const final_score( ( *fa_scorefxn )( pose ) );

			/// calc rmsd to pdb_pose
			Real rmsd( 0.0 );
			{
				/// trim some residues off either end, defaults are 2 and 5
				Size const ntrim( option[ my_options::ntrim_for_rmsd ]), ctrim( option[ my_options::ctrim_for_rmsd ]);
				runtime_assert( pose.total_residue() >= pdb_pose.total_residue() );
				using namespace core::id;
				AtomID_Map< AtomID > atom_map;
				initialize_atomid_map( atom_map, pose, id::GLOBAL_BOGUS_ATOM_ID );
				for ( Size i=ntrim+1; i<= pdb_pose.total_residue()-ctrim; ++i ) {
					atom_map[ AtomID( pose.residue(i).atom_index("CA"), i ) ] =
						AtomID( pdb_pose.residue(i).atom_index("CA"),i);
				}
				rmsd = rmsd_by_mapping( pose, pdb_pose, atom_map );
				superimpose_pose( pose, pdb_pose, atom_map ); // map pose into the same reference frame as pdb_pose
			}

			Real const maxsub( CA_maxsub( pose, pdb_pose, 4.0 ) ); // 4.0 rms threshold, is default

			string const outfilename( output_tag() + "refold_asym_L"+ string_of( line_index )+
				"_" + string_of( nrepeat ) +"_"+ string_of( repeatlen ) +
				"_S" + simtag +
				"_N"+ lead_zero_string_of( n, 4 )+".pdb" );

			bool const passed_score_filter
				( append_score_to_scorefile_and_filter( simtag, final_score, score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );


			//// mods 1/3/14 to compute over entire decoy...
			bools subset( pose.total_residue(), true );
			// Size const base_repeat( (nrepeat-1)/2 + 1 );
			// Size const base_repeat_offset( repeatlen*(base_repeat-1));
			// for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;

			Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa, buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score;
			compute_sasa_scores_for_subset_slow( 10 /*nrep*/, subset, pose, sasapack_score, norme_score, normsasa_score,
				exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa,
				buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score );


			Real unsatdonbb_per35aa, unsatdonsc_per35aa, unsataccbb_per35aa, unsataccsc_per35aa;
			get_buried_unsatisfied_counts_real_slow( subset, pose,
				unsatdonbb_per35aa, unsatdonsc_per35aa,
				unsataccbb_per35aa, unsataccsc_per35aa, 1.0, 10 );

			Real const normalizer( 35.0 / pose.total_residue() );
			unsatdonbb_per35aa *= normalizer;
			unsatdonsc_per35aa *= normalizer;
			unsataccbb_per35aa *= normalizer;
			unsataccsc_per35aa *= normalizer;

			string const buried_unsatisfied_string( get_buried_unsatisfied_string( pose ) );

			ostringstream approx_repeat_params_stats;
			// nrepeat will be 1 unless something like -nrepeats_to_fold was on the cmdline
			if ( nrepeat > 1 ) {
				// compute some approximate repeat params
				Size const ntrim( option[ my_options::ntrim_for_rmsd ]), ctrim( option[ my_options::ctrim_for_rmsd ]);
				Vectors fixcoords, movcoords;
				//Vector com(0,0,0);
				Size const nsup( (nrepeat-1) * repeatlen );
				for ( Size i=1; i<= nsup; ++i ) {
					if ( i<= ntrim || i+ctrim > nsup ) continue;
					Size const movpos( i ), fixpos( repeatlen+i );
					fixcoords.push_back( pose.residue(fixpos).xyz("CA"));
					movcoords.push_back( pose.residue(movpos).xyz("CA"));
				}

				Stub const stub1( movcoords[1], movcoords[2], movcoords[3] );
				superimpose_coords( fixcoords, movcoords );
				Stub const stub2( movcoords[1], movcoords[2], movcoords[3] );

				Real theta;
				Vector axis, center, t;

				get_stub_transform_data( stub1, stub2, center, axis, t, theta );
				//runtime_assert( fabs( axis.length_squared() - 1.0 )<1e-2 ); // confirm that axis is normal vector

				Real const rise( t.dot( axis ) );
				Real const twist( numeric::conversions::degrees( theta ) );
				approx_repeat_params_stats << " rise: " << F(9,3,rise) << " twist: " << F(9,3,twist) << ' ';
			}

			ostringstream out;
			out << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
				" nrepeat: " << nrepeat <<
				" repeatlen: " << repeatlen <<
				" repeatseq: " << repeatseq <<
				" rmsd: " << F(9,3,rmsd) <<
				" start_file: " << filename <<
				" maxsub: " << F(9,3,maxsub) <<
				" passed_score_filter: " << passed_score_filter <<
				" simtag: " << simtag <<
				" sasapack_score: " << F(9,3,sasapack_score) <<
				" norme_score: " << F(9,3,norme_score) <<
				" normsasa_score: " << F(9,3,normsasa_score) <<
				" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
				" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
				" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
				" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
				" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc)<<
				" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc)<<
				" packstat_score: " << F(9,3,packstat_score) <<
				" unsatdonbb_per35aa: " << F(9,3,unsatdonbb_per35aa) <<
				" unsatdonsc_per35aa: " << F(9,3,unsatdonsc_per35aa) <<
				" unsataccbb_per35aa: " << F(9,3,unsataccbb_per35aa) <<
				" unsataccsc_per35aa: " << F(9,3,unsataccsc_per35aa) <<
				" res_score: " << F(9,3,final_score/pose.total_residue()) <<
				' ' << approx_repeat_params_stats.str() <<
				' ' << buried_unsatisfied_string << '\n';

			string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

			if ( passed_score_filter && pdbfilename.size()>0 ) {
				run_command("gzip "+pdbfilename );
			}
			if ( option[ my_options::output_all_ca_coords ] ) { // dump the c-alpha coords
				cout << "CA_COORDS";
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					Vector const & xyz( pose.residue(i).xyz("CA"));
					cout << ' ' << i << ' ' << F(8,3,xyz.x()) << ' ' << F(8,3,xyz.y()) << ' ' << F(8,3,xyz.z());
				}
				cout << endl;
				fflush( stdout );
			}
			check_simtime();
		}
		signal_that_simfile_is_done( simfile );
	}
	signal_that_job_is_done();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
relax_designs_test()
{
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ) );

	strings files( start_files() );
	numeric::random::random_permutation( files, numeric::random::rg() );

	string const STAR("star"), IDL("idl");
	strings const known_protocols( make_vector1( STAR, IDL ) );
	map< string, Size > nstruct_map;
	//bool need_fragments( false );
	{
		strings protocols_nstruct( option[ my_options::protocols ] );
		for ( Size i=1; i<= protocols_nstruct.size(); ++i ) {
			string const p( protocols_nstruct[i] );
			strings const s( split_to_vector1( p, ":" ) );
			pbassert( s.size() == 2 && is_int( s[2] ) );
			nstruct_map[ s[1] ] = int_of( s[2] );
			cout << "nstruct_map: " << s[1] << ' ' << nstruct_map[ s[1] ] << endl;
			runtime_assert( has_element( known_protocols, s[1] ) );
		}
	}
	strings protocols( get_keys( nstruct_map ) );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );

		string const simfile( shared_output_tag() + filebase( filename) +".work" );
		if ( simfile_is_done( simfile ) ) continue;

		numeric::random::random_permutation( protocols, numeric::random::rg() );

		for ( Size ip=1; ip<= protocols.size(); ++ip ) {
			string const worktag( protocols[ip] );
			Size const nstruct( nstruct_map[ worktag ] );

			Size const first_n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct, simfile ) );

			if ( first_n > nstruct ) continue;


			Pose pdb_pose;
			pose_from_pdb( pdb_pose, filename );

			if ( pdb_pose.residue( pdb_pose.total_residue() ).name()=="VRT" ) {
				pdb_pose.conformation().delete_residue_slow( pdb_pose.total_residue() );
			}

			// look for cutpoints
			Real const devtot_threshold( 0.75 );
			Sizes cutpoints;
			for ( Size i=1; i< pdb_pose.total_residue(); ++i ) {
				Residue const & rsd1( pdb_pose.residue(i) ), &rsd2( pdb_pose.residue(i+1) );
				// check geometry at i->i+1 junction
				// C-N distance
				// Ca-C-N angle
				// C-N-CA angle
				Real const target_distance( 1.3285 ), target_angle1( 116.2 ), target_angle2( 121.7 );
				Real const distance( rsd1.xyz("C").distance( rsd2.xyz("N") ) ),
					angle1( numeric::angle_degrees( rsd1.xyz("CA"), rsd1.xyz("C"), rsd2.xyz("N" ) ) ),
					angle2( numeric::angle_degrees( rsd1.xyz("C" ), rsd2.xyz("N"), rsd2.xyz("CA") ) );
				Real const devtot( 50 * fabs( distance - target_distance ) +
					fabs( angle1 - target_angle1 ) + fabs( angle2 - target_angle2 ) );
				TR.Trace << "chainbreak: " << I(4,i) << F(9,3,devtot) << F(9,3,distance) << F(9,3,angle1) << F(9,3,angle2) <<
					' ' << rsd1.has_variant_type( CUTPOINT_LOWER ) <<
					' ' << rsd2.has_variant_type( CUTPOINT_UPPER ) <<
					' ' << filebase( filename ) << endl;
				if ( devtot > devtot_threshold || rsd1.has_variant_type( CUTPOINT_LOWER ) ||
						rsd2.has_variant_type( CUTPOINT_UPPER ) ) {
					cutpoints.push_back( i );
				}
			}

			for ( Size i=1; i<= pdb_pose.total_residue(); ++i ) { runtime_assert( pdb_pose.residue(i).is_protein() ); }

			if ( option[ my_options::test_resampling ] ) continue;

			bool first_time_through( true );

			while ( true ) {
				Size const n( first_time_through ?
					first_n : get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct, simfile ) );
				if ( n > nstruct ) break;
				first_time_through = false;

				Pose pose( pdb_pose );
				if ( worktag == IDL ) {
					// remove cutpoint variants
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						if ( pose.residue(i).has_variant_type( CUTPOINT_LOWER ) ) {
							remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, i );
						}
						if ( pose.residue(i).has_variant_type( CUTPOINT_UPPER ) ) {
							remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, i );
						}
					}

					protocols::idealize::IdealizeMover mover;
					mover.apply( pose );
					fa_scorefxn->set_weight( chainbreak, 0 );
				} else if ( worktag == STAR ) { // need to make sure that we have chainbreak energy on during relax...
					// jumps across all the chainbreaks
					FoldTree f( pose.total_residue() );
					// remove unexpected cutpoints
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						if ( pose.residue(i).has_variant_type( CUTPOINT_LOWER ) && !has_element( cutpoints, i   ) ) {
							remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, i );
						}
						if ( pose.residue(i).has_variant_type( CUTPOINT_UPPER ) && !has_element( cutpoints, i-1 ) ) {
							remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, i );
						}
					}
					for ( Sizes::const_iterator cut= cutpoints.begin(); cut!= cutpoints.end(); ++cut ) {
						f.new_jump( *cut, *cut+1, *cut );
						if ( !pose.residue(*cut).has_variant_type( CUTPOINT_LOWER ) ) {
							add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, *cut );
						}
						if ( !pose.residue(*cut+1).has_variant_type( CUTPOINT_UPPER ) ) {
							add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, *cut+1 );
						}
					}
					f.reorder(1);
					pose.fold_tree(f);
					fa_scorefxn->set_weight( chainbreak, 1 ); // too low?
				} else {
					utility_exit_with_message("unrecognized protocol "+worktag );
				}


				{
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					movemap->set_bb (true);
					movemap->set_chi(true);
					if ( worktag == STAR ) movemap->set_jump( true );
					//     movemap->set_bb ( pose.total_residue(), false );
					//     movemap->set_chi( pose.total_residue(), false );
					fastrelax.set_movemap( movemap );
					if ( !dry_run() ) fastrelax.apply( pose );
				}


				Real const final_score( ( *fa_scorefxn )( pose ) );

				/// calc rmsd to pdb_pose
				Real rmsd( 0.0 );
				{
					runtime_assert( pose.total_residue() == pdb_pose.total_residue() );
					using namespace core::id;
					AtomID_Map< AtomID > atom_map;
					initialize_atomid_map( atom_map, pose, id::GLOBAL_BOGUS_ATOM_ID );
					for ( Size i=1; i<= pdb_pose.total_residue(); ++i ) {
						atom_map[ AtomID( pose.residue(i).atom_index("CA"), i ) ] =
							AtomID( pdb_pose.residue(i).atom_index("CA"),i);
					}
					rmsd = rmsd_by_mapping( pose, pdb_pose, atom_map );
					superimpose_pose( pose, pdb_pose, atom_map ); // map pose into the same reference frame as pdb_pose
				}

				Real const maxsub( CA_maxsub( pose, pdb_pose, 4.0 ) ); // 4.0 rms threshold, is default

				string const outfilename( output_tag() + "relax_"+ filebase( filename )+"_S"+worktag+
					"_N"+ lead_zero_string_of( n, 4 )+".pdb" );

				bool const passed_score_filter
					( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
					score_filter_pass_early, simfile ) );


				//// mods 1/3/14 to compute over entire decoy...
				bools subset( pose.total_residue(), true );
				// Size const base_repeat( (nrepeat-1)/2 + 1 );
				// Size const base_repeat_offset( repeatlen*(base_repeat-1));
				// for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;

				Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
					buried_polar_sasa, buried_nonpolar_sasa, buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
					packstat_score;
				compute_sasa_scores_for_subset_slow( 10, subset, pose, sasapack_score, norme_score, normsasa_score,
					exposed_polar_sasa, exposed_nonpolar_sasa,
					buried_polar_sasa, buried_nonpolar_sasa,
					buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
					packstat_score );


				Real unsatdonbb_per35aa, unsatdonsc_per35aa, unsataccbb_per35aa, unsataccsc_per35aa;
				get_buried_unsatisfied_counts_real_slow( subset, pose,
					unsatdonbb_per35aa, unsatdonsc_per35aa,
					unsataccbb_per35aa, unsataccsc_per35aa, 1.0, 10 );

				Real const normalizer( 35.0 / pose.total_residue() );
				unsatdonbb_per35aa *= normalizer;
				unsatdonsc_per35aa *= normalizer;
				unsataccbb_per35aa *= normalizer;
				unsataccsc_per35aa *= normalizer;

				string const buried_unsatisfied_string( get_buried_unsatisfied_string( pose ) );


				cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
					//" nrepeat: " << nrepeat <<
					// " repeatlen: " << repeatlen <<
					// " repeatseq: " << repeatseq <<
					" rmsd: " << F(9,3,rmsd) <<
					" start_file: " << filename <<
					" maxsub: " << F(9,3,maxsub) <<
					" passed_score_filter: " << passed_score_filter <<
					" sasapack_score: " << F(9,3,sasapack_score) <<
					" norme_score: " << F(9,3,norme_score) <<
					" normsasa_score: " << F(9,3,normsasa_score) <<
					" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
					" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
					" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
					" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
					" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc)<<
					" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc)<<
					" packstat_score: " << F(9,3,packstat_score) <<
					" unsatdonbb_per35aa: " << F(9,3,unsatdonbb_per35aa) <<
					" unsatdonsc_per35aa: " << F(9,3,unsatdonsc_per35aa) <<
					" unsataccbb_per35aa: " << F(9,3,unsataccbb_per35aa) <<
					" unsataccsc_per35aa: " << F(9,3,unsataccsc_per35aa) <<
					" res_score: " << F(9,3,final_score/pose.total_residue()) <<
					' ' << buried_unsatisfied_string <<
					endl;

				if ( option[ my_options::output_all_ca_coords ] ) { // dump the c-alpha coords
					cout << "CA_COORDS";
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						Vector const & xyz( pose.residue(i).xyz("CA"));
						cout << ' ' << i << ' ' << F(8,3,xyz.x()) << ' ' << F(8,3,xyz.y()) << ' ' << F(8,3,xyz.z());
					}
					cout << endl;
				}

				fflush( stdout );

				if ( passed_score_filter ) {
					pose.dump_pdb( outfilename );
				}
				check_simtime();
			} // nstruct
		} // protocols
		signal_that_simfile_is_done( simfile );
	}
	signal_that_job_is_done();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
rsd_sasa_test()
{
	vector1< PoseOP > aa_poses;

	for ( Size i=1; i<= num_canonical_aas; ++i ) {
		AA const aa = AA(i);
		char const name1( oneletter_code_from_aa( aa ) );

		// create a 3 rsd pose
		Pose pose;
		pose.append_residue_by_bond( *get_vanilla_protein_residue( 'G' ), true );
		pose.append_residue_by_bond( *get_vanilla_protein_residue( name1 ), true );
		pose.append_residue_by_bond( *get_vanilla_protein_residue( 'G' ), true );

		aa_poses.push_back( PoseOP( new Pose( pose ) ) );

	}


	/// now read a long list of pdb files, for each non-terminal, good position, insert phi-psi-omegas for +/- 1 rsd
	/// window, and chi values for that rsd, into the aa_pose. Compute sasa (polar, nonpolar, total); write out.
	///

	strings const files( start_files() );

	//ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		set_ss_from_dssp( files[fi], pose );

		bools const goodrsd( identify_residues_with_all_single_occupancy( pose, files[fi] ) );

		add_termini_at_protein_chainbreaks( pose );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i) );
			/// aargh: should have checked for disulfide bonded CYS here
			if ( rsd.is_lower_terminus() || rsd.is_upper_terminus() || !rsd.is_protein() || !goodrsd[i] ) continue;

			if ( rsd.aa() == aa_tyr || rsd.aa() == aa_phe || rsd.aa() == aa_trp ) {
				cout << "AROCHI: " << rsd.name1() << ' ' << pose.secstruct(i) <<
					F(9,3,pose.phi(i)) <<F(9,3,pose.psi(i)) <<F(9,3,pose.omega(i))<<
					F(9,3,rsd.chi(1) ) << F(9,3,rsd.chi(2)) <<
					I(4,i) << ' ' << pose.pdb_info()->chain(i) << I(5,pose.pdb_info()->number(i)) << ' ' << files[fi] << endl;
			}

			//
			Pose & aapose( *aa_poses[ rsd.aa() ] );
			runtime_assert( aapose.total_residue() == 3 );
			runtime_assert( aapose.residue(2).aa() == rsd.aa() );

			for ( Size j=1; j<=3; ++j ) {
				aapose.set_phi  ( j, pose.phi  ( i+j-2 ) );
				aapose.set_psi  ( j, pose.psi  ( i+j-2 ) );
				aapose.set_omega( j, pose.omega( i+j-2 ) );
			}
			for ( Size j=1; j<= rsd.nchi(); ++j ) aapose.set_chi( j, 2, rsd.chi(j) );

			bools subset( 3, false ); subset[2] = true;
			Real polar_sasa, nonpolar_sasa;
			Reals rsd_polar_sasa, rsd_nonpolar_sasa, rsd_polar_sasa_sc, rsd_nonpolar_sasa_sc;
			get_rsd_sasas_by_atomtype( aapose, subset, polar_sasa, nonpolar_sasa, rsd_polar_sasa, rsd_nonpolar_sasa,
				rsd_polar_sasa_sc, rsd_nonpolar_sasa_sc );

			cout << "rsd_sasa: " << rsd.name1() <<
				" polar: " << F(9,3,polar_sasa) <<
				" nonpolar: " << F(9,3,nonpolar_sasa)<<
				" total: " << F(9,3,polar_sasa+nonpolar_sasa) <<
				" polar_sc: " << F(9,3,rsd_polar_sasa_sc[2] ) <<
				" nonpolar_sc: " << F(9,3,rsd_nonpolar_sasa_sc[2] ) <<
				" total_sc: " << F(9,3,rsd_polar_sasa_sc[2] + rsd_nonpolar_sasa_sc[2] ) <<
				I(6,i) << ' ' << files[fi] << endl;
			//    cout << "rsd_sasa: " << rsd.name1() <<
			//     " polar: " << F(9,3,polar_sasa) <<
			//     " nonpolar: " << F(9,3,nonpolar_sasa)<<
			//     I(6,i) << ' ' << files[fi] << endl;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





void
sasapack_test()
{

	strings const files( start_files() );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		bools const goodrsd( identify_residues_with_all_single_occupancy( pose, files[fi] ) );

		add_termini_at_protein_chainbreaks( pose );

		Reals rsd_sasapack_scores, rsd_sasa14_normalized;
		Real sasapack_score, average_normsasa;
		protocols::sasa_scores::compute_sasapack_scores( pose, rsd_sasapack_scores, rsd_sasa14_normalized,
			sasapack_score, average_normsasa );

		Real const big_probe_radius( 1.4 ), small_probe_radius( 0.5 );
		Reals rsd_sasa_big_probe, rsd_sasa_small_probe;
		protocols::sasa_scores::compute_residue_sasas_for_sasa_scores(   big_probe_radius, pose, rsd_sasa_big_probe );
		protocols::sasa_scores::compute_residue_sasas_for_sasa_scores( small_probe_radius, pose, rsd_sasa_small_probe );

		(*fa_scorefxn)( pose );

		cout << "total_sasapack_score: " << F(9,3,sasapack_score) << ' ' << files[fi] << endl;

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i));
			if ( !rsd.is_protein() ) continue;
			char n1( rsd.name1() );
			if ( rsd.aa() == aa_cys && rsd.has_variant_type( chemical::DISULFIDE ) ) n1 = 'c';

			EnergyMap const & emap( pose.energies().residue_total_energies(i) );

			cout << "SASAPACK " << I(4,i) << ' ' << n1 <<
				" termtag: " << rsd.is_lower_terminus() << rsd.is_upper_terminus() << ' '<<
				" good: " << goodrsd[i] <<
				I(5,pose.pdb_info()->number(i)) << A(2,pose.pdb_info()->chain(i)) <<
				" sasapack: " << F(9,3,rsd_sasapack_scores[i] ) <<
				" sasa14: " << F(9,3,rsd_sasa_big_probe[i] ) <<
				" sasa5: " << F(9,3,rsd_sasa_small_probe[i] ) <<
				" total_score: " << F(9,3,emap.dot( fa_scorefxn->weights() ) ) <<
				" weighted_energies: " << emap.weighted_string_of( fa_scorefxn->weights() ) << ' ' <<
				files[fi] << endl;

		}

		//   Reals rsd_norme_scores;
		//   Real const norme_score( compute_pdb_energies_scores( pose, rsd_norme_scores, rsd_sasa14_normalized, files[fi] ) );

		//   cout << "norme_score: " << F(9,3,norme_score) << ' ' << files[fi] << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef NOCOMPILE
void
sasapack_silent_test()
{
	strings silentfiles( start_files() );


	numeric::random::random_permutation( silentfiles, numeric::random::rg() );

	string const simfile( shared_output_tag() +"_silentfiles.work" );

	for ( Size fi=1; fi<= silentfiles.size(); ++fi ) {
		string const silentfile( silentfiles[fi] ), worktag( filebase( silentfile ) );

		Size const nstruct(1);
		Size const first_n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct, simfile ) );

		if ( first_n > nstruct ) continue;

		io::silent::SilentFileData sfd;
		sfd.read_file( silentfile );

		ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
		core::chemical::ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );


		string const outfilename( shared_output_tag()+"_rescore_"+filebase( silentfile )+".scores" );
		ofstream out( outfilename.c_str() );

		for ( io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {

			Pose pose;
			iter->fill_pose( pose, *rsd_set );
			string const filename( iter->decoy_tag() );

			Reals rsd_sasapack_scores, rsd_sasa14_normalized;
			Real sasapack_score, average_normsasa;
			protocols::sasa_scores::compute_sasapack_scores( pose, rsd_sasapack_scores, rsd_sasa14_normalized,
				sasapack_score, average_normsasa );
			Real const score( (*fa_scorefxn)( pose ) );
			Reals rsd_norme_scores;
			Real norme_score;
			protocols::sasa_scores::compute_avge_scores( pose, rsd_norme_scores, rsd_sasa14_normalized,
				norme_score, average_normsasa );

			Real total_sasa14_normalized( 0.0 );
			{
				Size count(0);
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					Residue const & rsd( pose.residue(i) );
					if ( !rsd.is_protein() ) continue;
					if ( rsd.is_lower_terminus() || rsd.is_upper_terminus() ) continue;
					if ( rsd.aa() == aa_cys && rsd.has_variant_type( chemical::DISULFIDE ) ) continue;
					total_sasa14_normalized += rsd_sasa14_normalized[i];
					++count;
				}
				if ( count ) total_sasa14_normalized /= count;
			}

			out << "final_scores " << F(9,3,score) << ' ' << filename <<
				" sasapack_score: " << F(9,3,sasapack_score) <<
				" norme_score: " << F(9,3,norme_score) <<
				" total_sasa14_normalized: " << F(9,3,total_sasa14_normalized);
			if ( pose.data().has( ( pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) ) {
				basic::datacache::CacheableStringFloatMap *data
					= dynamic_cast< basic::datacache::CacheableStringFloatMap* >
					( pose.data().get_raw_ptr( pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA) );

				runtime_assert( data != NULL );

				map< string, float > const & datamap( data->map() );
				out << " datamap_size: " << datamap.size();
				for ( map< string, float >::const_iterator it = datamap.begin(); it != datamap.end(); ++it ) {
					out << ' ' << it->first << ' ' << F(9,3,it->second);
				}
			}
			out << " weighted_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) <<
				'\n';
		}
		out.close();

	}
}
#endif

///////////////////////////////////////////////////////////////////////////////
void
old_sasa_test()
{

	strings const files( start_files() );

	ScoreFunctionOP fa_scorefxn;
	if ( option[ OK::dna::specificity::score_function ].user() ) {
		fa_scorefxn = get_score_function_from_command_line();
	} else {
		fa_scorefxn = ScoreFunctionFactory::create_score_function( "score12prime" );
	}

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		bools const goodrsd( identify_residues_with_all_single_occupancy( pose, files[fi] ) );

		add_termini_at_protein_chainbreaks( pose );

		//
		Real const big_probe_radius( 1.4 ), small_probe_radius( 0.5 );
		Reals rsd_sasa_big_probe, rsd_sasa_small_probe;
		protocols::sasa_scores::compute_residue_sasas_for_sasa_scores(   big_probe_radius, pose, rsd_sasa_big_probe );
		protocols::sasa_scores::compute_residue_sasas_for_sasa_scores( small_probe_radius, pose, rsd_sasa_small_probe );

		(*fa_scorefxn)( pose );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i));
			if ( !rsd.is_protein() ) continue;
			EnergyMap const & emap( pose.energies().residue_total_energies(i) );
			//bool is_lower_terminus( rsd.is_lower_terminus() ), is_upper_terminus( rsd.is_upper_terminus() );

			char n1( rsd.name1() );
			if ( rsd.aa() == aa_cys && rsd.has_variant_type( chemical::DISULFIDE ) ) n1 = 'c';

			cout << "SASA " << I(4,i) << ' ' << n1 << ' ' << rsd.is_lower_terminus() << rsd.is_upper_terminus() << ' '<<
				goodrsd[i] << ' ' <<
				F(9,3,rsd_sasa_big_probe[i] ) << ' ' <<
				F(9,3,rsd_sasa_small_probe[i] ) << ' ' <<
				F(9,3,emap.dot( fa_scorefxn->weights() ) ) <<
				" rsd_fa_atr_unweighted: " << F(9,3,emap[fa_atr]) <<
				' ' << emap.weighted_string_of( fa_scorefxn->weights() ) << endl;
		}


	}



}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
sasa_test()
{
	Size const nrepeat_for_sasa( option[ my_options::nrepeat ] );

	strings const files( start_files() );

	ScoreFunctionOP fa_scorefxn;
	if ( option[ OK::dna::specificity::score_function ].user() ) {
		fa_scorefxn = get_score_function_from_command_line();
	} else {
		fa_scorefxn = ScoreFunctionFactory::create_score_function( "score12prime" );
	}

	foreach_ ( string filename, files ) {

		try {
			Pose pose;
			pose_from_pdb( pose, filename );

			bools const goodrsd( identify_residues_with_all_single_occupancy( pose, filename ) );

			add_termini_at_protein_chainbreaks( pose );

			//
			Real const big_probe_radius( 1.4 ), small_probe_radius( 0.5 );
			Reals rsd_sasa_big_probe, rsd_sasa_big_probe_sdev, rsd_sasa_small_probe, rsd_sasa_small_probe_sdev;
			compute_residue_sasas_for_sasa_scores_slow(   big_probe_radius, pose,
				rsd_sasa_big_probe, rsd_sasa_big_probe_sdev, nrepeat_for_sasa );
			compute_residue_sasas_for_sasa_scores_slow( small_probe_radius, pose,
				rsd_sasa_small_probe, rsd_sasa_small_probe_sdev, nrepeat_for_sasa );

			(*fa_scorefxn)( pose );

			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const & rsd( pose.residue(i));
				if ( !rsd.is_protein() ) continue;

				EnergyMap const & emap( pose.energies().residue_total_energies(i) );
				//bool is_lower_terminus( rsd.is_lower_terminus() ), is_upper_terminus( rsd.is_upper_terminus() );

				char n1( rsd.name1() );
				if ( rsd.aa() == aa_cys && rsd.has_variant_type( chemical::DISULFIDE ) ) n1 = 'c';

				cout << "SASA " << I(4,i) << ' ' << n1 << ' ' <<
					" is_term: " << rsd.is_lower_terminus() << rsd.is_upper_terminus() <<
					" good: " << goodrsd[i] <<
					" sasa_big_probe: " << F(9,3,rsd_sasa_big_probe[i]) <<
					" sasa_big_probe_sdev: " << F(9,3,rsd_sasa_big_probe_sdev[i]) <<
					" sasa_small_probe: " << F(9,3,rsd_sasa_small_probe[i]) <<
					" sasa_small_probe_sdev: " << F(9,3,rsd_sasa_small_probe_sdev[i]) <<
					" rsdE: " << F(9,3,emap.dot( fa_scorefxn->weights() ) ) <<
					" rsd_fa_atr_unweighted: " << F(9,3,emap[fa_atr]) <<
					' ' << emap.weighted_string_of( fa_scorefxn->weights() ) << endl;
			}
		} catch ( utility::excn::Exception const & e ) {
			cout << "IO ERROR " << filename << endl;
			std::cout << "caught exception " << e.msg() << ' ' << filename << std::endl;
			//std::cerr << "caught exception " << e.msg() << ' ' << filename << std::endl;
		}
	}



}



///////////////////////////////////////////////////////////////////////////////
void
change_pose_repeatlen(
	Size const repeatlen,
	Size const nrepeat,
	Pose & pose
)
{
	runtime_assert( num_chains( pose ) == 3 );
	runtime_assert( pose.conformation().chain_end(1) % nrepeat == 0 );
	runtime_assert( nrepeat == ( pose.conformation().chain_end(3) - pose.conformation().chain_begin(3) + 1 ) );
	runtime_assert( nrepeat == ( pose.conformation().chain_end(2) - pose.conformation().chain_begin(2) + 1 ) );

	Size current_repeatlen( chain_end(1,pose) / nrepeat );

	while ( current_repeatlen < repeatlen ) {
		runtime_assert( chain_end(1,pose) == nrepeat * current_repeatlen );
		/// add residues
		remove_upper_terminus_type_from_pose_residue( pose, chain_end(1,pose ) );
		for ( Size i=nrepeat; i>= 1; --i ) {
			Size const repeatend( i*current_repeatlen );
			ResidueOP new_rsd( get_vanilla_protein_residue( 'A' ) );
			pose.append_polymer_residue_after_seqpos( *new_rsd, repeatend, true );
		}
		++current_repeatlen;
	}
	while ( current_repeatlen > repeatlen ) {
		runtime_assert( chain_end(1,pose) == nrepeat * current_repeatlen );
		/// delete residues
		for ( Size i=nrepeat; i>= 1; --i ) {
			Size const repeatend( i*current_repeatlen );
			pose.conformation().delete_residue_slow( repeatend );
		}
		--current_repeatlen;
	}
	//pose.dump_pdb("test_change_repeatlen_"+string_of( repeatlen )+".pdb");


}




///////////////////////////////////////////////////////////////////////////////
void
setup_repeatlen_and_anchorpos_from_commandline( Size & repeatlen, Size & anchorpos )
{

	repeatlen = random_element( option[ my_options::repeatlens ]() );

	// rescale the anchorpos distribution based on repeatlen
	anchorpos = 0;

	if ( option[ my_options::anchor_positions ].user() ) {
		anchorpos = random_element( option[ my_options::anchor_positions ]() );
		runtime_assert( repeatlen == 34 ); // otherwise should allow range to change as a function of repeatlen
	} else {
		Reals const anchor_position_range( option[ my_options::anchor_position_range ] );
		runtime_assert( anchor_position_range.size() == 2 );
		Size const min_anchorpos( int( anchor_position_range[1] * Real( repeatlen ) + 0.5 ) ),
			max_anchorpos( int( anchor_position_range[2] * Real( repeatlen ) + 0.5 ) );
		anchorpos = numeric::random::random_range( min_anchorpos, max_anchorpos );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
setup_repeatss_from_commandline(
	Size const repeatlen,
	Size const anchorpos,
	Size & anchorpos_loop_begin,
	Size & anchorpos_loop_end,
	Size & cutpoint_loop_begin,
	Size & cutpoint_loop_end,
	string & repeatss,
	string & repeatss_no_dashes,
	string & repeatbb, // empty unless we are forcing something
	string & repeatseq_for_fragpicking
)
{
	repeatss.clear(); repeatss.resize( repeatlen, 'H' ); repeatbb.clear(); repeatbb.resize( repeatlen, '-' );
	repeatseq_for_fragpicking.clear(); repeatseq_for_fragpicking.resize( repeatlen, '-' );

	/// pick anchorpos loop boundaries

	if ( option[ my_options::anchorpos_loop_bbs ].user() ) {
		string const bbl( random_element( option[ my_options::anchorpos_loop_bbs ]() ) ), // one position lowercased
			bbu( ObjexxFCL::uppercased( bbl ) ); // all uppercase

		anchorpos_loop_begin=0;
		for ( Size i=0; i<bbl.size(); ++i ) {
			if ( bbl[i] != bbu[i] ) {
				runtime_assert( !anchorpos_loop_begin );
				anchorpos_loop_begin = anchorpos - i;
				anchorpos_loop_end = anchorpos_loop_begin + bbl.size() - 1;
			}
		}
		runtime_assert( anchorpos_loop_begin );
		int const alpha_window_size( 3 );
		for ( int k=-1*alpha_window_size; k<int(bbu.size())+alpha_window_size; ++k ) {
			if ( k>= 0 && k < int(bbu.size()) ) repeatbb[ anchorpos_loop_begin+k-1 ] = bbu[k]; // 0-indexed into repeatbb!
			else repeatbb[ int(anchorpos_loop_begin)+k-1 ] = 'A';
		}
	} else if ( option[ my_options::anchorpos_loop_lengths ].user() ) {
		Size const anchorpos_loop_length( random_element( option[ my_options::anchorpos_loop_lengths ]() ) );
		Size const anchorpos_loop_index( int( uniform() * anchorpos_loop_length ) );
		anchorpos_loop_begin = anchorpos - anchorpos_loop_index;
		anchorpos_loop_end = anchorpos_loop_begin + anchorpos_loop_length-1;
	} else {
		Size tries(0);
		while ( true ) {
			++tries;
			runtime_assert( tries < 1000 );
			anchorpos_loop_begin = random_element( option[ my_options::anchorpos_loop_begins ]() );
			anchorpos_loop_end = random_element( option[ my_options::anchorpos_loop_ends ]() );
			if ( anchorpos_loop_begin <= anchorpos_loop_end &&
					anchorpos_loop_begin <= anchorpos &&
					anchorpos_loop_end   >= anchorpos ) break;
		}
	}

	/// pick cutpoint loop boundaries
	if ( option[ my_options::cutpoint_loop_lengths ].user() ) {
		Size const cutpoint_loop_length( random_element( option[ my_options::cutpoint_loop_lengths ]() ) );
		Size const cutpoint_loop_index( int( uniform() * cutpoint_loop_length ) );
		cutpoint_loop_begin = repeatlen - cutpoint_loop_index;
		cutpoint_loop_end = cutpoint_loop_length - cutpoint_loop_index - 1;
	} else {
		cutpoint_loop_begin = random_element( option[ my_options::cutpoint_loop_begins ]() );
		cutpoint_loop_end = random_element( option[ my_options::cutpoint_loop_ends ]() );
		runtime_assert( cutpoint_loop_end < cutpoint_loop_begin );
		runtime_assert( cutpoint_loop_begin <= repeatlen );
	}

	// set loop regions to L
	for ( Size i=1; i<= repeatlen; ++i ) {
		char ss('H');
		if ( i <= cutpoint_loop_end || i >= cutpoint_loop_begin ) ss = 'L';
		if ( i >= anchorpos_loop_begin && i <= anchorpos_loop_end ) ss = 'L';
		repeatss[i-1] = ss;
	}

	// add some wiggle room?
	repeatss_no_dashes = repeatss;

	if ( !option[ my_options::no_ss_wiggle ] ) { // add wiggle-room in ss?
		for ( Size i=0; i<repeatlen-3; ++i ) {
			if ( repeatss.substr(i,4) == "HHLL" ||
					repeatss.substr(i,4) == "LLHH" ) {
				if ( uniform() < 0.5 ) repeatss[i+1] = '-'; /// mod 8/3/2012: used to make both a '-'
				else                   repeatss[i+2] = '-';
			}
		}
		// e.g. at edges:
		for ( Size i=0; i<repeatlen-2; ++i ) {
			if ( repeatss.substr(i,3) == "HHL" ||
					repeatss.substr(i,3) == "HLL" ||
					repeatss.substr(i,3) == "LHH" ||
					repeatss.substr(i,3) == "LLH" ) {
				repeatss[i+1] = '-';
			}
		}
	}

	// add a break in the longer helix?
	if ( option[ my_options::break_longer_helix ] ||
			( uniform() < option[ my_options::sometimes_break_longer_helix ]() ) ) {
		Sizes break1_poslist, break2_poslist;
		for ( Size i=1; i<= repeatlen; ++i ) {
			if ( i > cutpoint_loop_end   && i < anchorpos_loop_begin ) {
				// if we convert this position to L, how long would the shorter helix be?
				Size const shorter_helix_len( min( i - cutpoint_loop_end-1, anchorpos_loop_begin - i - 1 ) );
				if ( shorter_helix_len >= 3 ) break1_poslist.push_back(i);
			}
			if ( i < cutpoint_loop_begin && i > anchorpos_loop_end   ) {
				Size const shorter_helix_len( min( i - anchorpos_loop_end-1, cutpoint_loop_begin - i - 1 ) );
				if ( shorter_helix_len >= 3 ) break2_poslist.push_back(i);
			}
		}
		if ( break1_poslist.empty() && break2_poslist.empty() ) {
			TR.Trace << "no helix breaks possible!" << endl;
		} else {
			Size breakpos;
			if ( break1_poslist.size() > break2_poslist.size() ||
					( break1_poslist.size() == break2_poslist.size() && uniform()<0.5 ) ) {
				breakpos = random_element( break1_poslist );
			} else breakpos = random_element( break2_poslist );
			TR.Trace << "helixbreak at pos: " << breakpos << endl;
			if ( option[ my_options::break_helix_using_proline_frags ] ) {
				repeatseq_for_fragpicking[ breakpos-1] = 'P';
			} else {
				// no smearing this one out
				repeatss[ breakpos-1] = repeatss_no_dashes[ breakpos-1 ] = 'L';
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
// void
// adjust_ref_weights_from_command_line(
//                    ScoreFunction & scorefxn
//                    )
// {
//  if ( !option[ my_options::adjust_ref_weights ].user() ) return;
//  scoring::methods::EnergyMethodOptions options( scorefxn.energy_method_options() );
//  Reals ref_weights( options.method_weights( ref ) );
//  runtime_assert( ref_weights.size() == 20 );
//  strings const l( option[ my_options::adjust_ref_weights ]() );
//  runtime_assert( l.size()%2 == 0 );
//  Size const nwts( l.size()/2 );
//  for ( Size ii=0; ii< nwts; ++ii ) {
//   AA const aa( aa_from_oneletter_code( l[2*ii+1 ][0] ) );
//   Real const adjustment( float_of( l[2*ii+2] ) );
//   cout << "Updating reference energy for " << aa << " oldval: " << F(9,3,ref_weights[aa]) <<
//    " newval: " << F(9,3,ref_weights[aa] + adjustment) << endl;
//   ref_weights[ aa ] += adjustment;
//  }
//  options.set_method_weights( ref, ref_weights );
//  scorefxn.set_energy_method_options( options );
// }

///////////////////////////////////////////////////////////////////////////////
//
// try symmetric fragment rebuilding from a starting, symmetrical structure, keeping position 13 frozen as the
// anchor
//
// dont forget to do the keep stub in residue thingy to maintain the orientation of the ca-cb vector of 13
//

void
frag_test()
{
	protocols::moves::PyMOLMoverOP pymol(0);
	if ( option[ my_options::show_pymol ] ) {
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol->keep_history(true);
	}

	option[ my_options::sasapack_datafile ]; // confirm that its present
	option[ my_options::pdb_energies_datafile ]; // confirm that its present

	runtime_assert( !option[ my_options::insert_frag_sequence ] ); // need to set up fragseq_poslist

	bool const flex_symdofs( false );
	Size const centroid_cycles( dry_run() ? 50 : n_outer() ); // -n_outer from cmdline
	Size const design_cycles( option[ my_options::design_cycles ] ); // was hard-coded to 5 until 10/21/12
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );
	Size const anchorpos_in( 13 );
	bool const fastrelax_before_design( option[ my_options::fastrelax_before_design ] );

	bool const make_interface_moves( option[ my_options::make_interface_moves ] );
	bool const randomize_aana( option[ my_options::randomize_aana ] );

	bool const single_sim( option[ my_options::single_sim ] );

	AllRTs all_rts;
	if ( make_interface_moves ) {
		load_rts_library( option[ my_options::interface_transforms_file ](), all_rts );
	}


	//
	strings protocols( option[ my_options::protocols] );

	// starting models
	strings files( start_files() );

	numeric::random::random_permutation( files, numeric::random::rg() );

	// scorefxns
	ScoreFunctionOP fa_scorefxn(0),
		cen_scorefxn( ScoreFunctionFactory::create_score_function( "score3" ) );
	{ // since we are using CENTROID_DNA rsd type set
		scoring::methods::EnergyMethodOptions options( cen_scorefxn->energy_method_options() );
		options.atom_vdw_atom_type_set_name( CENTROID_DNA );
		cen_scorefxn->set_energy_method_options( options );
	}

	if ( option[ OK::dna::specificity::score_function ].user() ) {
		fa_scorefxn = get_score_function_from_command_line();

		{ // check
			scoring::methods::EnergyMethodOptions emoptions( fa_scorefxn->energy_method_options() );
			emoptions.exclude_DNA_DNA( true );
			fa_scorefxn->set_energy_method_options( emoptions );
		}
	} else {
		fa_scorefxn = ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ); //"score12prime" );
	}

	// dont want intra-dna energy for different templates affecting the ranking of designs...
	runtime_assert( fa_scorefxn->energy_method_options().exclude_DNA_DNA() == true );

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );

	runtime_assert( pose::symmetry::is_symmetric( *cen_scorefxn ) );

	// initialize the fragment library
	FragLib tal_fraglib;
	if ( option[ my_options::use_tal_fragments ] ) {
		string const fragfile3( option[ my_options::fragfile3 ] ), fragfile9( option[ my_options::fragfile9 ] );
		tal_fraglib.library( 3 ).read_file( fragfile3, 3, 3 );
		tal_fraglib.library( 9 ).read_file( fragfile9, 9, 3 );
	}

	for ( Size fi=1; fi<= files.size(); ++fi ) { // loop over starting model
		string const filename( files[fi] );

		string const simfile( shared_output_tag() + filebase( filename ) +".work" );

		numeric::random::random_permutation( protocols, numeric::random::rg() );

		for ( Size pp=1; pp<= protocols.size(); ++pp ) {
			string const simtag( protocols[pp] );

			Size const first_n( get_next_decoy_number_and_reserve_if_not_done( simtag, nstruct(), simfile ) );

			if ( first_n > nstruct() ) continue;

			/// all parsing of simtag occurs here:
			bool const reversed( simtag.find("rev") != string::npos );
			bool const centroid_rebuild( simtag.find("cen") != string::npos );
			bool const design( simtag.find("design") != string::npos ); // will also be true if fastdesign
			bool const fastdesign( simtag.find("fastdesign") != string::npos );
			/// END


			// load model, create symmetric pdb

			Pose pose;
			pose_from_pdb( pose, filename );
			set_base_partner( pose );

			Size const nrepeat( pose.conformation().chain_end(2) - pose.conformation().chain_begin(2) + 1 );
			//Size const repeatlen( 34 ); // TAL canonical repeatlen

			runtime_assert( nrepeat == 8 );
			Size const base_repeat( 3 );

			{ // since these guys are scaled to the full pose
				cen_scorefxn->set_weight( dna_env, 1.0 / nrepeat );
				cen_scorefxn->set_weight( dna_pair, 1.0 / nrepeat );
				/// little worried about whether the derivs are correct for these guys...
				/// well, the energies are messed up, nevermind the derivs, at least as of the most recent trunk update
				// cen_scorefxn->set_weight( chainbreak, 1.0 / nrepeat );
				// fa_scorefxn->set_weight( chainbreak, 3.0 / nrepeat );
				// now using atom pair constraints instead
				cen_scorefxn->set_weight( atom_pair_constraint, 1.0 );
				fa_scorefxn->set_weight( atom_pair_constraint, 1.0 );
			}

			Pose const start_pose( pose );
			bool first_time_through( true );

			protocols::viewer::add_conformation_viewer( pose.conformation(), "test", 1200,1200 );

			while ( true ) {
				Size const n( first_time_through ?
					first_n : get_next_decoy_number_and_reserve_if_not_done( simtag, nstruct(), simfile ) );
				if ( n > nstruct() ) break;
				first_time_through = false;

				pose = start_pose;

				Size repeatlen, anchorpos;
				setup_repeatlen_and_anchorpos_from_commandline( repeatlen, anchorpos );

				change_pose_repeatlen( repeatlen, nrepeat, pose );
				set_base_partner( pose );

				//runtime_assert( anchorpos == 13 );

				FragLib fraglib;

				Size anchorpos_loop_begin(0), anchorpos_loop_end(0), cutpoint_loop_begin(0), cutpoint_loop_end(0);
				string repeatss_for_layer_design;
				if ( !centroid_rebuild ) {
					/// dont pick frags
				} else if ( option[ my_options::use_tal_fragments ] ) {
					fraglib.copy_fragments( tal_fraglib );
					runtime_assert( !option[ my_options::randomize_centroid_sequence ] );
				} else {

					string repeatss, repeatss_no_dashes, repeatbb, repeatseq_for_fragpicking;
					setup_repeatss_from_commandline( repeatlen, anchorpos, anchorpos_loop_begin, anchorpos_loop_end,
						cutpoint_loop_begin, cutpoint_loop_end, repeatss, repeatss_no_dashes,
						repeatbb, repeatseq_for_fragpicking );
					repeatss_for_layer_design = repeatss_no_dashes;

					string fullss, fullbb, fullseq_for_fragpicking;
					for ( Size i=1; i<= nrepeat; ++i ) {
						fullss += repeatss;
						fullbb += repeatbb;
						fullseq_for_fragpicking += repeatseq_for_fragpicking;
					}
					runtime_assert( fullss.size() == chain_end( 1, pose ) );

					TR.Trace << "repeatss: " << repeatss_no_dashes << ' ' << repeatss << ' ' << repeatbb << ' ' <<
						repeatseq_for_fragpicking <<
						" anchorpos " << anchorpos <<
						" anchorpos_loop: " << anchorpos_loop_begin << ' ' << anchorpos_loop_end <<
						" cutpoint_loop: " << cutpoint_loop_begin << ' ' << cutpoint_loop_end <<
						endl;

					Real const seq_weight( 1 ), ss_weight( 10.0 ), bb_weight( 100.0 );
					Sizes const fragsizes( make_vector1( 3, 9 ) );
					Size const nfrags( 200 );
					kinematics::MoveMap mm;
					for ( Size i=1; i<= chain_end( 1, pose ); ++i ) mm.set_bb( i, true );

					Sizes const homs_to_exclude; // empty list
					devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, pose, mm, fullss, seq_weight, ss_weight,
						fraglib, homs_to_exclude, bb_weight, fullbb,
						fullseq_for_fragpicking );

					/// note that we only do this sequence randomization if we are not using the tal fraglib...
					if ( centroid_rebuild && option[ my_options::randomize_centroid_sequence ] ) {
						///
						for ( Size i=1; i<= repeatlen; ++i ) {
							AA aa;
							if ( i == anchorpos && option[ my_options::freeze_anchorseq ] ) {
								aa = start_pose.residue( anchorpos_in ).aa();
							} else {
								if ( repeatss_no_dashes[i-1] == 'L' ) {
									aa = aa_from_oneletter_code( random_element( make_vector1( 'G','S','N','A' ) ) );
								} else {
									runtime_assert( repeatss_no_dashes[i-1] == 'H' );
									// helix position
									aa = aa_from_oneletter_code( random_element( make_vector1( 'A','L','V','I','D','S','K' ) ) );
								}
							}

							TR.Trace << "randomseq: " << I(3,i) << ' ' << repeatss_no_dashes[i-1] << ' ' << aa << endl;

							for ( Size j=0; j< nrepeat; ++j ) {
								Size const seqpos( j*repeatlen + i );

								make_sequence_change( seqpos, aa, pose );
							}
						}
					} // randomize sequence
				} // pick vall fragments based on ss string


				if ( randomize_aana ) {
					// randomize the dna base and matching
					vector1< AA_Pair > all_aanas( make_vector1( make_pair( aa_asn, na_ade ), // N:a
						make_pair( aa_asn, na_ade ), // N:a
						make_pair( aa_asn, na_ade ), // N:a
						make_pair( aa_gln, na_ade ), // N:g
						make_pair( aa_asn, na_gua ), // N:g
						make_pair( aa_his, na_gua ), // H:g
						make_pair( aa_arg, na_gua ), // R:g
						make_pair( aa_asp, na_cyt ), // D:c
						make_pair( aa_asp, na_cyt ), // D:c
						make_pair( aa_asp, na_cyt ), // D:c
						make_pair( aa_gly, na_thy ),  // G:t
						make_pair( aa_gly, na_thy )  // G:t
						) );
					AA_Pair const aana( random_element( all_aanas ) );
					for ( Size i=1; i<= nrepeat; ++i ) {
						Size const protpos( anchorpos + ( i-1 ) * repeatlen ), dnapos( chain_end(1,pose)+i );
						TR.Trace << "randomize_aana: " << aana.first << ' ' << aana.second <<
							I(4,i) << I(4,protpos ) << I(4,dnapos ) << endl;
						make_sequence_change( protpos, aana.first, pose );
						devel::dna::make_base_pair_mutation( pose, dnapos, aana.second );
					}
				}


				Size jump1, jump2, jump3, jump4;
				setup_symmetric_pose( nrepeat, repeatlen, base_repeat, anchorpos, flex_symdofs, jump1, jump2, jump3, jump4, pose,
					reversed, anchorpos_in );

				Pose const rmsd_pose( pose ); // possible that the pose is already wonky at this point...

				runtime_assert( pose::symmetry::is_symmetric( pose ) );


				// if ( true || dry_run() ) {
				//  //  do some debugging
				//  pose.dump_pdb("beforepack.pdb");

				//  // great! now let's try packing...
				//  pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
				//  task->initialize_from_command_line();
				//  task->or_include_current( true );
				//  if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

				//  bools is_protein( pose.total_residue(), false );
				//  for ( Size i=1; i<= chain_end(1,pose); ++i ) is_protein[i] = true;
				//  task->restrict_to_residues( is_protein );

				//  Size const nloop( option[ my_options::n_refold ] );
				//  protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( fa_scorefxn, task, nloop );
				//  packmover.apply( pose );
				//  pose.dump_pdb("afterpack.pdb");

				//  exit(0);
				// }


				Real centroid_score(0);
				clock_t starttime( clock() );
				if ( centroid_rebuild ) { // centroid fragment rebuild
					conformation::symmetry::SymmetryInfo const symminfo( *pose::symmetry::symmetry_info( pose ) );
					pose::symmetry::make_asymmetric_pose( pose );
					devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID_DNA );
					pose::symmetry::make_symmetric_pose( pose, symminfo );

					{ // add some chainbreak constraints
						pose.constraint_set(0);
						Size const base_protoffset( reversed ? (nrepeat-base_repeat)*repeatlen : (base_repeat-1)*repeatlen );
						Sizes const poslist( make_vector1( base_protoffset, base_protoffset + repeatlen ) );
						add_hacky_chainbreak_constraints_at_positions( true /*centroid*/, poslist, pose );
					}

					// randomize the pose conformation using fragments
					//bool const insert_frag_sequence_during_randomization( false ); // keep initial sequence at start (?)
					Sizes fragseq_poslist_tmp;
					//bool const insert_frag_sequence( option[ my_options::insert_frag_sequence ] );
					//pose.dump_pdb("test1.pdb");
					for ( Size i=1; i<= 100; ++i ) insert_protein_fragment( 9, fragseq_poslist_tmp, fraglib, pose );
					//pose.dump_pdb("test2.pdb");
					if ( make_interface_moves ) insert_interface_fragment_symmetrically( all_rts, pose );
					//pose.dump_pdb("test3.pdb");

					Sizes fragseq_poslist; // need to set this up, use freeze_helixseq, freeze_anchorseq

					cout << "checkpoint1" << endl; cerr << "checkpoint1" << endl; fflush( stdout );
					simple_frag_sim( centroid_cycles, *cen_scorefxn, fraglib, all_rts, make_interface_moves, fragseq_poslist,
						pose);
					cout << "checkpoint2" << endl; cerr << "checkpoint2" << endl; fflush( stdout );

					centroid_score = (*cen_scorefxn )( pose );
					//pose.dump_pdb("test4.pdb");

					pose::symmetry::make_asymmetric_pose( pose );
					devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
					pose::symmetry::make_symmetric_pose( pose, symminfo );
					//devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
				}

				clock_t stoptime( clock() );
				Real const centroid_simtime( ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

				Real const centroid_rmsd( scoring::nbr_atom_rmsd( pose, rmsd_pose ) );

				starttime = clock();

				//pose.dump_pdb("test5.pdb");

				{ // add some FULLATOM chainbreak constraints
					pose.constraint_set(0);
					Size const base_protoffset( reversed ? (nrepeat-base_repeat)*repeatlen : (base_repeat-1)*repeatlen );
					Sizes const poslist( make_vector1( base_protoffset, base_protoffset + repeatlen ) );
					add_hacky_chainbreak_constraints_at_positions( false /*centroid*/, poslist, pose );
				}

				if ( dry_run() ) {
					// do nothing
				} else if ( fastdesign ) { /// try a fast-relax w/ designing
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );
					utility_exit_with_message("no fastdesign");
					//fastrelax.allow_design( true );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					movemap->set_bb(true);
					movemap->set_chi(true);
					movemap->set_jump(true);
					movemap->set_jump(jump1,false); // turn off base-pair optimization
					fastrelax.set_movemap( movemap );

					/// also set a taskfactory
					pack::task::TaskFactoryOP taskfactory( new pack::task::TaskFactory() );
					{ using namespace pack::task::operation;
						taskfactory->push_back( pack::task::operation::TaskOperationCOP( new InitializeFromCommandline() ) );
						taskfactory->push_back( pack::task::operation::TaskOperationCOP( new IncludeCurrent() ) );
						// restrict to protein
						PreventRepackingOP preventop( new PreventRepacking() );
						for ( Size i=chain_begin(2,pose); i<= pose.total_residue(); ++i ) preventop->include_residue(i);
						taskfactory->push_back( preventop );
						// optionally freeze anchorseq
						if ( option[ my_options::freeze_anchorseq ] ) {
							RestrictResidueToRepackingOP restrictop( new RestrictResidueToRepacking() );
							for ( Size i=0; i< nrepeat; ++i ) restrictop->include_residue( i*repeatlen + anchorpos );
							taskfactory->push_back( restrictop );
						}
					}
					fastrelax.set_task_factory( taskfactory );
					fastrelax.apply( pose );

				} else if ( !design ) { /// now try a fast-relax
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					movemap->set_bb(true);
					movemap->set_chi(true);
					movemap->set_jump(true);
					movemap->set_jump(jump1,false); // turn off base-pair optimization
					fastrelax.set_movemap( movemap );
					fastrelax.apply( pose );
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////
				} else { // iterative design and relax ////////////////////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////
					if ( pymol ) {
						pymol->pymol_name("before_design_"+lead_zero_string_of(n,2));
						pymol->apply( pose );
					}
					SequenceConstraints sequence_constraints;
					if ( option[ my_options::layer_design ] ) {
						runtime_assert( repeatss_for_layer_design.size() == repeatlen );
						string fullss( pose.total_residue(), 'L' );
						for ( Size i=0; i<nrepeat; ++i ) {
							for ( Size j=0; j<repeatlen; ++j ) {
								fullss[ i*repeatlen + j ] = repeatss_for_layer_design[ j ];
							}
						}
						add_layer_design_constraints( pose, fullss, sequence_constraints );
					}
					bools is_protein( pose.total_residue(), false );
					for ( Size i=1; i<= pose.total_residue(); ++i ) is_protein[i] = pose.residue(i).is_protein();
					if ( pymol ) pymol->pymol_name("relaxdesign_"+lead_zero_string_of(n,2));
					for ( Size m=1; m<= design_cycles+1; ++m ) {
						bool const skip_relax( m == 1 && centroid_rebuild && option[ my_options::randomize_centroid_sequence ] &&
							!fastrelax_before_design );

						if ( !skip_relax ) { /// try fast relax
							protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

							MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
							movemap->set_jump(true);
							for ( Size i=1; i<= pose.total_residue(); ++i ) {
								movemap->set_bb ( i, is_protein[i] );
								movemap->set_chi( i, is_protein[i] );
							}
							movemap->set_jump(jump1,false); // turn off base-pair optimization
							fastrelax.set_movemap( movemap );
							cout << "checkpoint3" << endl; cerr << "checkpoint3" << endl; fflush( stdout );
							fastrelax.apply( pose );
							cout << "checkpoint4" << endl; cerr << "checkpoint4" << endl; fflush( stdout );
							if ( pymol ) pymol->apply( pose );
						}

						if ( m == design_cycles+1 ) break; // add an extra fastrelax at the end

						{
							// great! now let's try packing...
							pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
							task->initialize_from_command_line();
							task->or_include_current( true );
							if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

							task->restrict_to_residues( is_protein );
							if ( option[ my_options::freeze_anchorseq ] ) {
								for ( Size i=0; i< nrepeat; ++i ) {
									task->nonconst_residue_task( i*repeatlen + anchorpos ).restrict_to_repacking();
								}
							}
							{ /// add sequence forcing here?
								conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
								for ( SequenceConstraints::const_iterator
										it= sequence_constraints.begin(); it!= sequence_constraints.end(); ++it ) {
									Size const base_pos( it->first );
									runtime_assert( symminfo->bb_is_independent( base_pos ) );
									string const name1s( it->second );
									bools allowed_aas( num_canonical_aas, false );
									for ( Size k=0; k<name1s.size(); ++k ) {
										AA const aa( aa_from_oneletter_code( name1s[k] ) );
										allowed_aas[aa] = true;
									}
									if ( !allowed_aas[ pose.residue(base_pos).aa() ] ) {
										// current aa fails!
										cout << "sequence_constraints dont match current aa: " << base_pos << ' ' <<
											pose.residue( base_pos ).aa() << ' ' << name1s << endl;
										make_sequence_change( base_pos, aa_from_oneletter_code( name1s[0] ), pose );
									}
									for ( Size seqpos=1; seqpos<= pose.total_residue(); ++seqpos ) {
										if ( seqpos == base_pos || symminfo->bb_follows( seqpos ) == base_pos ) {
											TR_SYMDES_HH.Trace << "sequence_constraints: " << base_pos << ' ' << seqpos <<' '<< name1s << endl;
											task->nonconst_residue_task( seqpos ).restrict_absent_canonical_aas( allowed_aas );
										}
									}
								}
							}

							if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
								protocols::task_operations::LimitAromaChi2Operation lp_op;
								lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negs)
								lp_op.apply( pose, *task );
							}

							Size const nloop( 25 );
							ScoreFunctionOP design_scorefxn(0);
							if ( ( option[ my_options::use_softrep_for_early_design ] && m<= design_cycles/2 ) ||
									option[ my_options::use_softrep_for_design ] ) {
								design_scorefxn = ScoreFunctionFactory::create_score_function( soft_design_scorefxn_weights_tag );
								adjust_ref_weights_from_command_line( *design_scorefxn );
							} else {
								design_scorefxn = fa_scorefxn;
							}
							protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
							cout << "checkpoint5" << endl; cerr << "checkpoint5" << endl; fflush( stdout );
							packmover.apply( pose );
							if ( pymol ) pymol->apply( pose );
							cout << "checkpoint6" << endl; cerr << "checkpoint6" << endl; fflush( stdout );

							if ( option[ my_options::debug_derivs ] ) {
								TR.Trace << "debug_derivs:" << endl;
								Pose const save_pose( pose );
								MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
								movemap->set_jump(true);
								for ( Size i=1; i<= pose.total_residue(); ++i ) {
									movemap->set_bb ( i, is_protein[i] );
									movemap->set_chi( i, is_protein[i] );
								}
								movemap->set_jump(jump1,false); // turn off base-pair optimization
								protocols::minimization_packing::symmetry::SymMinMoverOP min_mover
									( new protocols::minimization_packing::symmetry::SymMinMover(movemap, fa_scorefxn, "dfpmin_armijo_nonmonotone",
									0.00001, true, true, true ));
								min_mover->apply( pose );
								pose = save_pose;
							}



						}
					} // cycles
				}
				stoptime = clock();

				Real const relax_simtime( ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

				Real const relax_rmsd( nbr_atom_rmsd( pose, rmsd_pose ) );


				/// compute some more stats
				///
				// compute sasapack scores
				//     Real sasapack_score_single_repeat(0.0), norme_score_single_repeat( 0.0 ), sasa14_normalized_single_repeat(0.0);
				Size base_repeat_begin(0), base_repeat_end(0);
				bools subset( pose.total_residue(), false );
				{
					conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
					for ( Size i=1; i<= chain_end(1,pose); ++i ) {
						if ( symminfo->chi_follows(i) == 0 ) {
							//TR.Trace << "repeatpos_sasapack: " << I(4,i) << F(9,3,rsd_sasapack_scores[i]) << endl;
							if ( !base_repeat_begin ) base_repeat_begin = i;
							subset[i] = true;
							base_repeat_end = i;
						}
					}
				}
				runtime_assert( base_repeat_begin && base_repeat_end );
				runtime_assert( base_repeat_end - base_repeat_begin+1 == repeatlen );
				//     Size const base_repeat_offset( repeatlen*(base_repeat-1));
				//     for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;

				Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
					buried_polar_sasa, buried_nonpolar_sasa, buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
					packstat_score;
				compute_sasa_scores_for_subset_slow( 5, subset, pose, sasapack_score, norme_score, normsasa_score,
					exposed_polar_sasa, exposed_nonpolar_sasa, buried_polar_sasa,
					buried_nonpolar_sasa,
					buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
					packstat_score );

				string const buried_unsatisfied_string( get_buried_unsatisfied_string( pose ) );
				Real n_buried_unsatisfied_donors_bb_per_repeat, n_buried_unsatisfied_acceptors_bb_per_repeat,
					n_buried_unsatisfied_donors_sc_per_repeat, n_buried_unsatisfied_acceptors_sc_per_repeat;
				{ // compute unsats over full pose
					Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
						n_buried_unsatisfied_acceptors_sc;
					bools subset_protein( pose.total_residue(), false );
					for ( Size i=1; i<= repeatlen*nrepeat; ++i ) subset_protein[i] = true;
					get_buried_unsatisfied_counts_real_slow( subset_protein, pose,
						n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc,
						n_buried_unsatisfied_acceptors_bb, n_buried_unsatisfied_acceptors_sc);

					n_buried_unsatisfied_donors_bb_per_repeat = Real( n_buried_unsatisfied_donors_bb )/nrepeat;
					n_buried_unsatisfied_donors_sc_per_repeat = Real( n_buried_unsatisfied_donors_sc )/nrepeat;
					n_buried_unsatisfied_acceptors_bb_per_repeat = Real( n_buried_unsatisfied_acceptors_bb )/nrepeat;
					n_buried_unsatisfied_acceptors_sc_per_repeat = Real( n_buried_unsatisfied_acceptors_sc )/nrepeat;
				}

				///
				Real const final_score( (*fa_scorefxn)(pose ) );

				Real unbound_score( 0 );
				{ // try constructing a pose in which the protein and dna are separated
					Pose unboundpose( pose );
					kinematics::Stub const jump3_stub( pose.conformation().upstream_jump_stub( jump3 ) ),
						jump2_stub( pose.conformation().upstream_jump_stub( jump2 ) );

					kinematics::Jump j( unboundpose.jump( jump2 ) );
					j.translation_along_axis( jump2_stub, jump3_stub.M.col_z(), 100.0 );
					//unboundpose.dump_pdb("before_ztrans.pdb");
					unboundpose.set_jump( jump2, j );
					//unboundpose.dump_pdb("after_ztrans.pdb");
					unbound_score = (*fa_scorefxn)( unboundpose );
				}


				/// score filter
				string const scorefiltertag( simtag + "_" + string_of( repeatlen ) ); // energies scale with repeatlen
				bool const passed_score_filter
					( append_score_to_scorefile_and_filter( scorefiltertag, final_score, score_filter_acceptance_rate,
					score_filter_pass_early, simfile ) );

				string const outfilename( output_tag() + "symfragrebuild_"+ filebase( filename )+
					"_S"+simtag+
					"_N"+lead_zero_string_of(n,4)+
					".pdb" );

				string const repeatseq( pose.sequence().substr(0,repeatlen) );
				string const repeatss( pose.secstruct().substr(0,repeatlen) );
				char const anchorseq( pose.residue( anchorpos ).name1() );

				EnergyMap interface_emap;
				{ // get prot-dna energies
					// the neighbor/energy links
					EnergyGraph const & energy_graph( pose.energies().energy_graph() );

					for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
						conformation::Residue const & resl( pose.residue( i ) );
						for ( utility::graph::Graph::EdgeListConstIter
								iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
								irue = energy_graph.get_node(i)->const_upper_edge_list_end();
								iru != irue; ++iru ) {
							EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
							Size const j( edge->get_second_node_ind() );
							conformation::Residue const & resu( pose.residue( j ) );

							if ( ( resl.is_protein() && resu.is_DNA() ) || ( resu.is_protein() && resl.is_DNA() ) ) {
								edge->add_to_energy_map( interface_emap );
							}
						}
					}
				}


				// abego string, use internal repeat to get well-defined torsions at cutpoints
				string const repeatbb( torsion2big_bin_string( base_repeat_begin, base_repeat_end, pose ) );
				runtime_assert( repeatbb.size() == repeatlen );


				Real sspred_match( 0.0 );
				string repeatsspred( repeatss );
				string fullsspred;
				string const sequence( pose.sequence().substr(0,chain_end(1,pose) ) );
				{ // run psipred to re-predict the secondary structure
					vector1< Reals > pred_HEL;
					string const secstruct( pose.secstruct().substr(0,chain_end(1,pose) ) );
					string sspred;

					run_psipred( sequence, sspred, pred_HEL );

					if ( sspred.size() >= 5*repeatlen ) {
						for ( Size i=1; i<= sspred.size(); ++i ) {
							char const ss( secstruct[i-1] );
							Real const prediction( ss == 'H' ? pred_HEL[i][1] : ( ss == 'E' ? pred_HEL[i][2] : pred_HEL[i][3] ) );
							sspred_match += prediction;
						}
						sspred_match /= sspred.size();
						repeatsspred = sspred.substr( 2* repeatlen, repeatlen );
						fullsspred = sspred;
					}
				}

				/// optionally do some de novo refolding
				ostringstream refolding_status;

				/// fragment quality:
				devel::blab::classic_frags::FragLib fraglib_refold;

				{
					pick_nnmake_fragments_from_single_sequence( sequence, fraglib_refold );
					vector1<Reals > frag_rmsds;
					Reals bb1_recovery, bb3_recovery, sub1_fraction;
					string fragss, fragbb;
					get_fragment_quality( nrepeat, repeatlen, fraglib_refold, pose, fragss, fragbb,
						bb1_recovery, bb3_recovery, sub1_fraction, frag_rmsds );

					runtime_assert( sub1_fraction.size() == repeatlen );
					Real min_sub1_fraction(1), avg_sub1_fraction(0);
					for ( Size i=1; i<= repeatlen; ++i ) {
						min_sub1_fraction = min( min_sub1_fraction, sub1_fraction[i] );
						avg_sub1_fraction += sub1_fraction[i];
					}
					avg_sub1_fraction /= repeatlen;

					refolding_status <<
						" min_sub1_fraction: " << F(9,6,min_sub1_fraction) <<
						" avg_sub1_fraction: " << F(9,3,avg_sub1_fraction);

				}


				Size const n_refold( option[ my_options::n_refold ] ), n_refold_sspred( option[ my_options::n_refold_sspred ] );
				if ( ( n_refold || n_refold_sspred ) && ( passed_score_filter || dry_run() ) ) {
					starttime = clock();
					for ( Size RR=2; RR<= 2; ++RR ) { // dont do the plain old refolding any more
						string const refoldtag( RR == 1 ? "refold" : "refold_sspred" );
						Size const n_refold_thistime( RR==1 ? n_refold : n_refold_sspred );
						refolding_status << " n_" << refoldtag << ": " << I(3,n_refold_thistime );
						Real min_rmsd( 1e6 );
						Size n_under_4(0);
						for ( Size rr=1; rr<= n_refold_thistime; ++rr ) {
							Pose protpose;

							bool const abrelax( false ), use_asymmetric_centroid_scoring( true );
							refold_repeat_pose( repeatseq, nrepeat, base_repeat, fraglib_refold, protpose,
								abrelax, 0, 0, use_asymmetric_centroid_scoring );


							// string const repss( RR == 1 ? repeatss : fullsspred );
							// refold_repeat_pose( repeatseq, repss, 3, 2, protpose );
							//        bool const use_asymmetric_centroid_scoring( true );
							//        refold_repeat_pose( repeatseq, repss, 5, 3, protpose, false, 0, 0, use_asymmetric_centroid_scoring );
							/// calc rmsd
							Real rmsd( 0.0 ), rmsd_single_repeat( 0.0 );
							{
								using namespace core::id;
								AtomID_Map< AtomID > atom_map;
								initialize_atomid_map( atom_map, protpose, id::GLOBAL_BOGUS_ATOM_ID );
								for ( Size i=1; i<= repeatlen; ++i ) {
									runtime_assert( protpose.residue(i).is_protein() );
									runtime_assert( pose.residue(i).is_protein() );
									atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
										AtomID( pose.residue(i).atom_index("CA"),i);
								}
								rmsd_single_repeat = rmsd_by_mapping( protpose, pose, atom_map );
								Size const nres_for_rmsd( min( chain_end( 1, pose ), chain_end( 1, protpose ) ) );
								for ( Size i=1; i<= nres_for_rmsd; ++i ) {
									runtime_assert( protpose.residue(i).is_protein() );
									runtime_assert( pose.residue(i).is_protein() );
									atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
										AtomID( pose.residue(i).atom_index("CA"),i);
								}
								rmsd = rmsd_by_mapping( protpose, pose, atom_map );
							}

							Real refold_handedness( 0.0 );
							{
								Vectors coords;
								for ( Size i=1; i<= chain_end(1,protpose); ++i ) coords.push_back( protpose.residue(i).xyz("CA") );
								refold_handedness = get_chirality( repeatlen, coords );
							}
							refolding_status << " R " << rr << F(9,3,rmsd) << F(9,3,rmsd_single_repeat) << F(9,3,refold_handedness);
							min_rmsd = min( min_rmsd, rmsd );
							if ( rmsd < 4.0 ) ++n_under_4;
						}
						refolding_status << " min_rmsd_" << refoldtag << ": " << F(9,3,min_rmsd) << " n_under_4_" << refoldtag <<
							": " << n_under_4 << ' ';
					} // RR=1,2
					stoptime = clock();
					Real const refolding_simtime = ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
					refolding_status << " refolding_simtime: " << F(9,3,refolding_simtime);
				}

				/// compute chirality of final structure
				Real handedness( 0.0 );
				{
					Vectors coords;
					for ( Size i=1; i<= chain_end(1,pose); ++i ) coords.push_back( pose.residue(i).xyz("CA") );
					handedness = get_chirality( repeatlen, coords );
				}

				cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
					" relax_rmsd: " << F(9,3,relax_rmsd) <<
					" passed_score_filter: " << passed_score_filter <<
					" centroid_score: "<< F(9,3,centroid_score) <<
					" centroid_rmsd: "<< F(9,3,centroid_rmsd) <<
					" centroid_simtime: " << F(9,3,centroid_simtime) <<
					" relax_simtime: " << F(9,3,relax_simtime) <<
					" anchorpos_loop_begin: " << anchorpos_loop_begin <<
					" anchorpos_loop_end: " << anchorpos_loop_end <<
					" cutpoint_loop_begin: " << cutpoint_loop_begin <<
					" cutpoint_loop_end: " << cutpoint_loop_end <<
					" anchorseq: " << anchorseq <<
					" repeatseq: " << repeatseq <<
					" repeatss: " << repeatss <<
					" repeatbb: " << repeatbb <<
					" na: " << pose.residue( chain_begin(2,pose) ).name1() <<
					" unbound_score: " << F(9,3,unbound_score ) <<
					" interface_score: " << F(9,3,interface_emap.dot( fa_scorefxn->weights() ) ) <<
					//      " sasapack_score_single_repeat: " << F(9,3,sasapack_score_single_repeat) <<
					//      " norme_score_single_repeat: " << F(9,3,norme_score_single_repeat) <<
					//      " sasa14_normalized_single_repeat: " << F(9,3,sasa14_normalized_single_repeat) <<
					" anchorpos: " << anchorpos <<
					" repeatlen: " << repeatlen <<
					" repeatsspred: " << repeatsspred <<
					" sasapack_score: " << F(9,3,sasapack_score) <<
					" norme_score: " << F(9,3,norme_score) <<
					" normsasa_score: " << F(9,3,normsasa_score) <<
					" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
					" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
					" exposed_polar_fraction: " << F(9,3,exposed_polar_sasa/(exposed_polar_sasa+exposed_nonpolar_sasa))<<
					" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
					" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
					" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc)<<
					" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc)<<
					" packstat_score: " << F(9,3,packstat_score) <<
					" n_buried_unsatisfied_donors_bb_per_repeat: " << F(9,3,n_buried_unsatisfied_donors_bb_per_repeat) <<
					" n_buried_unsatisfied_donors_sc_per_repeat: " << F(9,3,n_buried_unsatisfied_donors_sc_per_repeat) <<
					" n_buried_unsatisfied_acceptors_bb_per_repeat: " <<F(9,3,n_buried_unsatisfied_acceptors_bb_per_repeat)<<
					" n_buried_unsatisfied_acceptors_sc_per_repeat: " <<F(9,3,n_buried_unsatisfied_acceptors_sc_per_repeat)<<
					" handedness: " << F(9,3,handedness) <<
					" sspred_match: " << F(9,3,sspred_match) <<
					' ' << buried_unsatisfied_string <<
					' ' << refolding_status.str() <<
					" weighted_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) <<
					" weighted_interface_energies: "<< interface_emap.weighted_string_of( fa_scorefxn->weights() ) << endl;
				fflush( stdout );

				if ( passed_score_filter ) {
					pose.dump_pdb( outfilename );
					run_command("gzip "+outfilename );
				}
				if ( single_sim ) exit(0);
				fflush( stdout );
				check_simtime();
			}
		}
	}
}



///////////////////////////////////////////////////////////////////////////////
///
/// 06/26/12: results agree with mini branch for vdw, dna_env and dna_pair
///
///
void
centroid_test()
{

	Pose pose;
	pose_from_pdb( pose, start_file() );

	devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID_DNA );

	ScoreFunctionOP scorefxn( new ScoreFunction() );

	scorefxn->set_weight( vdw, 1.0 );
	scorefxn->set_weight( dna_env, 1.0 );
	scorefxn->set_weight( dna_pair, 1.0 );

	scorefxn->show( cout );

	Real const score( (*scorefxn)( pose ) );

	///
	cout << "centroid score: " << F(9,3,score)<< ' ' <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;


}


///////////////////////////////////////////////////////////////////////////////
Size
get_int_from_resample_bb_tag( string const & tag )
{
	if ( tag.find('-') != string::npos ) {
		strings const l( split_to_vector1( tag, "-" ) );
		runtime_assert( l.size() == 2 );
		runtime_assert( is_int( l[1] ) );
		runtime_assert( is_int( l[2] ) );
		return numeric::random::random_range( int_of( l[1] ), int_of( l[2] ) );

	} else {
		runtime_assert( is_int( tag ) );
		return int_of( tag );

	}


}
///////////////////////////////////////////////////////////////////////////////
// dashes or commas
//
Sizes
get_poslist_from_string( string const & tag )
{
	Sizes poslist;
	if ( tag.find('-') != string::npos ) {
		strings const l( split_to_vector1( tag, "-" ) );
		runtime_assert( l.size() == 2 );
		runtime_assert( is_int( l[1] ) );
		runtime_assert( is_int( l[2] ) );
		for ( int i= int_of( l[1] ); i<= int_of( l[2] ); ++i ) {
			poslist.push_back( Size(i) );
		}
	} else {
		strings const l( split_to_vector1( tag, "," ) );
		for ( Size i=1; i<= l.size(); ++i ) {
			runtime_assert( is_int( l[i] ) );
			poslist.push_back( int_of( l[i] ) );
		}
	}
	return poslist;
}


///////////////////////////////////////////////////////////////////////////////
void
centroid_build_test()
{
	Size const base_repeat( 3 );

	bool const helix2_is_a_strand( option[ my_options::helix2_is_a_strand ] );

	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
	runtime_assert( pose::symmetry::is_symmetric( *score3 ) );

	for ( Size n=1; n<= 10000; ++n ) {
		clock_t const starttime = clock();
		Pose pose;

		// create a symmetric pose
		Size const nrepeat( random_element( option[ my_options::nrepeats ]() )+2 ); // for overhang/base_repeat intxns

		Size const
			helix1_len( random_element( option[ my_options::helix1_lens ]() ) ),
			helix2_len( random_element( option[ my_options::helix2_lens ]() ) );
		string const
			turn1( random_element( option[ my_options::turn1s ]() ) ),
			turn2( random_element( option[ my_options::turn2s ]() ) );

		Size const repeatlen( helix1_len + helix2_len + turn1.size() + turn2.size() );

		TR.Trace << "repeatlen: " << repeatlen << ' ' << turn1 << ' ' << turn2 << ' ' <<
			helix1_len << ' ' << helix2_len << endl;

		char const helix2_bb( helix2_is_a_strand ? 'B' : 'A' ), helix2_ss( helix2_is_a_strand ? 'E' : 'H' );
		string repeatbb( string( helix1_len, 'A' ) + turn1 +
			string( helix2_len, helix2_bb ) + turn2 ),
			repeatss( string( helix1_len, 'H' ) + string( turn1.size(), 'L' ) +
			string( helix2_len, helix2_ss ) + string( turn2.size(), 'L' ) );

		string const target_repeatbb( repeatbb ), target_repeatss( repeatss );

		runtime_assert( repeatbb.size() == repeatss.size() && repeatbb.size() == repeatlen );

		Size nres_protein( nrepeat * repeatlen );
		for ( Size i=1; i<= nres_protein; ++i ) {
			pose.append_residue_by_bond( *get_vanilla_protein_residue( 'L' ), true ); // build_ideal_geometry
		}


		/// setup the centroid repeat sequence
		string const bbtag( string_of( helix1_len ) + turn1 + string_of( helix2_len ) + turn2 );

		string centroid_repeatseq;
		{
			for ( Size i=1; i<= repeatlen; ++i ) {
				if ( repeatss[i-1] != 'L' ) { // stupid backwards compatibility...
					centroid_repeatseq.push_back( oneletter_code_from_aa( aa_from_name( random_element(
						option[ my_options::centroid_helixseq ]() ) ) ) );
				} else centroid_repeatseq.push_back( 'G' );
			}
		}

		/// mutate the sequence in the helix region
		for ( Size i=1; i<= repeatlen; ++i ) {
			AA const aa( aa_from_oneletter_code( centroid_repeatseq[i-1] ) );
			for ( Size j=0; j< nrepeat; ++j ) {
				Size const seqpos( j*repeatlen + i );
				make_sequence_change( seqpos, aa, pose );
			}
		}

		/// bit of a hack: add a virtual residue at the end, fold starting from there...
		{
			//Residue const & lastrsd( pose.residue( nres_protein ) );
			ResidueOP vrtrsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
			pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...

			kinematics::FoldTree f( pose.total_residue() );
			f.reorder( pose.total_residue() );
			pose.fold_tree( f );
		}

		/// I think that these terminus types may be lost during replace residue calls?
		add_lower_terminus_type_to_pose_residue( pose, 1 );
		//add_upper_terminus_type_to_pose_residue( pose, nres_protein );
		pose.conformation().insert_chain_ending( nres_protein );

		Sizes fragseq_poslist;
		{
			Size const turn1_begin( helix1_len+1 ), turn1_end( helix1_len+turn1.size() ),
				turn2_begin( helix1_len+turn1.size()+helix2_len+1 );
			for ( Size i=1; i<= nres_protein; ++i ) {
				pose.set_phi  ( i, init_phi   );
				pose.set_psi  ( i, init_psi   );
				pose.set_omega( i, init_omega );
				pose.set_secstruct( i, 'L' );

				Size const rpos( (i-1)%repeatlen + 1 );

				if ( ( !option[ my_options::freeze_helixseq ] ) ||
						( rpos == 1 ||
						( rpos >= turn1_begin-1 && rpos <= turn1_end+1 ) ||
						( rpos >= turn2_begin-1 ) ) ) fragseq_poslist.push_back( i );
			}
		}

		//pose.dump_pdb( "test.pdb" );

		/// pick some fragments
		FragLib fraglib;
		string fullss, fullbb;
		for ( Size i=1; i<= nrepeat; ++i ) {
			fullss += repeatss;
			fullbb += repeatbb;
		}
		TR.Trace << "pick_frags: " << fullss << ' ' << fullbb << endl;
		// default for resample_centroid_seqs_fragweight is 1.0
		//
		Real const seq_weight( 0.0 ), ss_weight( 10.0 ), bb_weight( 100.0 );
		Sizes const fragsizes( make_vector1( 3, 9 ) );
		Size const nfrags( 200 );
		kinematics::MoveMap mm;
		for ( Size i=1; i<= chain_end( 1, pose ); ++i ) mm.set_bb( i, true );

		Sizes const homs_to_exclude;
		devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, pose, mm, fullss, seq_weight, ss_weight,
			fraglib, homs_to_exclude, bb_weight, fullbb );


		{ // get rid of frags that violate the bb constraints
			Sizes const fragsizes( make_vector1( 3, 9 ) );
			Sizes const min_nns( make_vector1( 50, 25 ) );

			for ( Size si=1; si<= 2; ++si ) {
				Size const fragsize( fragsizes[si] );
				Size const min_nn( min_nns[si] );
				devel::blab::classic_frags::TorsionFragmentLibrary & lib( fraglib.library( fragsize ) );

				for ( Size fragpos=1; fragpos<= lib.size(); ++fragpos ) {
					for ( Size nn=lib[ fragpos ].size(); nn> min_nn; --nn ) {
						bool badfrag( false );

						devel::blab::classic_frags::TorsionFragment const frag( lib[ fragpos ][ nn ] );

						for ( Size k=1; k<= fragsize; ++k ) {
							//char const ss( frag.get_secstruct(k) );
							Real const phi( frag.get_torsion(k,1 ) ), psi( frag.get_torsion( k, 2 ) ), omega( frag.get_torsion(k,3 ) );
							char const bigbin( torsion2big_bin( phi, psi, omega ) );
							Size const repeatpos( ( ( fragpos +k-1 ) -1 )%repeatlen + 1 );
							char const desired_bigbin( repeatbb[ repeatpos-1 ] );
							if ( desired_bigbin != 'X' && desired_bigbin != bigbin ) {
								runtime_assert( string("ABGEO").find( desired_bigbin ) != string::npos );
								badfrag = true;
								break;
							}
						}

						if ( badfrag ) {
							TR.Trace << "delete_fragment: " << I(4,fragsize) << I(4,fragpos) << I(4,nn) << endl;
							lib[ fragpos ].erase( nn );
						}
					} // nn
					TR.Trace << "final_nn: " << I(4,fragsize) << I(4,fragpos) << I(4,(fragpos-1)%repeatlen+1) <<
						I(4,lib[ fragpos ].size()) << endl;
				} // fragpos
			} // fragsize
		} // scope
		Real const frag_simtime( ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes


		/// setup symminfo?
		conformation::symmetry::SymmetryInfo symminfo;

		setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo );

		/// switch to centroid
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID ); // was CENTROID_DNA

		/// now make symmetric
		pose::symmetry::make_symmetric_pose( pose, symminfo );

		Real target_twist(0), target_rise( 0.0 );
		// if ( barrel_mode ) {
		{
			runtime_assert( base_repeat == 3 );
			target_twist = 360.0 / ( nrepeat-2 );
			target_rise = 0.0;
		}
		// } else if ( tpr_mode ) {
		//  target_twist = target_twist_tpr;
		//  target_rise = target_rise_tpr;
		// } else {
		//  runtime_assert( tal_mode );
		//  target_twist = 36.0;
		//  target_rise = 3.5; // plus or minus
		// }

		bool const use_twistrise_energy( true );
		use_linear_twist_rise_energy = true; /// HACKY GLOBAL VAR
		clock_t const centroid_starttime = clock();
		//Real const cycle_multiplier( 0.25 ); // reduce the cycles; note that they're already doubled due to fragseq_poslist
		fragseq_poslist.clear(); // no sequence changes during simulation !!!!!!!!

		//simple_fold_abinitio( fraglib, fragseq_poslist, pose, false /*use_twistrise_energy*/, target_twist, target_rise );
		simple_fold_abinitio( fraglib, fragseq_poslist, pose, use_twistrise_energy, target_twist, target_rise );
		//cycle_multiplier );
		Real const centroid_simtime( ( (double) clock() - centroid_starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

		//Real const final_centroid_score( (*score3)( pose ) );
		use_linear_twist_rise_energy = false; /// HACKY GLOBAL VAR

		/// get handedness
		Real handedness(0.0), helix_avg_rmsd(0), turn_avg_rmsd(0);
		if ( !helix2_is_a_strand ) { // will also fail if helices are shorter than length 7
			Vectors coords;
			for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
			handedness = get_chirality( repeatlen, coords );
			compute_bb_strain_rmsd( coords, nrepeat, repeatlen, helix1_len, helix2_len, turn1, turn2,
				helix_avg_rmsd, turn_avg_rmsd );
		}

		/// compute transform: twist, axis,
		Real twist, rise, min_radius, max_radius, com_radius;
		compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius );

		Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);

		compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
			helix1_dist, helix1_twist, helix2_dist, helix2_twist );

		Real const simtime( ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

		string const outfilename( output_tag() + "centroid_build_N"+lead_zero_string_of(n,4)+".pdb");

		Real const final_score( (*score3)( pose ) );

		cout << "final_scores: " << F(9,3,final_score) << ' ' << outfilename <<
			" nrepeat: " << nrepeat-2 <<
			" repeatlen: " << repeatlen <<
			" helix1_len: " << helix1_len <<
			" helix2_len: " << helix2_len <<
			" turn1: " << turn1 <<
			" turn2: " << turn2 <<
			" handedness: " << F(9,3,handedness) <<
			" twist: " << F(9,3,twist) <<
			" rise: " << F(9,3,rise) <<
			" helix1_dist: " << F(9,3,helix1_dist) <<
			" helix1_twist: " << F(9,3,helix1_twist) <<
			" helix2_dist: " << F(9,3,helix2_dist) <<
			" helix2_twist: " << F(9,3,helix2_twist) <<
			" min_radius: " << F(9,3,min_radius) <<
			" max_radius: " << F(9,3,max_radius) <<
			" com_radius: " << F(9,3,com_radius) <<
			" helix_avg_rmsd: " << F(9,3,helix_avg_rmsd) <<
			" turn_avg_rmsd: " << F(9,3,turn_avg_rmsd) <<
			" simtime: " << F(9,3,simtime) <<
			" centroid_simtime: " << F(9,3,centroid_simtime) <<
			" frag_simtime: " << F(9,3,frag_simtime) <<
			endl;

		if ( option[ my_options::output_pdb_files ] ) pose.dump_pdb( outfilename );

		fflush( stdout );

		check_simtime();
	}



}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
old_centroid_build_test()
{
	Size const nrepeat( 3 ), base_repeat( 2 );

	Sizes const helix_lens( option[ my_options::helix_lens ]() );
	strings const turns( option[ my_options::turns ] );

	for ( Size n=1; n<= 10000; ++n ) {

		Pose pose;

		// create a symmetric pose
		Size helix1_len( random_element( helix_lens ) ), helix2_len( random_element( helix_lens ) );
		string turn1( random_element( turns ) ), turn2( random_element( turns ) );
		if ( helix1_len > helix2_len ) {
			Size tmp( helix2_len );
			helix2_len = helix1_len;
			helix1_len = tmp;
		}
		Size const repeatlen( helix1_len + helix2_len + turn1.size() + turn2.size() );

		TR.Trace << "repeatlen: " << repeatlen << ' ' << turn1 << ' ' << turn2 << ' ' <<
			helix1_len << ' ' << helix2_len << endl;

		string repeatbb( string( helix1_len, 'A' ) + turn1 +
			string( helix2_len, 'A' ) + turn2 ),
			repeatss( string( helix1_len, 'H' ) + string( turn1.size(), 'L' ) +
			string( helix2_len, 'H' ) + string( turn2.size(), 'L' ) );

		string const target_repeatbb( repeatbb ), target_repeatss( repeatss );

		runtime_assert( repeatbb.size() == repeatss.size() && repeatbb.size() == repeatlen );

		Size nres_protein( nrepeat * repeatlen );
		for ( Size i=1; i<= nres_protein; ++i ) {
			pose.append_residue_by_bond( *get_vanilla_protein_residue( 'L' ), true ); // build_ideal_geometry
		}


		/// setup the centroid repeat sequence
		string const bbtag( string_of( helix1_len ) + turn1 + string_of( helix2_len ) + turn2 );

		string centroid_repeatseq;
		{
			for ( Size i=1; i<= repeatlen; ++i ) {
				if ( repeatss[i-1] == 'H' ) { // stupid backwards compatibility...
					centroid_repeatseq.push_back( oneletter_code_from_aa( aa_from_name( random_element(
						option[ my_options::centroid_helixseq ]() ) ) ) );
				} else centroid_repeatseq.push_back( 'G' );
			}
		}

		/// mutate the sequence in the helix region
		for ( Size i=1; i<= repeatlen; ++i ) {
			AA const aa( aa_from_oneletter_code( centroid_repeatseq[i-1] ) );
			for ( Size j=0; j< nrepeat; ++j ) {
				Size const seqpos( j*repeatlen + i );
				make_sequence_change( seqpos, aa, pose );
			}
		}

		/// bit of a hack: add a virtual residue at the end, fold starting from there...
		{
			//Residue const & lastrsd( pose.residue( nres_protein ) );
			ResidueOP vrtrsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
			pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...

			kinematics::FoldTree f( pose.total_residue() );
			f.reorder( pose.total_residue() );
			pose.fold_tree( f );
		}

		/// I think that these terminus types may be lost during replace residue calls?
		add_lower_terminus_type_to_pose_residue( pose, 1 );
		//add_upper_terminus_type_to_pose_residue( pose, nres_protein );
		pose.conformation().insert_chain_ending( nres_protein );

		Sizes fragseq_poslist;
		{
			Size const turn1_begin( helix1_len+1 ), turn1_end( helix1_len+turn1.size() ),
				turn2_begin( helix1_len+turn1.size()+helix2_len+1 );
			for ( Size i=1; i<= nres_protein; ++i ) {
				pose.set_phi  ( i, init_phi   );
				pose.set_psi  ( i, init_psi   );
				pose.set_omega( i, init_omega );
				pose.set_secstruct( i, 'L' );

				Size const rpos( (i-1)%repeatlen + 1 );

				if ( ( !option[ my_options::freeze_helixseq ] ) ||
						( rpos == 1 ||
						( rpos >= turn1_begin-1 && rpos <= turn1_end+1 ) ||
						( rpos >= turn2_begin-1 ) ) ) fragseq_poslist.push_back( i );
			}
		}

		//pose.dump_pdb( "test.pdb" );

		/// pick some fragments
		FragLib fraglib;
		string fullss, fullbb;
		for ( Size i=1; i<= nrepeat; ++i ) {
			fullss += repeatss;
			fullbb += repeatbb;
		}
		TR.Trace << "pick_frags: " << fullss << ' ' << fullbb << endl;
		// default for resample_centroid_seqs_fragweight is 1.0
		//
		Real const seq_weight( 0.0 ), ss_weight( 10.0 ), bb_weight( 100.0 );
		Sizes const fragsizes( make_vector1( 3, 9 ) );
		Size const nfrags( 200 );
		kinematics::MoveMap mm;
		for ( Size i=1; i<= chain_end( 1, pose ); ++i ) mm.set_bb( i, true );

		Sizes const homs_to_exclude;
		devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, pose, mm, fullss, seq_weight, ss_weight,
			fraglib, homs_to_exclude, bb_weight, fullbb );


		{ // get rid of frags that violate the bb constraints
			Sizes const fragsizes( make_vector1( 3, 9 ) );
			Sizes const min_nns( make_vector1( 50, 25 ) );

			for ( Size si=1; si<= 2; ++si ) {
				Size const fragsize( fragsizes[si] );
				Size const min_nn( min_nns[si] );
				devel::blab::classic_frags::TorsionFragmentLibrary & lib( fraglib.library( fragsize ) );

				for ( Size fragpos=1; fragpos<= lib.size(); ++fragpos ) {
					for ( Size nn=lib[ fragpos ].size(); nn> min_nn; --nn ) {
						bool badfrag( false );

						devel::blab::classic_frags::TorsionFragment const frag( lib[ fragpos ][ nn ] );

						for ( Size k=1; k<= fragsize; ++k ) {
							//char const ss( frag.get_secstruct(k) );
							Real const phi( frag.get_torsion(k,1 ) ), psi( frag.get_torsion( k, 2 ) ), omega( frag.get_torsion(k,3 ) );
							char const bigbin( torsion2big_bin( phi, psi, omega ) );
							Size const repeatpos( ( ( fragpos +k-1 ) -1 )%repeatlen + 1 );
							char const desired_bigbin( repeatbb[ repeatpos-1 ] );
							if ( desired_bigbin != 'X' && desired_bigbin != bigbin ) {
								runtime_assert( string("ABGEO").find( desired_bigbin ) != string::npos );
								badfrag = true;
								break;
							}
						}

						if ( badfrag ) {
							TR.Trace << "delete_fragment: " << I(4,fragsize) << I(4,fragpos) << I(4,nn) << endl;
							lib[ fragpos ].erase( nn );
						}
					} // nn
					TR.Trace << "final_nn: " << I(4,fragsize) << I(4,fragpos) << I(4,(fragpos-1)%repeatlen+1) <<
						I(4,lib[ fragpos ].size()) << endl;
				} // fragpos
			} // fragsize
		} // scope


		/// setup symminfo?
		conformation::symmetry::SymmetryInfo symminfo;

		setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo );

		/// switch to centroid
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID ); // was CENTROID_DNA

		/// now make symmetric
		pose::symmetry::make_symmetric_pose( pose, symminfo );

		Real target_twist(0), target_rise( 0.0 );
		// if ( barrel_mode ) {
		//  runtime_assert( base_repeat == 3 );
		//  target_twist = 360.0 / ( nrepeat-2 );
		//  target_rise = 0.0;
		// } else if ( tpr_mode ) {
		//  target_twist = target_twist_tpr;
		//  target_rise = target_rise_tpr;
		// } else {
		//  runtime_assert( tal_mode );
		//  target_twist = 36.0;
		//  target_rise = 3.5; // plus or minus
		// }

		bool const use_twistrise_energy( false );
		use_linear_twist_rise_energy = false; /// HACKY GLOBAL VAR
		//clock_t starttime = clock();
		Real const cycle_multiplier( 0.25 ); // reduce the cycles; note that they're already doubled due to fragseq_poslist
		simple_fold_abinitio( fraglib, fragseq_poslist, pose, use_twistrise_energy, target_twist, target_rise,
			cycle_multiplier );
		//Real const centroid_simtime( ( (double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

		//Real const final_centroid_score( (*score3)( pose ) );
		use_linear_twist_rise_energy = false; /// HACKY GLOBAL VAR

		/// get handedness
		Real handedness(0.0);
		{
			Vectors coords;
			for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
			handedness = get_chirality( repeatlen, coords );
		}

		/// compute transform: twist, axis,
		Real twist, rise, min_radius, max_radius, com_radius;
		compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius );

		Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);

		compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
			helix1_dist, helix1_twist, helix2_dist, helix2_twist );


		cout << "final_scores: " <<
			" repeatlen: " << repeatlen <<
			" helix1_len: " << helix1_len <<
			" helix2_len: " << helix2_len <<
			" turn1: " << turn1 <<
			" turn2: " << turn2 <<
			" handedness: " << F(9,3,handedness) <<
			" twist: " << F(9,3,twist) <<
			" rise: " << F(9,3,rise) <<
			" helix1_dist: " << F(9,3,helix1_dist) <<
			" helix1_twist: " << F(9,3,helix1_twist) <<
			" helix2_dist: " << F(9,3,helix2_dist) <<
			" helix2_twist: " << F(9,3,helix2_twist) <<
			" min_radius: " << F(9,3,min_radius) <<
			" max_radius: " << F(9,3,max_radius) <<
			" com_radius: " << F(9,3,com_radius) << endl;
		fflush( stdout );
	}



}




///////////////////////////////////////////////////////////////////////////////
void
unbound_frag_test()
{
	runtime_assert( option[ my_options::unfolded_sasas ].user() ); // make sure have datafile for buried sasa
	runtime_assert( option[ my_options::resample_centroid_seqs ].user() ||
		option[ my_options::centroid_helixseq ].user() );

	Real const min_helix2_len( 5 );
	bool const helix2_is_a_strand( option[ my_options::helix2_is_a_strand ] );

	bool const barrel_mode( option[ my_options::barrel_mode ] );
	bool const centroid_only( option[ my_options::centroid_only ] );
	bool const trefoil_mode( option[ my_options::trefoil_mode ] );
	if ( trefoil_mode ) { runtime_assert( barrel_mode ); }
	bool const pentafoil_mode( option[ my_options::pentafoil_mode ] );
	if ( pentafoil_mode ) { runtime_assert( barrel_mode ); }
	bool const tpr_mode( option[ my_options::tpr_mode ] );
	bool const free_mode( option[ my_options::free_mode ] ); // no constraints
	bool const tal_mode( option[ my_options::tal_mode ] );
	runtime_assert( barrel_mode || tpr_mode || free_mode || tal_mode );
	bool const force_triangles( option[ my_options::force_triangles ] );
	bool const force_biangles( option[ my_options::force_biangles ] );

	bool const force_handedness( option[ my_options::target_hand ].user() );

	Real const target_twist_tpr( option[ my_options::target_twist_tpr ] );
	Real const target_rise_tpr( option[ my_options::target_rise_tpr ] );
	//  Real const target_twist_tpr( 48.58 ); // mean values for 1n0a (regan tpr)
	//  Real const target_rise_tpr( 9.74 );

	bool const depth_filter( option[ my_options::min_depth ].user() );
	Real const min_depth( depth_filter ? option[ my_options::min_depth ]() : 0.0 );

	//bool const helix_angle_filter( option[ my_options::max_helix_angle ].user() );
	//Real const max_helix_angle( helix_angle_filter ? option[ my_options::max_helix_angle ]() : 0.0 );

	bool const resample_bbs( option[ my_options::resample_bbs ].user() );
	bool dont_force_first_helix_in( option[ my_options::dont_force_first_helix_in ] );
	if ( helix2_is_a_strand ) dont_force_first_helix_in = true;
	bool const force_second_helix_in( helix2_is_a_strand );

	bool const resample_centroid_seqs( option[ my_options::resample_centroid_seqs ].user() );
	map< string, strings > bbtag2centroid_seq;
	if ( resample_centroid_seqs ) {
		strings l( option[ my_options::resample_centroid_seqs ] );
		string bbtag;
		while ( !l.empty() ) {
			if ( is_int( l.front().substr(0,1) ) ) {
				bbtag = l.front();
			} else {
				runtime_assert( bbtag.size() );
				TR.Trace << "bbtag2centroid_seq: " << bbtag << ' ' << l.front() << endl;
				bbtag2centroid_seq[ bbtag ].push_back( l.front() );
			}
			l.erase( l.begin() );
		}
	}

	map< string, vector1< std::pair< Size, string > > > bbtag2force_seq;
	map< string, vector1< std::pair< Size, Sizes > > > bbtag2force_seq_relaxes;
	bool const force_sequence_positions( option[ my_options::force_sequence_positions ].user() );
	if ( force_sequence_positions ) {
		strings l( option[ my_options::force_sequence_positions ] );
		string bbtag;
		while ( !l.empty() ) {
			if ( is_int( l.front().substr(0,1) ) ) {
				bbtag = l.front();
			} else {
				runtime_assert( bbtag.size() );
				strings const ll( split_to_vector1( l.front(), ":" ) );
				if ( ll.size() == 2 ) {
					/// "STDNEQKRH:9"
					/// "p:9,13,16
					/// "h:15,17" // "h" = hydrophobic
					/// "p:4-11" //  "p" = polar
					Sizes poslist( get_poslist_from_string( ll[2] ) ); // uses "-" and also ","
					string name1s( ll[1] );
					if ( name1s == "h" ) name1s = "VAILMFYWGP";
					else if ( name1s == "p" ) name1s = "DENQKRHSTY"; // not sure about Y...
					for ( Sizes::const_iterator pos= poslist.begin(); pos != poslist.end(); ++pos ) {
						TR.Trace << "bbtag2force_seq: " << bbtag << ' ' << *pos << ' ' << name1s << endl;
						bbtag2force_seq[ bbtag ].push_back( make_pair( *pos, name1s ) );
					}
				} else {
					runtime_assert( ll.size() == 3 && ll[1] == "relax" );
					runtime_assert( is_int( ll[2] ) );
					Size const num2relax( int_of( ll[2] ) );
					strings const lll( split_to_vector1( ll[3], ",") );
					runtime_assert( num2relax <= lll.size() );
					Sizes relaxes;
					for ( strings::const_iterator s= lll.begin(); s != lll.end(); ++s ) {
						runtime_assert( is_int( *s ) );
						relaxes.push_back( int_of( *s ) );
					}
					bbtag2force_seq_relaxes[ bbtag ].push_back( make_pair( num2relax, relaxes ) );
				}
			}
			l.erase( l.begin() );
		}
	}


	bool const nodesign( option[ my_options::nodesign ] );
	if ( nodesign ) { runtime_assert( resample_centroid_seqs ); }

	/// for fullatom
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	//Real const centroid_score_filter_acceptance_rate( 0.05 );
	//bool const centroid_score_filter_pass_early( false );

	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
 	ScoreFunctionOP fa_scorefxn(0);
	if ( option[ OptionKeys::dna::specificity::score_function ].user() ) {
		fa_scorefxn = get_score_function_from_command_line();
	} else {
		fa_scorefxn = ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag );
	}
	//ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) ); // 09/16/15
	runtime_assert( pose::symmetry::is_symmetric( *score3 ) );

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	Real const donut_energy_weight( option[ my_options::donut_energy_weight ] );
	bool const use_donut_energy( donut_energy_weight>1e-3 );
	if ( free_mode ) { runtime_assert( !use_donut_energy ); }
	DonutWholeEnergy donut_whole_energy;
	Donut1B_Energy donut_1b_energy;
	if ( use_donut_energy ) {
		runtime_assert( barrel_mode );
		fa_scorefxn->add_extra_method( donut_whole, donut_energy_weight, donut_whole_energy );
		fa_scorefxn->add_extra_method( donut_1b   , donut_energy_weight, donut_1b_energy    );
	}


	TwistRiseEnergy    twist_rise_energy_tpr   ( target_twist_tpr, target_rise_tpr );
	TwistRise1B_Energy twist_rise_1b_energy_tpr( target_twist_tpr, target_rise_tpr );

	if ( tpr_mode ) {
		fa_scorefxn->add_extra_method( twistrise   , 0.5, twist_rise_energy_tpr );
		fa_scorefxn->add_extra_method( twistrise_1b, 0.5, twist_rise_1b_energy_tpr );
	}

	//  Size const nrepeat( option[ my_options::nrepeat ] ), base_repeat( option[ my_options::base_repeat ] );
	//Size const design_cycles( option[ my_options::design_cycles ] );

	clock_t search_starttime( clock() );

	Size const nstruct( option[ OptionKeys::out::nstruct ].user() ? option[ OptionKeys::out::nstruct ]() : 100000 );

	for ( Size n=1; n<= nstruct; ++n ) {
		clock_t sim_starttime = clock();

		Size base_repeat( option[ my_options::base_repeat ] ); // since we may fiddle with it later


		Pose pose;

		// create a symmetric pose
		Size nrepeat(0), helix1_len(0), helix2_len(0);
		string turn1, turn2;
		Size repeatlen(0 );
		char expected_hand( '-' );
		string resample_bb_tag;

		string tal_twist_type;

		if ( resample_bbs ) { // NOTE: convention here is that helix1 is the inner helix
			// EXAMPLE: " -resample_bbs  8-9,13,GB,13,GB,L  9,14-17,GBB,11,GBB,L "
			string const bb( random_element( option[ my_options::resample_bbs ]() ) );
			strings const l( split_to_vector1( bb, "," ) );
			runtime_assert( l.size() == 6 );
			//    runtime_assert( is_int( l[1] ) );
			//    runtime_assert( is_int( l[2] ) );
			//    runtime_assert( is_int( l[4] ) );
			runtime_assert( l[6].size() == 1 );
			nrepeat    = get_int_from_resample_bb_tag( l[1] );
			helix1_len = get_int_from_resample_bb_tag( l[2] );
			helix2_len = get_int_from_resample_bb_tag( l[4] );
			turn1 = l[3];
			turn2 = l[5];
			expected_hand = l[6][0];
			runtime_assert( expected_hand == 'R' || expected_hand == 'L' || expected_hand == 'U' );
			repeatlen = helix1_len + turn1.size() + helix2_len + turn2.size();
			resample_bb_tag = l[1];
			for ( Size j=2; j<= l.size(); ++j ) resample_bb_tag += "." + l[j];

		} else if ( helix2_is_a_strand ) {
			nrepeat = random_element( option[ my_options::nrepeats ]() ); // may be changed prior to relax
			helix1_len = random_element( option[ my_options::helix1_lens ]() );
			helix2_len = random_element( option[ my_options::helix2_lens ]() );
			turn1 = random_element( option[ my_options::turn1s ]() );
			turn2 = random_element( option[ my_options::turn2s ]() );
			runtime_assert( turn1 != "-" ); //temporary
			runtime_assert( turn2 != "-" );
			repeatlen = helix1_len + turn1.size() + helix2_len + turn2.size();

		} else if ( free_mode ) { // hacking
			nrepeat = random_element( option[ my_options::nrepeats ]() ); // may be changed prior to relax
			helix1_len = random_element( option[ my_options::helix1_lens ]() );
			helix2_len = random_element( option[ my_options::helix2_lens ]() );
			turn1 = random_element( option[ my_options::turn1s ]() );
			turn2 = random_element( option[ my_options::turn2s ]() );
			//runtime_assert( turn1 != "-" ); //temporary
			//runtime_assert( turn2 != "-" );
			repeatlen = helix1_len + get_turnlen( turn1 ) + helix2_len + get_turnlen( turn2 );

		} else if ( tal_mode ) { // hacking
			nrepeat = random_element( option[ my_options::nrepeats ]() ); // may be changed prior to relax
			helix1_len = random_element( option[ my_options::helix1_lens ]() );
			helix2_len = random_element( option[ my_options::helix2_lens ]() );
			turn1 = random_element( option[ my_options::turn1s ]() );
			turn2 = random_element( option[ my_options::turn2s ]() );
			//runtime_assert( turn1 != "-" ); //temporary
			//runtime_assert( turn2 != "-" );
			repeatlen = helix1_len + get_turnlen( turn1 ) + helix2_len + get_turnlen( turn2 );

			if ( repeatlen > Size( option[ my_options::max_repeatlen ] ) || helix2_len < helix1_len ) {
				--n;
				continue;
			}

			// what twist/rise to use?
			tal_twist_type = ( random_element( make_vector1( string("A"), string("B"), string("Z") ) ) );

		} else {
			nrepeat = random_element( option[ my_options::nrepeats ]() ); // may be changed prior to relax

			if ( option[ my_options::turn_pairs ].user() ) {
				strings const tp( split_to_vector1( random_element( option[ my_options::turn_pairs ]() ), "." ) );
				runtime_assert( tp.size() == 2 );
				turn1 = tp[1];
				turn2 = tp[2];
				helix1_len = random_element( option[ my_options::helix1_lens ]() );
				helix2_len = helix1_len + random_element( option[ my_options::helixlen_deltas ]() );
				repeatlen = helix1_len + helix2_len + turn1.size() + turn2.size();
			} else {
				repeatlen = random_element( option[ my_options::repeatlens ]() );
				Reals const helix1_len_range( option[ my_options::helix1_len_range ] );
				Size const min_helix1_len( int( helix1_len_range[1] * Real( repeatlen ) + 0.5 ) ),
					max_helix1_len( int( helix1_len_range[2] * Real( repeatlen ) + 0.5 ) );
				helix1_len = numeric::random::random_range( min_helix1_len, max_helix1_len );

				turn1 = random_element( option[ my_options::turns ]() );
				turn2 = random_element( option[ my_options::turns ]() );

				if ( turn1.size() + turn2.size() + helix1_len + min_helix2_len > repeatlen ) { --n; continue; }

				helix2_len = repeatlen - helix1_len - turn1.size() - turn2.size();

				if ( helix2_len < helix1_len ) { /// NEW CONVENTION: HELIX1 IS THE SHORTER HELIX
					Size const tmp( helix1_len );
					helix1_len = helix2_len;
					helix2_len = tmp;
				}
			}

		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// end setup
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		TR.Trace << "nrepeat: " << nrepeat << " repeatlen: " << repeatlen << ' ' << turn1 << ' ' << turn2 << ' ' <<
			helix1_len << ' ' << helix2_len << endl;


		char const helix2_bb( helix2_is_a_strand ? 'B' : 'A' ), helix2_ss( helix2_is_a_strand ? 'E' : 'H' );
		Size const turn1_len( get_turnlen( turn1 ) ), turn2_len( get_turnlen( turn2 ) );
		string repeatbb( string( helix1_len,       'A' ) + ( turn1_len==0 ? string() : turn1 ) +
			string( helix2_len, helix2_bb ) + ( turn2_len==0 ? string() : turn2 ) ),
			repeatss( string( helix1_len,       'H' ) + string( turn1_len, 'L' ) +
			string( helix2_len, helix2_ss ) + string( turn2_len, 'L' ) );

		string const target_repeatbb( repeatbb ), target_repeatss( repeatss );

		runtime_assert( repeatbb.size() == repeatss.size() && repeatbb.size() == repeatlen );

		Size nres_protein( nrepeat * repeatlen );
		for ( Size i=1; i<= nres_protein; ++i ) {
			pose.append_residue_by_bond( *get_vanilla_protein_residue( 'L' ), true ); // build_ideal_geometry
		}


		/// setup the centroid repeat sequence
		string const bbtag( string_of( helix1_len ) + turn1 + string_of( helix2_len ) + turn2 );

		string centroid_repeatseq;
		if ( resample_centroid_seqs ) {
			runtime_assert( bbtag2centroid_seq.count( bbtag ) );
			centroid_repeatseq = random_element( bbtag2centroid_seq[ bbtag ] );
			runtime_assert( centroid_repeatseq.size() == repeatlen );

		} else {
			for ( Size i=1; i<= repeatlen; ++i ) {
				if ( repeatss[i-1] != 'L' ) { // stupid backwards compatibility...
					centroid_repeatseq.push_back( oneletter_code_from_aa( aa_from_name( random_element(
						option[ my_options::centroid_helixseq ]() ) ) ) );
				} else centroid_repeatseq.push_back( 'G' );
			}
		}

		Sizes force_seq_relaxes;
		if ( force_sequence_positions ) {

			if ( bbtag2force_seq_relaxes.count( bbtag ) ) {
				vector1< std::pair< Size, Sizes > > const & relax_map( bbtag2force_seq_relaxes[ bbtag ] );
				for ( Size i=1; i<= relax_map.size(); ++i ) {
					Size const num2relax( relax_map[i].first );
					Sizes const & relaxes( relax_map[i].second );
					/// choosing WITH REPLACEMENT !!!
					for ( Size j=1; j<= num2relax; ++j ) force_seq_relaxes.push_back( random_element( relaxes ) );
				}
			}

			runtime_assert( bbtag2force_seq.count( bbtag ) );
			vector1< std::pair< Size, string > > const & poslist( bbtag2force_seq[ bbtag ] );

			for ( Size i=1; i<= poslist.size(); ++i ) {
				Size const repeatpos( poslist[i].first );
				runtime_assert( repeatpos <= repeatlen );
				string const name1s( poslist[i].second );
				if ( has_element( force_seq_relaxes, repeatpos ) ) {
					/// use whatever the sequence already is
					TR.Trace << "Relaxing sequence constraint at repeatpos: " << repeatpos << " using " <<
						centroid_repeatseq[ repeatpos-1 ] << " instead of " << name1s << endl;
					continue;
				}
				centroid_repeatseq[ repeatpos-1 ] = name1s[ int( uniform()*name1s.size() ) ];
			}
		}

		/// mutate the sequence in the helix region
		for ( Size i=1; i<= repeatlen; ++i ) {
			AA const aa( aa_from_oneletter_code( centroid_repeatseq[i-1] ) );
			for ( Size j=0; j< nrepeat; ++j ) {
				Size const seqpos( j*repeatlen + i );
				make_sequence_change( seqpos, aa, pose );
			}
		}

		/// bit of a hack: add a virtual residue at the end, fold starting from there...
		{
			//Residue const & lastrsd( pose.residue( nres_protein ) );
			ResidueOP vrtrsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
			pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...

			kinematics::FoldTree f( pose.total_residue() );
			f.reorder( pose.total_residue() );
			pose.fold_tree( f );
		}

		/// I think that these terminus types may be lost during replace residue calls?
		add_lower_terminus_type_to_pose_residue( pose, 1 );
		//add_upper_terminus_type_to_pose_residue( pose, nres_protein );
		pose.conformation().insert_chain_ending( nres_protein );

		Sizes fragseq_poslist;
		{
			Size const turn1_begin( helix1_len+1 ), turn1_end( helix1_len+turn1_len ),
				turn2_begin( helix1_len+turn1_len+helix2_len+1 );
			for ( Size i=1; i<= nres_protein; ++i ) {
				pose.set_phi  ( i, init_phi   );
				pose.set_psi  ( i, init_psi   );
				pose.set_omega( i, init_omega );
				pose.set_secstruct( i, 'L' );

				if ( !resample_centroid_seqs ) { // allow fragment sequence moves

					Size const rpos( (i-1)%repeatlen + 1 );

					if ( ( !option[ my_options::freeze_helixseq ] ) ||
							( rpos == 1 ||
							( rpos >= turn1_begin-1 && rpos <= turn1_end+1 ) ||
							( rpos >= turn2_begin-1 ) ) ) fragseq_poslist.push_back( i );
				}
			}
		}

		if ( option[ my_options::freeze_centroidseq ] ) fragseq_poslist.clear();

		//pose.dump_pdb( "test.pdb" );

		/// pick some fragments
		FragLib fraglib;
		string fullss, fullbb;
		for ( Size i=1; i<= nrepeat; ++i ) {
			fullss += repeatss;
			fullbb += repeatbb;
		}
		TR.Trace << "pick_frags: " << nrepeat << ' ' << repeatlen << ' ' << repeatss << ' ' << repeatbb << ' ' <<
			fullss << ' ' << fullbb << endl;
		// default for resample_centroid_seqs_fragweight is 1.0
		//
		Real const seq_weight( resample_centroid_seqs ? option[ my_options::resample_centroid_seqs_fragweight ]() : 0.0 ),
			ss_weight( 10.0 ), bb_weight( 100.0 );
		Sizes const fragsizes( make_vector1( 3, 9 ) );
		Size const nfrags( 200 );
		kinematics::MoveMap mm;
		for ( Size i=1; i<= chain_end( 1, pose ); ++i ) mm.set_bb( i, true );

		Sizes const homs_to_exclude;
		devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, pose, mm, fullss, seq_weight, ss_weight,
			fraglib, homs_to_exclude, bb_weight, fullbb );


		{ // get rid of frags that violate the bb constraints
			Sizes const fragsizes( make_vector1( 3, 9 ) );
			Sizes const min_nns( make_vector1( 50, 25 ) );

			for ( Size si=1; si<= 2; ++si ) {
				Size const fragsize( fragsizes[si] );
				Size const min_nn( min_nns[si] );
				devel::blab::classic_frags::TorsionFragmentLibrary & lib( fraglib.library( fragsize ) );

				for ( Size fragpos=1; fragpos<= lib.size(); ++fragpos ) {
					for ( Size nn=lib[ fragpos ].size(); nn> min_nn; --nn ) {
						bool badfrag( false );

						devel::blab::classic_frags::TorsionFragment const frag( lib[ fragpos ][ nn ] );

						for ( Size k=1; k<= fragsize; ++k ) {
							//char const ss( frag.get_secstruct(k) );
							Real const phi( frag.get_torsion(k,1 ) ), psi( frag.get_torsion( k, 2 ) ), omega( frag.get_torsion(k,3 ) );
							char const bigbin( torsion2big_bin( phi, psi, omega ) );
							Size const repeatpos( ( ( fragpos +k-1 ) -1 )%repeatlen + 1 );
							char const desired_bigbin( repeatbb[ repeatpos-1 ] );
							if ( desired_bigbin != 'X' && desired_bigbin != bigbin ) {
								runtime_assert( string("ABGEO").find( desired_bigbin ) != string::npos );
								badfrag = true;
								break;
							}
						}

						if ( badfrag ) {
							TR.Trace << "delete_fragment: " << I(4,fragsize) << I(4,fragpos) << I(4,nn) << endl;
							lib[ fragpos ].erase( nn );
						}
					} // nn
					TR.Trace << "final_nn: " << I(4,fragsize) << I(4,fragpos) << I(4,(fragpos-1)%repeatlen+1) <<
						I(4,lib[ fragpos ].size()) << endl;
				} // fragpos
			} // fragsize
		} // scope


		/// setup symminfo?
		conformation::symmetry::SymmetryInfo symminfo;

		setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo );

		/// switch to centroid
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID ); // was CENTROID_DNA

		/// now make symmetric
		pose::symmetry::make_symmetric_pose( pose, symminfo );


		bool use_twistrise_energy( true ); // default
		Real target_twist(0), target_rise( 0.0 );
		if ( barrel_mode ) {
			runtime_assert( base_repeat == 3 );
			if ( pentafoil_mode ) {
				runtime_assert( nrepeat-2==5 );
				target_twist = 720./5;
			} else {
				target_twist = 360.0 / ( nrepeat-2 );
			}
			target_rise = 0.0;
		} else if ( tpr_mode ) {
			target_twist = target_twist_tpr;
			target_rise = target_rise_tpr;
		} else if ( tal_mode ) {
			if ( tal_twist_type == "A" ) {
				target_twist = 32.7;
				target_rise = 2.9; // plus or minus
			} else if ( tal_twist_type == "B" ) {
				target_twist = 36.0;
				target_rise = 3.4; // plus or minus
			} else if ( tal_twist_type == "Z" ) {
				target_twist = 30.0;
				target_rise = -3.7; // plus or minus
			} else {
				utility_exit_with_message( "unrecognized tal_twist_type "+tal_twist_type );
			}

		} else if ( free_mode ) {
			if ( force_triangles ) {
				target_twist = 120.0;
				target_rise = 1000.0;
			} else if ( force_biangles ) {
				target_twist = 180.0;
				target_rise = 1000.0;
			} else {
				target_twist = target_rise = 0.0;
				use_twistrise_energy = false;
			}
		} else {
			utility_exit_with_message("bad unbound_frag_test mode");
		}

		use_linear_twist_rise_energy = true; /// HACKY GLOBAL VAR
		Real twistrise_energy_weight( 1.0 );
		if ( use_twistrise_energy && barrel_mode ) { // hack for really big barrels
			twistrise_energy_weight = max( 1.0, Real( nrepeat - 2 ) / 6.0 ); // multiplier of 1 for 6-repeat guys
		}

		//basic::prof_reset();
		cout << "checkpoint1" << endl; fflush( stdout );
		simple_fold_abinitio( fraglib, fragseq_poslist, pose, use_twistrise_energy, target_twist, target_rise, 0.0, twistrise_energy_weight );
		cout << "checkpoint2" << endl; fflush( stdout );

		if ( option[ my_options::make_trajectory ] ) {
			if ( n >= 10 ) break;
			continue; // talking about centroid trajectory I guess...
		}

		//basic::prof_show_oneliner();
		Real const centroid_simtime( ( (double) clock() - sim_starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

		Real const final_centroid_score( (*score3)( pose ) );
		use_linear_twist_rise_energy = false; /// HACKY GLOBAL VAR

		/// get handedness
		Real handedness(0.0);
		{
			Vectors coords;
			for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
			handedness = get_chirality( repeatlen, coords );
		}
		Real const centroid_handedness( handedness );

		/// compute transform: twist, axis,
		Real twist, rise, min_radius, max_radius, com_radius, depth;
		compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius, depth );
		Real const centroid_twist( twist ), centroid_rise( rise ), centroid_radius( min_radius );


		Real final_score( final_centroid_score ), relax_simtime( 0.0 );

		bool relax_this_decoy( false );

		Size const n2cdistpos( repeatlen*(nrepeat-2) ); // hacking

		Real n2cdist( pose.residue(1).xyz("N").distance( pose.residue(n2cdistpos).xyz("C") ) ),
			centroid_n2cdist( n2cdist );


		if ( barrel_mode ) {
			if ( option[ my_options::ignore_n2cdist ] ) {
				/// old way: ( rise > -2 && rise < 2 && twist > 0.9 * target_twist && twist < 1.1 * target_twist );
				Real const max_twist_dev(10.0);
				relax_this_decoy = ( rise > -2 && rise < 2 && fabs( twist - target_twist ) < max_twist_dev );
			} else {
				relax_this_decoy = n2cdist < 12.0;
			}
			cout << "relax_this_decoy: barrel_mode " << relax_this_decoy << ' ' << n2cdist << endl;

		} else if ( tpr_mode ) {
			relax_this_decoy = ( fabs( rise - target_rise ) <2 && fabs( twist - target_twist ) < 5 );

		} else if ( tal_mode ) {
			relax_this_decoy = ( fabs( rise - target_rise ) <2 && fabs( twist - target_twist ) < 5 );
			// this was the old code (pre-2/24/14)
			//
			// relax_this_decoy = ( fabs( 1.0 - fabs( handedness ) ) < 0.25 &&
			//            rise > -5 && rise < 5 &&
			//            twist < 50.0 ); // 1 or -1 for handedness, small theta
		} else if ( free_mode ) {
			if ( force_triangles ) {
				relax_this_decoy = ( twist > 115 && twist < 125 );
			} else if ( force_biangles ) {
				relax_this_decoy = ( twist > 160 );
			} else {
				relax_this_decoy = true;
			}
		}

		if ( trefoil_mode || pentafoil_mode ) {
			if ( trefoil_mode ) { runtime_assert( nrepeat == 5 ); } // silly
			else { runtime_assert( nrepeat == 7 ); }
			Size const nres_around( 1 + repeatlen*(nrepeat-2) ); // barrel mode nrepeat is real nrepeat+2  (!)
			Vectors coords;
			for ( Size i=1; i<= nres_around; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
			Stub const stub1( coords[1],coords[2],coords[3] ),
				stub2( coords[repeatlen+1],coords[repeatlen+2],coords[repeatlen+3] );
			Vector center, n, t;
			Real theta;
			get_stub_transform_data( stub1, stub2, center, n, t, theta );
			Real const total_rotation( get_total_rotation_degrees( coords, center, n ) ),
				winding( fabs( total_rotation)/360.0 ), desired_winding( 2.0 );
			// allow 20% fudge factor in either direction
			if (  winding <  desired_winding / 1.2 || winding > desired_winding * 1.2 ) relax_this_decoy = false;
			// hacking: write to stdout
			cout << "Xfoil_mode: total_rotation: " << F(9,3,total_rotation) <<
				" relax_this_decoy: " << relax_this_decoy <<
				" theta: " << F(9,3,degrees(theta)) <<
				" twist: " << F(9,3,twist) <<
				" target_twist: " << F(9,3,target_twist) <<
				" n2cdist: " << F(9,3,n2cdist) <<
				" rise: " << F(9,3,rise) <<
				endl;
			fflush(stdout);
			// if ( relax_this_decoy ) {
			// 	string outfilename( output_tag() + "unbound_frag_barrel_"+
			// 		string_of( helix1_len ) + turn1 + string_of( helix2_len ) + turn2 + "_"+
			// 		lead_zero_string_of( n,4 )+"_cen.pdb" );
			// 	pose.dump_pdb(outfilename);
			// }
		}

		if ( resample_bbs ) {

			Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);

			compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1_len, helix2_len,
				helix1_dist, helix1_twist, helix2_dist, helix2_twist );

			/// for resample_bbs, convention is that first helix is inner helix
			if ( !dont_force_first_helix_in && helix1_dist > helix2_dist ) relax_this_decoy = false;
			if ( force_second_helix_in && helix1_dist < helix2_dist ) relax_this_decoy = false;

			if ( ( expected_hand == 'R' && centroid_handedness < 0 ) ||
					( expected_hand == 'L' && centroid_handedness > 0 ) ) relax_this_decoy = false;

			cout << "relax_this_decoy: resample_bbs " << relax_this_decoy << ' '<< centroid_handedness << ' ' << helix1_dist << ' ' << helix2_dist << ' ' <<
				expected_hand << endl;

		} else if ( helix2_is_a_strand ) { // want the strand on the inside, I guess...

			Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);

			compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1_len, helix2_len,
				helix1_dist, helix1_twist, helix2_dist, helix2_twist );

			if ( force_second_helix_in && helix1_dist < helix2_dist ) relax_this_decoy = false;
		}

		if ( force_handedness ) {
			string const target_hand( option[ my_options::target_hand ] );
			if ( ( target_hand == "R" && centroid_handedness < 0 ) ||
					( target_hand == "L" && centroid_handedness > 0 ) ) {
				relax_this_decoy = false;
			}
		}

		if ( depth_filter && depth < min_depth ) {
			TR.Trace << "short: " << depth << ' ' << min_depth << endl;
			relax_this_decoy = false;
		}

		cout << "relax_this_decoy: final " << relax_this_decoy << endl;

		//Size nrepeat_this_decoy( nrepeat );

		ostringstream relax_status; // some scores only get written if we relaxed this decoy

		string fullsspred( "-" );
		Real searchtime(0);
		bool refold_this_decoy( true );
		if ( relax_this_decoy || dry_run() ) {
			searchtime = ((double) clock() - search_starttime )/( CLOCKS_PER_SEC*60 ); // records time other than in relax,refold

			// convert to fullatom
			conformation::symmetry::SymmetryInfo symminfo( *pose::symmetry::symmetry_info( pose ) );
			pose::symmetry::make_asymmetric_pose( pose );
			devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
			pose::symmetry::make_symmetric_pose( pose, symminfo );

			if ( barrel_mode ) {
				/// delete the last (base_repeat-1) repeats, actually just the last 2...
				nrepeat -= 2;
				if ( nrepeat <= base_repeat ) base_repeat = 2;
				runtime_assert( nrepeat >= base_repeat );
				nres_protein -= 2*repeatlen;
				conformation::symmetry::SymmetryInfo symminfo2;
				//Size const base_repeat_for_fullatom( nrepeat <= 3 ? 2 : 3 );
				setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo2 );

				pose::symmetry::make_asymmetric_pose( pose );
				{
					Size const delete_begin( nrepeat*repeatlen + 1 );
					pose.conformation().delete_residue_range_slow( delete_begin, pose.total_residue() );
					ResidueOP vrtrsd
						( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
					pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...
					pose.set_psi  ( nres_protein  , pose.psi  ( repeatlen   ) );
					pose.set_omega( nres_protein  , pose.omega( repeatlen   ) );
					//pose.set_phi  ( nres_protein+1, pose.phi  ( repeatlen+1 ) );
					pose.set_torsion( id::TorsionID( nres_protein+1, id::BB, 1 ), pose.phi  ( repeatlen+1 ) );

					kinematics::FoldTree f( pose.total_residue() );
					f.reorder( pose.total_residue() );
					pose.fold_tree( f );
					pose.conformation().insert_chain_ending( nres_protein );
				}
				pose::symmetry::make_symmetric_pose( pose, symminfo2 );
			}

			runtime_assert( pose.total_residue() == nres_protein+1 );

			clock_t const relax_starttime = clock();


			SequenceConstraints sequence_constraints;
			if ( force_sequence_positions ) {
				vector1< std::pair< Size, string > > const & poslist( bbtag2force_seq[ bbtag ] );
				for ( Size i=1; i<= poslist.size(); ++i ) {
					Size const repeatpos( poslist[i].first );
					if ( has_element( force_seq_relaxes, repeatpos ) ) {
						TR.Trace << "Relaxing design constraints at repeatpos " << repeatpos << endl;
						continue;
					}
					string const name1s( poslist[i].second );
					for ( Size k=0; k<nrepeat; ++k ) {
						Size const seqpos( k*repeatlen + repeatpos );
						sequence_constraints[ seqpos ] = name1s;
					}
				}
			}

			if ( nodesign ) {
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					if ( pose.residue(i).is_protein() ) {
						sequence_constraints[ i ] = string( 1, pose.residue(i).name1() );
					}
				}
			}


			if ( option[ my_options::force_helix_capping ] ) {
				Sizes ncaps( 1, repeatlen );
				Sizes ccaps( 1, helix1_len+1 );
				if ( !helix2_is_a_strand ) {
					ncaps.push_back( helix1_len + turn1_len );
					ccaps.push_back( helix1_len + turn1_len + helix2_len+1 );
				}
				for ( Size repeatpos=1; repeatpos<= repeatlen; ++repeatpos ) {
					string seqcst;
					if ( has_element( repeatpos, ccaps ) && repeatbb[ repeatpos-1 ] == 'G' ) seqcst = "G";
					if ( has_element( repeatpos, ncaps ) && repeatbb[ repeatpos-1 ] == 'B' ) seqcst = "STND";
					if ( seqcst.empty() || sequence_constraints.count( repeatpos ) ) continue; // earlier constraints have prio
					TR.Trace << "force_helix_capping: " << repeatpos << ' ' << seqcst << endl;
					for ( Size k=0; k<nrepeat; ++k ) {
						Size const seqpos( k*repeatlen + repeatpos );
						sequence_constraints[ seqpos ] = seqcst;
					}
				}
			}


			if ( option[ my_options::layer_design ] ) {
				///
				push_sequence_constraints_to_base_repeat( pose, sequence_constraints );
				string fullss;
				for ( Size i=1; i<= nrepeat; ++i ) fullss += target_repeatss;
				fullss.push_back( 'X' );
				runtime_assert( fullss.size() == pose.total_residue() );
				add_layer_design_constraints( pose, fullss, sequence_constraints );
			}

			if ( !centroid_only ) {
				bool const add_jump_flex( false ), use_atom_pair_constraints( false ), use_coordinate_constraints( false );
				bool const use_twistrise_energy( tal_mode && option[ my_options::use_fullatom_twistrise_energy ] );
				symmetric_design_and_relax( repeatlen, nrepeat, option[ my_options::use_softrep_for_design ],
					use_donut_energy, option[ my_options::design_cycles ], sequence_constraints, pose,
					add_jump_flex, use_atom_pair_constraints, use_coordinate_constraints,
					use_twistrise_energy, target_twist, target_rise );
			}

			relax_simtime = ((double) clock() - relax_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

			final_score = (*fa_scorefxn)( pose );
			Real const res_score( final_score / repeatlen );

			///
			bools subset( pose.total_residue(), false );
			Size const base_repeat_offset( repeatlen*(base_repeat-1));
			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;

			Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa, buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score;
			compute_sasa_scores_for_subset_slow( 5, subset, pose, sasapack_score, norme_score, normsasa_score,
				exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa,
				buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score );

			Real n_buried_unsatisfied_donors_bb_per_repeat, n_buried_unsatisfied_acceptors_bb_per_repeat,
				n_buried_unsatisfied_donors_sc_per_repeat, n_buried_unsatisfied_acceptors_sc_per_repeat;
			{ // compute unsats over full pose
				Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
					n_buried_unsatisfied_acceptors_sc;
				bools subset_full( pose.total_residue(), true ); subset_full[ pose.total_residue() ] = false;
				get_buried_unsatisfied_counts_real_slow( subset_full, pose,
					n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc,
					n_buried_unsatisfied_acceptors_bb, n_buried_unsatisfied_acceptors_sc );

				n_buried_unsatisfied_donors_bb_per_repeat = Real( n_buried_unsatisfied_donors_bb )/nrepeat;
				n_buried_unsatisfied_donors_sc_per_repeat = Real( n_buried_unsatisfied_donors_sc )/nrepeat;
				n_buried_unsatisfied_acceptors_bb_per_repeat = Real( n_buried_unsatisfied_acceptors_bb )/nrepeat;
				n_buried_unsatisfied_acceptors_sc_per_repeat = Real( n_buried_unsatisfied_acceptors_sc )/nrepeat;
			}

			if ( option[ my_options::pick_decoys ] ) {
				Real const exposed_polar_fraction( exposed_polar_sasa / ( exposed_polar_sasa + exposed_nonpolar_sasa ) );
				if ( exposed_polar_fraction > 0.55 || exposed_polar_fraction < 0.35 ) refold_this_decoy = false;
				if ( sasapack_score > 0 ) refold_this_decoy = false;

				if ( refold_this_decoy ) {
					Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
						n_buried_unsatisfied_acceptors_sc;
					get_buried_unsatisfied_counts_real_slow( subset, pose,
						n_buried_unsatisfied_donors_bb,
						n_buried_unsatisfied_donors_sc,
						n_buried_unsatisfied_acceptors_bb,
						n_buried_unsatisfied_acceptors_sc );
					if ( ( n_buried_unsatisfied_donors_bb + n_buried_unsatisfied_donors_sc > 1 ) ||
							( n_buried_unsatisfied_acceptors_bb + n_buried_unsatisfied_acceptors_sc > 0 ) ) refold_this_decoy = false;
				}
			}

			n2cdist = pose.residue(1).xyz("N").distance( pose.residue(chain_end(1,pose)).xyz("C") );

			/// compute psipred ss using single sequence
			string repeatsspred( "-" );
			Real sspred_match(0.0), sspred_turn1(0.0), sspred_turn2( 0.0 ), sspred_helix1(0.0), sspred_helix2(0.0);
			{ // run psipred to re-predict the secondary structure
				vector1< Reals > pred_HEL;
				string const sequence( pose.sequence().substr(0,chain_end(1,pose) ) );
				//string const secstruct( pose.secstruct().substr(0,chain_end(1,pose) ) );
				string sspred;
				string const repeatss( pose.secstruct().substr(0,repeatlen) ); // current, not the target

				run_psipred( sequence, sspred, pred_HEL );

				if ( sspred.size() == sequence.size() ) { // success
					fullsspred = sspred;
					runtime_assert( sspred.size() == nrepeat * repeatlen );
					repeatsspred = sspred.substr( (base_repeat-1)*repeatlen, repeatlen );
					Size const turn1_begin( helix1_len+1 ), turn1_end( helix1_len + turn1_len ),
						turn2_begin( helix1_len+turn1_len+helix2_len+1 );
					for ( Size i=1; i<= sspred.size(); ++i ) {
						Size const rpos( (i-1)%repeatlen+1 );
						char const ss( repeatss[ rpos-1 ] );
						Real const prediction( ss == 'H' ? pred_HEL[i][1] : ( ss == 'E' ? pred_HEL[i][2] : pred_HEL[i][3] ) );
						sspred_match += prediction;
						if ( rpos < turn1_begin                     ) sspred_helix1 += prediction;
						if ( rpos > turn1_end && rpos < turn2_begin ) sspred_helix2 += prediction;
						if ( rpos >= turn1_begin-1 && rpos <= turn1_end+1 ) sspred_turn1 += prediction;
						if ( rpos == 1 || rpos >= turn2_begin-1           ) sspred_turn2 += prediction;
					}
					sspred_match  /= (  repeatlen * nrepeat );
					sspred_helix1 /= ( helix1_len * nrepeat );
					sspred_helix2 /= ( helix2_len * nrepeat );
					sspred_turn1  /= ( ( turn1_len+2) * nrepeat );
					sspred_turn2  /= ( ( turn2_len+2) * nrepeat );
				}
			}

			Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);
			if ( helix1_len ) {
				runtime_assert( helix1_len + turn1_len + helix2_len + turn2_len == repeatlen );

				compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1_len, helix2_len,
					helix1_dist, helix1_twist, helix2_dist, helix2_twist );
			}

			Real const
				inner_helix_dist( min( helix1_dist, helix2_dist ) ),
				outer_helix_dist( max( helix1_dist, helix2_dist ) );
			Real const inner_helix_twist( helix1_dist < helix2_dist ? helix1_twist : helix2_twist );
			Real const outer_helix_twist( helix1_dist < helix2_dist ? helix2_twist : helix1_twist );

			Size const inner_helix_len( helix1_dist < helix2_dist ? helix1_len : helix2_len );
			Size const outer_helix_len( helix1_dist < helix2_dist ? helix2_len : helix1_len );
			string const inner_helix_turn( helix1_dist < helix2_dist ? turn1 : turn2 );
			string const outer_helix_turn( helix1_dist < helix2_dist ? turn2 : turn1 );

			string const buried_unsatisfied_string( get_buried_unsatisfied_string( pose ) );

			relax_status <<
				" sasapack_score: " << F(9,3,sasapack_score) <<
				" norme_score: " << F(9,3,norme_score) <<
				" normsasa_score: " << F(9,3,normsasa_score) <<
				" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
				" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
				" exposed_polar_fraction: " << F(9,3,exposed_polar_sasa/(exposed_polar_sasa+exposed_nonpolar_sasa))<<
				" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
				" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
				" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc)<<
				" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc)<<
				" packstat_score: " << F(9,3,packstat_score) <<
				" res_score: " << F(9,3,res_score ) <<
				" n_buried_unsatisfied_donors_bb_per_repeat: " << F(9,3,n_buried_unsatisfied_donors_bb_per_repeat) <<
				" n_buried_unsatisfied_donors_sc_per_repeat: " << F(9,3,n_buried_unsatisfied_donors_sc_per_repeat) <<
				" n_buried_unsatisfied_acceptors_bb_per_repeat: " <<F(9,3,n_buried_unsatisfied_acceptors_bb_per_repeat)<<
				" n_buried_unsatisfied_acceptors_sc_per_repeat: " <<F(9,3,n_buried_unsatisfied_acceptors_sc_per_repeat)<<
				" repeatsspred: " << repeatsspred <<
				" sspred_match: " << F(9,3,sspred_match) <<
				" sspred_helix1: " << F(9,3,sspred_helix1) <<
				" sspred_helix2: " << F(9,3,sspred_helix2) <<
				" sspred_turn1: " << F(9,3,sspred_turn1) <<
				" sspred_turn2: " << F(9,3,sspred_turn2) <<
				" helix1_dist: " << F(9,3,helix1_dist) <<
				" helix1_twist: " << F(9,3,helix1_twist) <<
				" helix2_dist: " << F(9,3,helix2_dist) <<
				" helix2_twist: " << F(9,3,helix2_twist) <<
				" inner_helix_dist: " << F(9,3,inner_helix_dist) <<
				" inner_helix_twist: " << F(9,3,inner_helix_twist) <<
				" outer_helix_dist: " << F(9,3,outer_helix_dist) <<
				" outer_helix_twist: " << F(9,3,outer_helix_twist) <<
				" inner_helix_len: " << I(4,inner_helix_len) <<
				" inner_helix_turn: " << inner_helix_turn <<
				" outer_helix_len: " << I(4,outer_helix_len) <<
				" outer_helix_turn: " << outer_helix_turn <<
				' ' << buried_unsatisfied_string <<
				' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() );


			/// RE-compute transform: twist, axis, after relax
			compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius );
			{
				Vectors coords;
				for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
				handedness = get_chirality( repeatlen, coords );
			}
		}

		string worktag( resample_bbs ? "resample_bb_"+resample_bb_tag : ( barrel_mode ? "barrel" : "tal" ) ),
			simfile( shared_output_tag()+"unbound_frag.work" );

		if ( option[ my_options::local_simfile ] ) {
			runtime_assert( shared_output_tag().substr(0,3) == "../" );
			simfile = shared_output_tag().substr(1)+"_unbound_frag.work";
		}

		worktag += "_" + string_of( repeatlen ) + "_" + string_of( nrepeat );
		if ( tal_mode ) worktag += "_TW" + tal_twist_type;
		if ( force_sequence_positions && force_seq_relaxes.size() ) {
			Size count(0);
			for ( Size i=1; i<= repeatlen; ++i ) if ( has_element( force_seq_relaxes, i ) ) ++count;
			worktag += "_rlxseq"+string_of(count);
		}

		//   bool const passed_centroid_score_filter
		//    ( append_score_to_scorefile_and_filter( worktag+"_cen", final_centroid_score,
		//                        centroid_score_filter_acceptance_rate,
		//                        centroid_score_filter_pass_early, simfile ) );

		bool const passed_score_filter( relax_this_decoy && // short-circuit evaluation
			append_score_to_scorefile_and_filter( worktag, final_score,
			score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );


		string const repeatseq( pose.sequence().substr(0,repeatlen) );
		Size const base_repeat_offset( repeatlen * ( base_repeat-1 ) );
		repeatbb = torsion2big_bin_string( base_repeat_offset+1, base_repeat_offset+repeatlen, pose );
		repeatss = pose.secstruct().substr(0,repeatlen);
		ostringstream refolding_status;

		/// try refolding de novo
		Size const n_refold( option[ my_options::n_refold ] ), n_refold_sspred( option[ my_options::n_refold_sspred ] );
		Real refolding_simtime(0);
		if ( refold_this_decoy && n_refold+n_refold_sspred > 0 && ( passed_score_filter || dry_run() ) ) {
			clock_t const refold_starttime = clock();
			for ( Size RR=1; RR<= 2; ++RR ) {
				Real min_refold_rmsd(1e6);
				Size n_under_4(0);
				if ( RR == 2 && fullsspred.size() != nrepeat * repeatlen ) continue;
				if ( RR==1 ) refolding_status << " n_refold: " << I(3,n_refold);
				else         refolding_status << " n_refold_sspred: " << I(3,n_refold_sspred);
				Size const rr_end( RR==1 ? n_refold : n_refold_sspred );
				for ( Size rr=1; rr<= rr_end; ++rr ) {
					Pose protpose;
					string const repss( RR == 1 ? repeatss : fullsspred );
					//Size const nrepeat_refold( barrel_mode ? nrepeat+2 : nrepeat );
					Size const nrepeat_refold( option[ my_options::nrepeat_refold ].user()?
						option[ my_options::nrepeat_refold ]() : ( barrel_mode ? nrepeat+2 : nrepeat ) );
					Size const base_repeat_refold( barrel_mode ? 3 : base_repeat );
					refold_repeat_pose( repeatseq, repss, nrepeat_refold, base_repeat_refold, protpose );
					/// calc rmsd
					Real rmsd( 0.0 ), rmsd_single_repeat( 0.0 );
					{
						using namespace core::id;
						AtomID_Map< AtomID > atom_map;
						initialize_atomid_map( atom_map, protpose, id::GLOBAL_BOGUS_ATOM_ID );
						for ( Size i=1; i<= repeatlen; ++i ) {
							atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
								AtomID( pose.residue(i).atom_index("CA"),i);
						}
						rmsd_single_repeat = rmsd_by_mapping( protpose, pose, atom_map );
						Size nrepeat_for_refolding_rmsd( min( nrepeat, nrepeat_refold ) );
						if ( option[ my_options::nrepeat_for_refolding_rmsd ].user() ) {
							nrepeat_for_refolding_rmsd = option[ my_options::nrepeat_for_refolding_rmsd ];
						}
						Size const nres_for_rmsd( min( chain_end( 1, protpose), nrepeat_for_refolding_rmsd * repeatlen ) );
						for ( Size i=1; i<= nres_for_rmsd; ++i ) {
							atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
								AtomID( pose.residue(i).atom_index("CA"),i);
						}
						rmsd = rmsd_by_mapping( protpose, pose, atom_map );
						min_refold_rmsd = min( min_refold_rmsd, rmsd );
						if ( rmsd < 4 ) ++n_under_4;
					}

					Real refold_handedness( 0.0 );
					{
						Vectors coords;
						for ( Size i=1; i<= chain_end( 1, protpose); ++i ) coords.push_back( protpose.residue(i).xyz("CA") );
						refold_handedness = get_chirality( repeatlen, coords );
					}
					Real tmptwist, tmprise, tmpmin_radius, tmpmax_radius, tmpcom_radius;
					compute_repeat_params( protpose, tmptwist, tmprise, tmpmin_radius, tmpmax_radius, tmpcom_radius );
					refolding_status << " R " << rr << F(9,3,rmsd) << F(9,3,rmsd_single_repeat) << F(9,3,refold_handedness) <<
						F(9,3,tmprise) << F(9,3,tmptwist);
				}
				refolding_status << " min_refold" << RR << "_rmsd: "<< F(9,3,min_refold_rmsd) <<
					" n_refold" << RR << "_under_4: " << n_under_4;

			}
			refolding_simtime = ((double) clock() - refold_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
		}


		string outfilename( output_tag() + "unbound_frag_"+worktag+"_"+
			string_of( helix1_len ) + turn1 + string_of( helix2_len ) + turn2 + "_"+
			lead_zero_string_of( n,4 )+".pdb" );

		Real const simtime = ((double) clock() - sim_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

		cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" handedness: " << F(9,3,handedness) <<
			" min_radius: " << F(9,3,min_radius ) <<
			" rise: " << F(9,3,rise ) <<
			" twist: " << F(9,3,twist) <<
			" final_centroid_score: "<< F(9,3,final_centroid_score) <<
			" passed_score_filter: " << passed_score_filter <<
			" relaxed: " << relax_this_decoy <<
			" min_radius: " << F(9,3,min_radius )<<
			" max_radius: " << F(9,3,max_radius )<<
			" com_radius: " << F(9,3,com_radius )<<
			" centroid_handedness: " << F(9,3,centroid_handedness) <<
			" centroid_radius: " << F(9,3,centroid_radius ) <<
			" centroid_rise: " << F(9,3,centroid_rise ) <<
			" centroid_twist: " << F(9,3,centroid_twist ) <<
			" repeatseq: " << repeatseq <<
			" repeatbb: " << repeatbb <<
			" repeatss: " << repeatss <<
			" turn1: " << turn1 <<
			" turn2: " << turn2 <<
			" helix1_len: " << helix1_len <<
			" helix2_len: " << helix2_len <<
			" repeatlen: " << repeatlen <<
			" nrepeat: " << nrepeat <<
			" n2cdist: " << F(9,3,n2cdist) <<
			" centroid_n2cdist: " << F(9,3,centroid_n2cdist) <<
			" target_repeatbb: " << target_repeatbb <<
			" target_repeatss: " << target_repeatss <<
			" " << relax_status.str() <<
			" " << refolding_status.str() <<
			" centroid_simtime: " << F(9,3,centroid_simtime) <<
			" relax_simtime: " << F(9,3,relax_simtime) <<
			" refolding_simtime: " << F(9,3,refolding_simtime) <<
			" searchtime: " << F(9,3,searchtime) <<
			" simtime: " << F(9,3,simtime) <<
			endl;

		if ( passed_score_filter || dry_run() ) {
			pose.dump_pdb( outfilename );
			run_command("gzip "+outfilename ); // NEW NEW NEW

			if ( option[ my_options::phase_fragment ].user() ) { // hacking
				Size const n_frag_repeats( option[ my_options::phase_fragment ] );
				outfilename += "_frag.pdb";

				pose::symmetry::make_asymmetric_pose( pose );

				Size const nres_frag( n_frag_repeats * repeatlen );

				pose.conformation().delete_residue_range_slow( nres_frag+1, pose.total_residue() );

				pose.dump_pdb( outfilename );
				string const fasta_file( outfilename+"_tmp.fasta" );
				{
					ofstream out( fasta_file.c_str() );
					out << ">tmp\n" << pose.sequence() << endl;
					out.close();
				}
				Real llg(0.0), tfz(0.0), rwork(1.0), rfree(1.0);
				string actual_space_group("X1"), refined_pdbfile("none");
				Size actual_num_models(0);
				run_phaser_and_phenix_refinement_generic( outfilename,
					option[ my_options::mtz_file ],
					fasta_file,
					option[ my_options::space_group ],
					option[ my_options::num_models ],
					option[ my_options::composition ],
					option[ my_options::llg_threshold_for_refinement ],
					actual_num_models, llg, tfz, rwork, rfree,
					actual_space_group, refined_pdbfile );
				run_command("gzip "+outfilename+"*");

				cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
					" passed_score_filter: " << passed_score_filter <<
					" llg: " << F(9,3,llg) <<
					" tfz: " << F(9,3,tfz) <<
					" rwork: " << F(9,3,rwork) <<
					" rfree: " << F(9,3,rfree) <<
					" actual_num_models: " << actual_num_models <<
					" target_num_models: " << option[ my_options::num_models ] <<
					" space_group: " << actual_space_group <<
					" refined_pdbfile: " << refined_pdbfile <<
					endl;
			} // phase_fragment
		} // passed_score_filter

		fflush( stdout );
		check_simtime();

		search_starttime = clock(); // restart the search timer
	}




}


///////////////////////////////////////////////////////////////////////////////
//


///////////////////////////////////////////////////////////////////////////////
void
refold_test()
{
	Size const nrepeat( 5 ), base_repeat( 3 ); // make these configurable?

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );


	string const datafile( start_file() );

	vector1< strings > targets;
	{
		ifstream data( datafile.c_str() );
		string line;
		while ( getline( data, line ) ) {
			istringstream l( line );
			string repeatseq, repeatss, filename;
			l >> repeatseq >> repeatss >> filename;
			runtime_assert( !l.fail() );
			runtime_assert( repeatseq.size() == repeatss.size() );
			targets.push_back( make_vector1( repeatseq, repeatss, filename ) );
		}
	}

	numeric::random::random_permutation( targets, numeric::random::rg() );


	for ( Size ti=1; ti<= targets.size(); ++ti ) {
		string const repeatseq( targets[ti][1] );
		string const repeatss( targets[ti][2] );
		string const filename( targets[ti][3] );
		Size const repeatlen( repeatseq.size() );

		TR.Trace << "target: " << I(4,repeatlen) << ' ' << repeatseq << ' ' << repeatss << ' ' << filename << endl;

		string const simfile( shared_output_tag() + filebase( filename )+".work" ), worktag("tmp" );
		Size const first_n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );

		if ( first_n > nstruct() ) continue;

		/// read pdb file
		Pose bound_pose;
		pose_from_pdb( bound_pose, filename );

		runtime_assert( bound_pose.sequence().substr(0,repeatlen) == repeatseq );

		///
		bool first_time_through( true );

		while ( true ) {
			Size const n( first_time_through ?
				first_n : get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( n > nstruct() ) break;
			first_time_through = false;


			Pose pose;
			clock_t starttime( clock() );
			Real final_score( refold_repeat_pose( repeatseq, repeatss, nrepeat, base_repeat, pose ) );
			clock_t stoptime( clock() );
			Real const centroid_simtime( ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

			/// do we want to fastrelax?
			Real relax_simtime( 0.0 );

			if ( option[ my_options::fastrelax_after_refold ] ) {
				ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
				runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );
				protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

				MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
				movemap->set_bb(true);
				movemap->set_chi(true);
				fastrelax.set_movemap( movemap );
				starttime = clock();
				if ( true || !dry_run() ) fastrelax.apply( pose );
				stoptime = clock();
				relax_simtime = ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

				final_score = (*fa_scorefxn)( pose );
			}


			Real rmsd( 0.0 ), rmsd_single_repeat( 0.0 );
			{
				using namespace core::id;
				AtomID_Map< AtomID > atom_map;
				initialize_atomid_map( atom_map, pose, id::GLOBAL_BOGUS_ATOM_ID );
				for ( Size i=1; i<= repeatlen; ++i ) {
					atom_map[ AtomID( pose.residue(i).atom_index("CA"), i ) ] = AtomID( bound_pose.residue(i).atom_index("CA"),i);
				}
				rmsd_single_repeat = rmsd_by_mapping( pose, bound_pose, atom_map );
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					atom_map[ AtomID( pose.residue(i).atom_index("CA"), i ) ] = AtomID( bound_pose.residue(i).atom_index("CA"),i);
				}
				rmsd = rmsd_by_mapping( pose, bound_pose, atom_map );
				//superimpose_pose( pose, bound_pose, atom_map );
			}

			Real handedness( 0.0 ),bound_pose_handedness( 0.0 );
			{
				Vectors coords;
				for ( Size i=1; i<= pose.total_residue(); ++i ) coords.push_back( pose.residue(i).xyz("CA") );
				handedness = get_chirality( repeatlen, coords );

				coords.clear();
				for ( Size i=1; i<= chain_end( 1, bound_pose ); ++i ) coords.push_back( bound_pose.residue(i).xyz("CA") );
				bound_pose_handedness = get_chirality( repeatlen, coords );
			}

			bool const passed_score_filter
				( append_score_to_scorefile_and_filter( worktag, rmsd, score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );

			string const outfilename( output_tag() + "refold_"+ filebase( filename )+"_N"+lead_zero_string_of(n,4)+".pdb" );

			cout << "final_scores " << F(9,3,final_score ) << ' ' << outfilename <<
				" rmsd: " << F(9,3,rmsd ) <<
				" rmsd_single_repeat: " << F(9,3,rmsd_single_repeat ) <<
				" handedness: " << F(9,3,handedness) <<
				" bound_pose_handedness: " << F(9,3,bound_pose_handedness) <<
				" repeatlen: " << repeatlen <<
				" repeatseq: " << repeatseq <<
				" repeatss: " << repeatss <<
				" centroid_simtime: " << F(9,3,centroid_simtime) <<
				" relax_simtime: " << F(9,3,relax_simtime) <<
				" passed_score_filter: " << passed_score_filter << endl;
			fflush( stdout );

			// dump pdb file?
			if ( passed_score_filter ) {
				pose.dump_pdb( outfilename );
			}
		} // nstruct
	} // targets


}

///////////////////////////////////////////////////////////////////////////////////////
void
simple_repack_and_minimize(
	bools const & is_flexible,
	ScoreFunctionOP scorefxn,
	Pose & pose
)
{
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	task->or_include_current( true );

	task->restrict_to_repacking();
	task->restrict_to_residues( is_flexible );

	Size const nloop( 25 );
	protocols::minimization_packing::PackRotamersMover packmover( scorefxn, task, nloop );

	packmover.apply( pose );


	MoveMapOP mm( new MoveMap );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( is_flexible[i] ) mm->set_chi( i, true );
	}

	protocols::minimization_packing::MinMoverOP min_mover
		( new protocols::minimization_packing::MinMover( mm, scorefxn, "dfpmin_armijo", 0.001, true ) );

	min_mover->apply( pose );

}

///////////////////////////////////////////////////////////////////////////////////////
void
talspecsym_test()
{
	option[ my_options::score_filter_acceptance_rate]; // check to see if present

	//ScoreFunctionOP fa_scorefxn( get_score_function_from_command_line() );
	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );

	/// specificity calculations using symmpose pdbs with nrepeat = 8
	///
	Size const repeatlen( 34 ), nrepeat( 3 ), // 3 copies of the RVD
		dna_context( 2 ), protein_context( 1 ),
		nbp( nrepeat + 2 * dna_context ), nrepeat_protein( nrepeat + 2*protein_context ),
		anchorpos( 13 ); // rvd pos 2

	strings files( start_files() );

	numeric::random::random_permutation( files, numeric::random::rg() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {

		Pose pose;
		pose_from_pdb( pose, files[fi] );

		pbassert( num_chains( pose ) == 3 );

		{ // trim if necessary
			Size const nbp_current( chain_end( pose,2) - chain_begin( pose,2) + 1 );
			if ( nbp_current > nbp ) {
				Size const ndel( nbp_current - nbp );
				pose.conformation().delete_residue_range_slow( chain_begin( pose, 3 ), chain_begin( pose,3)+ndel-1 );
				pose.conformation().delete_residue_range_slow( chain_end( pose, 2 )-ndel+1, chain_end( pose,2) );
				//set_DNA_terminus_variants_from_chains( pose );
			}

			pbassert( chain_end( pose,1)%repeatlen == 0 );
			pbassert( dna_context >= protein_context );
			Size const nrepeat_protein_current( chain_end( pose,1)/repeatlen );
			if ( nrepeat_protein_current > nrepeat_protein ) {
				Size const ndel( repeatlen * ( nrepeat_protein_current - nrepeat_protein ) );
				Size const ndel_nterm( repeatlen * ( dna_context - protein_context ) );
				Size const ndel_cterm( ndel - ndel_nterm );
				pose.conformation().delete_residue_range_slow( chain_end( pose,1) - ndel_cterm + 1, chain_end( pose,1) );
				pose.conformation().delete_residue_range_slow( chain_begin( pose,1), chain_begin( pose,1) + ndel_nterm -1 );
			}
			set_base_partner( pose );
		}
		pbassert( chain_end( pose,1)%repeatlen == 0 );
		pbassert( chain_end( pose,1)/repeatlen == nrepeat_protein );
		pbassert( chain_end( pose,2)-chain_begin( pose,2)+1 == nbp );

		vector1< std::pair< AA, AA > > workpairs;
		for ( Size aa_i=1; aa_i<= 20; ++aa_i ) {
			AA const aa = AA( aa_i );

			for ( Size na_i=first_DNA_aa; na_i<= last_DNA_aa; ++na_i ) {
				AA const na = AA( na_i );
				workpairs.push_back( make_pair( aa, na ) );
			}
		}
		pbassert( workpairs.size() == 80 );

		numeric::random::random_permutation( workpairs, numeric::random::rg() );

		for ( Size ii=1; ii<= workpairs.size(); ++ii ) {
			AA const aa( workpairs[ii].first ), na( workpairs[ii].second );

			string const simfile( shared_output_tag() + filebase( files[fi] ) +".work" );
			string const worktag( string_of( aa )+"_"+ string_of( na ) );


			Size const first_n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( first_n > nstruct() ) continue;

			/// mutate the rvd positions to aa; dna target site to na
			/// also setup flexibility bools
			Size const nres( pose.total_residue() ), nres_protein( chain_end( pose,1) );
			bools is_flexible( nres, false );

			char wtaa('X'), wtna('x');

			for ( Size i=0; i<= nrepeat+1; ++i ) {
				Size const ppos( (protein_context+i-1)*repeatlen + anchorpos );
				Size const dpos( nres_protein + dna_context + i ), dposp( retrieve_base_partner_from_pose( pose )[ dpos ] );

				is_flexible[ ppos ] = is_flexible[ ppos-1 ] = true;
				is_flexible[ dpos ] = is_flexible[ dposp ] = true;

				if ( i>=1 && i <= nrepeat ) { // make the mutations
					wtaa = pose.residue(ppos).name1();
					wtna = pose.residue(dpos).name1();
					make_sequence_change( ppos, aa, pose );
					devel::dna::make_base_pair_mutation( pose, dpos, na );
				}
			}

			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				TR.Trace << "is_flexible: " << I(4,i) << " " << is_flexible[i] << endl;
			}

			bool first_time_through( true );
			Pose const start_pose( pose );

			while ( true ) {
				Size const n( first_time_through ?
					first_n : get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
				if ( n > nstruct() ) break;
				first_time_through = false;

				pose = start_pose;

				clock_t starttime = clock();
				simple_repack_and_minimize( is_flexible, fa_scorefxn, pose );
				clock_t stoptime = clock();
				Real const simtime( ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ) );

				///
				Real const final_score( (*fa_scorefxn)( pose ) );

				Real unbound_rescored_protein_energy, unbound_relaxed_protein_energy, unbound_rescored_dna_energy(0.0),
					unbound_relaxed_dna_energy(0.0);
				//Pose unbound_relaxed_dna_pose;
				Pose protpose( pose );
				starttime = clock();
				{ // compute unbound energies



					protpose.conformation().delete_residue_range_slow( chain_begin( protpose,2), chain_end( protpose,3) );
					simple_repack_and_minimize( is_flexible, fa_scorefxn, pose );

					unbound_relaxed_protein_energy = (*fa_scorefxn)( protpose );
				}
				stoptime = clock();
				Real const unbound_simtime( ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ) );

				string const outfilename( output_tag() + "talsym_" + filebase( files[fi]) +"_" +
					worktag + "_N"+ lead_zero_string_of( n,4)+".pdb" );

				bool const passed_score_filter
					( append_score_to_scorefile_and_filter( worktag, final_score,
					option[ my_options::score_filter_acceptance_rate],
					option[ my_options::score_filter_pass_early], simfile ) );

				cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
					" passed_score_filter: " << passed_score_filter <<
					" aa: " << aa <<
					" na: " << na <<
					//      " inner_cycles: " << I(4,inner_cycles) <<
					//      " outer_cycles: " << I(4,outer_cycles) <<
					" unbound_rescored_protein_energy: " << F(9,3,unbound_rescored_protein_energy) <<
					" unbound_relaxed_protein_energy: " << F(9,3,unbound_relaxed_protein_energy) <<
					" unbound_rescored_dna_energy: " << F(9,3,unbound_rescored_dna_energy) <<
					" unbound_relaxed_dna_energy: " << F(9,3,unbound_relaxed_dna_energy) <<
					" wtaa: " << wtaa <<
					" wtna: " << wtna <<
					" simtime: " << F(9,3,simtime) <<
					" unbound_simtime: " << F(9,3,unbound_simtime) <<
					endl;

				fflush( stdout );

				if ( passed_score_filter ) {
					pose.dump_pdb( outfilename );
					//append_hbond_info_to_pdb_file( pose, *fa_scorefxn, outfilename );
					//unbound_relaxed_dna_pose.dump_pdb( outfilename+".unbound_relaxed_dna.pdb" );
					protpose.dump_pdb( outfilename+".unbound_relaxed_protein.pdb" );
				}

			} // nstruct
		} // aa
	} // files


}

///////////////////////////////////////////////////////////////////////////////
void
donut_deriv_test()
{
	runtime_assert( option[ basic::options::OptionKeys::out::file::output_virtual ] );
	Size const nrepeat( 7 ), repeatlen( 35 ), base_repeat( 3 );

	Pose pose;

	if ( false ) {

		pose_from_pdb( pose, start_file() );

	} else {
		for ( Size i=1; i<= nrepeat * repeatlen; ++i ) {
			pose.append_residue_by_bond( *get_vanilla_protein_residue( 'A' ), true ); // build_ideal_geometry
		}

	}


	Size const nres_protein( pose.total_residue() );
	runtime_assert( nres_protein == nrepeat * repeatlen );

	/// ADD A VIRTUAL RESIDUE AT THE END
	remove_upper_terminus_type_from_pose_residue( pose, pose.total_residue() );
	ResidueOP vrtrsd
		( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
	pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...
	pose.conformation().insert_chain_ending( nres_protein );

	// copy the torsions over from an internal repeat
	pose.set_psi  ( nres_protein  , pose.psi  ( repeatlen   ) );
	pose.set_omega( nres_protein  , pose.omega( repeatlen   ) );
	//pose.set_phi  ( nres_protein+1, pose.phi  ( repeatlen+1 ) );
	pose.set_torsion( id::TorsionID( nres_protein+1, id::BB, 1 ), pose.phi  ( repeatlen+1 ) );

	pose.dump_pdb( "start.pdb" );

	FoldTree f( pose.total_residue() );
	f.reorder( pose.total_residue() );
	pose.fold_tree( f );

	/// now need to setup some symm info

	conformation::symmetry::SymmetryInfo symminfo;
	setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo );


	pose::symmetry::make_symmetric_pose( pose, symminfo );

	Size const base_repeat_offset( ( base_repeat-1 )*repeatlen );

	for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) {
		pose.set_phi( i, pose.phi(i) );
		pose.set_psi( i, pose.psi(i) );
		pose.set_omega( i, pose.omega(i) );
	}


	pose.dump_pdb("symmetrized.pdb");


	/// now scoring
	scoring::symmetry::SymmetricScoreFunctionOP scorefxn( new scoring::symmetry::SymmetricScoreFunction() );
	//ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );

	runtime_assert( pose::symmetry::is_symmetric( *scorefxn ) );


	if ( false ) { // testing donut energy here

		DonutWholeEnergy donut_whole_energy;
		Donut1B_Energy donut_1b_energy;

		scorefxn->add_extra_method( donut_whole, 10.0, donut_whole_energy );
		scorefxn->add_extra_method( donut_1b, 10.0, donut_1b_energy );

	} else { // testing twist rise energy here

		Real const target_twist( 360.0 / nrepeat ), target_rise( 0.0 );
		TwistRiseEnergy    twist_rise_energy   ( target_twist, target_rise );
		TwistRise1B_Energy twist_rise_1b_energy( target_twist, target_rise );
		scorefxn->add_extra_method( twistrise   , 2.0, twist_rise_energy );
		scorefxn->add_extra_method( twistrise_1b, 2.0, twist_rise_1b_energy );

	}

	Real const score_premin( (*scorefxn)( pose ) );

	cout << "score_premin: " << F(9,3,score_premin) << ' ' <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	kinematics::MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
	movemap->set_bb(true);
	movemap->set_chi(false);
	movemap->set_jump(false);
	movemap->set_bb ( pose.total_residue(), false ); // silly that this is necessary

	protocols::minimization_packing::symmetry::SymMinMoverOP min_mover
		( new protocols::minimization_packing::symmetry::SymMinMover(movemap, scorefxn, "dfpmin", 0.00001, true,
		true, true ));

	min_mover->apply( pose );

	Real const score_postmin( (*scorefxn)( pose ) );

	cout << "score_postmin: " << F(9,3,score_postmin) << ' ' <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

	pose.dump_pdb("aftermin.pdb");

}



///////////////////////////////////////////////////////////////////////////////
void
unsatisfied_test()
{
	strings const files( start_files() );

	for ( Size fi=1; fi <= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		vector1< id::AtomID > buried_unsatisfied_acceptors, buried_unsatisfied_donors;
		bools subset( pose.total_residue(), true );
		subset[1] = subset[ chain_end( 1,pose ) ] = false;

		find_buried_unsatisfied_polars( subset, pose, buried_unsatisfied_donors, buried_unsatisfied_acceptors, files[fi] );

		{ // hacking
			Size repeatlen(5);
			string const seq( pose.sequence() );
			while ( seq.substr( 0, repeatlen ) != seq.substr( repeatlen, repeatlen ) ) {
				++repeatlen;
				if ( repeatlen > 50 ) utility_exit_with_message("repeatlen problem");
			}
			runtime_assert( chain_end( 1, pose )%repeatlen == 0 );
			Size const nrepeat( chain_end( 1,pose )/repeatlen );

			Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
				n_buried_unsatisfied_acceptors_sc;

			get_buried_unsatisfied_counts_real_slow( subset, pose,
				n_buried_unsatisfied_donors_bb,
				n_buried_unsatisfied_donors_sc,
				n_buried_unsatisfied_acceptors_bb,
				n_buried_unsatisfied_acceptors_sc, 1.0 );

			cout << "n_buried_unsatisfied10: " <<
				" donbb: " << n_buried_unsatisfied_donors_bb << F(6,2,Real(n_buried_unsatisfied_donors_bb)/nrepeat) <<
				" donsc: " << n_buried_unsatisfied_donors_sc << F(6,2,Real(n_buried_unsatisfied_donors_sc)/nrepeat) <<
				" accbb: " << n_buried_unsatisfied_acceptors_bb << F(6,2,Real(n_buried_unsatisfied_acceptors_bb)/nrepeat) <<
				" accsc: " << n_buried_unsatisfied_acceptors_sc << F(6,2,Real(n_buried_unsatisfied_acceptors_sc)/nrepeat) <<
				' ' << nrepeat << ' ' << repeatlen << ' ' << filebase( files[fi] ) << endl;

			get_buried_unsatisfied_counts_real_slow( subset, pose,
				n_buried_unsatisfied_donors_bb,
				n_buried_unsatisfied_donors_sc,
				n_buried_unsatisfied_acceptors_bb,
				n_buried_unsatisfied_acceptors_sc, 1.2 );

			cout << "n_buried_unsatisfied12: " <<
				" donbb: " << n_buried_unsatisfied_donors_bb << F(6,2,Real(n_buried_unsatisfied_donors_bb)/nrepeat) <<
				" donsc: " << n_buried_unsatisfied_donors_sc << F(6,2,Real(n_buried_unsatisfied_donors_sc)/nrepeat) <<
				" accbb: " << n_buried_unsatisfied_acceptors_bb << F(6,2,Real(n_buried_unsatisfied_acceptors_bb)/nrepeat) <<
				" accsc: " << n_buried_unsatisfied_acceptors_sc << F(6,2,Real(n_buried_unsatisfied_acceptors_sc)/nrepeat) <<
				' ' << nrepeat << ' ' << repeatlen << ' ' << filebase( files[fi] ) << endl;

			get_buried_unsatisfied_counts_real_slow( subset, pose,
				n_buried_unsatisfied_donors_bb,
				n_buried_unsatisfied_donors_sc,
				n_buried_unsatisfied_acceptors_bb,
				n_buried_unsatisfied_acceptors_sc, 1.4 );

			cout << "n_buried_unsatisfied14: " <<
				" donbb: " << n_buried_unsatisfied_donors_bb << F(6,2,Real(n_buried_unsatisfied_donors_bb)/nrepeat) <<
				" donsc: " << n_buried_unsatisfied_donors_sc << F(6,2,Real(n_buried_unsatisfied_donors_sc)/nrepeat) <<
				" accbb: " << n_buried_unsatisfied_acceptors_bb << F(6,2,Real(n_buried_unsatisfied_acceptors_bb)/nrepeat) <<
				" accsc: " << n_buried_unsatisfied_acceptors_sc << F(6,2,Real(n_buried_unsatisfied_acceptors_sc)/nrepeat) <<
				' ' << nrepeat << ' ' << repeatlen << ' ' << filebase( files[fi] ) << endl;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
symdock_dna_test()
{
	using kinematics::Stub;

	ScoreFunction censcorefxn;
	censcorefxn.set_weight( vdw, 1.0 );

	Size const base_repeat(3);

	/// read the logfile
	string const logfile( start_file() );

	ifstream data( logfile.c_str() );

	Pose dnapose;
	cenpose_from_pdb( dnapose, "symmpose_01.pdb" );
	dnapose.conformation().delete_residue_range_slow( chain_begin( 1, dnapose ), chain_end( 1, dnapose ) );

	Vector dna_axis, dna_center;
	Real dna_rise;
	{
		Size nb( dnapose.total_residue()/2 );
		Stub const stub1( scoring::dna::get_base_pair_stub_slow( dnapose.residue(1), dnapose.residue(2*nb) ) );
		Stub const stub2( scoring::dna::get_base_pair_stub_slow( dnapose.residue(2), dnapose.residue(2*nb-1) ) );
		Vector t;
		Real theta;
		get_stub_transform_data( stub1, stub2, dna_center, dna_axis, t, theta );
		dna_rise = t.dot( dna_axis );
		cout << "dna_twist: " << numeric::conversions::degrees( theta ) << " dna_rise: " << dna_rise << endl;

		/// try out adding some more residues
		Pose movpose( dnapose );
		numeric::xyzMatrix< Real > R( numeric::rotation_matrix( dna_axis, theta*nb ) );
		Vector const v( dna_center + t * nb - R*dna_center );
		movpose.apply_transform_Rx_plus_v( R, v );

		for ( Size i=1; i<= nb; ++i ) {
			dnapose.append_polymer_residue_after_seqpos( movpose.residue(i), chain_end( 1, dnapose ), false );
		}
		for ( Size i=0; i< nb; ++i ) {
			dnapose.prepend_polymer_residue_before_seqpos( movpose.residue(movpose.total_residue()-i),
				chain_begin( 2, dnapose ), false );
		}
		nb*=2;


		Size const midbp( (nb-1)/2+1 );
		Stub const midstub( scoring::dna::get_base_pair_stub_slow( dnapose.residue(midbp), dnapose.residue(2*nb-midbp+1) ) );
		dna_center = dna_center + ( midstub.v - dna_center ).dot( dna_axis ) * dna_axis;
		runtime_assert( is_small( dna_axis.dot( dna_center - midstub.v ) ) );
	}



	string line;
	Size counter(0),filecounter(0);
	while ( getline( data, line ) ) {
		++counter;
		if ( counter%10000 == 0 ) cerr << counter << ' ' << filecounter << endl;

		//istringstream l( line );
		strings const l( split_to_vector1( line ) );
		string const dirtag( l[1] ), filename( dirtag.substr( 0, dirtag.find_last_of( "/" ) ) +"/" + l[3] );

		bool passed_score_filter( false ), relaxed( true );
		Size repeatlen(0), nrepeat(0);
		string repeatss, turn1, turn2, helix1_len, helix2_len;


		for ( Size i=1; i< l.size(); ++i ) {
			if      ( l[i] == "passed_score_filter:" ) passed_score_filter = ( l[i+1] == "1" );
			else if ( l[i] == "passed_fullatom_score_filter:" ) passed_score_filter = ( l[i+1] == "1" );
			else if ( l[i] == "repeatlen:" ) repeatlen = int_of( l[i+1] );
			else if ( l[i] == "nrepeat:" ) nrepeat = int_of( l[i+1] );
			else if ( l[i] == "relaxed:" ) relaxed = ( l[i+1] == "1" );
			else if ( l[i] == "repeatss:" ) repeatss = l[i+1];
			else if ( l[i] == "helix1_len:" ) helix1_len = int_of( l[i+1] );
			else if ( l[i] == "helix2_len:" ) helix2_len = int_of( l[i+1] );
			else if ( l[i] == "turn1:" ) turn1 = l[i+1];
			else if ( l[i] == "turn2:" ) turn2 = l[i+1];
		}

		if ( !relaxed ) passed_score_filter = false;

		if ( !passed_score_filter ) continue;

		runtime_assert( repeatlen * nrepeat > 0 );

		if ( !utility::file::file_exists( filename ) ) {
			cout << "MISSING file: " << filename << endl;
			continue;
		}

		Size const base_repeat_offset( (base_repeat-1)*repeatlen );

		cerr << "read " << filename << endl;
		++filecounter;

		Pose pose;
		cenpose_from_pdb( pose, filename );

		runtime_assert( num_chains( pose ) == 1 );
		runtime_assert( chain_end( 1, pose ) == nrepeat * repeatlen );


		//// try docking symmetrical tal dna into this guy
		/// need to figure out the symmetry axis of both the protein and the dna, and align them

		Vector symmetry_center(0,0,0), symmetry_axis(0,0,0);
		Real twist_radians;
		{
			Size const rpos1( 1 ), rpos2( repeatlen/3 ), rpos3( (2*repeatlen)/3 );
			for ( Size i=1; i<nrepeat; ++i ) {
				Size const off1( (i-1)*repeatlen ), off2( i*repeatlen );
				Stub const stub1( pose.residue( off1+rpos1 ).xyz("CA"),
					pose.residue( off1+rpos2 ).xyz("CA"),
					pose.residue( off1+rpos3 ).xyz("CA") );
				Stub const stub2( pose.residue( off2+rpos1 ).xyz("CA"),
					pose.residue( off2+rpos2 ).xyz("CA"),
					pose.residue( off2+rpos3 ).xyz("CA") );
				Vector center,n,t;
				get_stub_transform_data( stub1, stub2, center, n, t, twist_radians );
				if ( i == 1 ) {
					symmetry_center = center;
					symmetry_axis = n;

					/// shift symmetry_center so that it aligns with closest rsd in base repeat
					Real minradius( 1e6 );
					Size minpos(0);
					for ( Size i=1; i<= repeatlen; ++i ) {
						Vector const caxyz( pose.residue(base_repeat_offset+i).xyz("CA") );
						Vector radius( caxyz - center );
						radius -= n.dot( radius ) * n;
						runtime_assert( is_small( radius.dot(n) ) );
						if ( radius.length() < minradius ) {
							minradius = radius.length();
							minpos = base_repeat_offset+i;
						}
					}

					Vector const caxyz( pose.residue(minpos).xyz("CA") );
					symmetry_center = symmetry_center + ( caxyz - symmetry_center ).dot(n) * n;
					runtime_assert( is_small( n.dot( symmetry_center - caxyz ) ) );

				} else {
					cout << "symdev: "<< F(9,3, ( center - symmetry_center ).normalized().cross( symmetry_axis ).length() ) <<
						F(9,3,numeric::conversions::degrees( acos( symmetry_axis.dot(n) )) ) << endl;
				}
			}
		}

		//pose.dump_pdb("withoutdna.pdb");


		// translate dna so that dna axis aligns with symmetry_axis, dna center aligns with symmetry center
		Size const n_angle_steps( 12 ), n_rise_steps( 3 );

		//Real minvdw( 1e6 );

		Reals dna_vdws, dna_dis2s;

		for ( Size i_angle=0; i_angle< n_angle_steps; ++i_angle ) {
			for ( Size i_rise=0; i_rise< n_rise_steps; ++i_rise ) {


				Pose movpose( dnapose );

				/// twist/rise about original dna axes
				{
					Real const twist_angle( i_angle * numeric::constants::d::pi * 2.0 / Real ( n_angle_steps ) );
					Real const rise_step( i_rise * dna_rise / Real( n_rise_steps ) );
					numeric::xyzMatrix< Real > R( numeric::rotation_matrix( dna_axis, twist_angle ) );
					Vector const v( dna_center + rise_step * dna_axis - R*dna_center );
					movpose.apply_transform_Rx_plus_v( R, v );
				}


				Real const rotangle( acos( dna_axis.dot( symmetry_axis ) ) );
				Vector const rotaxis( dna_axis.cross( symmetry_axis ) );
				numeric::xyzMatrix< Real > R( numeric::rotation_matrix( rotaxis, rotangle ) );
				runtime_assert( is_small( symmetry_axis.distance_squared( R*dna_axis ) ) );
				Vector const v( symmetry_center - R*dna_center );
				runtime_assert( is_small( symmetry_center.distance_squared( R*dna_center + v ) ) );

				movpose.apply_transform_Rx_plus_v( R, v );

				Pose bigpose( pose );

				for ( Size i=1; i<= movpose.total_residue(); ++i ) {
					Residue const & rsd( movpose.residue(i) );
					if ( rsd.is_lower_terminus() ) {
						bigpose.append_residue_by_jump( rsd, 1 );
						bigpose.conformation().insert_chain_ending( bigpose.total_residue()-1 );
					} else {
						bigpose.append_residue_by_bond( rsd );
					}
				}

				Real const vdw( censcorefxn( bigpose ) );
				TR.Trace << "vdw: " << I(3,i_angle) << I(3,i_rise) << F(9,3,vdw) << endl;


				// check distances to dna
				Real mindis2( 1e6 );
				for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) {
					for ( Size j=chain_begin( 2, bigpose ); j<= chain_end( 3, bigpose ); ++j ) {
						Residue const & irsd( bigpose.residue(i)), &jrsd( bigpose.residue(j) );
						Vector jxyz;
						if ( jrsd.aa() == na_cyt ) jxyz = jrsd.xyz("C5");
						else if ( jrsd.aa() == na_gua ) jxyz = jrsd.xyz("N7");
						else utility_exit_with_message("bad jrsd "+jrsd.name());
						mindis2 = min( mindis2, irsd.xyz("CA").distance_squared( jxyz ) );
					}
				}

				/// add to lists
				dna_vdws.push_back( vdw );
				dna_dis2s.push_back( mindis2 );

				//bigpose.dump_pdb("withdna_"+lead_zero_string_of( i_angle,2 )+"_"+lead_zero_string_of( i_rise,2 )+".pdb");
			} // i_rise
		} // i_angle

		Real const minvdw( min( dna_vdws ) );
		Real mindis( 1e6 );
		for ( Size i=1; i<= dna_vdws.size(); ++i ) {
			if ( dna_vdws[i] < minvdw+0.2 ) {
				mindis = min( mindis, sqrt( dna_dis2s[i] ) );
			}
		}

		Real const protein_vdw( censcorefxn( pose ) ), dna_only_vdw( censcorefxn( dnapose ) );

		cout << line <<
			" min_dna_vdw: " << F(9,3,minvdw-protein_vdw-dna_only_vdw) <<
			" min_dna_dis: " << F(9,3,mindis ) <<
			endl;

		// superimpose pose onto itself, shifted by 1 repeat -- this is more robust than using
		//   {
		//    Pose fixpose( pose ), movpose( pose );
		//    using namespace core::id;
		//    AtomID_Map< AtomID > atom_map;
		//    initialize_atomid_map( atom_map, movpose, id::GLOBAL_BOGUS_ATOM_ID );

		//    for ( Size i=1; i<nrepeat; ++i ) {
		//     for ( Size j=1; j<= repeatlen; ++j ) {
		//      Size const pos1( (i-1)*repeatlen + j ), pos2( i*repeatlen + j );
		//      atom_map[ AtomID( movpose.residue( pos1 ).atom_index("CA"), pos1 ) ] =
		//       AtomID( fixpose.residue( pos2 ).atom_index("CA"), pos2 );
		//     }
		//    }

		//    superimpose_pose( movpose, fixpose, atom_map );
	}
	data.close();


}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
shift_test()
{
	Size const repeatlen( 35 );
	strings const files( start_files() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		runtime_assert( pose.total_residue()%repeatlen == 0 );

		remove_lower_terminus_type_from_pose_residue( pose, 1 );
		remove_upper_terminus_type_from_pose_residue( pose, pose.total_residue() );
		pose.prepend_polymer_residue_before_seqpos( pose.residue( pose.total_residue()   ), 1, true );
		pose.prepend_polymer_residue_before_seqpos( pose.residue( pose.total_residue()-1 ), 1, true );
		for ( Size i=1; i<= 3; ++i ) { // need to include phi of rsd 3, at the very least
			pose.set_phi  ( i, pose.phi  ( i+repeatlen ) );
			pose.set_psi  ( i, pose.psi  ( i+repeatlen ) );
			pose.set_omega( i, pose.omega( i+repeatlen ) );
		}

		pose.conformation().delete_residue_range_slow( pose.total_residue()-1, pose.total_residue() );
		add_lower_terminus_type_to_pose_residue( pose, 1);
		add_upper_terminus_type_to_pose_residue( pose,  pose.total_residue() );
		pose.dump_pdb(files[fi]+"_shift2.pdb");
	}

}



///////////////////////////////////////////////////////////////////////////////

void
spin_test()
{
	Size const nrepeat(6), repeatlen(35);

	Pose pose;
	pose_from_pdb( pose, start_file() );

	Pose const start_pose( pose );

	for ( Size shift=1; shift<nrepeat; ++shift ) {
		using namespace id;
		pose = start_pose;
		// create mapping from pose to start_pose
		AtomID_Map< AtomID > atom_map;
		initialize_atomid_map( atom_map, pose, id::GLOBAL_BOGUS_ATOM_ID );
		for ( Size i=1; i<= nrepeat*repeatlen; ++i ) {
			Size const irep( (i-1)/repeatlen+1 ), pos( (i-1)%repeatlen+1 ), jrep( (irep+shift-1)%nrepeat+1 ),
				j( (jrep-1)*repeatlen+pos);
			runtime_assert( (irep-1)*repeatlen+pos == i );

			atom_map[ AtomID( pose.residue(i).atom_index("CA"), i ) ] = AtomID( start_pose.residue(j).atom_index("CA"),j);
		}
		Real const rmsd( rmsd_by_mapping( pose, start_pose, atom_map ) );
		string const outfilename( filebase( start_file())+"_spin"+string_of(shift)+".pdb");
		cout << "spin_rmsd: " << F(9,3,rmsd) << ' ' << outfilename << endl;
		superimpose_pose( pose, start_pose, atom_map );
		pose.dump_pdb( outfilename );

	}

}

#ifdef ABBA
///////////////////////////////////////////////////////////////////////////////
// try out various a-b-b-a structures
//
void
abba_test()
{
	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );
	Size const design_cycles( option[ my_options::design_cycles ] );

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	string const sses( "HLELELH");
	Size const num_sses( sses.size() );

	vector1< Sizes > const sse_size_ranges
		( make_vector1( make_vector1( 7, 8 ),
										make_vector1( 2, 3, 4 ),
										make_vector1( 4, 5, 6 ),
										make_vector1( 2 ),
										make_vector1( 4, 5, 6 ),
										make_vector1( 2, 3, 4 ),
										make_vector1( 7, 8 ) ) );
	runtime_assert( num_sses == sse_size_ranges.size() );


	string const simfile( shared_output_tag()+"abba.work" );

	while ( true ) {
		// setup the pose
		Sizes sse_sizes;
		string ss;
		for ( Size i=1; i<= num_sses; ++i ) {
			sse_sizes.push_back( random_element( sse_size_ranges[ i ] ) );
			ss += string( sses[i-1], sse_sizes.back() );
		}

		// create the pose
		Pose pose;
		Size const nres( ss.size() );
		Sizes fragseq_poslist;// positions that can vary during centroid rebuild, using fragment sequence
		for ( Size i=1; i<= nres; ++i ) {
			char const iss( ss[i-1] );
			char name1;
			if ( iss == 'H' ) name1 = random_element( make_vector1( 'L','L','L','I','I','V' ) );
			else if ( iss == 'E' ) name1 = random_element( make_vector1( 'I','I','V','V' ) );
			else name1 = 'G';
			pose.append_residue_by_bond( *get_vanilla_protein_residue( name1 ), true ); // build_ideal_geometry
			pose.set_secstruct( i, iss );
		}


		// pick frags


		// try frag-building simulation

		// fastrelax/design

		// dump pdb


	}


	signal_that_job_is_done();

}

#endif

///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Matrix
// random_reorientation_matrix()
// {
//  // I got this method from http://www.mech.utah.edu/~brannon/public/rotation.pdf
//  // They dont have a proof, but they did some simulation tests...
//  Vector x( random_unit_vector() ), y( random_unit_vector() );
//  // if x and y are parallel this wont work, but that seems unlikely
//  y = ( y - y.dot(x) * x ).normalized();
//  Vector const z( x.cross(y) );

//  return Matrix::cols( x, y, z );
// }

void
get_interface_between_positions(
	Sizes const & partner1,
	Sizes const & partner2,
	Pose const & pose,
	Sizes & partner1_interface,
	Sizes & partner2_interface
)
{
	partner1_interface.clear();
	partner2_interface.clear();

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( Sizes::const_iterator pos1= partner1.begin(); pos1 != partner1.end(); ++pos1 ) {
		for ( utility::graph::Graph::EdgeListConstIter
				ir  = energy_graph.get_node( *pos1 )->const_edge_list_begin(),
				ire = energy_graph.get_node( *pos1 )->const_edge_list_end();
				ir != ire; ++ir ) {
			Size const pos2( (*ir)->get_other_ind( *pos1 ) );
			if ( has_element( partner2, pos2 ) ) {
				if ( !has_element( partner1_interface, *pos1 ) ) partner1_interface.push_back( *pos1 );
				if ( !has_element( partner2_interface,  pos2 ) ) partner2_interface.push_back(  pos2 );
			}
		}
	}
	TR.Trace << "get_interface_between_positions: " << partner1.size() << ' ' << partner2.size() << ' ' <<
		partner1_interface.size() << ' ' << partner2_interface.size() << endl;

}


Vector
get_centroid_of_positions(
	Pose const & pose,
	Sizes const & poslist
)
{
	Vector centroid(0.0);
	for ( Sizes::const_iterator pos=poslist.begin(); pos != poslist.end(); ++pos ) {
		centroid += pose.residue( *pos ).nbr_atom_xyz();
	}
	centroid /= poslist.size();
	return centroid;
}


///////////////////////////////////////////////////////////////////////////////////////

void
set_jump_rb_centers( Pose & pose )
{
	runtime_assert( num_chains( pose ) == 2 );
	runtime_assert( pose.fold_tree().num_jump() == 1 );

	/// sets centers for rotation in minimization and gaussian moves
	/// figure out the interface residues on both sides
	/// take centroid of the interface residues
	/// if no interface residues on either side, take centroid of entire protein
	/// needs a prior scorefxn call to ensure that the nbr graph is up to date


	Sizes partner1, partner2, partner1_interface, partner2_interface;
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.chain(i) == 1 ) partner1.push_back(i);
		else if ( pose.chain(i) == 2 ) partner2.push_back(i);
	}

	get_interface_between_positions( partner1, partner2, pose, partner1_interface, partner2_interface );

	Vector const centroid1( partner1_interface.empty() ? get_centroid_of_positions( pose, partner1 ):
		get_centroid_of_positions( pose, partner1_interface ) );
	Vector const centroid2( partner2_interface.empty() ? get_centroid_of_positions( pose, partner2 ):
		get_centroid_of_positions( pose, partner2_interface ) );


	Size const jump_number( 1 );
	pbassert( has_element( partner1, Size( pose.fold_tree().upstream_jump_residue( jump_number ) ) ) &&
		has_element( partner2, Size( pose.fold_tree().downstream_jump_residue( jump_number ) ) ) );

	kinematics::Jump j( pose.jump(jump_number) );
	j.set_rb_center(  1, pose.conformation().downstream_jump_stub( jump_number ), centroid2 );
	j.set_rb_center( -1, pose.conformation().  upstream_jump_stub( jump_number ), centroid1 );
	pose.set_jump( jump_number, j );

}

///////////////////////////////////////////////////////////////////////////////

void
setup_interface_movemap(
	Pose const & pose,
	kinematics::MoveMap & mm
)
{

	bools is_partner1( pose.total_residue(), false ), is_partner2( pose.total_residue(), false );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		is_partner1[i] = ( pose.chain(i) == 1 );
		is_partner2[i] = ( pose.chain(i) == 2 );
	}
	bools nbrs1, nbrs2;
	find_neighbors( is_partner1, pose, nbrs1 );
	find_neighbors( is_partner2, pose, nbrs2 );


	mm.set_bb( false );
	mm.set_chi( false );
	Size int1_count(0), int2_count(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( ( is_partner1[i] && nbrs2[i] ) || ( is_partner2[i] && nbrs1[i] ) ) {
			TR.Trace << "flexpos: " << I(4,i) << endl;
			mm.set_chi(i, true );
			int1_count += is_partner1[i];
			int2_count += is_partner2[i];
		}
	}
	TR.Trace << "numflex: " << int1_count << ' ' << int2_count << endl;
	mm.set_jump( 1, true );

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
setup_cage_design_movemap(
	Pose const & pose,
	kinematics::MoveMap & mm
)
{

	bools is_partner1( pose.total_residue(), false ), is_partner2( pose.total_residue(), false );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		is_partner1[i] = ( pose.chain(i) == 1 );
		is_partner2[i] = ( pose.chain(i) == 2 );
	}
	bools nbrs1, nbrs2;
	find_neighbors( is_partner1, pose, nbrs1 );
	find_neighbors( is_partner2, pose, nbrs2 );


	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );
	Size int1_count(0), int2_count(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( is_partner1[i] && nbrs2[i] ) { // || ( is_partner2[i] && nbrs1[i] ) ) {
			TR.Trace << "flexpos: " << I(4,i) << endl;
			mm.set_chi(i, true );
			int1_count += is_partner1[i];
			int2_count += is_partner2[i];
		}
	}
	TR.Trace << "numflex: " << int1_count << ' ' << int2_count << endl;
	runtime_assert( pose.fold_tree().num_jump() == 2 + 2 + 1 );
	runtime_assert( num_chains( pose ) == 3 );
	mm.set_jump(3,true); // jump from central to axial vrtrsd for base repeat


}

///////////////////////////////////////////////////////////////////////////////

void
get_toroid_axis(
	Pose const & pose,
	Size const nrepeat,
	Size const repeatlen,
	Vector & axis,
	Vector & center
)
{
	using namespace id;
	Pose fix_pose( pose ), mov_pose( pose );

	AtomID_Map< AtomID > atom_map;

	initialize_atomid_map( atom_map, mov_pose, id::GLOBAL_BOGUS_ATOM_ID );

	Size const shift( (nrepeat/2) * repeatlen ), nres_monomer( nrepeat * repeatlen );

	Vector centroid(0,0,0);
	for ( Size i=1; i<= nres_monomer; ++i ) {
		Size const j( (i+shift-1)%nres_monomer+1 );
		atom_map[ AtomID( mov_pose.residue(i).atom_index("CA"), i ) ] =
			AtomID( fix_pose.residue(j).atom_index("CA"),j);
		centroid += fix_pose.residue(i).xyz("CA");
	}
	centroid /= nres_monomer;
	Real const rmsd( rmsd_by_mapping( mov_pose, fix_pose, atom_map ) );
	TR.Trace << "get_toroid_axis: spin_rmsd= " << F(9,3,rmsd) << endl;
	superimpose_pose( mov_pose, fix_pose, atom_map );

	Size const rpos( 5 );
	Stub const stub1( fix_pose.residue( rpos     ).xyz("CA"),
		fix_pose.residue( rpos + 1 ).xyz("CA"),
		fix_pose.residue( rpos + 2 ).xyz("CA") );
	Stub const stub2( mov_pose.residue( rpos     ).xyz("CA"),
		mov_pose.residue( rpos + 1 ).xyz("CA"),
		mov_pose.residue( rpos + 2 ).xyz("CA") );
	Real theta;
	Vector t;
	get_stub_transform_data( stub1, stub2, center, axis, t, theta );
	runtime_assert( fabs( axis.length_squared() - 1.0 )<1e-2 ); // confirm that axis is normal vector

	/// place center at projection of toroid centroid onto axis
	Vector old_center( center ); // for debugging
	center += axis * ( centroid-center ).dot(axis);
	runtime_assert( fabs( (centroid - center).dot(axis) )<1e-2 );
	runtime_assert( ( center - old_center ).cross( axis ).length_squared() < 1e-2 );

	TR.Trace << "get_toroid_axis: theta: " << numeric::conversions::degrees( theta ) << endl;
}


void
reconstruct_planar_2c_test()
{

	Pose dimer_pose;
	pose_from_pdb( dimer_pose, start_file() );
	devel::blab::loops::switch_pose_to_residue_type_set( dimer_pose, CENTROID );

	runtime_assert( num_chains( dimer_pose ) == 2 );
	Pose pose2( dimer_pose );
	pose2.conformation().delete_residue_range_slow( 1, chain_end(1,dimer_pose) );
	Size nrepeat2( 6 ), repeatlen2( pose2.total_residue()/nrepeat2 );
	runtime_assert( pose2.total_residue()%nrepeat2 == 0 );
	Vector axis2, center2;
	get_toroid_axis( pose2, nrepeat2, repeatlen2, axis2, center2 );

	Real const axis_distance( fabs( center2.y() ) );
	cout << "small? " << axis2.x() << ' ' << axis2.y() << ' ' << center2.x() << endl;

	Pose pose( dimer_pose );

	for ( int k=-2; k<= 2; ++k ) { // layer in the x direction
		for ( int kk=0; kk<= 1; ++kk ) { // shift in the y-direction
			Vector const translation( axis_distance * Vector( Real(k) * sqrt(3.0)/2,
				3.0 * (Real(kk) + fabs(Real(abs(k)%2))/2),
				0.0 ) );
			if ( translation.length()<1 ) continue;
			numeric::xyzMatrix_double R( numeric::xyzMatrix_double::identity() );
			Pose dpose( dimer_pose );
			dpose.apply_transform_Rx_plus_v( R, translation );
			for ( Size i=1; i<= dpose.total_residue(); ++i ) {
				if ( dpose.residue(i).is_lower_terminus() ) {
					pose.append_residue_by_jump( dpose.residue(i), 1 );
					pose.conformation().insert_chain_ending( pose.total_residue()-1 );
				} else {
					pose.append_residue_by_bond( dpose.residue(i) );
				}
			}
		}
	}
	pose.dump_pdb("planar_2c.pdb");
}


///////////////////////////////////////////////////////////////////////////////
//
// in two_component_test, we use two central virtual residues, both with x-axis along the x-axis (i),
// one with y=j and z=k axes also standard, and one with a rotation of axis_angle about the x-axis
// so y =  cos(axis_angle) j + sin(axis_angle) k
//    z = -sin(axis_angle) j + cos(axis_angle) k
//
void
reconstruct_2c_test()
{

	string const symm_type( option[ my_options::symm_type ] );//"I2" );

	Vectors vertices;
	vector1< Sizes > nbr_vertices;
	get_symm_type_vertices( symm_type[0], vertices, nbr_vertices );
	//axis_angle = 0.5 * std::acos( vertices[1].normalized().dot( nbr_vertices[1].normalized() ) );

	Pose dimer_pose;
	pose_from_pdb( dimer_pose, start_file() );
	devel::blab::loops::switch_pose_to_residue_type_set( dimer_pose, CENTROID );

	// let's make #(vertices) copies of pose
	//

	Pose pose;
	for ( Size i=1; i<= vertices.size(); ++i ) {
		Vector const vertex( vertices[i] );
		// rotate dimer pose so that its z-axis is parallel to the vector from the origin to the vertex
		for ( Size j=1; j<= nbr_vertices[i].size(); ++j ) {
			Size const nbr_vertex_index( nbr_vertices[i][j] );
			runtime_assert( nbr_vertex_index<= vertices.size() );
			Vector const nbr_vertex( vertices[ nbr_vertex_index ] );

			Vector z( vertex.normalized() ), y( vertex - nbr_vertex ); // y( nbr_vertices[i] - vertices[i] );
			y = ( y - z.dot( y ) * z ).normalized();
			Vector const x( y.cross(z) );

			numeric::xyzMatrix_double R( numeric::xyzMatrix_double::cols( x, y, z ) );
			Pose dpose( dimer_pose );
			dpose.apply_transform_Rx_plus_v( R, Vector(0,0,0) );
			for ( Size k=1; k<= dpose.total_residue(); ++k ) {
				if ( dpose.chain(k) == 1 && j > 1 ) continue;
				if ( dpose.chain(k) == 2 && nbr_vertex_index < i ) continue;
				if ( dpose.residue(k).is_lower_terminus() && pose.total_residue() ) {
					pose.append_residue_by_jump( dpose.residue(k), 1 );
					pose.conformation().insert_chain_ending( pose.total_residue()-1 );
				} else {
					pose.append_residue_by_bond( dpose.residue(k) );
				}
			}
		}
	}
	pose.dump_pdb(filebase( start_file() )+"_rec2c.pdb" );


}




///////////////////////////////////////////////////////////////////////////////
void
cage_test()
{
	using namespace kinematics;

	string const symm_type( option[ my_options::symm_type ] );

	/// try making some symmetrical structures

	///
	Pose monomer_pose;
	pose_from_pdb( monomer_pose, start_file() );

	// figure out where the symmetry axis is, this assumes perfect symmetry
	Vector center, axis;
	{
		Size const rpos(5), repeatlen( 35 ), nrepeat( 6 ); // rpos is arbitrary
		runtime_assert( monomer_pose.total_residue() == nrepeat * repeatlen );
		Pose const & pose( monomer_pose ); // doh
		Stub const stub1( pose.residue( rpos     ).xyz("CA"),
			pose.residue( rpos + 1 ).xyz("CA"),
			pose.residue( rpos + 2 ).xyz("CA") );
		Stub const stub2( pose.residue( rpos + repeatlen    ).xyz("CA"),
			pose.residue( rpos + repeatlen + 1 ).xyz("CA"),
			pose.residue( rpos + repeatlen + 2 ).xyz("CA") );
		Real theta;
		Vector t;
		get_stub_transform_data( stub1, stub2, center, axis, t, theta );
		TR.Trace << "theta: " << numeric::conversions::degrees( theta ) << endl;
	}

	Size const nres_monomer( monomer_pose.total_residue() );

	Pose pose;
	{


		//// get the vertices:
		Vectors vertices;

		if ( symm_type == "O" ) {
			for ( int i=-1; i<= 1; i+= 2 ) {
				for ( int j=-1; j<= 1; j+= 2 ) {
					for ( int k=-1; k<= 1; k+= 2 ) {
						vertices.push_back( Vector( i, j, k ) );
					}
				}
			}

		} else if ( symm_type == "T" ) {
			for ( int i=-1; i<= 1; i+= 2 ) {
				vertices.push_back( Vector( i, 0, -1.0 / sqrt(2) ) );
				vertices.push_back( Vector( 0, i,  1.0 / sqrt(2) ) );
			}

		} else if ( symm_type == "I" ) {
			/// vertices of cube:
			for ( int i=-1; i<= 1; i+= 2 ) {
				for ( int j=-1; j<= 1; j+= 2 ) {
					for ( int k=-1; k<= 1; k+= 2 ) {
						vertices.push_back( Vector( i, j, k ) );
					}
				}
			}
			Real const golden( 0.5 * ( 1 + sqrt(5) ) ), inv_golden( 1.0/ golden );
			for ( int i=-1; i<= 1; i+= 2 ) {
				for ( int j=-1; j<= 1; j+= 2 ) {
					vertices.push_back( Vector( 0, i * inv_golden, j * golden ) );
					vertices.push_back( Vector( i * inv_golden, j * golden, 0 ) );
					vertices.push_back( Vector( i * golden, 0, j * inv_golden ) );
				}
			}
		} else {
			utility_exit_with_message("unrecognized symm_type: "+symm_type );
		}

		Size const n_monomers( vertices.size() );

		/// get nbr relationships among vertices
		Vectors nbr_vertices( n_monomers );
		Size nbr_count1(0);
		for ( Size i=1; i<= n_monomers; ++i ) {
			Vector const & vi( vertices[i] );
			Real mindis2( 1e6 ), epsilon( 1e-2 );
			Size nbr_count(0);
			for ( Size j=1; j<= n_monomers; ++j ) {
				if ( j==i ) continue;
				Vector const & vj( vertices[j] );
				Real const dis2( vi.distance_squared(vj) );
				if ( dis2<mindis2-epsilon ) {
					mindis2 = dis2;
					nbr_vertices[i] = vj;
					nbr_count = 1;
				} else if ( dis2<mindis2+epsilon ) {
					++nbr_count;
				}
			}
			TR.Trace << "nbr_count: " << i << ' ' << nbr_count << F(9,3,sqrt(mindis2)) << endl;
			if ( i==1 ) {
				nbr_count1 = nbr_count;
			} else {
				runtime_assert( nbr_count == nbr_count1 );
			}
		}

		runtime_assert( nbr_count1 == 3 ); // 3-fold symmetries right now


		ResidueOPs vrt_rsds;

		// create a virtual residue in the middle of the monomer, with z-axis parallel to monomer symmetry axis
		for ( Size i=1; i<= n_monomers; ++i ) {
			ResidueOP rsd
				( conformation::ResidueFactory::create_residue( monomer_pose.residue_type_set_for_pose()->name_map("VRT")));
			rsd->set_xyz("ORIG", center);
			Vector z( axis ), x( monomer_pose.residue(1).xyz("CA") - center );
			x = ( x - z * z.dot(x) ).normalized();
			Vector y( z.cross( x ) );
			rsd->set_xyz("X", center+x);
			rsd->set_xyz("Y", center+y);
			vrt_rsds.push_back( rsd->clone() );
		}


		for ( Size i=1; i<= n_monomers; ++i ) {
			Vector const & vertex( vertices[i] ), &nbr_vertex( nbr_vertices[i] );
			Vector z( vertex.normalized() ), y( ( nbr_vertex - z * nbr_vertex.dot(z) ).normalized() );
			Vector x( y.cross(z) );
			ResidueOP rsd
				( conformation::ResidueFactory::create_residue( monomer_pose.residue_type_set_for_pose()->name_map("VRT")));
			rsd->set_xyz( "ORIG", Vector(0,0,0) );
			rsd->set_xyz( "X", x );
			rsd->set_xyz( "Y", y );
			vrt_rsds.push_back( rsd );
		}

		for ( Size i=1; i<= n_monomers; ++i ) {
			for ( Size j=1; j<= monomer_pose.total_residue(); ++j ) {
				if ( j==1 && i> 1 ) {
					pose.append_residue_by_jump( monomer_pose.residue(j), 1 );
				} else {
					pose.append_residue_by_bond( monomer_pose.residue(j) );
				}
			}
			if ( i>1 ) pose.conformation().insert_chain_ending( nres_monomer*(i-1) );
		}


		runtime_assert( vrt_rsds.size() == 2 * n_monomers );
		for ( Size i=1; i<= 2*n_monomers; ++i ) pose.append_residue_by_jump( *vrt_rsds[i],1 );

		/// now set the proper foldtree
		Size const nres_protein( n_monomers * nres_monomer );
		pose.conformation().insert_chain_ending( nres_protein );

		FoldTree f( pose.total_residue() );

		/// first n vrt rsds are located on the symmetry axes of corresponding monomers
		for ( Size i=1; i<= n_monomers; ++i ) {
			Size const protpos( (i-1)*nres_monomer + 1 ), vrtpos( nres_protein+i ), cutpoint( nres_monomer*i );
			f.new_jump( protpos, vrtpos, cutpoint );
		}

		/// next n vrt rsds are at the origin, defining the symmetry
		for ( Size i=1; i<=n_monomers; ++i ) {
			Size const vrtpos1( nres_protein+i ), vrtpos2( nres_protein + n_monomers + i );
			f.new_jump( vrtpos1, vrtpos2, vrtpos1 );
		}
		for ( Size i=1; i<n_monomers; ++i ) {
			Size const vrtpos1( nres_protein+n_monomers+i ), vrtpos2( nres_protein + n_monomers + i+1 );
			f.new_jump( vrtpos1, vrtpos2, vrtpos1 );
		}
		Size const root( nres_protein + n_monomers + 1 ); //the first of the origin vrt rsds
		f.reorder( root );
		pose.fold_tree( f );

		/// set an initial jump
		numeric::xyzMatrix_double R( numeric::xyzMatrix_double::identity() );
		Real distance( 100 );
		Jump jump( RT( R, Vector(0,0, distance) ) );

		/// 1 -> n_monomers: jumps from the axial vrtrsds to the protein
		/// n_monomers+1 --> 2*n_monomers: jumps from central vrtrsds
		runtime_assert( pose.fold_tree().num_jump() == 3 * n_monomers-1 );
		for ( Size i=1; i<= n_monomers; ++i ) {
			Size const jumpno( i+n_monomers );
			runtime_assert( pose.fold_tree().upstream_jump_residue(jumpno) == nres_protein + n_monomers + i );
			runtime_assert( pose.fold_tree().downstream_jump_residue(jumpno) == nres_protein + i );
			pose.set_jump( jumpno, jump );
		}


		{ /// try slide into contact
			devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

			ScoreFunction scorefxn;
			scorefxn.set_weight( vdw, 1.0 );

			Real start_score( scorefxn( pose ) );
			while ( true ) {
				distance -= 1.0;


				for ( Size i=1; i<= n_monomers; ++i ) {
					Size const jumpno( i+n_monomers );
					pose.set_jump( jumpno, Jump( RT( numeric::xyzMatrix_double::identity(), Vector(0,0, distance) ) ) );
				}

				Real const score( scorefxn( pose ) );
				TR.Trace << "slide_into_contact: " << F(9,3,distance) << F(9,3,score) << endl;
				if ( score - start_score>5 ) break;
			}
		}


		pose.dump_pdb("cage_"+symm_type+".pdb" );
	}


}


Real
my_calc_total_sasa(
	pose::Pose const & pose,
	Real const probe_radius
)
{
	using namespace id;

	AtomID_Map< Real > atom_sasa;
	Reals rsd_sasa;

	bool const use_big_polar_H( false ), // default
		use_naccess_sasa_radii( false ), // default
		expand_polar_radii( false ), // default
		include_probe_radius_in_atom_radii( true ), // default
		use_lj_radii( true ); // NON-default !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Real const polar_expansion_radius( 1.0 ); // default, unused

	AtomID_Map< bool > atom_subset;
	atom_subset.clear();
	core::pose::initialize_atomid_map( atom_subset, pose, true ); // use all atoms

	return calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius,
		use_big_polar_H, atom_subset,
		use_naccess_sasa_radii, expand_polar_radii, polar_expansion_radius,
		include_probe_radius_in_atom_radii, use_lj_radii );

}

///////////////////////////////////////////////////////////////////////////////
void
count_unsats(
	Pose const & pose_in,
	Size & n_unsat_donors,
	Size & n_unsat_donors_bb,
	Size & n_unsat_acceptors,
	Size & n_unsat_acceptors_bb,
	string const tag = "tmp",
	Real const probe_radius_for_unsat_calcs = 1.0
)
{
	n_unsat_donors = n_unsat_donors_bb = n_unsat_acceptors = n_unsat_acceptors_bb = 0;

	pose::Pose pose( pose_in );

	id::AtomID_Map< Real > atom_sasa, total_hbond_energy;
	utility::vector1< Real > rsd_sasa;

	scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius_for_unsat_calcs,
		true ); // use_big_polar_H

	// HACK!!!!!!!!!!!!!! TO TRIGGER 10A nbr graph update
	scoring::hbonds::HBondSet hbond_set;
	{
		using namespace scoring;
		ScoreFunction sf;
		sf.set_weight( hbond_sc, 1.0 );
		sf(pose);
		pose.update_residue_neighbors();
		scoring::hbonds::fill_hbond_set( pose, false, hbond_set );
	}

	core::pose::initialize_atomid_map( total_hbond_energy, pose, 0.0 );

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		if ( true ) { // hbond_set.allow_hbond(i) ) { ////////////////////// makes more lenient
			scoring::hbonds::HBond const & hb( hbond_set.hbond(i) );
			id::AtomID const hatm( hb.don_hatm(), hb.don_res() );
			id::AtomID const aatm( hb.acc_atm() , hb.acc_res() );
			total_hbond_energy[ hatm ] += hb.energy(); // unweighted
			total_hbond_energy[ aatm ] += hb.energy(); // unweighted
		}
	}

	// now find unsat+buried
	Real const hbond_energy_threshold( -0.01 );
	Real const burial_threshold( 0.01 );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );

		// donors
		for ( chemical::AtomIndices::const_iterator
				hnum  = rsd.Hpos_polar().begin(),
				hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			id::AtomID const hatm( *hnum, i );
			if ( total_hbond_energy[ hatm ] >= hbond_energy_threshold &&
					atom_sasa[ hatm ] <= burial_threshold ) {
				++n_unsat_donors;
				if ( rsd.atom_is_backbone( *hnum ) ) {
					TR.Trace << "unsat_bb_don " << I(4,i) << ' ' << rsd.name1() << A(5,rsd.atom_name(*hnum) ) << ' '<< tag << endl;
					++n_unsat_donors_bb;
				}
			}
		}

		// acceptors
		for ( chemical::AtomIndices::const_iterator
				anum  = rsd.accpt_pos().begin(),
				anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
			id::AtomID const aatm( *anum, i );
			if ( total_hbond_energy[ aatm ] >= hbond_energy_threshold &&
					atom_sasa[ aatm ] <= burial_threshold ) {
				++n_unsat_acceptors;
				if ( rsd.atom_is_backbone( *anum ) ) {
					TR.Trace << "unsat_bb_acc " << I(4,i) << ' ' << rsd.name1() << A(5,rsd.atom_name(*anum) ) << ' '<< tag << endl;
					++n_unsat_acceptors_bb;
				}
			}
		}
	} // i



}


///////////////////////////////////////////////////////////////////////////////

Size non_neg_delta( Size const a, Size const b )
{
	if ( b<a ) return a-b;
	else return 0;
}

///////////////////////////////////////////////////////////////////////////////
void
analyze_interface(
	Pose const & complex_pose,
	Real & bsasa5,
	Real & bsasa14,
	/// these are at the interface (caused by complexation)
	Size & n_unsat_donors,
	Size & n_unsat_donors_bb,
	Size & n_unsat_acceptors,
	Size & n_unsat_acceptors_bb,
	string const tag = "tmp"
	// vector1< id::AtomID > & buried_unsatisfied_donors,
	// vector1< id::AtomID > & buried_unsatisfied_acceptors
)
{
	// let's look at buried sasa

	// also unsatisfied donors and acceptors at interface
	Pose partner1, partner2;
	for ( Size i=1; i<= complex_pose.total_residue(); ++i ) {
		if ( complex_pose.chain(i) == 1 ) partner1.append_residue_by_bond( complex_pose.residue(i) );
		if ( complex_pose.chain(i) == 2 ) partner2.append_residue_by_bond( complex_pose.residue(i) );
	}



	//Real const probe_radius_for_unsat_calcs( 1.0 ); // less than 1.4 --> fewer things "buried" --> more conservative

	//////////////
	// buried sasa
	Real const total_sasa_complex14( my_calc_total_sasa( complex_pose, 1.4 ) );
	Real const total_sasa_complex5 ( my_calc_total_sasa( complex_pose, 0.5 ) );

	Real const total_sasa_partner1_14( my_calc_total_sasa( partner1, 1.4 ) );
	Real const total_sasa_partner1_5 ( my_calc_total_sasa( partner1, 0.5 ) );

	Real const total_sasa_partner2_14( my_calc_total_sasa( partner2, 1.4 ) );
	Real const total_sasa_partner2_5 ( my_calc_total_sasa( partner2, 0.5 ) );

	bsasa14 = ( total_sasa_partner1_14 + total_sasa_partner2_14 - total_sasa_complex14 )/2;
	bsasa5  = ( total_sasa_partner1_5  + total_sasa_partner2_5  - total_sasa_complex5  )/2;

	////////////////////////////////////
	// // now for buried unsatisfied hbonds

	Size n_unsat_donors_partner1, n_unsat_donors_bb_partner1, n_unsat_acceptors_partner1, n_unsat_acceptors_bb_partner1,
		n_unsat_donors_partner2, n_unsat_donors_bb_partner2, n_unsat_acceptors_partner2, n_unsat_acceptors_bb_partner2,
		n_unsat_donors_complex, n_unsat_donors_bb_complex, n_unsat_acceptors_complex, n_unsat_acceptors_bb_complex;

	count_unsats( partner1, n_unsat_donors_partner1, n_unsat_donors_bb_partner1,
		n_unsat_acceptors_partner1, n_unsat_acceptors_bb_partner1, tag+"_partner1" );

	count_unsats( partner2, n_unsat_donors_partner2, n_unsat_donors_bb_partner2,
		n_unsat_acceptors_partner2, n_unsat_acceptors_bb_partner2, tag+"_partner2" );

	count_unsats( complex_pose, n_unsat_donors_complex, n_unsat_donors_bb_complex,
		n_unsat_acceptors_complex, n_unsat_acceptors_bb_complex, tag+"_complex" );

	/// silly:
	n_unsat_donors = non_neg_delta( n_unsat_donors_complex, n_unsat_donors_partner1 + n_unsat_donors_partner2 );
	n_unsat_donors_bb = non_neg_delta( n_unsat_donors_bb_complex, n_unsat_donors_bb_partner1 + n_unsat_donors_bb_partner2);

	n_unsat_acceptors = non_neg_delta( n_unsat_acceptors_complex, n_unsat_acceptors_partner1 + n_unsat_acceptors_partner2);
	n_unsat_acceptors_bb = non_neg_delta( n_unsat_acceptors_bb_complex, n_unsat_acceptors_bb_partner1 +
		n_unsat_acceptors_bb_partner2 );



	// pose::Pose pose;
	// pose = complex_pose; // since we have to perform a non-const operation, updating the nbrs

	// id::AtomID_Map< Real > atom_sasa_complex, total_hbond_energy;
	// utility::vector1< Real > rsd_sasa_complex;
	// scoring::calc_per_atom_sasa( pose, atom_sasa_complex, rsd_sasa_complex, probe_radius_for_unsat_calcs,
	//                true ); // use_big_polar_H

	// // HACK!!!!!!!!!!!!!! TO TRIGGER 10A nbr graph update
	// scoring::hbonds::HBondSet hbond_set;
	// {
	//  using namespace scoring;
	//  ScoreFunction sf;
	//  sf.set_weight( hbond_sc, 1.0 );
	//  sf(pose);
	//  pose.update_residue_neighbors();
	//  scoring::hbonds::fill_hbond_set( pose, false, hbond_set );
	// }

	// core::pose::initialize_atomid_map( total_hbond_energy, pose, 0.0 );

	// for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
	//  if ( true ) { // hbond_set.allow_hbond(i) ) { ////////////////////// makes more lenient
	//   scoring::hbonds::HBond const & hb( hbond_set.hbond(i) );
	//   id::AtomID const hatm( hb.don_hatm(), hb.don_res() );
	//   id::AtomID const aatm( hb.acc_atm() , hb.acc_res() );
	//   total_hbond_energy[ hatm ] += hb.energy(); // unweighted
	//   total_hbond_energy[ aatm ] += hb.energy(); // unweighted
	//  }
	// }

	// // now find unsat+buried
	// buried_unsatisfied_donors.clear();
	// buried_unsatisfied_acceptors.clear();

	// for ( Size i=1; i<= pose.total_residue(); ++i ) {
	//  conformation::Residue const & rsd( pose.residue(i) );

	//  // donors
	//  for ( chemical::AtomIndices::const_iterator
	//      hnum  = rsd.Hpos_polar().begin(),
	//      hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
	//   id::AtomID const hatm( *hnum, i );
	//   if ( total_hbond_energy[ hatm ] >= hbond_energy_threshold &&
	//      atom_sasa_complex[ hatm ] <= burial_threshold ) {
	//    TR << "UNSAT_DONOR " << I(4,i) << ' ' << rsd.name3() << ' ' << rsd.atom_name(*hnum) << std::endl;
	//    buried_unsatisfied_donors.push_back( hatm );
	//    //++buried_unsatisfied_donors;
	//   }
	//  }

	//  // acceptors
	//  for ( chemical::AtomIndices::const_iterator
	//      anum  = rsd.accpt_pos().begin(),
	//      anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
	//   id::AtomID const aatm( *anum, i );
	//   if ( total_hbond_energy[ aatm ] >= hbond_energy_threshold &&
	//      atom_sasa_complex[ aatm ] <= burial_threshold ) {
	//    TR << "UNSAT_ACCEPTOR " << I(4,i) << ' ' << rsd.name3() << ' ' << rsd.atom_name(*anum) << std::endl;
	//    buried_unsatisfied_acceptors.push_back( aatm );
	//    //++buried_unsatisfied_acceptors;
	//   }
	//  }
	// } // i



}


///////////////////////////////////////////////////////////////////////////////
void
analyze_interface_test()
{
	strings const files( start_files() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		if ( num_chains( pose ) > 2 ) {
			pose.conformation().delete_residue_range_slow( chain_begin( 3, pose ), pose.total_residue() );
		}
		runtime_assert( num_chains( pose ) == 2 );

		Real bsasa5, bsasa14;
		Size n_unsat_donors, n_unsat_donors_bb, n_unsat_acceptors, n_unsat_acceptors_bb;
		analyze_interface( pose, bsasa5, bsasa14,
			n_unsat_donors, n_unsat_donors_bb,
			n_unsat_acceptors, n_unsat_acceptors_bb, files[fi] );

		Real bsasa_ratio(0.0);
		if ( bsasa14 >1e-3 ) bsasa_ratio = bsasa5 / bsasa14;

		cout << "analyze_interface: " <<
			" bsasa5: " << F(9,3,bsasa5) <<
			" bsasa14: " << F(9,3,bsasa14) <<
			" bsasa_ratio: " << F(9,3,bsasa_ratio) <<
			" n_unsat_donors: " << I(4,n_unsat_donors) <<
			" n_unsat_donors_bb: " << I(4,n_unsat_donors_bb) <<
			" n_unsat_acceptors: " << I(4,n_unsat_acceptors) <<
			" n_unsat_acceptors_bb: " << I(4,n_unsat_acceptors_bb) <<
			' '<< files[fi] << endl;
	} // files loop

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
cage_design_test()
{
	using namespace kinematics;
	// string const repeatseq("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLG");
	// string const design_me("--------------------++--++-++-++++-"); // no positions in helix1

	Real const cys_dist_ca( 8.2 ); // distance between cysteine c-alphas for gold bridging... approx
	Real const cys_dist_cb( 5.3 ); // distance between cysteine c-alphas for gold bridging... approx

	// redo this for 12x, first the rosetta model
	// lindsey's pdb (5byo)
	string const repeatseq("ISVEELLKLAKAAYYSGTTVEEAYKLALKLG");
	string const design_me("--------------------++--+---+--");
	// rosetta
	// string const repeatseq("VEELLKLAKAAYYSGTTVEEAYKLALKLGIS");
	// string const design_me("------------------++--+---+----");
	Size const nrepeat( 12 );
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	Size const design_cycles( option[ my_options::design_cycles ] );

	bool const nodesign( option[ my_options::nodesign ] );

	// ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	ScoreFunctionOP fa_scorefxn( get_score_function_from_command_line() );
	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );
	adjust_ref_weights_from_command_line( *fa_scorefxn );
	Real const cst_weight( 2.0 );
	fa_scorefxn->set_weight( atom_pair_constraint, cst_weight );



	//string const symm_type( option[ my_options::symm_type ] );
	strings symm_types( make_vector1( string("W0") ) );
	// strings symm_types( make_vector1( string("T0"), string("T1"),
	// 	string("O0"), string("O1"),
	// 	string("I0"), string("I1") ) );

	/// try making some symmetrical structures


	// expand to allow multiple templates
	strings files( start_files() );
	numeric::random::random_permutation( files, numeric::random::rg() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {

		string const monomer_filename( files[fi] );

		///
		Size const repeatlen( repeatseq.size() ), nres_monomer( nrepeat * repeatlen ); // rpos is arbitrary
		Pose monomer_pose;
		pose_from_pdb( monomer_pose, monomer_filename );

		runtime_assert( monomer_pose.sequence().substr(0,repeatlen) == repeatseq );

		// trim off C-terminal stuff
		if ( monomer_pose.total_residue() > nres_monomer ) {
			monomer_pose.conformation().delete_residue_range_slow( nres_monomer+1, monomer_pose.total_residue() );
		}

		// figure out where the symmetry axis is, this assumes perfect symmetry
		Vector center, axis;
		{
			using namespace id;
			Pose fix_pose( monomer_pose ), mov_pose( monomer_pose );

			AtomID_Map< AtomID > atom_map;

			initialize_atomid_map( atom_map, mov_pose, id::GLOBAL_BOGUS_ATOM_ID );

			Size const shift( 3*repeatlen );
			for ( Size i=1; i<= nres_monomer; ++i ) {
				Size const j( (i+shift-1)%nres_monomer+1 );
				atom_map[ AtomID( mov_pose.residue(i).atom_index("CA"), i ) ] =
					AtomID( fix_pose.residue(j).atom_index("CA"),j);
			}
			Real const rmsd( rmsd_by_mapping( mov_pose, fix_pose, atom_map ) );
			TR.Trace << "spin_180_rmsd: " << F(9,3,rmsd) << endl;
			superimpose_pose( mov_pose, fix_pose, atom_map );

			Size const rpos( 5 );
			Stub const stub1( fix_pose.residue( rpos     ).xyz("CA"),
				fix_pose.residue( rpos + 1 ).xyz("CA"),
				fix_pose.residue( rpos + 2 ).xyz("CA") );
			Stub const stub2( mov_pose.residue( rpos     ).xyz("CA"),
				mov_pose.residue( rpos + 1 ).xyz("CA"),
				mov_pose.residue( rpos + 2 ).xyz("CA") );
			Real theta;
			Vector t;
			get_stub_transform_data( stub1, stub2, center, axis, t, theta );
			TR.Trace << "theta: " << numeric::conversions::degrees( theta ) << endl;
		}
		// the old way, perfectly symmetrical template
		// {
		//  Size const rpos(5), repeatlen( 35 ), nrepeat( 6 ); // rpos is arbitrary
		//  runtime_assert( monomer_pose.total_residue() == nrepeat * repeatlen );
		//  Pose const & pose( monomer_pose ); // doh
		//  Stub const stub1( pose.residue( rpos     ).xyz("CA"),
		//           pose.residue( rpos + 1 ).xyz("CA"),
		//           pose.residue( rpos + 2 ).xyz("CA") );
		//  Stub const stub2( pose.residue( rpos + repeatlen    ).xyz("CA"),
		//           pose.residue( rpos + repeatlen + 1 ).xyz("CA"),
		//           pose.residue( rpos + repeatlen + 2 ).xyz("CA") );
		//  Real theta;
		//  Vector t;
		//  Real const dev( get_stub_transform_data( stub1, stub2, center, axis, t, theta ) );
		//  TR.Trace << "theta: " << numeric::conversions::degrees( theta ) << endl;
		// }

		//Size const nres_monomer( monomer_pose.total_residue() );


		numeric::random::random_permutation( symm_types, numeric::random::rg() );

		//Vector const start_axis( axis );
		for ( Size sti=1; sti<= symm_types.size(); ++sti ) {
			string const symm_type( symm_types[ sti ] ), simfile( shared_output_tag()+"_cage_design_"+symm_type+".work" );
			if ( simfile_is_done( simfile ) ) continue;


			//// get the vertices:
			Vectors vertices, nbr_vertices;
			get_symm_type_vertices( symm_type[0], vertices, nbr_vertices );


			/// just a pair of nbring vertices
			vertices = make_vector1( vertices[1], nbr_vertices[1] );
			nbr_vertices = make_vector1( vertices[2], vertices[1] );
			Size const n_monomers( 2 );// = vertices.size();

			ResidueOPs vrt_rsds;

			// create a virtual residue in the middle of the monomer, with z-axis parallel to monomer symmetry axis
			Real const axis_factor( symm_type[1] == '0' ? -1.0 : 1.0 );
			for ( Size i=1; i<= n_monomers; ++i ) {
				ResidueOP rsd
					( conformation::ResidueFactory::create_residue( monomer_pose.residue_type_set_for_pose()->name_map("VRT")));
				rsd->set_xyz("ORIG", center);
				Vector z( axis_factor * axis ), x( monomer_pose.residue(1).xyz("CA") - center );
				x = ( x - z * z.dot(x) ).normalized();
				Vector y( z.cross( x ) );
				rsd->set_xyz("X", center+x);
				rsd->set_xyz("Y", center+y);
				vrt_rsds.push_back( rsd->clone() );
			}


			for ( Size i=1; i<= n_monomers; ++i ) {
				Vector const & vertex( vertices[i] ), &nbr_vertex( nbr_vertices[i] );
				Vector z( vertex.normalized() ), y( ( nbr_vertex - z * nbr_vertex.dot(z) ).normalized() );
				Vector x( y.cross(z) );
				ResidueOP rsd
					( conformation::ResidueFactory::create_residue( monomer_pose.residue_type_set_for_pose()->name_map("VRT")));
				rsd->set_xyz( "ORIG", Vector(0,0,0) );
				rsd->set_xyz( "X", x );
				rsd->set_xyz( "Y", y );
				vrt_rsds.push_back( rsd );
			}

			Pose pose;
			for ( Size i=1; i<= n_monomers; ++i ) {
				for ( Size j=1; j<= monomer_pose.total_residue(); ++j ) {
					if ( j==1 && i> 1 ) {
						pose.append_residue_by_jump( monomer_pose.residue(j), 1 );
					} else {
						pose.append_residue_by_bond( monomer_pose.residue(j) );
					}
				}
				if ( i>1 ) pose.conformation().insert_chain_ending( nres_monomer*(i-1) );
			}

			runtime_assert( vrt_rsds.size() == 2 * n_monomers );
			for ( Size i=1; i<= 2*n_monomers; ++i ) pose.append_residue_by_jump( *vrt_rsds[i],1 );

			/// now set the proper foldtree
			Size const nres_protein( n_monomers * nres_monomer );
			pose.conformation().insert_chain_ending( nres_protein );

			FoldTree f( pose.total_residue() );

			/// first n vrt rsds are located on the symmetry axes of corresponding monomers
			for ( Size i=1; i<= n_monomers; ++i ) {
				Size const protpos( (i-1)*nres_monomer + 1 ), vrtpos( nres_protein+i ), cutpoint( nres_monomer*i );
				f.new_jump( protpos, vrtpos, cutpoint );
			}

			/// next n vrt rsds are at the origin, defining the symmetry
			for ( Size i=1; i<=n_monomers; ++i ) {
				Size const vrtpos1( nres_protein+i ), vrtpos2( nres_protein + n_monomers + i );
				f.new_jump( vrtpos1, vrtpos2, vrtpos1 );
			}
			for ( Size i=1; i<n_monomers; ++i ) {
				Size const vrtpos1( nres_protein+n_monomers+i ), vrtpos2( nres_protein + n_monomers + i+1 );
				f.new_jump( vrtpos1, vrtpos2, vrtpos1 );
			}
			Size const root( nres_protein + n_monomers + 1 ); //the first of the origin vrt rsds
			f.reorder( root );
			pose.fold_tree( f );

			/// set an initial jump
			numeric::xyzMatrix_double R( numeric::xyzMatrix_double::identity() );
			Jump jump( RT( R, Vector(0,0, 30.0) ) );

			/// 1 -> n_monomers: jumps from the axial vrtrsds to the protein
			/// n_monomers+1 --> 2*n_monomers: jumps from central vrtrsds
			runtime_assert( pose.fold_tree().num_jump() == 3 * n_monomers-1 );
			for ( Size i=1; i<= n_monomers; ++i ) {
				Size const jumpno( i+n_monomers );
				runtime_assert( pose.fold_tree().upstream_jump_residue(jumpno) == nres_protein + n_monomers + i );
				runtime_assert( pose.fold_tree().downstream_jump_residue(jumpno) == nres_protein + i );
				pose.set_jump( jumpno, jump );
			}


			conformation::symmetry::SymmetryInfo symminfo;

			{ /// let's try making it a symmetric pose
				runtime_assert( n_monomers == 2 );
				runtime_assert( pose.fold_tree().num_jump() == 2 + 2 + 1 ); //
				symminfo.num_virtuals( 2 * n_monomers );
				symminfo.set_use_symmetry( true );
				for ( Size i=1; i<= nres_monomer; ++i ) {
					symminfo.add_bb_clone( i, nres_monomer+i );
					symminfo.add_chi_clone( i, nres_monomer+i );
				}
				// the jumps from the axial vrt to the monomers (these will be fixed)
				symminfo.add_jump_clone( 1, 2, 0.0 ); // 0.0 is the wt, dont think its used here
				// the jumps from the central vrt rsds to the axial vrt rsds (these will have z rotation and translation flex)
				symminfo.add_jump_clone( 3, 4, 0.0 ); // 0.0 is the wt, dont think its used here

				{ // the SymDof's -- not clear what these do, exactly...
					using core::conformation::symmetry::SymDof;
					map< Size, SymDof > symdofs;
					if ( true ) {
						SymDof symdof;
						symdof.read( "z angle_z" );
						symdofs[ 3 ] = symdof;
					}
					symminfo.set_dofs( symdofs );
				}
				symminfo.set_flat_score_multiply( pose.total_residue(), 1 ); // takes Size not Real !
				for ( Size i=1; i<= nres_monomer; ++i ) symminfo.set_score_multiply( i, 2 );
				symminfo.update_score_multiply_factor();
			}


			Pose const start_pose( pose ); // start_pose is not a symmetric pose, but it is geometrically symmetric


			while ( true ) {
				string const worktag( filebase( monomer_filename ) );
				Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
				if ( n > nstruct() ) break;





				{ /// try slide into contact, centroid simulation
					devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );
					pose::symmetry::make_symmetric_pose( pose, symminfo );

					{ /// randomly rotate
						Size const jumpno( n_monomers+1 );
						Jump jump( pose.jump( jumpno ) );
						/// only 60 degree spin b/c of symmetry...
						jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, numeric::random::uniform()*360.0 ); // full spin for asym mons
						//jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, numeric::random::uniform()*60.0 ); // ie 6
						jump.fold_in_rb_deltas();
						pose.set_jump( jumpno, jump );
					}

					{ /// move far out
						Size const jumpno( n_monomers+1 );
						Jump jump( pose.jump( jumpno ) );
						jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, 100.0 + numeric::random::uniform() ); // since we step integers
						jump.fold_in_rb_deltas();
						pose.set_jump( jumpno, jump );
					}

					scoring::symmetry::SymmetricScoreFunction scorefxn;
					scorefxn.set_weight( vdw, 1.0 );

					Real start_score( scorefxn( pose ) );
					while ( true ) {
						{ /// move back in 1A
							Size const jumpno( n_monomers+1 );
							Jump jump( pose.jump( jumpno ) );
							jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, -1.0 );
							jump.fold_in_rb_deltas();
							pose.set_jump( jumpno, jump );
						}

						Real const score( scorefxn( pose ) );
						TR.Trace << "slide_into_contact: " << F(9,3,score) << endl;
						if ( score - start_score>5 ) {
							{ /// move back out a bit
								Size const jumpno( n_monomers+1 );
								Jump jump( pose.jump( jumpno ) );
								jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, numeric::random::uniform()*1.5 ); // not sure about this
								jump.fold_in_rb_deltas();
								pose.set_jump( jumpno, jump );
							}
							break;
						}
					}
				} // slide into contact (first converts to centroid symmetric)


				// store good poses for each design position
				Sizes design_positions;
				for ( Size i=1; i<= repeatlen; ++i ) {
					if ( design_me[ i-1 ] == '+' ) {
						design_positions.push_back(i);
						TR.Trace << "is_designable: " << i << endl;
					}
				}

				vector1< std::pair< Size, PoseOP > > good_cenposes;
				//map< Size, PoseOP > good_cenposes;


				{ // centroid simulation
					scoring::symmetry::SymmetricScoreFunction scorefxn;
					scorefxn.set_weight( vdw, 1.0 );

					/// search through a grid of slides and spins
					Pose const start_pose( pose );
					Real const start_vdw( scorefxn( pose ) );
					Real const max_vdw_increase( 4.0 );

					Real const z_stepsize( 0.25 ), a_stepsize( 0.5 ); // Angstroms, degrees

					for ( int z_step=-5; z_step<=20; ++z_step ) { // start in, move out
						for ( int a_step=-30; a_step<= 30; ++a_step ) {

							pose = start_pose;


							{ /// rotate by a_step
								Size const jumpno( n_monomers+1 );
								Jump jump( pose.jump( jumpno ) );
								jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, a_step * a_stepsize );
								jump.fold_in_rb_deltas();
								pose.set_jump( jumpno, jump );
							}

							{ /// slide by z_step
								Size const jumpno( n_monomers+1 );
								Jump jump( pose.jump( jumpno ) );
								jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, z_step*z_stepsize );
								jump.fold_in_rb_deltas();
								pose.set_jump( jumpno, jump );
							}

							// don't want clashes...
							if ( scorefxn( pose ) - start_vdw > max_vdw_increase ) continue;

							// look for positions that could be bridged
							for ( Size dpos : design_positions ) {
								// we want a copy of this position with a close approach to the
								Real min_disdev(1e6);
								Size min_rep1(0);
								for ( Size rep1=0; rep1<nrepeat; ++rep1 ) {
									Size const rep2( (rep1+1)%nrepeat );
									Size const pos1( dpos + rep1*repeatlen ), pos2( dpos + rep2*repeatlen + nres_monomer );
									Real const disdev(
										fabs( pose.residue(pos1).xyz("CA").distance( pose.residue(pos2).xyz("CA")) - cys_dist_ca ) +
										fabs( pose.residue(pos1).xyz("CB").distance( pose.residue(pos2).xyz("CB")) - cys_dist_cb ) );
									if ( disdev < min_disdev ) {
										min_disdev = disdev;
										min_rep1 = rep1;
									}
								}
								TR.Trace << "min_disdev: " << min_disdev <<
									" min_rep1: " << min_rep1 <<
									" dpos: " << dpos <<
									" z_step: " << z_step <<
									" a_step: " << a_step << endl;
								if ( min_disdev<0.7 || ( (dpos == 29 || dpos==22 ) && min_disdev<1.4 ) ) {
									// dump pdb
									// string const outfilename( "trap_cage_"+string_of(dpos)+"_"+string_of(z_step)+"_"+string_of(a_step)+
									// 	"_"+string_of(min_disdev)+".pdb" );
									// pose.dump_pdb(outfilename);
									PoseOP pose_ptr( new Pose() );
									*pose_ptr = pose;
									good_cenposes.push_back( make_pair( dpos, pose_ptr ) );
								}
							} // dpos
						} // a_step
					} // z_step

				}



				if ( good_cenposes.empty() ) { // failure!
					cerr << "no good_cenposes" << endl;
					continue;
				}

				// now we choose a random one of the good_cenposes
				std::pair< Size, PoseOP > pp( random_element( good_cenposes ) );
				Size const cys_pos( pp.first ); // between 1...repeatlen
				pose = *(pp.second);

				//// fullatom
				pose::symmetry::make_asymmetric_pose( pose );
				devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );

				for ( Size i=1; i<= nres_protein; ++i ) {
					pose.replace_residue( i, start_pose.residue(i), true ); // orient_backbone
				}

				// mutate the cysteine positions
				for ( Size i=1; i<= nres_protein; ++i ) {
					if ( (i-1)%repeatlen == cys_pos-1 ) { // mutate to cysteine
						make_sequence_change( i, aa_cys, pose );
					}
				}

				pose::symmetry::make_symmetric_pose( pose, symminfo );

				Size num_mutations(0);
				if ( nodesign ) {
					// make random polar aa choices at the non-cys design positions
					// then fastrelax all those positions
					string const polars( "DEHKRSTNQ" ); // too many tyrosines
					//string const polars( "DEHKRSTNQY" );
					MoveMapOP mm( new MoveMap );
					for ( Size dpos : design_positions ) {
						if ( dpos == cys_pos ) continue;
						AA const new_aa( aa_from_oneletter_code( polars[ random_range( 0, polars.size()-1 ) ] ) );
						AA const old_aa( pose.residue(dpos).aa() );
						TR.Trace << "new_aa: " << new_aa << " dpos: " << dpos << " cys_pos: " << cys_pos <<
							" old_aa: " << old_aa << endl;
						num_mutations += ( new_aa != old_aa );
						for ( Size i=0; i< nrepeat; ++i ) {
							make_sequence_change( i*repeatlen+dpos, new_aa, pose );
							mm->set_chi( i*repeatlen+dpos, true );
							mm->set_chi( i*repeatlen+cys_pos, true );
						}
					}
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );
					fastrelax.set_movemap( mm );
					TR.Trace << "relax the mutated positions" << endl;
					fastrelax.apply( pose );

				}


				// setup rotamer linking
				// pack::rotamer_set::RotamerLinksOP rotamer_links( new pack::rotamer_set::RotamerLinks() );
				// rotamer_links->resize( pose.total_residue());
				// for ( Size i=1; i<= pose.total_residue(); ++i ) {
				// 	rotamer_links->set_equiv(i,i);
				// }
				// for ( Size dpos : design_positions ) {
				// 	if ( dpos == cys_pos ) continue;
				// 	vector1<int> posl;
				// 	for ( Size i=0; i<2*nrepeat; ++i ) {
				// 		posl.push_back( i*repeatlen + dpos );
				// 	}
				// 	for ( Size pos: posl ) {
				// 		rotamer_links->set_equiv( pos, posl );
				// 	}
				// }


				// add constraints between the cysteines and their partners
				// find the repeat number
				Size rep1(0), rep2(0); // NOTE these are 0-indexed
				{
					Real min_disdev(1e6);
					Size min_rep1(0);
					for ( Size rep1=0; rep1<nrepeat; ++rep1 ) {
						Size const rep2( (rep1+1)%nrepeat );
						Size const pos1( cys_pos + rep1*repeatlen ), pos2( cys_pos + rep2*repeatlen + nres_monomer );
						Real const disdev(
							fabs( pose.residue(pos1).xyz("CA").distance( pose.residue(pos2).xyz("CA")) - cys_dist_ca ) +
							fabs( pose.residue(pos1).xyz("CB").distance( pose.residue(pos2).xyz("CB")) - cys_dist_cb ) );
						if ( disdev < min_disdev ) {
							min_disdev = disdev;
							min_rep1 = rep1;
						}
					}
					TR.Trace << "min_disdev: " << min_disdev <<
						" min_rep1: " << min_rep1 <<
						" cys_pos: " << cys_pos << endl;

					rep1 = min_rep1;
					rep2 = (rep1+1)%nrepeat;

					// add constraints
					using namespace scoring::constraints;
					using namespace scoring::func;
					using namespace id;

					utility::vector1< ConstraintCOP > new_constraints;

					for ( Size r=1; r<=2; ++r ) {
						// in monomer1
						Size const pos1( r == 1 ? rep1*repeatlen + cys_pos : rep2*repeatlen + cys_pos );
						// in monomer2
						Size const pos2( r == 1 ? rep2*repeatlen + cys_pos + nres_monomer : rep1*repeatlen + cys_pos + nres_monomer);
						AtomID const pos1_sg( pose.residue( pos1 ).atom_index( "SG"), pos1 );
						AtomID const pos1_cb( pose.residue( pos1 ).atom_index( "CB"), pos1 );
						AtomID const pos2_sg( pose.residue( pos2 ).atom_index( "SG"), pos2 );
						AtomID const pos2_cb( pose.residue( pos2 ).atom_index( "CB"), pos2 );
						FuncOP sg_sg_func( new HarmonicFunc( 4.9, 0.75 ) );
						FuncOP sg_cb_func( new MinMaxFunc( 5.3, 1000.0, 2.0, 2.0 ) );
						FuncOP cb_cb_func( new MinMaxFunc( 5.3, 1000.0, 2.0, 2.0 ) );
						new_constraints.push_back( ConstraintCOP( new AtomPairConstraint( pos1_sg, pos2_sg, sg_sg_func->clone() )));
						new_constraints.push_back( ConstraintCOP( new AtomPairConstraint( pos1_sg, pos2_cb, sg_cb_func->clone() )));
						new_constraints.push_back( ConstraintCOP( new AtomPairConstraint( pos1_cb, pos2_sg, sg_cb_func->clone() )));
						new_constraints.push_back( ConstraintCOP( new AtomPairConstraint( pos1_cb, pos2_cb, cb_cb_func->clone() )));
						Real const start_dist( pose.xyz(pos1_sg).distance( pose.xyz(pos2_sg) ) );
						TR.Trace << "start_dist: "<< start_dist << " start_cst_score: " << sg_sg_func->func( start_dist ) << endl;
					}
					pose.add_constraints( new_constraints );
				}






				if ( false ) { // repack all sidechains, not necessary if starting pose is perfect
					TR.Trace << "optimize template sidechains " << monomer_filename << endl;
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );
					MoveMapOP mm( new MoveMap );
					// the old way was to optimize *ALL* the side chains, but that takes too long for the 12x
					//mm->set_chi( true );
					setup_cage_design_movemap( pose, *mm );
					mm->set_jump(false);
					fastrelax.set_movemap( mm );
					fastrelax.apply( pose );
				}


				bools is_designable( pose.total_residue(), false );
				if ( !nodesign ) {
					Size const repeatlen( design_me.size() );
					runtime_assert( nres_monomer % repeatlen == 0 );
					for ( Size i=1; i<= nres_monomer; ++i ) {
						is_designable[i] = ( design_me[ (i-1)%repeatlen ] == '+' && ( (i-1)%repeatlen != cys_pos-1 ) );
						if ( is_designable[i] ) {
							TR.Trace << "is_designable: " << i << ' ' << (i-1)%repeatlen+1 << endl;
						}
					}
				}


				//Real starttime = clock();
				for ( Size m=1; m<= design_cycles+1; ++m ) {


					// { // deriv test
					//  MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					//  setup_cage_design_movemap( pose, *movemap );
					//  //movemap->set_chi( false );

					//  protocols::minimization_packing::symmetry::SymMinMoverOP min_mover
					//   ( new protocols::minimization_packing::symmetry::SymMinMover(movemap, fa_scorefxn, "dfpmin_armijo_nonmonotone", 0.00001,
					//                              true, true, true ) );

					//  min_mover->apply( pose );


					//  {
					//   pose::symmetry::make_asymmetric_pose( pose );
					//   ScoreFunctionOP scorefxn( new ScoreFunction() );
					//   scorefxn->set_weight( fa_atr, 1.0 );
					//   scorefxn->set_weight( fa_rep, 1.0 );

					//   protocols::minimization_packing::MinMoverOP min_mover
					//    ( new protocols::minimization_packing::MinMover(movemap, scorefxn, "dfpmin", 0.00001, true,
					//                        true, true ));
					//   min_mover->apply( pose );
					//  }

					//  exit(0);
					// }



					bool const skip_relax( m == 1 && !option[ my_options::fastrelax_before_design ] );

					if ( !skip_relax ) { // relax
						protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

						MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
						setup_cage_design_movemap( pose, *movemap );
						//movemap->set_chi( false );
						fastrelax.set_movemap( movemap );
						// static Size counter(0);
						// ++counter;
						// pose.dump_pdb("before_fastrelax_"+lead_zero_string_of(counter,4)+".pdb" );
						if ( !dry_run() ) fastrelax.apply( pose );
						// pose.dump_pdb("after_fastrelax_"+lead_zero_string_of(counter,4)+".pdb" );
					}

					if ( m>design_cycles ) break; // relax is cheap, so do one extra one at the end...

					if ( !nodesign ) { // design



						pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
						task->initialize_from_command_line();
						task->or_include_current( true );
						if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

						bools is_flexible( pose.total_residue(), false );
						MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
						setup_cage_design_movemap( pose, *movemap );
						for ( Size i=1; i<= nres_monomer; ++i ) is_flexible[i] = ( movemap->get_chi(i) );
						task->restrict_to_residues( is_flexible );
						Size n_repack(0), n_design(0);
						for ( Size i=1; i<= pose.total_residue(); ++i ) {
							if ( is_flexible[i] ) {
								if ( is_designable[i] ) {
									++n_design;
								} else {
									++n_repack;
									task->nonconst_residue_task(i).restrict_to_repacking();
								}
							}
						}
						//task->rotamer_links( rotamer_links );
						TR.Trace << "n_design: " << n_design << " n_repack: " << n_repack << endl;

						if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
							protocols::task_operations::LimitAromaChi2Operation lp_op;
							lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negatives)
							lp_op.apply( pose, *task );
						}
						Size const nloop( 25 );
						ScoreFunctionOP design_scorefxn(0);

						if ( m == design_cycles && option[ my_options::final_round_design_score_function ].user() ) {
							design_scorefxn = setup_score_function( option[ my_options::final_round_design_score_function ] );
							adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
						} else if ( option[ my_options::design_score_function ].user() ) {
							design_scorefxn = setup_score_function( option[ my_options::design_score_function ] );
							adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
						} else {
							design_scorefxn = fa_scorefxn; // already adjusted refwts
						}
						design_scorefxn->set_weight( atom_pair_constraint, cst_weight );
						protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
						if ( !dry_run() ) packmover.apply( pose );
					}

				} // cycles

				//Real const relax_simtime( ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes
				Real const final_score( (*fa_scorefxn)( pose ) );

				// compute interface energies
				EnergyMap interface_emap;
				{
					EnergyGraph const & energy_graph( pose.energies().energy_graph() );
					for ( Size pos1 = 1; pos1 <= nres_monomer; ++pos1 ) {
						runtime_assert( pose.chain(pos1) == 1 );
						for ( utility::graph::Graph::EdgeListConstIter
								ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
								ire = energy_graph.get_node( pos1 )->const_edge_list_end();
								ir != ire; ++ir ) {
							EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
							Size const pos2( edge->get_other_ind( pos1 ) );
							if ( pose.chain(pos2) == 2 ) {
								edge->add_to_energy_map( interface_emap );
							}
						}
					}
				}
				Real interface_energy( fa_scorefxn->weights().dot( interface_emap ) );
				Real const cst_energy( cst_weight * pose.energies().total_energies()[ atom_pair_constraint ] );
				interface_energy += cst_energy;

				Real bsasa5, bsasa14;
				Size n_unsat_donors, n_unsat_donors_bb, n_unsat_acceptors, n_unsat_acceptors_bb;
				{
					Pose interface_pose;
					for ( Size i=1; i<= nres_monomer; ++i ) interface_pose.append_residue_by_bond( pose.residue( i ) );
					interface_pose.append_residue_by_jump( pose.residue( nres_monomer+1 ), 1 );
					for ( Size i=nres_monomer+2; i<= 2*nres_monomer; ++i ) {
						interface_pose.append_residue_by_bond( pose.residue( i ) );
					}
					interface_pose.conformation().insert_chain_ending( nres_monomer );
					//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;

					analyze_interface( interface_pose, bsasa5, bsasa14,
						n_unsat_donors, n_unsat_donors_bb,
						n_unsat_acceptors, n_unsat_acceptors_bb );
				}


				bool passed_score_filter
					( append_score_to_scorefile_and_filter( worktag+string_of(cys_pos), final_score, score_filter_acceptance_rate,
					score_filter_pass_early, simfile ) );

				// now also filter on interface energy
				passed_score_filter = passed_score_filter ||
					( append_score_to_scorefile_and_filter( worktag+string_of(cys_pos)+"_int", interface_energy,
						score_filter_acceptance_rate, score_filter_pass_early, simfile ) );

				if ( num_mutations==0 ) passed_score_filter = true;

				string const outfilename( output_tag() + "cage_design_"+symm_type+"_cys_"+string_of(cys_pos)+
					"_T"+ filebase( monomer_filename) +
					"_N"+ lead_zero_string_of( n,4)+".pdb" );

				Real bsasa_ratio(0.0);
				if ( bsasa14 >1e-3 ) bsasa_ratio = bsasa5 / bsasa14;

				string design_seqs;
				for ( Size rep : make_vector1(rep1,rep2) ) {
					design_seqs += ".";
					for( Size pos:design_positions ) {
						design_seqs.push_back( pose.sequence()[ rep*repeatlen + pos-1 ] );
					}
				}

				ostringstream out;
				out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
					//" simtime: " << F(9,3,simtime) <<
					" cys_pos: " << cys_pos <<
					" design_seqs: " << design_seqs <<
					" interface_energy: " << F(9,3,interface_energy) <<
					" num_mutations: " << num_mutations <<
					" cst_energy: " << F(9,3,cst_energy) <<
					" monomer_filename: " << monomer_filename <<
					" passed_score_filter: " << passed_score_filter <<
					" bsasa5: " << F(9,3,bsasa5) <<
					" bsasa14: " << F(9,3,bsasa14) <<
					" bsasa_ratio: " << F(9,3,bsasa_ratio) <<
					" n_unsat_donors: " << n_unsat_donors <<
					" n_unsat_donors_bb: " << n_unsat_donors_bb <<
					" n_unsat_acceptors: " << n_unsat_acceptors <<
					" n_unsat_acceptors_bb: " << n_unsat_acceptors_bb <<
					" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
					" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

				string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

				basic::prof_show_oneliner();

				if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
					append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
				}
				check_simtime();
			} // nstruct
			signal_that_simfile_is_done( simfile );
		} // symm_types loop
	} // monomer files
	signal_that_job_is_done();


	/// test symmetry machinery
	// Size const dir(1);
	// for ( Size i=1; i<= 10; ++i ) {
	//  Size const jumpno( n_monomers+1 );
	//  Jump jump( pose.jump( jumpno ) );
	//  jump.set_rb_delta( kinematics::Jump::TRANS_Z, dir, 1.0 ); // ie 3
	//  jump.fold_in_rb_deltas();
	//  pose.set_jump( jumpno, jump );
	//  pose.dump_pdb("test1_"+symm_type+lead_zero_string_of(i,2)+".pdb" );
	// }
	// for ( Size i=1; i<= 10; ++i ) {
	//  Size const jumpno( n_monomers+1 );
	//  Jump jump( pose.jump( jumpno ) );
	//  jump.set_rb_delta( kinematics::Jump::ROT_Z, dir, 10.0 ); // ie 6
	//  jump.fold_in_rb_deltas();
	//  pose.set_jump( jumpno, jump );
	//  pose.dump_pdb("test2_"+symm_type+lead_zero_string_of(i,2)+".pdb" );
	// }


}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
generic_cage_design_test()
{
	using namespace kinematics;
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	Size const design_cycles( option[ my_options::design_cycles ] );

	// ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	ScoreFunctionOP fa_scorefxn( core::scoring::get_score_function() );
	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );
	adjust_ref_weights_from_command_line( *fa_scorefxn );


	Size const nrepeat( option[ my_options::nrepeat ] );

	//string const symm_type( option[ my_options::symm_type ] );
	strings symm_types( make_vector1( string("T0"), string("T1"),
		string("O0"), string("O1"),
		string("I0"), string("I1") ) );

	/// try making some symmetrical structures


	// expand to allow multiple templates
	strings files( start_files() );
	numeric::random::random_permutation( files, numeric::random::rg() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {

		string const monomer_filename( files[fi] );

		///
		Pose monomer_pose;
		pose_from_pdb( monomer_pose, monomer_filename );

		runtime_assert( num_chains( monomer_pose ) == 1 );
		runtime_assert( monomer_pose.total_residue() % nrepeat == 0 );
		Size const repeatlen( monomer_pose.total_residue()/nrepeat ), nres_monomer( nrepeat * repeatlen );

		// figure out where the symmetry axis is, this assumes perfect symmetry
		Vector center, axis;
		{
			using namespace id;
			Pose fix_pose( monomer_pose ), mov_pose( monomer_pose );

			AtomID_Map< AtomID > atom_map;

			initialize_atomid_map( atom_map, mov_pose, id::GLOBAL_BOGUS_ATOM_ID );

			Size const shift( 3*repeatlen );
			for ( Size i=1; i<= nres_monomer; ++i ) {
				Size const j( (i+shift-1)%nres_monomer+1 );
				atom_map[ AtomID( mov_pose.residue(i).atom_index("CA"), i ) ] =
					AtomID( fix_pose.residue(j).atom_index("CA"),j);
			}
			Real const rmsd( rmsd_by_mapping( mov_pose, fix_pose, atom_map ) );
			TR.Trace << "spin_180_rmsd: " << F(9,3,rmsd) << endl;
			superimpose_pose( mov_pose, fix_pose, atom_map );

			Size const rpos( 5 );
			Stub const stub1( fix_pose.residue( rpos     ).xyz("CA"),
				fix_pose.residue( rpos + 1 ).xyz("CA"),
				fix_pose.residue( rpos + 2 ).xyz("CA") );
			Stub const stub2( mov_pose.residue( rpos     ).xyz("CA"),
				mov_pose.residue( rpos + 1 ).xyz("CA"),
				mov_pose.residue( rpos + 2 ).xyz("CA") );
			Real theta;
			Vector t;
			get_stub_transform_data( stub1, stub2, center, axis, t, theta );
			TR.Trace << "theta: " << numeric::conversions::degrees( theta ) << endl;
		}
		// the old way, perfectly symmetrical template
		// {
		//  Size const rpos(5), repeatlen( 35 ), nrepeat( 6 ); // rpos is arbitrary
		//  runtime_assert( monomer_pose.total_residue() == nrepeat * repeatlen );
		//  Pose const & pose( monomer_pose ); // doh
		//  Stub const stub1( pose.residue( rpos     ).xyz("CA"),
		//           pose.residue( rpos + 1 ).xyz("CA"),
		//           pose.residue( rpos + 2 ).xyz("CA") );
		//  Stub const stub2( pose.residue( rpos + repeatlen    ).xyz("CA"),
		//           pose.residue( rpos + repeatlen + 1 ).xyz("CA"),
		//           pose.residue( rpos + repeatlen + 2 ).xyz("CA") );
		//  Real theta;
		//  Vector t;
		//  Real const dev( get_stub_transform_data( stub1, stub2, center, axis, t, theta ) );
		//  TR.Trace << "theta: " << numeric::conversions::degrees( theta ) << endl;
		// }

		//Size const nres_monomer( monomer_pose.total_residue() );


		numeric::random::random_permutation( symm_types, numeric::random::rg() );

		//Vector const start_axis( axis );
		for ( Size sti=1; sti<= symm_types.size(); ++sti ) {
			string const symm_type( symm_types[ sti ] ), simfile( shared_output_tag()+"_cage_design_"+symm_type+".work" );
			if ( simfile_is_done( simfile ) ) continue;


			//// get the vertices:
			Vectors vertices, nbr_vertices;
			get_symm_type_vertices( symm_type[0], vertices, nbr_vertices );


			/// just a pair of nbring vertices
			vertices = make_vector1( vertices[1], nbr_vertices[1] );
			nbr_vertices = make_vector1( vertices[2], vertices[1] );
			Size const n_monomers( 2 );// = vertices.size();

			ResidueOPs vrt_rsds;

			// create a virtual residue in the middle of the monomer, with z-axis parallel to monomer symmetry axis
			Real const axis_factor( symm_type[1] == '0' ? -1.0 : 1.0 );
			for ( Size i=1; i<= n_monomers; ++i ) {
				ResidueOP rsd
					( conformation::ResidueFactory::create_residue( monomer_pose.residue_type_set_for_pose()->name_map("VRT")));
				rsd->set_xyz("ORIG", center);
				Vector z( axis_factor * axis ), x( monomer_pose.residue(1).xyz("CA") - center );
				x = ( x - z * z.dot(x) ).normalized();
				Vector y( z.cross( x ) );
				rsd->set_xyz("X", center+x);
				rsd->set_xyz("Y", center+y);
				vrt_rsds.push_back( rsd->clone() );
			}


			for ( Size i=1; i<= n_monomers; ++i ) {
				Vector const & vertex( vertices[i] ), &nbr_vertex( nbr_vertices[i] );
				Vector z( vertex.normalized() ), y( ( nbr_vertex - z * nbr_vertex.dot(z) ).normalized() );
				Vector x( y.cross(z) );
				ResidueOP rsd
					( conformation::ResidueFactory::create_residue( monomer_pose.residue_type_set_for_pose()->name_map("VRT")));
				rsd->set_xyz( "ORIG", Vector(0,0,0) );
				rsd->set_xyz( "X", x );
				rsd->set_xyz( "Y", y );
				vrt_rsds.push_back( rsd );
			}

			Pose pose;
			for ( Size i=1; i<= n_monomers; ++i ) {
				for ( Size j=1; j<= monomer_pose.total_residue(); ++j ) {
					if ( j==1 && i> 1 ) {
						pose.append_residue_by_jump( monomer_pose.residue(j), 1 );
					} else {
						pose.append_residue_by_bond( monomer_pose.residue(j) );
					}
				}
				if ( i>1 ) pose.conformation().insert_chain_ending( nres_monomer*(i-1) );
			}

			runtime_assert( vrt_rsds.size() == 2 * n_monomers );
			for ( Size i=1; i<= 2*n_monomers; ++i ) pose.append_residue_by_jump( *vrt_rsds[i],1 );

			/// now set the proper foldtree
			Size const nres_protein( n_monomers * nres_monomer );
			pose.conformation().insert_chain_ending( nres_protein );

			FoldTree f( pose.total_residue() );

			/// first n vrt rsds are located on the symmetry axes of corresponding monomers
			for ( Size i=1; i<= n_monomers; ++i ) {
				Size const protpos( (i-1)*nres_monomer + 1 ), vrtpos( nres_protein+i ), cutpoint( nres_monomer*i );
				f.new_jump( protpos, vrtpos, cutpoint );
			}

			/// next n vrt rsds are at the origin, defining the symmetry
			for ( Size i=1; i<=n_monomers; ++i ) {
				Size const vrtpos1( nres_protein+i ), vrtpos2( nres_protein + n_monomers + i );
				f.new_jump( vrtpos1, vrtpos2, vrtpos1 );
			}
			for ( Size i=1; i<n_monomers; ++i ) {
				Size const vrtpos1( nres_protein+n_monomers+i ), vrtpos2( nres_protein + n_monomers + i+1 );
				f.new_jump( vrtpos1, vrtpos2, vrtpos1 );
			}
			Size const root( nres_protein + n_monomers + 1 ); //the first of the origin vrt rsds
			f.reorder( root );
			pose.fold_tree( f );

			/// set an initial jump
			numeric::xyzMatrix_double R( numeric::xyzMatrix_double::identity() );
			Jump jump( RT( R, Vector(0,0, 30.0) ) );

			/// 1 -> n_monomers: jumps from the axial vrtrsds to the protein
			/// n_monomers+1 --> 2*n_monomers: jumps from central vrtrsds
			runtime_assert( pose.fold_tree().num_jump() == 3 * n_monomers-1 );
			for ( Size i=1; i<= n_monomers; ++i ) {
				Size const jumpno( i+n_monomers );
				runtime_assert( pose.fold_tree().upstream_jump_residue(jumpno) == nres_protein + n_monomers + i );
				runtime_assert( pose.fold_tree().downstream_jump_residue(jumpno) == nres_protein + i );
				pose.set_jump( jumpno, jump );
			}


			conformation::symmetry::SymmetryInfo symminfo;

			{ /// let's try making it a symmetric pose
				runtime_assert( n_monomers == 2 );
				runtime_assert( pose.fold_tree().num_jump() == 2 + 2 + 1 ); //
				symminfo.num_virtuals( 2 * n_monomers );
				symminfo.set_use_symmetry( true );
				for ( Size i=1; i<= nres_monomer; ++i ) {
					symminfo.add_bb_clone( i, nres_monomer+i );
					symminfo.add_chi_clone( i, nres_monomer+i );
				}
				// the jumps from the axial vrt to the monomers (these will be fixed)
				symminfo.add_jump_clone( 1, 2, 0.0 ); // 0.0 is the wt, dont think its used here
				// the jumps from the central vrt rsds to the axial vrt rsds (these will have z rotation and translation flex)
				symminfo.add_jump_clone( 3, 4, 0.0 ); // 0.0 is the wt, dont think its used here

				{ // the SymDof's -- not clear what these do, exactly...
					using core::conformation::symmetry::SymDof;
					map< Size, SymDof > symdofs;
					if ( true ) {
						SymDof symdof;
						symdof.read( "z angle_z" );
						symdofs[ 3 ] = symdof;
					}
					symminfo.set_dofs( symdofs );
				}
				symminfo.set_flat_score_multiply( pose.total_residue(), 1 ); // takes Size not Real !
				for ( Size i=1; i<= nres_monomer; ++i ) symminfo.set_score_multiply( i, 2 );
				symminfo.update_score_multiply_factor();
			}


			Pose const start_pose( pose ); // start_pose is not a symmetric pose, but it is geometrically symmetric


			while ( true ) {
				string const worktag( filebase( monomer_filename ) );
				Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
				if ( n > nstruct() ) break;





				{ /// try slide into contact, centroid simulation
					devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );
					pose::symmetry::make_symmetric_pose( pose, symminfo );

					{ /// randomly rotate
						Size const jumpno( n_monomers+1 );
						Jump jump( pose.jump( jumpno ) );
						/// only 60 degree spin b/c of symmetry...
						jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, numeric::random::uniform()*360.0 ); // full spin for asym mons
						//jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, numeric::random::uniform()*60.0 ); // ie 6
						jump.fold_in_rb_deltas();
						pose.set_jump( jumpno, jump );
					}

					{ /// move far out
						Size const jumpno( n_monomers+1 );
						Jump jump( pose.jump( jumpno ) );
						jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, 100.0 + numeric::random::uniform() ); // since we step integers
						jump.fold_in_rb_deltas();
						pose.set_jump( jumpno, jump );
					}

					scoring::symmetry::SymmetricScoreFunction scorefxn;
					scorefxn.set_weight( vdw, 1.0 );

					Real start_score( scorefxn( pose ) );
					while ( true ) {
						{ /// move back in 1A
							Size const jumpno( n_monomers+1 );
							Jump jump( pose.jump( jumpno ) );
							jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, -1.0 );
							jump.fold_in_rb_deltas();
							pose.set_jump( jumpno, jump );
						}

						Real const score( scorefxn( pose ) );
						TR.Trace << "slide_into_contact: " << F(9,3,score) << endl;
						if ( score - start_score>5 ) {
							{ /// move back out a bit
								Size const jumpno( n_monomers+1 );
								Jump jump( pose.jump( jumpno ) );
								jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, numeric::random::uniform()*1.5 ); // not sure about this
								jump.fold_in_rb_deltas();
								pose.set_jump( jumpno, jump );
							}
							break;
						}
					}
				}


				{ // centroid simulation
					// skip, for now


				}

				//// fullatom
				pose::symmetry::make_asymmetric_pose( pose );
				for ( Size i=1; i<= nres_protein; ++i ) {
					pose.replace_residue( i, start_pose.residue(i), true ); // orient_backbone
				}

				pose::symmetry::make_symmetric_pose( pose, symminfo );

				{ // repack all sidechains, not necessary if starting pose is perfect
					TR.Trace << "optimize template sidechains " << monomer_filename << endl;
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );
					MoveMapOP mm( new MoveMap );
					mm->set_chi( true );
					fastrelax.set_movemap( mm );
					fastrelax.apply( pose );
				}


				bools is_designable( pose.total_residue(), true );
				// {
				//  Size const repeatlen( design_me.size() );
				//  runtime_assert( nres_monomer % repeatlen == 0 );
				//  for ( Size i=1; i<= nres_monomer; ++i ) {
				//   is_designable[i] = ( design_me[ (i-1)%repeatlen ] == '+' );
				//   if ( is_designable[i] ) {
				//    TR.Trace << "is_designable: " << i << ' ' << (i-1)%repeatlen+1 << endl;
				//   }
				//  }
				// }


				//Real starttime = clock();
				for ( Size m=1; m<= design_cycles+1; ++m ) {

					bool const skip_relax( m == 1 && !option[ my_options::fastrelax_before_design ] );

					if ( !skip_relax ) { // relax
						protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

						MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
						setup_cage_design_movemap( pose, *movemap );
						//movemap->set_chi( false );
						fastrelax.set_movemap( movemap );
						// static Size counter(0);
						// ++counter;
						// pose.dump_pdb("before_fastrelax_"+lead_zero_string_of(counter,4)+".pdb" );
						if ( !dry_run() ) fastrelax.apply( pose );
						// pose.dump_pdb("after_fastrelax_"+lead_zero_string_of(counter,4)+".pdb" );
					}

					if ( m>design_cycles ) break; // relax is cheap, so do one extra one at the end...

					{ // design

						pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
						task->initialize_from_command_line();
						task->or_include_current( true );
						if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

						bools is_flexible( pose.total_residue(), false );
						MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
						setup_cage_design_movemap( pose, *movemap );
						for ( Size i=1; i<= nres_monomer; ++i ) is_flexible[i] = ( movemap->get_chi(i) );
						task->restrict_to_residues( is_flexible );
						Size n_repack(0), n_design(0);
						for ( Size i=1; i<= pose.total_residue(); ++i ) {
							if ( is_flexible[i] ) {
								if ( is_designable[i] ) {
									++n_design;
								} else {
									++n_repack;
									task->nonconst_residue_task(i).restrict_to_repacking();
								}
							}
						}
						TR.Trace << "n_design: " << n_design << " n_repack: " << n_repack << endl;

						if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
							protocols::task_operations::LimitAromaChi2Operation lp_op;
							lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negatives)
							lp_op.apply( pose, *task );
						}
						Size const nloop( 25 );
						ScoreFunctionOP design_scorefxn(0);
						if ( ( option[ my_options::use_softrep_for_early_design ] && m<= design_cycles/2 ) ||
								option[ my_options::use_softrep_for_design ] ) {
							design_scorefxn = ScoreFunctionFactory::create_score_function( "soft_rep_design.wts" );
							adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
						} else {
							design_scorefxn = fa_scorefxn; // already adjusted refwts
						}
						protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
						if ( !dry_run() ) packmover.apply( pose );
					}

				} // cycles

				//Real const relax_simtime( ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes
				Real const final_score( (*fa_scorefxn)( pose ) );

				// compute interface energies
				EnergyMap interface_emap;
				{
					EnergyGraph const & energy_graph( pose.energies().energy_graph() );
					for ( Size pos1 = 1; pos1 <= nres_monomer; ++pos1 ) {
						runtime_assert( pose.chain(pos1) == 1 );
						for ( utility::graph::Graph::EdgeListConstIter
								ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
								ire = energy_graph.get_node( pos1 )->const_edge_list_end();
								ir != ire; ++ir ) {
							EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
							Size const pos2( edge->get_other_ind( pos1 ) );
							if ( pose.chain(pos2) == 2 ) {
								edge->add_to_energy_map( interface_emap );
							}
						}
					}
				}
				Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );

				Real bsasa5, bsasa14;
				Size n_unsat_donors, n_unsat_donors_bb, n_unsat_acceptors, n_unsat_acceptors_bb;
				{
					Pose interface_pose;
					for ( Size i=1; i<= nres_monomer; ++i ) interface_pose.append_residue_by_bond( pose.residue( i ) );
					interface_pose.append_residue_by_jump( pose.residue( nres_monomer+1 ), 1 );
					for ( Size i=nres_monomer+2; i<= 2*nres_monomer; ++i ) {
						interface_pose.append_residue_by_bond( pose.residue( i ) );
					}
					interface_pose.conformation().insert_chain_ending( nres_monomer );
					//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;

					analyze_interface( interface_pose, bsasa5, bsasa14,
						n_unsat_donors, n_unsat_donors_bb,
						n_unsat_acceptors, n_unsat_acceptors_bb );
				}


				bool passed_score_filter
					( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
					score_filter_pass_early, simfile ) );

				// now also filter on interface energy
				passed_score_filter = passed_score_filter ||
					( append_score_to_scorefile_and_filter( worktag+"_int", interface_energy, score_filter_acceptance_rate,
					score_filter_pass_early, simfile ) );
				//string const outfilename( output_tag() + "_dock_N"+ lead_zero_string_of( n,4 )+".pdb");
				string const outfilename( output_tag() + "cage_design_"+symm_type+
					"_T"+ filebase( monomer_filename) +
					"_N"+ lead_zero_string_of( n,4)+".pdb" );

				Real bsasa_ratio(0.0);
				if ( bsasa14 >1e-3 ) bsasa_ratio = bsasa5 / bsasa14;
				Real const interface_quality( 100.0 * interface_energy / bsasa14 );

				ostringstream out;
				out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
					//" simtime: " << F(9,3,simtime) <<
					" interface_energy: " << F(9,3,interface_energy) <<
					" monomer_filename: " << monomer_filename <<
					" passed_score_filter: " << passed_score_filter <<
					" bsasa5: " << F(9,3,bsasa5) <<
					" bsasa14: " << F(9,3,bsasa14) <<
					" bsasa_ratio: " << F(9,3,bsasa_ratio) <<
					" interface_quality: " << F(9,3,interface_quality) <<
					" n_unsat_donors: " << n_unsat_donors <<
					" n_unsat_donors_bb: " << n_unsat_donors_bb <<
					" n_unsat_acceptors: " << n_unsat_acceptors <<
					" n_unsat_acceptors_bb: " << n_unsat_acceptors_bb <<
					" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
					" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

				string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

				basic::prof_show_oneliner();

				if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
					append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
				}
				check_simtime();
			} // nstruct
			signal_that_simfile_is_done( simfile );
		} // symm_types loop
	} // monomer files
	signal_that_job_is_done();


	/// test symmetry machinery
	// Size const dir(1);
	// for ( Size i=1; i<= 10; ++i ) {
	//  Size const jumpno( n_monomers+1 );
	//  Jump jump( pose.jump( jumpno ) );
	//  jump.set_rb_delta( kinematics::Jump::TRANS_Z, dir, 1.0 ); // ie 3
	//  jump.fold_in_rb_deltas();
	//  pose.set_jump( jumpno, jump );
	//  pose.dump_pdb("test1_"+symm_type+lead_zero_string_of(i,2)+".pdb" );
	// }
	// for ( Size i=1; i<= 10; ++i ) {
	//  Size const jumpno( n_monomers+1 );
	//  Jump jump( pose.jump( jumpno ) );
	//  jump.set_rb_delta( kinematics::Jump::ROT_Z, dir, 10.0 ); // ie 6
	//  jump.fold_in_rb_deltas();
	//  pose.set_jump( jumpno, jump );
	//  pose.dump_pdb("test2_"+symm_type+lead_zero_string_of(i,2)+".pdb" );
	// }


}






///////////////////////////////////////////////////////////////////////////////
//
// read in two lists of repeat-protein components
//
// read in set of potential interface angles (assume we don't have the parallel case right now
//
// given pair of components and target axis-intersection angle:
// -- add virtual residues for each monomer along the rotation axis
// -- add central virtual residues at origin, pointing along axes
// -- setup fold tree
//
// -- random setup:
// --- rotate each component randomly about axis
// --- slide into contact (a little trickier without symmetry)
// ----- choose distance ratio? slide with those ratios...
//

void
two_component_test()
{
	using namespace kinematics;
	// string const repeatseq("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLG");
	// string const design_me("--------------------++--++-++-++++-"); // no positions in helix1
	// string const design_me("---++--+------+-----++--++-++-++++-");

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );
	Real const centroid_bsasa14_threshold( option[ my_options::centroid_bsasa14_threshold ] );
	//Size const design_cycles( option[ my_options::design_cycles ] );

	// ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	ScoreFunctionOP fa_scorefxn( core::scoring::get_score_function() ),
		cen_scorefxn( new ScoreFunction() );

	// # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
	// env     1.0
	// pair    1.0
	// cbeta   1.0
	// vdw     1.0
	// rg      3.0
	// cenpack 1.0
	// hs_pair 1.0
	// ss_pair 1.0
	// rsigma  1.0
	// sheet   1.0

	cen_scorefxn->set_weight( env, 1.0 );
	cen_scorefxn->set_weight( scoring::pair, 1.0 );
	cen_scorefxn->set_weight( cbeta, 1.0 );
	cen_scorefxn->set_weight( vdw, 1.0 );
	cen_scorefxn->set_weight( cenpack, 1.0 ); // not sure about this guy...

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	Size const nrepeat1( 6 ), nrepeat2( 6 );



	strings symm_types( option[ my_options::symm_types ]() );

	numeric::random::random_permutation( symm_types, numeric::random::rg() );


	for ( Size sti=1; sti<= symm_types.size(); ++sti ) { // loop over symm-types
		string const symm_type( symm_types[sti] );

		// simfile
		string const simfile( shared_output_tag()+"_S"+symm_type+".work" );
		if ( simfile_is_done( simfile ) ) continue;

		Real axis_angle(0);
		{
			Vectors vertices, nbr_vertices;
			get_symm_type_vertices( symm_type[0], vertices, nbr_vertices );
			axis_angle = 0.5 * std::acos( vertices[1].normalized().dot( nbr_vertices[1].normalized() ) );
			if ( symm_type == "O" ) axis_angle *= 2; // all components at vertices
			TR.Trace << "axis_angle: " << F(9,3,numeric::conversions::degrees( axis_angle )) << ' ' << symm_type << endl;
		}

		while ( true ) {
			clock_t starttime( clock() );

			string const worktag( "tmp" );
			Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( n > nstruct() ) break;

			// choose comp1, comp2 randomly
			string const pose1_filename( random_element( option[ my_options::pose1_pdbs ]() ) );
			string pose2_filename( random_element( option[ my_options::pose2_pdbs ]() ) );
			if ( option[ my_options::self_dock_2c ] ) pose2_filename = pose1_filename;

			// choose their flips randomly too (inside/outside)

			// setup pose

			Pose pose1, pose2;
			pose_from_pdb( pose1, pose1_filename );
			pose_from_pdb( pose2, pose2_filename );

			Size const nres_monomer1( pose1.total_residue() ), nres_monomer2( pose2.total_residue() );
			runtime_assert( nres_monomer1%nrepeat1 == 0 );
			Size const repeatlen1( nres_monomer1/nrepeat1 );

			runtime_assert( nres_monomer2%nrepeat2 == 0 );
			Size const repeatlen2( nres_monomer2/nrepeat2 );

			Vector axis1, axis2, center1, center2;
			get_toroid_axis( pose1, nrepeat1, repeatlen1, axis1, center1 );

			get_toroid_axis( pose2, nrepeat2, repeatlen2, axis2, center2 );

			Real const axis1_factor( ( numeric::random::uniform()<0.5 ) ? -1 : 1 );
			Real const axis2_factor( ( numeric::random::uniform()<0.5 ) ? -1 : 1 );


			Pose pose( pose1 );
			pose.append_residue_by_jump( pose2.residue(1), 1 );
			for ( Size i=2; i<= nres_monomer2; ++i ) pose.append_residue_by_bond( pose2.residue(i) );

			// now virtual residues
			for ( Size r=1; r<= 2; ++r ) { //
				Vector const axis( r == 1 ? axis1 : axis2 ), center( r == 1 ? center1 : center2 ),
					orient( r == 1 ? pose1.residue(1).xyz("CA") : pose2.residue(1).xyz("CA") );
				Real const axis_factor( r == 1 ? axis1_factor : axis2_factor );
				ResidueOP rsd
					( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map("VRT")));
				rsd->set_xyz("ORIG", center);
				Vector z( axis_factor * axis ), x( orient - center );
				x = ( x - z * z.dot(x) ).normalized();
				Vector y( z.cross( x ) );
				rsd->set_xyz("X", center+x);
				rsd->set_xyz("Y", center+y);
				pose.append_residue_by_jump( *rsd, 1 );
			}

			// central virtual residues: rotate about the x-axis, so that the z-axes make an angle of axis_angle
			for ( Size r=1; r<= 2; ++r ) {
				Real const angle( r == 1 ? 0.0 : axis_angle );
				ResidueOP rsd
					( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map("VRT")));
				rsd->set_xyz( "ORIG", Vector(0,0,0) );
				rsd->set_xyz( "X", Vector( 1,0,0 ) );
				rsd->set_xyz( "Y", Vector( 0.0, cos( angle ), sin( angle ) ) ); // Y rotates in y-z plane
				pose.append_residue_by_jump( *rsd, 1 );
			}

			// now the fold_tree:
			/// now set the proper foldtree
			Size const nres_protein( nres_monomer1 + nres_monomer2 );
			pose.conformation().insert_chain_ending( nres_monomer1 );
			pose.conformation().insert_chain_ending( nres_protein );

			//pose.dump_pdb("built.pdb");

			FoldTree f( pose.total_residue() );

			/// first 2 (fixed) jumps are between monomers and their axial vrt rsds
			f.new_jump( 1, nres_protein+1, nres_monomer1 );
			f.new_jump( nres_monomer1+1, nres_protein+2, nres_protein );

			/// next 2 (Z-flexible) jumps are from the central to the axial vrt rsds
			f.new_jump( nres_protein+1, nres_protein+3, nres_protein+1 );
			f.new_jump( nres_protein+2, nres_protein+4, nres_protein+2 );

			/// then there's a (fixed) jump between the central vrt rsds
			f.new_jump( nres_protein+3, nres_protein+4, nres_protein+3 );
			Size const root( nres_protein + 3 ); //the first of the origin vrt rsds
			f.reorder( root );
			pose.fold_tree( f );

			/// set an initial jump
			Jump jump( RT( numeric::xyzMatrix_double::identity(), Vector(0,0,0) ) ); // axial vrt rsd also at origin

			pose.set_jump( 3, jump );
			pose.set_jump( 4, jump );

			//pose.dump_pdb("start.pdb");

			/// between 3/4 and 4/3, makes sense?
			Real const slide_scale( 0.75 + ( 1.333 - 0.75 ) * numeric::random::uniform() );

			/// randomly spin each monomer about its axis
			for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
				Jump jump( pose.jump( jumpno ) );
				jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, numeric::random::uniform()*360.0 ); // full spin for asym mons
				jump.fold_in_rb_deltas();
				pose.set_jump( jumpno, jump );
			}

			/// now move far out
			for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
				Jump jump( pose.jump( jumpno ) );
				Real const distance( jumpno==3 ? 100.0 : slide_scale * 100.0 );
				jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, distance );
				jump.fold_in_rb_deltas();
				pose.set_jump( jumpno, jump );
			}

			// convert to centroid for slide into contact and centroid simulation
			//
			devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );
			{ /// now we need to slide into contact

				ScoreFunction scorefxn;
				scorefxn.set_weight( vdw, 1.0 );

				Real start_score( scorefxn( pose ) );
				while ( true ) {
					for ( Size jumpno=3; jumpno<= 4; ++jumpno ) { /// move back in 1A
						Jump jump( pose.jump( jumpno ) );
						Real const step( jumpno==3 ? -1.0 : -1.0 * slide_scale );
						jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, step );
						jump.fold_in_rb_deltas();
						pose.set_jump( jumpno, jump );
					}

					Real const score( scorefxn( pose ) );
					TR.Trace << "slide_into_contact: " << F(9,3,score) << endl;
					if ( score - start_score>5 ) {
						// { /// move back out a bit
						//  Size const jumpno( n_monomers+1 );
						//  Jump jump( pose.jump( jumpno ) );
						//  jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, numeric::random::uniform()*1.5 ); // not sure about this
						//  jump.fold_in_rb_deltas();
						//  pose.set_jump( jumpno, jump );
						// }
						break;
					}
				}
			} // slide into contact

			//pose.dump_pdb("slid.pdb");

			for ( Size i=1; i<= nres_protein; ++i ) {
				// make sequence mutations everywhere
				AA const new_aa( aa_from_oneletter_code( random_element( make_vector1( 'L', 'L', 'L',
					'A', 'A', 'A',
					'I', 'I',
					'G', 'G',
					'V', 'V', 'V', 'V' ) ) ) );
				make_sequence_change( i, new_aa, pose );

			}

			//pose.dump_pdb("randomseq.pdb");


			Real centroid_bsasa14( 0.0 );
			{ // centroid simulation

				(*cen_scorefxn)( pose );

				Size const n_outer( 5 ), n_inner( 500 );

				Real const mc_hitemp( 4 ), mc_lotemp( 0.5 );
				Real const gamma = std::pow( (mc_lotemp/mc_hitemp), 1.0/Real(n_inner*n_outer));
				Real mc_temp( mc_hitemp / gamma );

				protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *cen_scorefxn, mc_temp ) );

				for ( Size n=1; n<= n_outer; ++n ) {
					for ( Size m=1; m<= n_inner; ++m ) {
						mc_temp *= gamma;
						mc->set_temperature( mc_temp );

						/// randomly perturb one of the partners
						Size const jumpno( numeric::random::uniform() < 0.5 ? 3 : 4 );

						Real const slide( numeric::random::gaussian() * 1.5 ),
							spin( numeric::random::gaussian() * 10 );

						Jump jump( pose.jump( jumpno ) );
						jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, slide );
						jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, spin );
						jump.fold_in_rb_deltas();
						pose.set_jump( jumpno, jump );

						//(*cen_scorefxn)( pose );
						bool const mc_accept( mc->boltzmann( pose ) );
						TR.Trace << "mc_accept " << I(4,n) << I(4,m) << I(3,mc_accept) << F(9,3,mc->last_accepted_score()) <<
							F(9,3,mc->lowest_score() ) << F(9,3,mc->temperature() ) << endl;
						//set_jump_rb_centers( pose );
					}
					mc->recover_low( pose );
				}

				// Real const final_score( (*cen_scorefxn)( pose ) );
				// static Size counter(0);
				// ++counter;
				// string const outfilename( "aftercen_"+lead_zero_string_of( counter,4)+".pdb");
				// cout << "final_score: " << F(9,3,final_score) << ' ' << outfilename << endl;
				//pose.dump_pdb( outfilename );


				devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );


				{ // for screening the models, lets check the sasa after very quickly optimizing the interface sidechains
					pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
					task->set_bump_check( false );
					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					setup_interface_movemap( pose, *movemap );
					bools is_flexible( pose.total_residue(), false ), is_designable( pose.total_residue(), true ); // at the moment
					for ( Size i=1; i<= nres_protein; ++i ) is_flexible[i] = ( movemap->get_chi(i) );
					task->restrict_to_repacking();
					task->restrict_to_residues( is_flexible );
					pack::pack_rotamers( pose, *fa_scorefxn, task );
					/// now compute buried sasa
					Real bsasa5(0), bsasa14(0);
					//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;
					Size n_unsat_donors(0), n_unsat_donors_bb(0), n_unsat_acceptors(0), n_unsat_acceptors_bb(0);
					analyze_interface( pose, bsasa5, bsasa14,
						n_unsat_donors, n_unsat_donors_bb,
						n_unsat_acceptors, n_unsat_acceptors_bb );
					centroid_bsasa14 = bsasa14;
					TR.Trace << "centroid_bsasa14: " << F(9,3,centroid_bsasa14) << endl;
				}


				// recover the monomer sidechains
				for ( Size i=1; i<= nres_monomer1; ++i ) {
					pose.replace_residue( i, pose1.residue(i), true ); // orient_backbone
				}
				for ( Size i=1; i<= nres_monomer2; ++i ) {
					pose.replace_residue( nres_monomer1+i, pose2.residue(i), true ); // orient_backbone
				}

			} // centroid simulation

			bool const do_design( centroid_bsasa14 >= centroid_bsasa14_threshold );

			if ( do_design ) { // design/relax
				Size const design_cycles( option[ my_options::design_cycles ] );

				for ( Size m=1; m<= design_cycles+1; ++m ) {
					bool const skip_relax( m == 1 );

					if ( !skip_relax ) { // relax
						protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

						MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
						setup_interface_movemap( pose, *movemap );
						movemap->set_jump( 1, false ); // was set by setup_interface_movemap
						for ( Size jumpno=3; jumpno<=4; ++jumpno ) {
							id::DOF_ID const & slide_z_id
								( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jumpno,id::JUMP,3)));
							id::DOF_ID const & spin_z_id
								( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jumpno,id::JUMP,6)));
							TR.Trace << "slide_z_id: " << slide_z_id << " spin_z_id: " << spin_z_id << endl;
							movemap->set( slide_z_id, true );
							movemap->set( spin_z_id, true );
						}
						fastrelax.set_movemap( movemap );
						if ( !dry_run() ) fastrelax.apply( pose );
					}

					if ( m>design_cycles ) break; // relax is cheap, so do one extra one at the end...

					{ // design
						pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
						task->initialize_from_command_line();
						task->or_include_current( true );
						// if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

						MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
						setup_interface_movemap( pose, *movemap );

						bools is_flexible( pose.total_residue(), false ), is_designable( pose.total_residue(), true ); // at the moment
						for ( Size i=1; i<= nres_protein; ++i ) is_flexible[i] = ( movemap->get_chi(i) );
						task->restrict_to_residues( is_flexible );
						Size n_repack(0), n_design(0);
						for ( Size i=1; i<= pose.total_residue(); ++i ) {
							if ( is_flexible[i] ) {
								if ( is_designable[i] ) {
									++n_design;
								} else {
									++n_repack;
									task->nonconst_residue_task(i).restrict_to_repacking();
								}
							}
						}
						TR.Trace << "n_design: " << n_design << " n_repack: " << n_repack << endl;
						//if ( nodesign ) task->restrict_to_repacking();

						if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
							protocols::task_operations::LimitAromaChi2Operation lp_op;
							lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negatives)
							lp_op.apply( pose, *task );
						}
						Size const nloop( 1 ); // in case linmem ig
						ScoreFunctionOP design_scorefxn(0);
						// if ( ( option[ my_options::use_softrep_for_early_design ] && m<= design_cycles/2 ) ||
						//    option[ my_options::use_softrep_for_design ] ) {
						design_scorefxn = ScoreFunctionFactory::create_score_function( "soft_rep_design.wts" );
						adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
						// } else {
						//  design_scorefxn = fa_scorefxn; // already adjusted refwts
						// }
						protocols::minimization_packing::PackRotamersMover packmover( design_scorefxn, task, nloop );
						// protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
						if ( !dry_run() ) packmover.apply( pose );
					}
				} // cycles
			} // fullatom simulation


			//Real const relax_simtime( ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes
			Real const final_score( (*fa_scorefxn)( pose ) );

			// compute interface energies
			EnergyMap interface_emap;
			{
				EnergyGraph const & energy_graph( pose.energies().energy_graph() );
				for ( Size pos1 = 1; pos1 <= nres_monomer1; ++pos1 ) {
					runtime_assert( pose.chain(pos1) == 1 );
					for ( utility::graph::Graph::EdgeListConstIter
							ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
							ire = energy_graph.get_node( pos1 )->const_edge_list_end();
							ir != ire; ++ir ) {
						EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
						Size const pos2( edge->get_other_ind( pos1 ) );
						if ( pose.chain(pos2) == 2 ) {
							edge->add_to_energy_map( interface_emap );
						}
					}
				}
			}
			Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );

			Real bsasa5, bsasa14;
			Size n_unsat_donors, n_unsat_donors_bb, n_unsat_acceptors, n_unsat_acceptors_bb;
			{
				Pose interface_pose;
				for ( Size i=1; i<= nres_monomer1; ++i ) interface_pose.append_residue_by_bond( pose.residue( i ) );
				interface_pose.append_residue_by_jump( pose.residue( nres_monomer1+1 ), 1 );
				for ( Size i=nres_monomer1+2; i<= nres_protein; ++i ) {
					interface_pose.append_residue_by_bond( pose.residue( i ) );
				}
				interface_pose.conformation().insert_chain_ending( nres_monomer1 );
				//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;

				analyze_interface( interface_pose, bsasa5, bsasa14,
					n_unsat_donors, n_unsat_donors_bb,
					n_unsat_acceptors, n_unsat_acceptors_bb );
			}

			if ( bsasa14<1 ) bsasa14 = 1; // prevent division by 0

			/// what do we filter on?
			///
			/// -- interface energy
			/// -- total energy per residue
			///
			Real const residue_energy( final_score / nres_protein );
			Real const interface_quality( 100.0 * interface_energy / bsasa14 );
			Real const bsasa_ratio( bsasa5 / bsasa14 );
			bool const small_interface( bsasa14 < 500 );


			bool passed_score_filter
				( do_design &&
				append_score_to_scorefile_and_filter( worktag, residue_energy, score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );

			// now also filter on interface energy
			passed_score_filter = passed_score_filter ||
				( do_design &&
				append_score_to_scorefile_and_filter( worktag+"_int", interface_energy, score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );

			passed_score_filter = passed_score_filter && !small_interface; // no small interfaces

			//string const outfilename( output_tag() + "_dock_N"+ lead_zero_string_of( n,4 )+".pdb");
			string const outfilename( output_tag() + "twocomp_S"+symm_type+
				"_N"+ lead_zero_string_of( n,4)+".pdb" );

			Real const simtime = ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

			ostringstream out;
			out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
				//" simtime: " << F(9,3,simtime) <<
				" residue_energy: " << F(9,3,residue_energy) <<
				" interface_energy: " << F(9,3,interface_energy) <<
				" interface_quality: " << F(9,3,interface_quality) <<
				" pose1_filename: " << pose1_filename <<
				" pose2_filename: " << pose2_filename <<
				" passed_score_filter: " << passed_score_filter <<
				" bsasa5: " << F(9,3,bsasa5) <<
				" bsasa14: " << F(9,3,bsasa14) <<
				" bsasa_ratio: " << F(9,3,bsasa_ratio) <<
				" do_design: " << do_design <<
				" centroid_bsasa14: " << F(9,3,centroid_bsasa14) <<
				" simtime: " << F(9,3,simtime) <<
				" n_unsat_donors: " << n_unsat_donors <<
				" n_unsat_donors_bb: " << n_unsat_donors_bb <<
				" n_unsat_acceptors: " << n_unsat_acceptors <<
				" n_unsat_acceptors_bb: " << n_unsat_acceptors_bb <<
				" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
				" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

			string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

			TR.Trace << out.str() << endl;

			basic::prof_show_oneliner();

			if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
				append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
			}
			check_simtime();
		} // nstruct
		signal_that_simfile_is_done( simfile );
	} // symm_types loop
	signal_that_job_is_done();


}



///////////////////////////////////////////////////////////////////////////////
void
planar_two_component_test()
{
	protocols::moves::PyMOLMoverOP pymol(0);
	if ( option[ my_options::show_pymol ] ) {
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol->keep_history(true);
	}
	TR.Trace << "show_pymol: " << option[ my_options::show_pymol ] << ' ' << pymol << endl;

	using namespace kinematics;
	// string const repeatseq("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLG");
	// string const design_me("--------------------++--++-++-++++-"); // no positions in helix1
	// string const design_me("---++--+------+-----++--++-++-++++-");

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );
	Real const centroid_bsasa14_threshold( option[ my_options::centroid_bsasa14_threshold ] );
	//Size const design_cycles( option[ my_options::design_cycles ] );

	// ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	ScoreFunctionOP fa_scorefxn( get_score_function_from_command_line() ), // was: core::scoring::get_score_function() ),
		cen_scorefxn( new ScoreFunction() );

	// # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
	// env     1.0
	// pair    1.0
	// cbeta   1.0
	// vdw     1.0
	// rg      3.0
	// cenpack 1.0
	// hs_pair 1.0
	// ss_pair 1.0
	// rsigma  1.0
	// sheet   1.0

	cen_scorefxn->set_weight( env, 1.0 );
	cen_scorefxn->set_weight( scoring::pair, 1.0 );
	cen_scorefxn->set_weight( cbeta, 1.0 );
	cen_scorefxn->set_weight( vdw, 1.0 );
	cen_scorefxn->set_weight( cenpack, 1.0 ); // not sure about this guy...

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	Size nrepeat1( 6 ), nrepeat2( 6 );

	if ( option[ my_options::nrepeats ].user() ) {
		runtime_assert( option[ my_options::nrepeats ].size() == 2 );
		nrepeat1 = option[ my_options::nrepeats ]()[1];
		nrepeat2 = option[ my_options::nrepeats ]()[2];
	}

	Sizes design_positions;
	if ( option[ my_options::design_positions ].user() ) {
		design_positions = option[ my_options::design_positions ]();
	}

	// strings symm_types( option[ my_options::symm_types ]() );
	// numeric::random::random_permutation( symm_types, numeric::random::rg() );
	string const symm_type( option[ my_options::symm_type ] );//"SP2" );
	runtime_assert( symm_type == "SP0" || symm_type == "SP90" ); // 0-degree or 90-degree angle between axes



	// simfile
	string const simfile( shared_output_tag()+"_S"+symm_type+".work" );
	if ( simfile_is_done( simfile ) ) return;

	// Real axis_angle(0);
	// {
	//  Vectors vertices, nbr_vertices;
	//  get_symm_type_vertices( symm_type[0], vertices, nbr_vertices );
	//  axis_angle = 0.5 * std::acos( vertices[1].normalized().dot( nbr_vertices[1].normalized() ) );
	//  if ( symm_type == "O" ) axis_angle *= 2; // all components at vertices
	//  TR.Trace << "axis_angle: " << F(9,3,numeric::conversions::degrees( axis_angle )) << ' ' << symm_type << endl;
	// }

	while ( true ) {
		clock_t starttime( clock() );

		string const worktag( "tmp" );
		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		// choose comp1, comp2 randomly
		string const pose1_filename( random_element( option[ my_options::pose1_pdbs ]() ) );
		string pose2_filename( random_element( option[ my_options::pose2_pdbs ]() ) );
		if ( option[ my_options::self_dock_2c ] ) pose2_filename = pose1_filename;

		// choose their flips randomly too (inside/outside)

		// setup pose

		Pose pose1, pose2;
		pose_from_pdb( pose1, pose1_filename );
		pose_from_pdb( pose2, pose2_filename );

		Size const nres_monomer1( pose1.total_residue() ), nres_monomer2( pose2.total_residue() );
		runtime_assert( nres_monomer1%nrepeat1 == 0 );
		Size const repeatlen1( nres_monomer1/nrepeat1 );

		runtime_assert( nres_monomer2%nrepeat2 == 0 );
		Size const repeatlen2( nres_monomer2/nrepeat2 );

		if ( !design_positions.empty() ) {
			runtime_assert( repeatlen1 == repeatlen2 && nres_monomer1 == nres_monomer2 );
		}

		Vector axis1, axis2, center1, center2;
		get_toroid_axis( pose1, nrepeat1, repeatlen1, axis1, center1 );

		get_toroid_axis( pose2, nrepeat2, repeatlen2, axis2, center2 );

		Real const axis1_factor( 1 ); // no need to flip both, right? //( numeric::random::uniform()<0.5 ) ? -1 : 1 );
		Real const axis2_factor( ( numeric::random::uniform()<0.5 ) ? -1 : 1 );


		Pose pose( pose1 );
		pose.append_residue_by_jump( pose2.residue(1), 1 );
		for ( Size i=2; i<= nres_monomer2; ++i ) pose.append_residue_by_bond( pose2.residue(i) );

		// now virtual residues
		for ( Size r=1; r<= 2; ++r ) { //
			Vector const axis( r == 1 ? axis1 : axis2 ), center( r == 1 ? center1 : center2 ),
				orient( r == 1 ? pose1.residue(1).xyz("CA") : pose2.residue(1).xyz("CA") );
			Real const axis_factor( r == 1 ? axis1_factor : axis2_factor );
			ResidueOP rsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map("VRT")));
			rsd->set_xyz("ORIG", center);
			Vector z( axis_factor * axis ), x( orient - center );
			x = ( x - z * z.dot(x) ).normalized();
			Vector y( z.cross( x ) );
			rsd->set_xyz("X", center+x);
			rsd->set_xyz("Y", center+y);
			pose.append_residue_by_jump( *rsd, 1 );
		}

		// central virtual residues:
		for ( Size r=1; r<= 2; ++r ) {
			ResidueOP rsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map("VRT")));
			rsd->set_xyz( "ORIG", Vector(0,0,0) );
			if ( r == 1 || symm_type == "SP0" ) {
				rsd->set_xyz( "X", Vector( 1,0,0 ) );
				rsd->set_xyz( "Y", Vector( 0,1,0 ) );
			} else if ( r == 2 && symm_type == "SP90" ) { // rotate the second VRT 90 degrees about the y-axis
				rsd->set_xyz( "X", Vector( 0,0,1 ) );
				rsd->set_xyz( "Y", Vector( 0,1,0 ) );
			} else {
				utility_exit_with_message( "unrecognized "+symm_type );
			}
			pose.append_residue_by_jump( *rsd, 1 );
		}

		// now the fold_tree:
		/// now set the proper foldtree
		Size const nres_protein( nres_monomer1 + nres_monomer2 );
		pose.conformation().insert_chain_ending( nres_monomer1 );
		pose.conformation().insert_chain_ending( nres_protein );

		//pose.dump_pdb("built.pdb");

		FoldTree f( pose.total_residue() );

		/// first 2 (fixed) jumps are between monomers and their axial vrt rsds
		f.new_jump( 1, nres_protein+1, nres_monomer1 );
		f.new_jump( nres_monomer1+1, nres_protein+2, nres_protein );

		/// next 2 jumps are from the central to the axial vrt rsds
		/// the  first of these (#3) will be Z-rotation-flexible
		/// the second of these (#4) will be Z-rotation and Z-translation flexible
		f.new_jump( nres_protein+1, nres_protein+3, nres_protein+1 );
		f.new_jump( nres_protein+2, nres_protein+4, nres_protein+2 );

		/// then there's a jump (#5) between the central vrt rsds
		/// this will be Y-translation flexible
		f.new_jump( nres_protein+3, nres_protein+4, nres_protein+3 );
		Size const root( nres_protein + 3 ); //the first of the origin vrt rsds
		f.reorder( root );
		pose.fold_tree( f );
		runtime_assert( pose.fold_tree().num_jump() == 5 );

		/// set an initial jump
		Jump jump( RT( numeric::xyzMatrix_double::identity(), Vector(0,0,0) ) ); // axial vrt rsd also at origin

		pose.set_jump( 3, jump );
		pose.set_jump( 4, jump );

		//
		if ( pymol ) { pymol->pymol_name("start"); pymol->apply( pose ); pose.dump_pdb("start.pdb"); }

		/// randomly spin each monomer about its axis
		for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
			Jump jump( pose.jump( jumpno ) );
			jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, numeric::random::uniform()*360.0 ); // full spin for asym mons
			jump.fold_in_rb_deltas();
			pose.set_jump( jumpno, jump );
		}
		//pose.dump_pdb("spin.pdb");
		if ( pymol ) { pymol->pymol_name("spin"); pymol->apply( pose ); pose.dump_pdb("spin.pdb"); }

		/// randomly shift monomer2 up or down a little bit
		{
			Jump jump( pose.jump( 4 ) );
			jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, numeric::random::gaussian()*5.0 );
			jump.fold_in_rb_deltas();
			pose.set_jump( 4, jump );
		}
		//pose.dump_pdb("shift.pdb");
		if ( pymol ) { pymol->pymol_name("shift"); pymol->apply( pose ); pose.dump_pdb("shift.pdb"); }

		/// now move far away
		{
			Jump jump( pose.jump( 5 ) );
			Real const distance( 200.0+uniform() );
			jump.set_rb_delta( kinematics::Jump::TRANS_Y, 1, distance );
			jump.fold_in_rb_deltas();
			pose.set_jump( 5, jump );
		}

		//pose.dump_pdb("far.pdb");
		if ( pymol ) { pymol->pymol_name("far"); pymol->apply( pose ); pose.dump_pdb("far.pdb");}

		// convert to centroid for slide into contact and centroid simulation
		//
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );
		{ /// now we need to slide into contact

			ScoreFunction scorefxn;
			scorefxn.set_weight( vdw, 1.0 );

			Real start_score( scorefxn( pose ) );
			while ( true ) {
				{
					Jump jump( pose.jump( 5 ) );
					Real const step( -0.75 );
					jump.set_rb_delta( kinematics::Jump::TRANS_Y, 1, step );
					jump.fold_in_rb_deltas();
					pose.set_jump( 5, jump );
				}

				Real const score( scorefxn( pose ) );
				TR.Trace << "slide_into_contact: " << F(9,3,score) << endl;
				if ( score - start_score>5 ) {
					break;
				}
			}

			// step back out a random amount, occasionally back in
			{
				Jump jump( pose.jump( 5 ) );
				Real const step( -0.25+uniform() ); // step in [-0.25,0.75]
				jump.set_rb_delta( kinematics::Jump::TRANS_Y, 1, step );
				jump.fold_in_rb_deltas();
				pose.set_jump( 5, jump );
			}

		} // slide into contact

		//pose.dump_pdb("slide.pdb");
		if ( pymol ) { pymol->pymol_name("slide"); pymol->apply( pose ); pose.dump_pdb("slide.pdb"); }

		for ( Size i=1; i<= nres_protein; ++i ) {
			// make sequence mutations everywhere
			AA const new_aa( aa_from_oneletter_code( random_element( make_vector1(
							'L', 'L', 'L',
							'A', 'A', 'A',
							'S', 'S', 'S',
							'I', 'I',
							'G', 'G',
							'D', 'D','D',
							'K', 'K','K',
							'V', 'V', 'V', 'V'
						) ) ) );
			if ( option[ my_options::self_dock_2c ] ) {
				runtime_assert( nres_monomer1 == nres_monomer2 );
				if ( i<= nres_monomer1 ) {
					make_sequence_change( i, new_aa, pose );
					make_sequence_change( i+nres_monomer1, new_aa, pose );
				}
			} else {
				make_sequence_change( i, new_aa, pose );
			}

		}

		//pose.dump_pdb("randomseq.pdb");


		Real centroid_bsasa14( 0.0 );
		{ // centroid simulation

			(*cen_scorefxn)( pose );

			Size const n_outer( 5 ), n_inner( 500 );

			Real const mc_hitemp( 4 ), mc_lotemp( 0.5 );
			Real const gamma = std::pow( (mc_lotemp/mc_hitemp), 1.0/Real(n_inner*n_outer));
			Real mc_temp( mc_hitemp / gamma );

			protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *cen_scorefxn, mc_temp ) );

			for ( Size n=1; n<= n_outer; ++n ) {
				for ( Size m=1; m<= n_inner; ++m ) {
					mc_temp *= gamma;
					mc->set_temperature( mc_temp );

					//using kinematics::Jump;
					{ // slide
						Real const slide( numeric::random::gaussian() * 0.5 );
						Jump jump( pose.jump( 5 ) );
						jump.set_rb_delta( kinematics::Jump::TRANS_Y, 1, slide );
						jump.fold_in_rb_deltas();
						pose.set_jump( 5, jump );
					}

					{ // shift
						Real const shift( numeric::random::gaussian() * 0.5 );
						Jump jump( pose.jump( 4 ) );
						jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, shift );
						jump.fold_in_rb_deltas();
						pose.set_jump( 4, jump );
					}

					// spins
					for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
						Jump jump( pose.jump( jumpno ) );
						Real const spin( numeric::random::gaussian() * 5.0 );
						jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, spin );
						jump.fold_in_rb_deltas();
						pose.set_jump( jumpno, jump );
					}

					//(*cen_scorefxn)( pose );
					bool const mc_accept( mc->boltzmann( pose ) );
					TR.Trace << "mc_accept " << I(4,n) << I(4,m) << I(3,mc_accept) << F(9,3,mc->last_accepted_score()) <<
						F(9,3,mc->lowest_score() ) << F(9,3,mc->temperature() ) << endl;
					//set_jump_rb_centers( pose );
				}
				mc->recover_low( pose );
			}

			// Real const final_score( (*cen_scorefxn)( pose ) );
			// static Size counter(0);
			// ++counter;
			// string const outfilename( "aftercen_"+lead_zero_string_of( counter,4)+".pdb");
			// cout << "final_score: " << F(9,3,final_score) << ' ' << outfilename << endl;
			//pose.dump_pdb( outfilename );
			//pose.dump_pdb( "after_cen.pdb" );


			devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );


			{ // for screening the models, lets check the sasa after very quickly optimizing the interface sidechains
				pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
				task->set_bump_check( false );
				MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
				setup_interface_movemap( pose, *movemap );
				bools is_flexible( pose.total_residue(), false ), is_designable( pose.total_residue(), true ); // at the moment

				for ( Size i=1; i<= nres_protein; ++i ) is_flexible[i] = ( movemap->get_chi(i) );
				task->restrict_to_repacking();
				task->restrict_to_residues( is_flexible );
				pack::pack_rotamers( pose, *fa_scorefxn, task );
				/// now compute buried sasa
				Real bsasa5(0), bsasa14(0);
				//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;
				Size n_unsat_donors(0), n_unsat_donors_bb(0), n_unsat_acceptors(0), n_unsat_acceptors_bb(0);
				analyze_interface( pose, bsasa5, bsasa14,
					n_unsat_donors, n_unsat_donors_bb,
					n_unsat_acceptors, n_unsat_acceptors_bb );
				centroid_bsasa14 = bsasa14;
				if ( pymol ) { pymol->pymol_name("aftercen"); pymol->apply( pose ); pose.dump_pdb("aftercen.pdb"); }
				TR.Trace << "centroid_bsasa14: " << F(9,3,centroid_bsasa14) << endl;
			}


			// recover the monomer sidechains
			for ( Size i=1; i<= nres_monomer1; ++i ) {
				pose.replace_residue( i, pose1.residue(i), true ); // orient_backbone
			}
			for ( Size i=1; i<= nres_monomer2; ++i ) {
				pose.replace_residue( nres_monomer1+i, pose2.residue(i), true ); // orient_backbone
			}

		} // centroid simulation

		bool const do_design( centroid_bsasa14 >= centroid_bsasa14_threshold );

		if ( do_design ) { // design/relax
			Size const design_cycles( option[ my_options::design_cycles ] );

			for ( Size m=1; m<= design_cycles+1; ++m ) {
				bool const skip_relax( m == 1 );

				if ( !skip_relax ) { // relax
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					setup_interface_movemap( pose, *movemap );
					movemap->set_jump( 1, false ); // was set by setup_interface_movemap
					id::DOF_ID const & spin_z3_id
						( pose.conformation().dof_id_from_torsion_id(id::TorsionID(3,id::JUMP,6)));
					id::DOF_ID const & slide_z4_id
						( pose.conformation().dof_id_from_torsion_id(id::TorsionID(4,id::JUMP,3)));
					id::DOF_ID const & spin_z4_id
						( pose.conformation().dof_id_from_torsion_id(id::TorsionID(4,id::JUMP,6)));
					id::DOF_ID const & slide_y5_id
						( pose.conformation().dof_id_from_torsion_id(id::TorsionID(5,id::JUMP,2)));
					movemap->set( spin_z3_id, true );
					movemap->set( spin_z4_id, true );
					movemap->set( slide_z4_id, true );
					movemap->set( slide_y5_id, true );
					fastrelax.set_movemap( movemap );
					if ( !dry_run() ) fastrelax.apply( pose );
				}

				if ( m>design_cycles ) break; // relax is cheap, so do one extra one at the end...

				{ // design
					pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
					task->initialize_from_command_line();
					task->or_include_current( true );
					// if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					setup_interface_movemap( pose, *movemap );

					bools is_flexible( pose.total_residue(), false ), is_designable( pose.total_residue(), true ); // at the moment
					if ( !design_positions.empty() ) {
						for ( Size i=1; i<= nres_protein; ++i ) {
							Size const rpos( (i-1)%repeatlen1 + 1 );
							if ( !has_element( rpos, design_positions ) ) is_designable[i] = false;
							TR.Trace << "is_designable: " << is_designable[i] << I(6,i) << I(3,rpos) << endl;
						}
					}
					for ( Size i=1; i<= nres_protein; ++i ) is_flexible[i] = ( movemap->get_chi(i) );
					task->restrict_to_residues( is_flexible );
					Size n_repack(0), n_design(0);
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						if ( is_flexible[i] ) {
							if ( is_designable[i] ) {
								++n_design;
							} else {
								++n_repack;
								task->nonconst_residue_task(i).restrict_to_repacking();
							}
						}
					}
					TR.Trace << "n_design: " << n_design << " n_repack: " << n_repack << endl;
					//if ( nodesign ) task->restrict_to_repacking();

					if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
						protocols::task_operations::LimitAromaChi2Operation lp_op;
						lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negatives)
						lp_op.apply( pose, *task );
					}
					Size const nloop( 1 ); // in case linmem ig
					ScoreFunctionOP design_scorefxn(0);
					// if ( ( option[ my_options::use_softrep_for_early_design ] && m<= design_cycles/2 ) ||
					//    option[ my_options::use_softrep_for_design ] ) {
					// before 2/10/2016:
					// design_scorefxn = ScoreFunctionFactory::create_score_function( "soft_rep_design.wts" );
					// adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
					runtime_assert( option[ my_options::design_score_function ].user() );
					design_scorefxn = setup_score_function( option[ my_options::design_score_function ] );
					adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
					// } else {
					//  design_scorefxn = fa_scorefxn; // already adjusted refwts
					// }
					protocols::minimization_packing::PackRotamersMover packmover( design_scorefxn, task, nloop );
					// protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
					if ( !dry_run() ) packmover.apply( pose );
				}
			} // cycles
		} // fullatom simulation


		//Real const relax_simtime( ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes
		Real const final_score( (*fa_scorefxn)( pose ) );

		// compute interface energies
		EnergyMap interface_emap;
		{
			EnergyGraph const & energy_graph( pose.energies().energy_graph() );
			for ( Size pos1 = 1; pos1 <= nres_monomer1; ++pos1 ) {
				runtime_assert( pose.chain(pos1) == 1 );
				for ( utility::graph::Graph::EdgeListConstIter
						ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
						ire = energy_graph.get_node( pos1 )->const_edge_list_end();
						ir != ire; ++ir ) {
					EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
					Size const pos2( edge->get_other_ind( pos1 ) );
					if ( pose.chain(pos2) == 2 ) {
						edge->add_to_energy_map( interface_emap );
					}
				}
			}
		}
		Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );

		Real bsasa5, bsasa14;
		Size n_unsat_donors, n_unsat_donors_bb, n_unsat_acceptors, n_unsat_acceptors_bb;
		{
			Pose interface_pose;
			for ( Size i=1; i<= nres_monomer1; ++i ) interface_pose.append_residue_by_bond( pose.residue( i ) );
			interface_pose.append_residue_by_jump( pose.residue( nres_monomer1+1 ), 1 );
			for ( Size i=nres_monomer1+2; i<= nres_protein; ++i ) {
				interface_pose.append_residue_by_bond( pose.residue( i ) );
			}
			interface_pose.conformation().insert_chain_ending( nres_monomer1 );
			//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;

			analyze_interface( interface_pose, bsasa5, bsasa14,
				n_unsat_donors, n_unsat_donors_bb,
				n_unsat_acceptors, n_unsat_acceptors_bb );
		}

		if ( bsasa14<1 ) bsasa14 = 1; // prevent division by 0

		/// what do we filter on?
		///
		/// -- interface energy
		/// -- total energy per residue
		///
		Real const residue_energy( final_score / nres_protein );
		Real const interface_quality( 100.0 * interface_energy / bsasa14 );
		Real const bsasa_ratio( bsasa5 / bsasa14 );
		bool const small_interface( bsasa14 < 500 );


		bool passed_score_filter
			( do_design &&
			append_score_to_scorefile_and_filter( worktag, residue_energy, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		// now also filter on interface energy
		passed_score_filter = passed_score_filter ||
			( do_design &&
			append_score_to_scorefile_and_filter( worktag+"_int", interface_energy, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );

		passed_score_filter = passed_score_filter && !small_interface; // no small interfaces

		//string const outfilename( output_tag() + "_dock_N"+ lead_zero_string_of( n,4 )+".pdb");
		string const outfilename( output_tag() + "twocomp_S"+symm_type+
			"_N"+ lead_zero_string_of( n,4)+".pdb" );

		Real const simtime = ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

		ostringstream out;
		out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
			//" simtime: " << F(9,3,simtime) <<
			" residue_energy: " << F(9,3,residue_energy) <<
			" interface_energy: " << F(9,3,interface_energy) <<
			" interface_quality: " << F(9,3,interface_quality) <<
			" pose1_filename: " << pose1_filename <<
			" pose2_filename: " << pose2_filename <<
			" passed_score_filter: " << passed_score_filter <<
			" bsasa5: " << F(9,3,bsasa5) <<
			" bsasa14: " << F(9,3,bsasa14) <<
			" bsasa_ratio: " << F(9,3,bsasa_ratio) <<
			" do_design: " << do_design <<
			" centroid_bsasa14: " << F(9,3,centroid_bsasa14) <<
			" simtime: " << F(9,3,simtime) <<
			" n_unsat_donors: " << n_unsat_donors <<
			" n_unsat_donors_bb: " << n_unsat_donors_bb <<
			" n_unsat_acceptors: " << n_unsat_acceptors <<
			" n_unsat_acceptors_bb: " << n_unsat_acceptors_bb <<
			" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
			" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

		string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

		TR.Trace << out.str() << endl;

		basic::prof_show_oneliner();

		if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
			append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
		}
		check_simtime();
	} // nstruct
	signal_that_simfile_is_done( simfile );
	signal_that_job_is_done();


}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
planar_one_component_test()
{
	protocols::moves::PyMOLMoverOP pymol(0);
	if ( option[ my_options::show_pymol ] ) {
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol->keep_history(true);
	}
	TR.Trace << "show_pymol: " << option[ my_options::show_pymol ] << ' ' << pymol << endl;

	using namespace kinematics;

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );
	Real const centroid_bsasa14_threshold( option[ my_options::centroid_bsasa14_threshold ] );
	//Size const design_cycles( option[ my_options::design_cycles ] );

	// ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	ScoreFunctionOP fa_scorefxn( get_score_function_from_command_line() ), // was: core::scoring::get_score_function() ),
		cen_scorefxn( new ScoreFunction() );

	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );
	runtime_assert( !pose::symmetry::is_symmetric( *cen_scorefxn ) );

	// # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
	// env     1.0
	// pair    1.0
	// cbeta   1.0
	// vdw     1.0
	// rg      3.0
	// cenpack 1.0
	// hs_pair 1.0
	// ss_pair 1.0
	// rsigma  1.0
	// sheet   1.0

	cen_scorefxn->set_weight( env, 1.0 );
	cen_scorefxn->set_weight( scoring::pair, 1.0 );
	cen_scorefxn->set_weight( cbeta, 1.0 );
	cen_scorefxn->set_weight( vdw, 1.0 );
	cen_scorefxn->set_weight( cenpack, 1.0 ); // not sure about this guy...

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	Size const nrepeat( option[ my_options::nrepeat ] );

	Sizes design_positions;
	if ( option[ my_options::design_positions ].user() ) {
		design_positions = option[ my_options::design_positions ]();
	}

	// strings symm_types( option[ my_options::symm_types ]() );
	// numeric::random::random_permutation( symm_types, numeric::random::rg() );
	// string const symm_type( option[ my_options::symm_type ] );//"SP2" );
	// runtime_assert( symm_type == "SP0" || symm_type == "SP90" ); // 0-degree or 90-degree angle between axes


	strings files( start_files() );
	numeric::random::random_permutation( files, numeric::random::rg() );

	foreach_( string const filename, files ) {

		// simfile
		string const simfile( shared_output_tag()+"_p1c_"+filebase(filename)+".work" );
		if ( simfile_is_done( simfile ) ) {
			check_if_job_is_done();
			continue;
		}

		Pose pdb_pose;
		pose_from_pdb( pdb_pose, filename );

		runtime_assert( pdb_pose.residue( pdb_pose.total_residue() ).is_protein() ); // no VRT residue

		Size const nres_monomer( pdb_pose.total_residue() );

		runtime_assert( nres_monomer%nrepeat == 0 );
		Size const repeatlen( nres_monomer/nrepeat );

		Vector center, axis;
		get_toroid_axis( pdb_pose, nrepeat, repeatlen, axis, center );

		for ( Size n=1; n<= nstruct(); ++n ) { // all independent
			clock_t starttime( clock() );

			string const worktag( "tmp" );

			string const symm_type( random_element( split_to_vector1("A P")));


			Pose pose( pdb_pose );
			pose.append_residue_by_jump( pdb_pose.residue(1), 1 );
			for ( Size i=2; i<= nres_monomer; ++i ) pose.append_residue_by_bond( pdb_pose.residue(i) );

			// now virtual residues
			// first one for each monomer which is in the center, with Z parallel to the axis
			// these will move along with the monomers
			//

			for ( Size r=1; r<= 2; ++r ) { //
				Vector const orient( pdb_pose.residue(nres_monomer/2).xyz("CA") ); // x will point toward the midpoint of the toroid
				ResidueOP rsd
					( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map("VRT")));
				rsd->set_xyz("ORIG", center);
				Vector z( axis ), x( orient - center );
				x = ( x - z * z.dot(x) ).normalized();
				Vector y( z.cross( x ) );
				rsd->set_xyz("X", center+x);
				rsd->set_xyz("Y", center+y);
				pose.append_residue_by_jump( *rsd, 1 );
			}


			// central virtual residues: position at the origin, with Y-axes in opposite directions (+y and -y direcions)
			// these will stay completely fixed
			//
			for ( Size r=1; r<= 2; ++r ) {
				ResidueOP rsd
					( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map("VRT")));
				rsd->set_xyz( "ORIG", Vector(0,0,0) );
				if ( r == 1 ) {
					rsd->set_xyz( "X", Vector( 1,0,0 ) );
					rsd->set_xyz( "Y", Vector( 0,1,0 ) );
				} else {
					if ( symm_type == "P" ) { // toroid axes are parallel to one another
						// rotate 180 degrees around z
						rsd->set_xyz( "X", Vector( -1,  0, 0 ) );
						rsd->set_xyz( "Y", Vector(  0, -1, 0 ) );
					} else {
						// rotate 180 degrees around x
						runtime_assert( symm_type == "A" );
						rsd->set_xyz( "X", Vector(  1,  0, 0 ) );
						rsd->set_xyz( "Y", Vector(  0, -1, 0 ) );
					}
				}
				pose.append_residue_by_jump( *rsd, 1 );
			}

			// now the fold_tree:
			/// now set the proper foldtree
			Size const nres_protein( 2*nres_monomer );
			pose.conformation().insert_chain_ending( nres_monomer );
			pose.conformation().insert_chain_ending( nres_protein );

			//pose.dump_pdb("built.pdb");

			FoldTree f( pose.total_residue() );

			/// first 2 (fixed) jumps are between monomers and their axial vrt rsds
			f.new_jump( 1, nres_protein+1, nres_monomer );
			f.new_jump( nres_monomer+1, nres_protein+2, nres_protein );

			/// next 2 jumps are from the central to the axial vrt rsds
			/// these will both be Y- and X- translation flexible (and Z-translation if symm_type=="A"), and Z-rotation flexible
			f.new_jump( nres_protein+1, nres_protein+3, nres_protein+1 );
			f.new_jump( nres_protein+2, nres_protein+4, nres_protein+2 );

			/// then there's a jump (#5) between the central vrt rsds
			/// this will be completely frozen
			f.new_jump( nres_protein+3, nres_protein+4, nres_protein+3 );
			Size const root( nres_protein + 3 ); //the first of the origin vrt rsds
			f.reorder( root );
			pose.fold_tree( f );
			runtime_assert( pose.fold_tree().num_jump() == 5 );

			/// set an initial jump
			Jump jump( RT( numeric::xyzMatrix_double::identity(), Vector(0,0,0) ) ); // axial vrt rsd also at origin

			pose.set_jump( 3, jump );
			pose.set_jump( 4, jump );

			//
			if ( pymol ) { pymol->pymol_name("start"); pymol->apply( pose ); pose.dump_pdb("start.pdb"); }

			/// randomly spin each monomer about its axis
			/// how much do we need to search over?
			Real const max_spin( 360.0 / nrepeat ), spin_angle( numeric::random::uniform()*max_spin );
			for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
				Jump jump( pose.jump( jumpno ) );
				jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, spin_angle ); // full spin for asym mons
				jump.fold_in_rb_deltas();
				pose.set_jump( jumpno, jump );
			}
			if ( pymol ) { pymol->pymol_name("spin"); pymol->apply( pose ); pose.dump_pdb("spin.pdb"); }

			/// randomly shift monomer2 up or down a little bit
			if ( symm_type == "A" ) {
				Real const z_shift( numeric::random::gaussian()*5.0 );
				for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
					Jump jump( pose.jump( jumpno ) );
					jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, z_shift );
					jump.fold_in_rb_deltas();
					pose.set_jump( jumpno, jump );
				}
			}
			if ( pymol ) { pymol->pymol_name("shift"); pymol->apply( pose ); pose.dump_pdb("shift.pdb"); }

			/// now move far away
			Real const far_distance( 100.0+uniform() );
			{
				for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
					Jump jump( pose.jump( jumpno ) );
					jump.set_rb_delta( kinematics::Jump::TRANS_Y, 1, far_distance );
					jump.fold_in_rb_deltas();
					pose.set_jump( jumpno, jump );
				}
			}

			if ( pymol ) { pymol->pymol_name("far"); pymol->apply( pose ); pose.dump_pdb("far.pdb");}

			// convert to centroid for slide into contact and centroid simulation
			//
			devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );
			{ /// now we need to slide into contact
				Real const slide_stepsize( -0.4 ); // times 2 for total distance change
				Size const max_slide_steps( int( (3*far_distance)/slide_stepsize ) ); // this is a big over-estimate

				ScoreFunction scorefxn;
				scorefxn.set_weight( vdw, 1.0 );

				Real start_score( scorefxn( pose ) );
				Size nsteps(0);
				while ( true ) {
					++nsteps;
					if ( nsteps>max_slide_steps ) break;
					for( Size jumpno=3; jumpno<= 4; ++jumpno ) {
						Jump jump( pose.jump( jumpno ) );
						jump.set_rb_delta( kinematics::Jump::TRANS_Y, 1, slide_stepsize );
						jump.fold_in_rb_deltas();
						pose.set_jump( jumpno, jump );
					}

					Real const score( scorefxn( pose ) );
					TR.Trace << "slide_into_contact: " << F(9,3,score) << endl;
					if ( score - start_score>5 ) {
						break;
					}
				}

				if ( nsteps > max_slide_steps ) {
					--n;
					continue; // go to another nstruct
				}

				// step back out a random amount, occasionally back in
				Real const step( -0.25+uniform() ); // step in [-0.25,0.75]
				for( Size jumpno=3; jumpno<= 4; ++jumpno ) {
					Jump jump( pose.jump( jumpno ) );
					jump.set_rb_delta( kinematics::Jump::TRANS_Y, 1, step );
					jump.fold_in_rb_deltas();
					pose.set_jump( jumpno, jump );
				}
			} // slide into contact

			//pose.dump_pdb("slide.pdb");
			if ( pymol ) { pymol->pymol_name("slide"); pymol->apply( pose ); pose.dump_pdb("slide.pdb"); }

			for ( Size i=1; i<= nres_monomer; ++i ) {
				// make sequence mutations everywhere
				AA const new_aa( aa_from_oneletter_code( random_element( make_vector1(
								'L', 'L', 'L',
								'A', 'A', 'A',
								'S', 'S', 'S',
								'I', 'I',
								'G', 'G',
								'D', 'D','D',
								'K', 'K','K',
								'V', 'V', 'V', 'V'
							) ) ) );
				make_sequence_change( i, new_aa, pose );
				make_sequence_change( i+nres_monomer, new_aa, pose );
			}

			//pose.dump_pdb("randomseq.pdb");


			Real centroid_bsasa14( 0.0 );
			{ // centroid simulation

				(*cen_scorefxn)( pose );

				Size const n_outer( dry_run() ? 1 : 5 ), n_inner( dry_run() ? 5 : 500 );

				Real const mc_hitemp( 4 ), mc_lotemp( 0.5 );
				Real const gamma = std::pow( (mc_lotemp/mc_hitemp), 1.0/Real(n_inner*n_outer));
				Real mc_temp( mc_hitemp / gamma );

				protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *cen_scorefxn, mc_temp ) );

				for ( Size n=1; n<= n_outer; ++n ) {
					for ( Size m=1; m<= n_inner; ++m ) {
						mc_temp *= gamma;
						mc->set_temperature( mc_temp );

						//using kinematics::Jump;
						{ // slide
							Real const slide( numeric::random::gaussian() * 0.25 );
							for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
								Jump jump( pose.jump( jumpno ) );
								jump.set_rb_delta( kinematics::Jump::TRANS_Y, 1, slide );
								jump.fold_in_rb_deltas();
								pose.set_jump( jumpno, jump );
							}
						}

						if ( symm_type == "A" ) { // shift
							Real const shift( numeric::random::gaussian() * 0.25 );
							for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
								Jump jump( pose.jump( jumpno ) );
								jump.set_rb_delta( kinematics::Jump::TRANS_Z, 1, shift );
								jump.fold_in_rb_deltas();
								pose.set_jump( jumpno, jump );
							}
						}

						// spins
						Real const spin( numeric::random::gaussian() * 2.5 );
						for ( Size jumpno=3; jumpno<= 4; ++jumpno ) {
							Jump jump( pose.jump( jumpno ) );
							jump.set_rb_delta( kinematics::Jump::ROT_Z, 1, spin );
							jump.fold_in_rb_deltas();
							pose.set_jump( jumpno, jump );
						}

						//(*cen_scorefxn)( pose );
						bool const mc_accept( mc->boltzmann( pose ) );
						TR.Trace << "mc_accept " << I(4,n) << I(4,m) << I(3,mc_accept) << F(9,3,mc->last_accepted_score()) <<
							F(9,3,mc->lowest_score() ) << F(9,3,mc->temperature() ) << endl;
						//set_jump_rb_centers( pose );
					}
					mc->recover_low( pose );
				}

				// Real const final_score( (*cen_scorefxn)( pose ) );
				// static Size counter(0);
				// ++counter;
				// string const outfilename( "aftercen_"+lead_zero_string_of( counter,4)+".pdb");
				// cout << "final_score: " << F(9,3,final_score) << ' ' << outfilename << endl;
				//pose.dump_pdb( outfilename );
				//pose.dump_pdb( "after_cen.pdb" );


				devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );


				{ // for screening the models, lets check the sasa after very quickly optimizing the interface sidechains
					// pose is still asym at this point...
					conformation::symmetry::turn_symmetry_off();
					ScoreFunctionOP pack_scorefxn( get_score_function_from_command_line() );
					runtime_assert( !pose::symmetry::is_symmetric( *pack_scorefxn ) );
					pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
					task->set_bump_check( false );
					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					setup_interface_movemap( pose, *movemap );
					bools is_flexible( pose.total_residue(), false ), is_designable( pose.total_residue(), true ); // at the moment

					for ( Size i=1; i<= nres_protein; ++i ) is_flexible[i] = ( movemap->get_chi(i) );
					task->restrict_to_repacking();
					task->restrict_to_residues( is_flexible );
					pack::pack_rotamers( pose, *pack_scorefxn, task );
					/// now compute buried sasa
					Real bsasa5(0), bsasa14(0);
					//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;
					Size n_unsat_donors(0), n_unsat_donors_bb(0), n_unsat_acceptors(0), n_unsat_acceptors_bb(0);
					analyze_interface( pose, bsasa5, bsasa14,
						n_unsat_donors, n_unsat_donors_bb,
						n_unsat_acceptors, n_unsat_acceptors_bb );
					centroid_bsasa14 = bsasa14;
					if ( pymol ) { pymol->pymol_name("aftercen"); pymol->apply( pose ); pose.dump_pdb("aftercen.pdb"); }
					TR.Trace << "centroid_bsasa14: " << F(9,3,centroid_bsasa14) << endl;
					conformation::symmetry::turn_symmetry_on();
				}


				// recover the monomer sidechains
				for ( Size i=1; i<= nres_monomer; ++i ) {
					pose.replace_residue( i, pdb_pose.residue(i), true ); // orient_backbone
					pose.replace_residue( nres_monomer+i, pdb_pose.residue(i), true ); // orient_backbone
				}

			} // centroid simulation

			bool const do_design( centroid_bsasa14 >= centroid_bsasa14_threshold );

			if ( !do_design ) {
				cout << "not designing small interface: " << centroid_bsasa14 << ' ' << centroid_bsasa14_threshold << endl;
				--n;
				continue;
			}

			// now going to make the pose symmetric
			conformation::symmetry::SymmetryInfo symminfo;
			for ( Size i=1; i<= nres_monomer; ++i ) {
				symminfo.add_bb_clone ( i, i+nres_monomer );
				symminfo.add_chi_clone( i, i+nres_monomer );
			}

			symminfo.add_jump_clone( 1, 2, 0.0 );
			symminfo.add_jump_clone( 3, 4, 0.0 );

			{
				using core::conformation::symmetry::SymDof;
				map< Size, SymDof > symdofs;
				SymDof symdof;
				if ( symm_type == "A" ) {
					symdof.read( "x y z angle_z" );
				} else {
					symdof.read( "x y angle_z" );
				}
				symdofs[ 3 ] = symdof;
				symminfo.set_dofs( symdofs );
			}

			symminfo.num_virtuals( 4 );
			symminfo.set_use_symmetry( true );
			symminfo.set_flat_score_multiply( pose.total_residue(), 1 );
			for ( Size i=1; i<= nres_monomer; ++i ) symminfo.set_score_multiply( i, 2 );
			symminfo.update_score_multiply_factor();

			pose::symmetry::make_symmetric_pose( pose, symminfo ); // now symmetric

			pose.set_jump( 1, pose.jump(1) );
			pose.set_jump( 3, pose.jump(3) );
			if ( pymol ) { pymol->pymol_name("aftersymm"); pymol->apply( pose ); pose.dump_pdb("aftersymm.pdb"); }


			if ( do_design ) { // design/relax
				Size const design_cycles( option[ my_options::design_cycles ] );

				for ( Size m=1; m<= design_cycles+1; ++m ) {
					bool const skip_relax( m == 1 );

					if ( !skip_relax ) { // relax
						protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

						MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
						setup_interface_movemap( pose, *movemap );
						movemap->set_jump( 1, false ); // was set by setup_interface_movemap
						movemap->set_jump( 2, false );
						movemap->set_jump( 3, true );
						movemap->set_jump( 4, false );
						movemap->set_jump( 5, false );
						// I think the symdofs take care of all this:
						// id::DOF_ID const & spin_z3_id
						// 	( pose.conformation().dof_id_from_torsion_id(id::TorsionID(3,id::JUMP,6)));
						// id::DOF_ID const & spin_z4_id
						// 	( pose.conformation().dof_id_from_torsion_id(id::TorsionID(4,id::JUMP,6)));
						// id::DOF_ID const & slide_z3_id
						// 	( pose.conformation().dof_id_from_torsion_id(id::TorsionID(3,id::JUMP,3)));
						// id::DOF_ID const & slide_z4_id
						// 	( pose.conformation().dof_id_from_torsion_id(id::TorsionID(4,id::JUMP,3)));
						// id::DOF_ID const & slide_y3_id
						// 	( pose.conformation().dof_id_from_torsion_id(id::TorsionID(3,id::JUMP,2)));
						// id::DOF_ID const & slide_y4_id
						// 	( pose.conformation().dof_id_from_torsion_id(id::TorsionID(4,id::JUMP,2)));
						// movemap->set( spin_z3_id, true );
						// movemap->set( spin_z4_id, true );
						// if ( symm_type=="A" ) {
						// 	movemap->set( slide_z3_id, true );
						// 	movemap->set( slide_z4_id, true );
						// }
						// movemap->set( slide_y3_id, true );
						// movemap->set( slide_y4_id, true );
						fastrelax.set_movemap( movemap );
						if ( !dry_run() ) fastrelax.apply( pose );
					}

					if ( m>design_cycles ) break; // relax is cheap, so do one extra one at the end...

					{ // design
						pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
						task->initialize_from_command_line();
						task->or_include_current( true );
						// if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

						MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
						setup_interface_movemap( pose, *movemap );

						bools is_flexible( pose.total_residue(), false ), is_designable( pose.total_residue(), true ); // at the moment
						if ( !design_positions.empty() ) {
							for ( Size i=1; i<= nres_protein; ++i ) {
								Size const rpos( (i-1)%repeatlen + 1 );
								if ( !has_element( rpos, design_positions ) ) is_designable[i] = false;
								TR.Trace << "is_designable: " << is_designable[i] << I(6,i) << I(3,rpos) << endl;
							}
						}
						for ( Size i=1; i<= nres_monomer; ++i ) is_flexible[i] = ( movemap->get_chi(i) );
						task->restrict_to_residues( is_flexible );
						Size n_repack(0), n_design(0);
						for ( Size i=1; i<= pose.total_residue(); ++i ) {
							if ( is_flexible[i] ) {
								if ( is_designable[i] ) {
									++n_design;
								} else {
									++n_repack;
									task->nonconst_residue_task(i).restrict_to_repacking();
								}
							}
						}
						TR.Trace << "n_design: " << n_design << " n_repack: " << n_repack << endl;
						//if ( nodesign ) task->restrict_to_repacking();

						if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
							protocols::task_operations::LimitAromaChi2Operation lp_op;
							lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negatives)
							lp_op.apply( pose, *task );
						}
						Size const nloop( 1 ); // in case linmem ig
						ScoreFunctionOP design_scorefxn(0);
						// if ( ( option[ my_options::use_softrep_for_early_design ] && m<= design_cycles/2 ) ||
						//    option[ my_options::use_softrep_for_design ] ) {
						// before 2/10/2016:
						// design_scorefxn = ScoreFunctionFactory::create_score_function( "soft_rep_design.wts" );
						// adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
						runtime_assert( option[ my_options::design_score_function ].user() );
						design_scorefxn = setup_score_function( option[ my_options::design_score_function ] );
						runtime_assert( pose::symmetry::is_symmetric( *design_scorefxn ) );
						adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
						// } else {
						//  design_scorefxn = fa_scorefxn; // already adjusted refwts
						// }
						//protocols::minimization_packing::PackRotamersMover packmover( design_scorefxn, task, nloop );
						protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
						if ( !dry_run() ) packmover.apply( pose );
					}
				} // cycles
			} // fullatom simulation


			//Real const relax_simtime( ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes
			Real const final_score( (*fa_scorefxn)( pose ) );

			// compute interface energies
			EnergyMap interface_emap;
			{
				EnergyGraph const & energy_graph( pose.energies().energy_graph() );
				for ( Size pos1 = 1; pos1 <= nres_monomer; ++pos1 ) {
					runtime_assert( pose.chain(pos1) == 1 );
					for ( utility::graph::Graph::EdgeListConstIter
									ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
									ire = energy_graph.get_node( pos1 )->const_edge_list_end();
								ir != ire; ++ir ) {
						EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
						Size const pos2( edge->get_other_ind( pos1 ) );
						if ( pose.chain(pos2) == 2 ) {
							edge->add_to_energy_map( interface_emap );
						}
					}
				}
			}
			Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );

			Real bsasa5, bsasa14;
			Size n_unsat_donors, n_unsat_donors_bb, n_unsat_acceptors, n_unsat_acceptors_bb;
			{
				Pose interface_pose;
				for ( Size i=1; i<= nres_monomer; ++i ) interface_pose.append_residue_by_bond( pose.residue( i ) );
				interface_pose.append_residue_by_jump( pose.residue( nres_monomer+1 ), 1 );
				for ( Size i=nres_monomer+2; i<= nres_protein; ++i ) {
					interface_pose.append_residue_by_bond( pose.residue( i ) );
				}
				interface_pose.conformation().insert_chain_ending( nres_monomer );
				//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;

				analyze_interface( interface_pose, bsasa5, bsasa14,
					n_unsat_donors, n_unsat_donors_bb,
					n_unsat_acceptors, n_unsat_acceptors_bb );
			}

			if ( bsasa14<1 ) bsasa14 = 1; // prevent division by 0

			/// what do we filter on?
			///
			/// -- interface energy
			/// -- total energy per residue
			///
			Real const residue_energy( final_score / nres_protein );
			Real const interface_quality( 100.0 * interface_energy / bsasa14 );
			Real const bsasa_ratio( bsasa5 / bsasa14 );
			bool const small_interface( bsasa14 < 500 );


			bool passed_score_filter
				( do_design &&
					append_score_to_scorefile_and_filter( worktag, residue_energy, score_filter_acceptance_rate,
						score_filter_pass_early, simfile ) );

			// now also filter on interface energy
			passed_score_filter = passed_score_filter ||
				( do_design &&
					append_score_to_scorefile_and_filter( worktag+"_int", interface_energy, score_filter_acceptance_rate,
						score_filter_pass_early, simfile ) );

			passed_score_filter = passed_score_filter && !small_interface; // no small interfaces

			//string const outfilename( output_tag() + "_dock_N"+ lead_zero_string_of( n,4 )+".pdb");
			string const outfilename( output_tag() + "onecomp_S"+symm_type+
				"_N"+ lead_zero_string_of( n,4)+".pdb" );

			Real const simtime = ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

			string mutstring;
			Size nmutations(0);
			{
				ostringstream out;
				for ( Size i=1; i<= nres_monomer; ++i ) {
					char const name1( pose.residue(i).name1() ), natname1( pdb_pose.residue(i).name1() );
					runtime_assert( pose.residue(i).name1() == pose.residue(i+nres_monomer).name1()); // confirm symm
					if ( name1 != natname1 ) {
						++nmutations;
						if ( out.str().size() ) out << ',';
						out << i << ':' << natname1 << name1;
					}
				}
				mutstring = out.str();
			}


			ostringstream out;
			out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
				//" simtime: " << F(9,3,simtime) <<
				" residue_energy: " << F(9,3,residue_energy) <<
				" interface_energy: " << F(9,3,interface_energy) <<
				" interface_quality: " << F(9,3,interface_quality) <<
				" start_filename: " << filename <<
				" passed_score_filter: " << passed_score_filter <<
				" bsasa5: " << F(9,3,bsasa5) <<
				" bsasa14: " << F(9,3,bsasa14) <<
				" bsasa_ratio: " << F(9,3,bsasa_ratio) <<
				" do_design: " << do_design <<
				" centroid_bsasa14: " << F(9,3,centroid_bsasa14) <<
				" nmutations: " << nmutations <<
				" mutstring: " << mutstring <<
				" simtime: " << F(9,3,simtime) <<
				" n_unsat_donors: " << n_unsat_donors <<
				" n_unsat_donors_bb: " << n_unsat_donors_bb <<
				" n_unsat_acceptors: " << n_unsat_acceptors <<
				" n_unsat_acceptors_bb: " << n_unsat_acceptors_bb <<
				" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
				" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

			string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

			TR.Trace << out.str() << endl;

			basic::prof_show_oneliner();

			if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
				append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
			}
			check_simtime();
		} // nstruct
		signal_that_simfile_is_done( simfile ); // kind of dumb since we are using our own nstruct
	} // files
	signal_that_job_is_done();

}


///////////////////////////////////////////////////////////////////////////////
void
spin_coords(
	Size const spin,
	Size const repeatlen,
	Vectors & coords
)
{
	runtime_assert( spin<repeatlen );
	runtime_assert( coords.size()%repeatlen == 0 );


	for ( Size i=0; i< spin; ++i ) {
		for ( Size j=1; j<= repeatlen; ++j ) {
			coords.push_back( coords.front() );
			coords.erase( coords.begin() );
		}
	}

}

///////////////////////////////////////////////////////////////////////////////

void
compare_docked_toroid_coords(
	Size const nrepeat,
	Size const repeatlen,
	Vectors const & coords1,
	Vectors const & coords2_in,
	Size & chain1_spin,
	Size & chain2_spin,
	bool & flip_chains,
	Real & min_rmsd
)
{
	Size const nres_monomer( nrepeat * repeatlen );

	min_rmsd = 1e6;
	chain1_spin = chain2_spin = 0;
	flip_chains = false;

	for ( Size spin1=0; spin1< nrepeat; ++spin1 ) {
		Vectors coords2_chain1( coords2_in );
		coords2_chain1.erase( coords2_chain1.begin() + nres_monomer, coords2_chain1.end() );
		runtime_assert( coords2_chain1.size() == nres_monomer );
		spin_coords( spin1, repeatlen, coords2_chain1 );
		runtime_assert( coords2_chain1.size() == nres_monomer );

		for ( Size spin2=0; spin2< nrepeat; ++spin2 ) {

			Vectors coords2_chain2( coords2_in );
			coords2_chain2.erase( coords2_chain2.begin(), coords2_chain2.begin()+nres_monomer );
			runtime_assert( coords2_chain2.size() == nres_monomer );
			spin_coords( spin2, repeatlen, coords2_chain2 );
			runtime_assert( coords2_chain2.size() == nres_monomer );

			for ( Size ch=1; ch<= 2; ++ch ) {
				Vectors coords2;
				if ( ch == 1 ) {
					for ( Size i=1; i<= nres_monomer; ++i ) coords2.push_back( coords2_chain1[i] );
					for ( Size i=1; i<= nres_monomer; ++i ) coords2.push_back( coords2_chain2[i] );
				} else {
					for ( Size i=1; i<= nres_monomer; ++i ) coords2.push_back( coords2_chain2[i] );
					for ( Size i=1; i<= nres_monomer; ++i ) coords2.push_back( coords2_chain1[i] );
				}
				Real const rmsd( numeric::model_quality::calc_rms( coords1, coords2 ) );
				if ( rmsd < min_rmsd ) {
					min_rmsd = rmsd;
					chain1_spin = spin1;
					chain2_spin = spin2;
					flip_chains = ( ch == 2 );
					TR.Trace << "min_rmsd: " << F(9,3,min_rmsd) << I(4,spin1) << I(4,spin2) << I(2,ch) << endl;
				}
			}
		}
	}
}



///////////////////////////////////////////////////////////////////////////////
class ToroidCoordsDistanceMetric {
public:
	Real
	operator()( Vectors const & coords1, Vectors const & coords2 ) const
	{
		Size chain1_spin, chain2_spin;
		bool flip_chains;
		Real rmsd;
		Size const nrepeat( 6 ); // dont hardcode this
		Size const repeatlen( 35 ); // this either
		compare_docked_toroid_coords( nrepeat, repeatlen, coords1, coords2, chain1_spin, chain2_spin, flip_chains, rmsd );
		return rmsd;
	}
};


///////////////////////////////////////////////////////////////////////////////
void
read_cage_designs_test()
{
	ScoreFunctionOP fa_scorefxn( get_score_function() );

	Real min_threshold, max_threshold;
	Size min_cluster_size, min_top_cluster_size, try_top_cluster_size, max_top_cluster_size,
		max_clusters, max_decoys_per_cluster_pdbfile;
	bool const do_clustering( option[ my_options::clustering_params ].user() );

	if ( do_clustering ) {
		parse_clustering_params( min_threshold, max_threshold, min_cluster_size, min_top_cluster_size, try_top_cluster_size,
			max_top_cluster_size, max_clusters, max_decoys_per_cluster_pdbfile );
	}

	string const repeatseg("ddaaaaaaaaaaaaaabbbcccccccccccccccd");
	//string const repeatseg("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGGW");

	Size const nrepeat( 6 ), repeatlen( 35 ), nres_monomer( nrepeat*repeatlen );

	strings const files( start_files() );

	Pose template_pose;
	pose_from_pdb( template_pose, option[ my_options::template_pdb ]() );

	vector1< Vectors > all_coords;
	//PoseOPs all_poses;
	strings all_filenames;
	// Reals all_int_energies, all_tot_energies, all_bsasa14s, all_bsasa_ratios;
	vector1< map< string, Real > > all_scores;

	bool const calc_sasa( true );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		runtime_assert( pose.total_residue() == 2 * nrepeat * repeatlen );

		// count sequence mutations
		Size n_mismatches( 0 ), n_trp( 0 );
		for ( Size i=1; i<= nres_monomer; ++i ) {
			Size const i_rpos( (i-1)%repeatlen+1 );
			Size const i_repeatno( (i-1)/repeatlen+1 );
			char const desaa( pose.sequence()[i-1] );
			if ( desaa == 'W' ) ++n_trp;
			string repeataas;
			for ( Size k=0; k< nrepeat; ++k ) {
				repeataas += template_pose.residue( k*repeatlen + i_rpos ).name1();
			}
			if ( repeataas.find( desaa ) == string::npos ) {
				// how close are we to partner 2
				Real mindis2( 1e6 );
				Residue const & irsd( pose.residue(i) );
				for ( Size j=nres_monomer+1; j<= 2*nres_monomer; ++j ) {
					Residue const & jrsd( pose.residue(j) );
					for ( Size ii=1; ii<= irsd.nheavyatoms(); ++ii ) {
						for ( Size jj=1; jj<= jrsd.nheavyatoms(); ++jj ) {
							mindis2 = min( mindis2, irsd.xyz(ii).distance_squared(jrsd.xyz(jj)));
						}
					}
				}
				cout << "mismatch: " << F(9,3,sqrt(mindis2)) << ' ' << i_repeatno << repeatseg[i_rpos-1] <<
					I(3,i_rpos) << I(4,i) << ' ' <<
					pose.residue(i).name1() << ' ' << repeataas << ' ' << files[fi] << endl;
				++n_mismatches;
			}
		}

		Real bsasa5(0), bsasa14(0);
		//vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_acceptors;
		Size n_unsat_donors(0), n_unsat_donors_bb(0), n_unsat_acceptors(0), n_unsat_acceptors_bb(0);
		if ( calc_sasa ) {
			analyze_interface( pose, bsasa5, bsasa14,
				n_unsat_donors, n_unsat_donors_bb,
				n_unsat_acceptors, n_unsat_acceptors_bb );
		}

		Real const final_score( (*fa_scorefxn)( pose ) );

		// compute interface energies
		EnergyMap interface_emap;
		{
			EnergyGraph const & energy_graph( pose.energies().energy_graph() );
			for ( Size pos1 = 1; pos1 <= nres_monomer; ++pos1 ) {
				runtime_assert( pose.chain(pos1) == 1 );
				for ( utility::graph::Graph::EdgeListConstIter
						ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
						ire = energy_graph.get_node( pos1 )->const_edge_list_end();
						ir != ire; ++ir ) {
					EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
					Size const pos2( edge->get_other_ind( pos1 ) );
					if ( pose.chain(pos2) == 2 ) {
						edge->add_to_energy_map( interface_emap );
					}
				}
			}
		}
		Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );

		Real bsasa_ratio(0.0);
		if ( bsasa14 >1e-3 ) bsasa_ratio = bsasa5 / bsasa14;

		cout << "final_scores: " << F(9,3,final_score) << ' ' << files[fi] <<
			" interface_energy: " << F(9,3,interface_energy) <<
			" n_mismatches: " << n_mismatches <<
			" bsasa5: " << F(9,3,bsasa5) <<
			" bsasa14: " << F(9,3,bsasa14) <<
			" bsasa_ratio: " << F(9,3,bsasa_ratio) <<
			" n_unsat_donors: " << I(4,n_unsat_donors) <<
			" n_unsat_donors_bb: " << I(4,n_unsat_donors_bb) <<
			" n_unsat_acceptors: " << I(4,n_unsat_acceptors) <<
			" n_unsat_acceptors_bb: " << I(4,n_unsat_acceptors_bb) <<
			' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << endl;

		if ( do_clustering ) {
			Vectors coords;
			for ( Size i=1; i<= pose.total_residue(); ++i ) coords.push_back( pose.residue(i).xyz("CA"));
			all_coords.push_back( coords );
			all_filenames.push_back( files[fi] );
			map< string, Real > scores;
			scores["totE"] = final_score;
			scores["intE"] = interface_energy;
			scores["bsasa14"] = bsasa14;
			scores["bsasa_ratio"] = bsasa_ratio;
			scores["n_unsat_donors_bb"] = n_unsat_donors_bb;
			scores["n_unsat_acceptors_bb"] = n_unsat_acceptors_bb;
			scores["n_mismatches"] = n_mismatches;
			scores["n_trp"] = n_trp;
			all_scores.push_back( scores );
		}
	} // files loop

	if ( !do_clustering ) exit(0);

	ToroidCoordsDistanceMetric metric;
	vector1< Sizes > cluster_members;
	Real threshold;
	cout << "clustering " << all_coords.size() << " decoys" << endl;
	fflush( stdout );
	devel::blab::cluster::simple_cluster( all_coords, metric,
		min_cluster_size,
		min_top_cluster_size, try_top_cluster_size, max_top_cluster_size,
		min_threshold, max_threshold,
		threshold, cluster_members );

	cout << "threshold: " << F(9,3,threshold) << " top-cluster-size: " << cluster_members[1].size() << endl;

	string const prefix( output_tag() + "cluster_N"+string_of( all_filenames.size())+
		"_T"+string_of( int( 100*threshold)) );

	string const distfile( prefix+"_distances.txt" );
	ofstream distout( distfile.c_str() );

	Sizes cluster_centers;
	for ( Size i=1; i<= cluster_members.size(); ++i ) {
		Size const center( cluster_members[i][1] );
		cluster_centers.push_back( center );

		/// write out distances to other cluster centers
		Vectors const & center_cooords( all_coords[ center ] );
		for ( Size j=1; j<= i; ++j ) {
			distout << ' ' << F(9,3,metric( center_cooords, all_coords[ cluster_centers[j] ] ) );
		}
		distout << '\n';

		/// log msg to stdout
		cout << "cluster " << I(4,i) << ' ' << I(4,cluster_members[i].size() ) << ' ' <<
			all_filenames[ center ] << endl;

		/// make a superposition pdb file
		string const cluster_file( prefix + "_"+lead_zero_string_of( i, 3 )+
			"_"+lead_zero_string_of( cluster_members[i].size(),3)+ ".pdb" );
		string cmd("models.py "+all_filenames[ center ] );
		Sizes members_shuffled( cluster_members[i] );
		members_shuffled.erase( members_shuffled.begin() ); // remove center
		random_permutation( members_shuffled, numeric::random::rg() );
		for ( Size j=1; j<= members_shuffled.size() && j < max_decoys_per_cluster_pdbfile; ++j ) {
			cmd += " "+all_filenames[ members_shuffled[j ] ];
		}
		cmd += " > "+cluster_file;
		run_command( cmd );

	}

	distout.close();


	strings const scorekeys( get_keys( all_scores[1] ) );
	for ( strings::const_iterator sk = scorekeys.begin(); sk!= scorekeys.end(); ++sk ) {

		/// write out energies, etc
		string const scoresfile( prefix+"_"+(*sk)+".txt" );
		ofstream out( scoresfile.c_str() );
		for ( Size i=1; i<= cluster_members.size(); ++i ) {
			for ( Size j=1; j<= cluster_members[i].size(); ++j ) {
				out << ' ' << F(9,3,all_scores[ cluster_members[i][j] ][ *sk ]);
			}
			out << '\n';
		}
		out.close();

		string outfile_prefix(prefix+"_tree_"+(*sk) );
		strings const suffixlist( make_vector1( string(".ps"), string(".png")));
		string const extra_args( (*sk).find( "unsat" ) == string::npos ? "  " : " -color_by_average_score " );
		for ( strings::const_iterator suf = suffixlist.begin(); suf != suffixlist.end(); ++suf ) {
			run_command( "python /home/pbradley/python/make_color_trees_simple.py "+distfile+
				" -o "+outfile_prefix+(*suf)+ extra_args +
				" -scores_file "+scoresfile+" 0 " ); // scores start at col# 0
		}
	}// types of score-colored trees
}


///////////////////////////////////////////////////////////////////////////////
Real
min_helix_rmsd( Vectors const & coords1_in, Vectors const & coords2_in )
{
	Vectors const & coords1( coords1_in.size() >= coords2_in.size() ? coords1_in : coords2_in );
	Vectors const & coords2( coords1_in.size() >= coords2_in.size() ? coords2_in : coords1_in );
	runtime_assert( coords1.size() >= coords2.size() );
	Real min_rmsd(1e6 );
	for ( Size shift=0; shift<= coords1.size()-coords2.size(); ++shift ) {
		Real rmsd(0);
		for ( Size i=1; i<= coords2.size(); ++i ) {
			rmsd += ( coords1[ shift+i].distance_squared( coords2[i] ) );
		}
		min_rmsd = min( min_rmsd, rmsd );
	}
	return min_rmsd/coords2.size(); // still squared, not divided by #coords

}

///////////////////////////////////////////////////////////////////////////////
class TwoComponentCoordsDistanceMetric {
public:
	Real
	operator()( vector1< Vectors > const & coords1,
		vector1< Vectors > const & coords2 ) const
	{
		runtime_assert( coords1.size() == 6 );
		runtime_assert( coords2.size() == 6 );

		// find best matches
		vector1< Reals > helix_pair_rmsds( 5 );
		for ( Size i=1; i<= 5; ++i ) helix_pair_rmsds[i].resize( 5, 1e6 );

		for ( Size i=1; i<= 5; ++i ) {
			if ( i == 3 ) continue;
			for ( Size j=1; j<= 5; ++j ) {
				if ( j == 3 ) continue;
				// this could be made more efficient
				helix_pair_rmsds[i][j] = 0.5 * ( min_helix_rmsd( coords1[i  ], coords2[j  ] ) +
					min_helix_rmsd( coords1[i+1], coords2[j+1] ) );
			}
		}
		return sqrt( 0.5 * ( min( min( helix_pair_rmsds[1] ),  min( helix_pair_rmsds[2] ) ) +
			min( min( helix_pair_rmsds[4] ),  min( helix_pair_rmsds[5] ) ) ) );
	}

};


///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
void
cluster_2c_test()
{
	//ScoreFunctionOP fa_scorefxn( get_score_function() );

	Real min_threshold, max_threshold;
	Size min_cluster_size, min_top_cluster_size, try_top_cluster_size, max_top_cluster_size,
		max_clusters, max_decoys_per_cluster_pdbfile;

	parse_clustering_params( min_threshold, max_threshold, min_cluster_size, min_top_cluster_size, try_top_cluster_size,
		max_top_cluster_size, max_clusters, max_decoys_per_cluster_pdbfile );

	//string const repeatseg("ddaaaaaaaaaaaaaabbbcccccccccccccccd");
	//string const repeatseg("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGVSLEQALKILKVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAELGTTVEEAVKRALKLKTKLGVSLEQALKILEVAAKLGTTVEEAVKRALKLKTKLGGW");

	//Size const nrepeat( 6 ), repeatlen( 35 ), nres_monomer( nrepeat*repeatlen );
	Size const nrepeat( 6 );

	//strings const files( start_files() ); // old way
	string const logfile( start_file() );

	bool const planar( option[ my_options::planar ] );

	// Pose template_pose;
	// pose_from_pdb( template_pose, option[ my_options::template_pdb ]() );

	vector1< vector1< Vectors > > all_coords;
	//PoseOPs all_poses;
	strings all_filenames, all_scorelines;
	// Reals all_int_energies, all_tot_energies, all_bsasa14s, all_bsasa_ratios;
	vector1< map< string, Real > > all_scores;
	bools all_flips;

	//bool const calc_sasa( true );
	TwoComponentCoordsDistanceMetric metric;

	ifstream data( logfile.c_str() );
	string line;

	strings const good_tags( make_vector1( string("interface_quality"),
		string("bsasa14"),
		string("bsasa_ratio"),
		string("n_unsat_donors_bb"),
		string("n_unsat_acceptors_bb"),
		string("outer_helix_twist"),
		string("min_refold_rmsd") ) );

	while ( getline( data, line ) ) {
		if ( line.find("final_score") == string::npos ) continue;
		strings const l( split_to_vector1( line ) );
		runtime_assert( l[1].substr( l[1].size()-17 ) == ".out:final_scores" );
		string const pdbdir( l[1].substr( 0, l[1].size() - 17 )+".pdbs/" );
		string const filename( pdbdir + l[3] );

		if ( !utility::file::file_exists( filename ) ) {
			cerr << "missing: " << filename << endl;
			continue;
		}

		map< string, Real > scores;
		for ( Size i=1; i< l.size(); ++i ) {
			if ( l[i][ l[i].size()-1 ]!= ':' ) continue;
			string const tag( l[i].substr(0,l[i].size()-1) );
			if ( has_element(good_tags,tag) && is_float( l[i+1] ) ) {
				scores[ tag ] = float_of( l[i+1] );
			}
		}

		vector1< Vectors > coords( read_CA_coords_from_complex_file( filename ));

		/// should we flip coords?
		bool flipped( false );
		if ( planar && coords[1][1].z() <0 ) {
			// rotate about the y-axis
			for ( Size r=1; r<= 2; ++r ) {
				for ( Size i=1; i<= coords[r].size(); ++i ) {
					coords[r][i].x( -1 * coords[r][i].x() );
					coords[r][i].z( -1 * coords[r][i].z() );
				}
			}
			flipped = true;
		}

		runtime_assert( coords[1].size()%nrepeat == 0 );
		Size const repeatlen1( coords[1].size()/nrepeat );

		runtime_assert( coords[2].size()%nrepeat == 0 );
		Size const repeatlen2( coords[2].size()/nrepeat );

		TR.Trace << "read_coords " << repeatlen1 << ' ' << repeatlen2 << ' ' << filename << endl;

		// find the interface helices...
		vector1< Vectors > repeat_coords; // we dont know where the helices are!

		// count inter-repeat contacts
		vector1< Sizes > contact_counts(2);
		contact_counts[1].resize( nrepeat, 0 );
		contact_counts[2].resize( nrepeat, 0 );
		Real const dis2_threshold( 9 * 9 ); // C-alpha dis2 threshold
		for ( Size i=1; i<= coords[1].size(); ++i ) {
			for ( Size j=1; j<= coords[2].size(); ++j ) {
				if ( coords[1][i].distance_squared( coords[2][j] ) < dis2_threshold ) {
					++contact_counts[1][ (i-1)/repeatlen1+1 ];
					++contact_counts[2][ (j-1)/repeatlen2+1 ];
				}
			}
		}

		// now add the contacting repeats to the list
		for ( Size r=1; r<= 2; ++r ) {
			Size const repeatlen( r==1 ? repeatlen1 : repeatlen2  );
			Size const rep( arg_max( contact_counts[r] ) );
			TR.Trace << "core_repeat " << r << ' '<< rep << I(4,contact_counts[r][rep] ) << ' ' << filename << endl;
			Sizes contact_repeats;
			if ( rep == 1 ) contact_repeats.push_back( nrepeat );
			else contact_repeats.push_back( rep-1 );
			contact_repeats.push_back( rep );
			if ( rep == nrepeat ) contact_repeats.push_back( 1 );
			else contact_repeats.push_back( rep+1 );
			for ( Sizes::const_iterator rr= contact_repeats.begin(); rr != contact_repeats.end(); ++rr ) {
				Vectors rr_coords;
				for ( Size i=1; i<= repeatlen; ++i ) {
					rr_coords.push_back( coords[r][ ((*rr)-1)*repeatlen + i ] );
				}
				repeat_coords.push_back( rr_coords );
			}
		}

		runtime_assert( repeat_coords.size() == 6 );
		all_coords.push_back( repeat_coords );
		all_filenames.push_back( filename );
		all_flips.push_back( flipped );
		//map< string, Real > scores;
		//scores["lenAVG"] = coords[1].size() + coords[2].size();
		all_scores.push_back( scores );
		all_scorelines.push_back( line );
		// if ( all_coords.size()>1 ) {
		//  TR.Trace << "metric " << metric( all_coords.front(), all_coords.back() ) << endl;
		// }
	}


	//ToroidCoordsDistanceMetric metric;
	vector1< Sizes > cluster_members;
	Real threshold;
	cout << "clustering " << all_coords.size() << " decoys" << endl;
	fflush( stdout );
	devel::blab::cluster::simple_cluster( all_coords, metric,
		min_cluster_size,
		min_top_cluster_size, try_top_cluster_size, max_top_cluster_size,
		min_threshold, max_threshold,
		threshold, cluster_members );

	cout << "threshold: " << F(9,3,threshold) << " top-cluster-size: " << cluster_members[1].size() << endl;

	string const prefix( output_tag() + "cluster_N"+string_of( all_filenames.size())+
		"_T"+string_of( int( 100*threshold)) );

	string const distfile( prefix+"_distances.txt" );
	ofstream distout( distfile.c_str() );

	Sizes cluster_centers;
	for ( Size i=1; i<= cluster_members.size(); ++i ) {
		Size const center( cluster_members[i][1] );
		cluster_centers.push_back( center );

		/// write out distances to other cluster centers
		vector1< Vectors > const & center_cooords( all_coords[ center ] );
		for ( Size j=1; j<= i; ++j ) {
			distout << ' ' << F(9,3,metric( center_cooords, all_coords[ cluster_centers[j] ] ) );
		}
		distout << '\n';

		/// log msg to stdout
		cout << "cluster " << I(4,i) << ' ' << I(4,cluster_members[i].size() ) << ' ' <<
			all_filenames[ center ] << endl;

		/// make a superposition pdb file
		string const cluster_file( prefix + "_"+lead_zero_string_of( i, 3 )+
			"_"+lead_zero_string_of( cluster_members[i].size(),3)+ ".pdb" );
		string cmd("models.py ");//+all_filenames[ center ] );
		Sizes members_shuffled( cluster_members[i] );
		members_shuffled.erase( members_shuffled.begin() ); // remove center
		random_permutation( members_shuffled, numeric::random::rg() );
		members_shuffled.insert( members_shuffled.begin(), center );
		Size tmpcounter(0);
		for ( Size j=1; j<= members_shuffled.size() && j <= max_decoys_per_cluster_pdbfile; ++j ) {
			Size const member( members_shuffled[j] );
			string filename( all_filenames[ member ] );
			if ( all_flips[ member ] ) {
				// make a temporary file
				++tmpcounter;
				string const tmpfile("tmpc2c"+lead_zero_string_of(tmpcounter,4)+".pdb");
				ofstream out( tmpfile.c_str() );
				ifstream data( filename.c_str() );
				string line;
				while ( getline( data, line ) ) {
					if ( line.substr(0,4) == "ATOM" ) {
						Real x( -1 * float_of( line.substr(30,8) ) ), y ( float_of( line.substr(38,8))),
							z( -1 * float_of( line.substr(46,8) ) );
						out << line.substr(0,30) << F(8,3,x) << F(8,3,y) << F(8,3,z) << line.substr(54) << '\n';
					}
				}
				out.close();
				data.close();
				filename = tmpfile;
			}
			cmd += " "+filename;
		}
		cmd += " > "+cluster_file;
		run_command( cmd );
		ofstream out( cluster_file.c_str(), std::ios_base::app );
		for ( Size j=1; j<= cluster_members[i].size(); ++j ) {
			out << "REMARK SCORELINE " << I(4,i) << I(4,j) << ' ' << all_scorelines[ cluster_members[i][j] ] << '\n';
		}
		out.close();

	}

	distout.close();

	if ( !all_scores.empty() ) {
		strings const scorekeys( get_keys( all_scores[1] ) );
		for ( strings::const_iterator sk = scorekeys.begin(); sk!= scorekeys.end(); ++sk ) {

			/// write out energies, etc
			string const scoresfile( prefix+"_"+(*sk)+".txt" );
			ofstream out( scoresfile.c_str() );
			for ( Size i=1; i<= cluster_members.size(); ++i ) {
				for ( Size j=1; j<= cluster_members[i].size(); ++j ) {
					out << ' ' << F(9,3,all_scores[ cluster_members[i][j] ][ *sk ]);
				}
				out << '\n';
			}
			out.close();

			string outfile_prefix(prefix+"_tree_"+(*sk) );
			strings const suffixlist( make_vector1( string(".ps"), string(".png")));
			string const extra_args( (*sk).find( "interface_quality" ) != string::npos ? "  " : " -color_by_average_score " );
			for ( strings::const_iterator suf = suffixlist.begin(); suf != suffixlist.end(); ++suf ) {
				run_command( "python /home/pbradley/python/make_color_trees_simple.py "+distfile+
					" -o "+outfile_prefix+(*suf)+ extra_args +
					" -scores_file "+scoresfile+" 0 " ); // scores start at col# 0
			}
		}// types of score-colored trees
	}
}


///////////////////////////////////////////////////////////////////////////////
void
read_int_designs_test()
{
	Real const rmsd_threshold( 2.0 );
	Size const nrepeat( 6 ), repeatlen( 35 ), nres_monomer( nrepeat*repeatlen );

	strings const files( start_files() );

	vector1< Vectors > all_coords;
	PoseOPs all_poses;
	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		runtime_assert( pose.total_residue() == 2 * nrepeat * repeatlen );


		Vector center1, axis1, center2, axis2;
		{
			Size const rpos(5), repeatlen( 35 ); //, nrepeat( 6 ); // rpos is arbitrary
			Stub const stub1( pose.residue( rpos     ).xyz("CA"),
				pose.residue( rpos + 1 ).xyz("CA"),
				pose.residue( rpos + 2 ).xyz("CA") );
			Stub const stub2( pose.residue( rpos + repeatlen    ).xyz("CA"),
				pose.residue( rpos + repeatlen + 1 ).xyz("CA"),
				pose.residue( rpos + repeatlen + 2 ).xyz("CA") );
			Real theta;
			Vector t;
			get_stub_transform_data( stub1, stub2, center1, axis1, t, theta );
			TR.Trace << "theta: " << numeric::conversions::degrees( theta ) << endl;
		}
		{
			Size const rpos(nrepeat*repeatlen+5), repeatlen( 35 );//, nrepeat( 6 ); // rpos is arbitrary
			Stub const stub1( pose.residue( rpos     ).xyz("CA"),
				pose.residue( rpos + 1 ).xyz("CA"),
				pose.residue( rpos + 2 ).xyz("CA") );
			Stub const stub2( pose.residue( rpos + repeatlen    ).xyz("CA"),
				pose.residue( rpos + repeatlen + 1 ).xyz("CA"),
				pose.residue( rpos + repeatlen + 2 ).xyz("CA") );
			Real theta;
			Vector t;
			get_stub_transform_data( stub1, stub2, center2, axis2, t, theta );
			TR.Trace << "theta: " << numeric::conversions::degrees( theta ) << endl;
		}


		/// relationship between axes...
		Real const angle( numeric::conversions::degrees( std::acos( axis1.dot( axis2 ) ) ) );
		Real const v12( axis1.dot( axis2 ) ), v11( axis1.dot(axis1)), v22( axis2.dot(axis2)),
			c11( center1.dot( axis1 ) ), c12( center1.dot( axis2) ), c21( center2.dot( axis1 ) ), c22( center2.dot( axis2) ),
			a( c21 - c11 ), b( c22 - c12 ),
			l1( ( a * v22 - b * v12 )/( v11 * v22 - v12 * v12 ) ),
			l2( ( a * v12 - b * v11 )/( -1*( v12 * v12 - v22 * v11 ) ) );
		// check work
		Vector const p1( center1 + l1 * axis1 ), p2( center2 + l2 * axis2 );
		cout << "test " << (p1-p2).dot( axis1) << ' ' << (p1-p2).dot( axis2 ) << endl;


		cout << "axes: angle: " << F(9,3,angle) << " distance: " << F(9,3,p1.distance(p2)) << ' ' << files[fi] << endl;


		Vectors coords;
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			coords.push_back( pose.residue(i).xyz("CA") );
		}
		all_coords.push_back( coords );
		all_poses.push_back( PoseOP( new Pose( pose ) ) );
	}

	for ( Size ii=1; ii<= all_coords.size(); ++ii ) {
		for ( Size jj=1; jj<= all_coords.size(); ++jj ) {
			if ( ii==jj ) continue;
			Size chain1_spin, chain2_spin;
			bool flip_chains;
			Real min_rmsd;
			compare_docked_toroid_coords( nrepeat, repeatlen, all_coords[ii], all_coords[jj],
				chain1_spin, chain2_spin, flip_chains,
				min_rmsd );

			cout << "compare: " << I(5,ii) << I(5,jj) << F(9,3,min_rmsd) <<
				I(2,chain1_spin) << I(2,chain2_spin) << I(2,flip_chains) <<
				' ' << files[ii] << ' ' << files[jj] << endl;

			if ( min_rmsd < rmsd_threshold ) {
				using namespace id;

				Pose const & pose1( *all_poses[ii] );
				Pose pose2( *all_poses[jj] ); // make a copy

				/// superimpose pose2 onto pose1
				AtomID_Map< AtomID > atom_map;

				initialize_atomid_map( atom_map, pose2, id::GLOBAL_BOGUS_ATOM_ID );

				for ( Size chain2=1; chain2<= 2; ++chain2 ) { // chain number in pose2
					Size const spin( chain2 == 1 ? chain1_spin : chain2_spin );
					Size const chain1( flip_chains ? chain2%2+1 : chain2 ); // chain number in pose1

					for ( Size rep1=1; rep1<= nrepeat; ++rep1 ) {
						Size const rep2( (rep1+spin-1)%nrepeat+1 );

						for ( Size pos=1; pos<= repeatlen; ++pos ) {
							Size const pos1( (chain1-1)*nres_monomer + (rep1-1)*repeatlen + pos );
							Size const pos2( (chain2-1)*nres_monomer + (rep2-1)*repeatlen + pos );

							atom_map[ AtomID( pose2.residue(pos2).atom_index("CA"), pos2 ) ] =
								AtomID( pose1.residue(pos1).atom_index("CA"),pos1);
						}
					}
				}
				Real const rmsd( rmsd_by_mapping( pose2, pose1, atom_map ) );
				TR.Trace << "recompute_rmsd: " << F(9,3,min_rmsd) << F(9,3,rmsd) << endl;

				superimpose_pose( pose2, pose1, atom_map );

				pose2.dump_pdb( "sup_"+lead_zero_string_of( ii, 2 )+"_"+lead_zero_string_of( jj, 2 )+".pdb" );
				// copy of pose1, unshifted
				pose1.dump_pdb( "sup_"+lead_zero_string_of( ii, 2 )+"_"+lead_zero_string_of( ii, 2 )+".pdb" );

			}
		}
	}


}



///////////////////////////////////////////////////////////////////////////////
void
interface_test()
{
	string const repeatseq("VSLEQALKILKVAAELGTTVEEAVKRALKLKTKLG");
	string const design_me("---++--+------+-----++--++-++-++++-");

	Real const trans_mag( 1.5 ), rot_mag( 15 );

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) ),
		cen_scorefxn( new ScoreFunction() );

	// # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
	// env     1.0
	// pair    1.0
	// cbeta   1.0
	// vdw     1.0
	// rg      3.0
	// cenpack 1.0
	// hs_pair 1.0
	// ss_pair 1.0
	// rsigma  1.0
	// sheet   1.0

	cen_scorefxn->set_weight( env, 1.0 );
	cen_scorefxn->set_weight( scoring::pair, 1.0 );
	cen_scorefxn->set_weight( cbeta, 1.0 );
	cen_scorefxn->set_weight( vdw, 1.0 );
	cen_scorefxn->set_weight( cenpack, 1.0 ); // not sure about this guy...

	/// start out, just do some self-docking, see what comes out...
	Pose monomer_pose;
	pose_from_pdb( monomer_pose, start_file() );

	Real const monomer_score( (*fa_scorefxn)( monomer_pose ) );
	TR.Trace << "monomer_score: " << F(9,3,monomer_score) << ' ' << start_file() <<
		' ' << monomer_pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << endl;

	Size const nres_monomer( monomer_pose.total_residue() );

	Pose pose( monomer_pose );

	pose.append_residue_by_jump( monomer_pose.residue(1), 1 );
	for ( Size i=2; i<= monomer_pose.total_residue(); ++i ) {
		pose.append_residue_by_bond( monomer_pose.residue(i) );
	}

	pose.conformation().insert_chain_ending( monomer_pose.total_residue() );

	bools is_designable( pose.total_residue(), false );
	{
		Size const repeatlen( design_me.size() );
		runtime_assert( pose.total_residue() % repeatlen == 0 );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			is_designable[i] = ( design_me[ (i-1)%repeatlen ] == '+' );
			if ( is_designable[i] ) {
				TR.Trace << "is_designable: " << i << ' ' << (i-1)%repeatlen+1 << endl;
			}
		}
	}

	string const simfile( shared_output_tag() +"dock.work" ), worktag("tmp");

	Pose const start_pose( pose );

	while ( true ) {

		Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
		if ( n > nstruct() ) break;

		pose = start_pose;

		// randomly mutate the designable positions to hydrophobic residues to encourage contacts
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( is_designable[i] ) {
				AA const new_aa( aa_from_oneletter_code( random_element( make_vector1( 'L', 'L', 'L',
					'A', 'A', 'A',
					'I', 'I',
					'G', 'G',
					'V', 'V', 'V', 'V' ) ) ) );
				make_sequence_change( i, new_aa, pose );
			}
		}


		{ // randomly re-orient the two poses
			using namespace kinematics;
			using namespace numeric;
			xyzMatrix_double R( protocols::geometry::random_reorientation_matrix( 360.0, 360.0 ) );

			/// first pre-rotate the second monomer by our random rotation
			pose.set_jump( 1, Jump( RT( R, Vector(0,0,0) ) ) );

			/// figure out where the COM of mon2 is
			Vector com1(0,0,0), com2(0,0,0);
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( i<= nres_monomer ) com1 += pose.residue(i).xyz("CA");
				else com2 += pose.residue(i).xyz("CA");
			}
			com1/= nres_monomer;
			com2/= nres_monomer;

			/// so we get a center of mass of "com" with a translation of 0,0,0
			Vector const com_offset_with_zero_translation( com2-com1 );


			/// now translate

			// slide into contact
			Real const contact_dis2_threshold( 7.0 * 7.0 );
			Real distance( 100.0 );
			Vector const translation_unit(  random_unit_vector() );
			while ( true ) {
				distance -= 1.;
				Vector const desired_com_offset( distance * translation_unit );
				Stub const upstub( pose.conformation().upstream_jump_stub(1) );

				pose.set_jump( 1, Jump( RT( R, upstub.global2local( desired_com_offset-com_offset_with_zero_translation+upstub.v))));

				{ // debugging
					Vector com1(0,0,0), com2(0,0,0);
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						if ( i<= nres_monomer ) com1 += pose.residue(i).xyz("CA");
						else com2 += pose.residue(i).xyz("CA");
					}
					com1/= nres_monomer;
					com2/= nres_monomer;

					/// so we get a center of mass of "com" with a translation of 0,0,0
					Vector const actual_com_offset( com2-com1 );
					Real const com_error( actual_com_offset.distance( desired_com_offset ) );
					TR.Trace << "com_error " << F(9,3,com_error ) << endl;
				}

				Real mindis2( 1e6 );
				for ( Size i=1; i<= nres_monomer; ++i ) {
					for ( Size j=nres_monomer+1; j<= pose.total_residue(); ++j ) {
						Real const dis2( pose.residue(i).nbr_atom_xyz().distance_squared( pose.residue(j).nbr_atom_xyz() ) );
						mindis2 = std::min( dis2, mindis2 );
					}
				}
				TR.Trace << "mindis: " << F(9,3,distance) << F(9,3,sqrt( mindis2 ) ) << endl;
				if ( mindis2 < contact_dis2_threshold ) break;
			}

			//pose.dump_pdb("test_"+string_of(n)+".pdb");
		}


		{ /// centroid perturbations
			Size const n_outer( 5 ), n_inner( 500 );

			devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

			(*cen_scorefxn)( pose );
			set_jump_rb_centers( pose );


			Real const mc_hitemp( 4 ), mc_lotemp( 0.5 );
			Real const gamma = std::pow( (mc_lotemp/mc_hitemp), 1.0/Real(n_inner*n_outer));
			Real mc_temp( mc_hitemp / gamma );

			protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *cen_scorefxn, mc_temp ) );

			for ( Size n=1; n<= n_outer; ++n ) {
				for ( Size m=1; m<= n_inner; ++m ) {
					mc_temp *= gamma;
					mc->set_temperature( mc_temp );

					protocols::rigid::gaussian_jump_move( pose, 1, trans_mag, rot_mag, 0 ); // dir=0 --> choose randomly
					(*cen_scorefxn)( pose );
					bool const mc_accept( mc->boltzmann( pose ) );
					TR.Trace << "mc_accept " << I(4,n) << I(4,m) << I(3,mc_accept) << F(9,3,mc->last_accepted_score()) <<
						F(9,3,mc->lowest_score() ) << F(9,3,mc->temperature() ) << endl;
					set_jump_rb_centers( pose );
				}
				mc->recover_low( pose );
			}



			devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );

			// recover the sidechains
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				pose.replace_residue( i, start_pose.residue(i), true ); // orient_backbone
			}
		}


		if ( true ) { // design
			Size const design_cycles( option[ my_options::design_cycles ] );

			for ( Size m=1; m<= design_cycles+1; ++m ) {
				bool const skip_relax( m == 1 );

				if ( !skip_relax ) { // relax
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					setup_interface_movemap( pose, *movemap );
					fastrelax.set_movemap( movemap );
					if ( !dry_run() ) fastrelax.apply( pose );
				}

				if ( m>design_cycles ) break; // relax is cheap, so do one extra one at the end...

				{ // design
					pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
					task->initialize_from_command_line();
					task->or_include_current( true );
					// if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					setup_interface_movemap( pose, *movemap );

					bools is_flexible( pose.total_residue(), false );
					for ( Size i=1; i<= pose.total_residue(); ++i ) is_flexible[i] = ( movemap->get_chi(i) );
					task->restrict_to_residues( is_flexible );
					Size n_repack(0), n_design(0);
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						if ( is_flexible[i] ) {
							if ( is_designable[i] ) {
								++n_design;
							} else {
								++n_repack;
								task->nonconst_residue_task(i).restrict_to_repacking();
							}
						}
					}
					TR.Trace << "n_design: " << n_design << " n_repack: " << n_repack << endl;
					//if ( nodesign ) task->restrict_to_repacking();

					if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
						protocols::task_operations::LimitAromaChi2Operation lp_op;
						lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negatives)
						lp_op.apply( pose, *task );
					}
					Size const nloop( 1 ); // in case linmem ig 25 );
					ScoreFunctionOP design_scorefxn(0);
					// if ( ( option[ my_options::use_softrep_for_early_design ] && m<= design_cycles/2 ) ||
					//    option[ my_options::use_softrep_for_design ] ) {
					design_scorefxn = ScoreFunctionFactory::create_score_function( "soft_rep_design.wts" );
					adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
					// } else {
					//  design_scorefxn = fa_scorefxn; // already adjusted refwts
					// }
					protocols::minimization_packing::PackRotamersMover packmover( design_scorefxn, task, nloop );
					// protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
					if ( !dry_run() ) packmover.apply( pose );
				}

			} // cycles



		} else {
			/// fullatom "refinement" -- just chi-angles/rots at interface and rb-offset
			// bools is_partner1( pose.total_residue(), false ), is_partner2( pose.total_residue(), false );
			// for ( Size i=1; i<= pose.total_residue(); ++i ) {
			//  is_partner1[i] = ( pose.chain(i) == 1 );
			//  is_partner2[i] = ( pose.chain(i) == 2 );
			// }
			// bools nbrs1, nbrs2;
			// find_neighbors( is_partner1, pose, nbrs1 );
			// find_neighbors( is_partner2, pose, nbrs2 );


			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP mm( new MoveMap );
			setup_interface_movemap( pose, *mm );
			// mm->set_bb( false );
			// mm->set_chi( false );
			// Size int1_count(0), int2_count(0);
			// for ( Size i=1; i<= pose.total_residue(); ++i ) {
			//  if ( ( is_partner1[i] && nbrs2[i] ) || ( is_partner2[i] && nbrs1[i] ) ) {
			//   TR.Trace << "flexpos: " << I(4,i) << endl;
			//   mm->set_chi(i, true );
			//   int1_count += is_partner1[i];
			//   int2_count += is_partner2[i];
			//  }
			// }
			// TR.Trace << "numflex: " << int1_count << ' ' << int2_count << endl;
			// mm->set_jump(1,true);

			fastrelax.set_movemap( mm );
			fastrelax.apply( pose );
		}


		Real const final_score( (*fa_scorefxn )( pose ) );

		// compute interface energies
		EnergyMap interface_emap;
		{
			EnergyGraph const & energy_graph( pose.energies().energy_graph() );
			for ( Size pos1 = 1; pos1 <= nres_monomer; ++pos1 ) {
				runtime_assert( pose.chain(pos1) == 1 );
				for ( utility::graph::Graph::EdgeListConstIter
						ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
						ire = energy_graph.get_node( pos1 )->const_edge_list_end();
						ir != ire; ++ir ) {
					EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
					Size const pos2( edge->get_other_ind( pos1 ) );
					if ( pose.chain(pos2) == 2 ) {
						edge->add_to_energy_map( interface_emap );
					}
				}
			}
		}
		Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );


		bool const passed_score_filter
			( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
			score_filter_pass_early, simfile ) );
		string const outfilename( output_tag() + "_dock_N"+ lead_zero_string_of( n,4 )+".pdb");

		ostringstream out;
		out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
			//" simtime: " << F(9,3,simtime) <<
			" interface_energy: " << F(9,3,interface_energy) <<
			" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
			" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

		string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

		basic::prof_show_oneliner();

		if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
			append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
		}
		check_simtime();

	}

	signal_that_job_is_done();

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
redock_test()
{

	Real const trans_mag( 1.5 ), rot_mag( 15 );

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( core::scoring::get_score_function() ),
		cen_scorefxn( new ScoreFunction() );

	// # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
	// env     1.0
	// pair    1.0
	// cbeta   1.0
	// vdw     1.0
	// rg      3.0
	// cenpack 1.0
	// hs_pair 1.0
	// ss_pair 1.0
	// rsigma  1.0
	// sheet   1.0

	cen_scorefxn->set_weight( env, 1.0 );
	cen_scorefxn->set_weight( scoring::pair, 1.0 );
	cen_scorefxn->set_weight( cbeta, 1.0 );
	cen_scorefxn->set_weight( vdw, 1.0 );
	cen_scorefxn->set_weight( cenpack, 1.0 ); // not sure about this guy...



	strings files( start_files() );
	numeric::random::random_permutation( files, numeric::random::rg() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );
		string const simfile( shared_output_tag()+"_redock_"+filebase(filename)+".work" );

		if ( simfile_is_done( simfile ) ) continue;


		Pose pose;
		pose_from_pdb( pose, filename );

		runtime_assert( num_chains( pose ) == 2 );
		runtime_assert( pose.total_residue() == 2 * chain_end( 1, pose ) );
		Size const nres_monomer( chain_end( 1, pose ) );

		{
			FoldTree f( pose.total_residue() );
			f.new_jump( 5, nres_monomer+5, nres_monomer );
			f.reorder(5);
			pose.fold_tree( f );
		}

		Pose const start_pose( pose );

		while ( true ) {
			string const worktag("tmp");
			Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( n > nstruct() ) break;

			pose = start_pose;


			{ // randomly re-orient the two poses
				using namespace kinematics;
				using namespace numeric;
				xyzMatrix_double R( protocols::geometry::random_reorientation_matrix( 360.0, 360.0 ) );

				/// first pre-rotate the second monomer by our random rotation
				pose.set_jump( 1, Jump( RT( R, Vector(0,0,0) ) ) );

				/// figure out where the COM of mon2 is
				Vector com1(0,0,0), com2(0,0,0);
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					if ( i<= nres_monomer ) com1 += pose.residue(i).xyz("CA");
					else com2 += pose.residue(i).xyz("CA");
				}
				com1/= nres_monomer;
				com2/= nres_monomer;

				/// so we get a center of mass of "com" with a translation of 0,0,0
				Vector const com_offset_with_zero_translation( com2-com1 );


				/// now translate

				// slide into contact
				Real const contact_dis2_threshold( 7.0 * 7.0 );
				Real distance( 100.0 );
				Vector const translation_unit(  random_unit_vector() );
				while ( true ) {
					distance -= 1.;
					Vector const desired_com_offset( distance * translation_unit );
					Stub const upstub( pose.conformation().upstream_jump_stub(1) );

					pose.set_jump( 1, Jump( RT( R, upstub.global2local( desired_com_offset-com_offset_with_zero_translation+upstub.v))));

					{ // debugging
						Vector com1(0,0,0), com2(0,0,0);
						for ( Size i=1; i<= pose.total_residue(); ++i ) {
							if ( i<= nres_monomer ) com1 += pose.residue(i).xyz("CA");
							else com2 += pose.residue(i).xyz("CA");
						}
						com1/= nres_monomer;
						com2/= nres_monomer;

						/// so we get a center of mass of "com" with a translation of 0,0,0
						Vector const actual_com_offset( com2-com1 );
						Real const com_error( actual_com_offset.distance( desired_com_offset ) );
						TR.Trace << "com_error " << F(9,3,com_error ) << endl;
					}

					Real mindis2( 1e6 );
					for ( Size i=1; i<= nres_monomer; ++i ) {
						for ( Size j=nres_monomer+1; j<= pose.total_residue(); ++j ) {
							Real const dis2( pose.residue(i).nbr_atom_xyz().distance_squared( pose.residue(j).nbr_atom_xyz() ) );
							mindis2 = std::min( dis2, mindis2 );
						}
					}
					TR.Trace << "mindis: " << F(9,3,distance) << F(9,3,sqrt( mindis2 ) ) << endl;
					if ( mindis2 < contact_dis2_threshold ) break;
				}

				//pose.dump_pdb("test_"+string_of(n)+".pdb");
			}


			{ /// centroid perturbations
				Size const n_outer( 5 ), n_inner( 500 );

				devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

				(*cen_scorefxn)( pose );
				set_jump_rb_centers( pose );


				Real const mc_hitemp( 4 ), mc_lotemp( 0.5 );
				Real const gamma = std::pow( (mc_lotemp/mc_hitemp), 1.0/Real(n_inner*n_outer));
				Real mc_temp( mc_hitemp / gamma );

				protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *cen_scorefxn, mc_temp ) );

				for ( Size n=1; n<= n_outer; ++n ) {
					for ( Size m=1; m<= n_inner; ++m ) {
						mc_temp *= gamma;
						mc->set_temperature( mc_temp );

						protocols::rigid::gaussian_jump_move( pose, 1, trans_mag, rot_mag, 0 ); // dir=0 --> choose randomly
						(*cen_scorefxn)( pose );
						bool const mc_accept( mc->boltzmann( pose ) );
						TR.Trace << "mc_accept " << I(4,n) << I(4,m) << I(3,mc_accept) << F(9,3,mc->last_accepted_score()) <<
							F(9,3,mc->lowest_score() ) << F(9,3,mc->temperature() ) << endl;
						set_jump_rb_centers( pose );
					}
					mc->recover_low( pose );
				}



				devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );

				// recover the sidechains
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					pose.replace_residue( i, start_pose.residue(i), true ); // orient_backbone
				}
			}


			{
				/// fullatom "refinement" -- just chi-angles/rots at interface and rb-offset

				protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

				MoveMapOP mm( new MoveMap );
				setup_interface_movemap( pose, *mm );

				fastrelax.set_movemap( mm );
				fastrelax.apply( pose );
			}


			Real const final_score( (*fa_scorefxn )( pose ) );

			// compute interface energies
			EnergyMap interface_emap;
			{
				EnergyGraph const & energy_graph( pose.energies().energy_graph() );
				for ( Size pos1 = 1; pos1 <= nres_monomer; ++pos1 ) {
					runtime_assert( pose.chain(pos1) == 1 );
					for ( utility::graph::Graph::EdgeListConstIter
							ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
							ire = energy_graph.get_node( pos1 )->const_edge_list_end();
							ir != ire; ++ir ) {
						EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
						Size const pos2( edge->get_other_ind( pos1 ) );
						if ( pose.chain(pos2) == 2 ) {
							edge->add_to_energy_map( interface_emap );
						}
					}
				}
			}
			Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );

			// figure out the rmsd between the starting model and the final model
			// probably also worth calculating interface rmsd as well
			Real const rmsd( scoring::nbr_atom_rmsd( pose, start_pose ) );


			bool const passed_score_filter
				( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );
			string const outfilename( output_tag() + "_dock_N"+ lead_zero_string_of( n,4 )+".pdb");

			ostringstream out;
			out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
				//" simtime: " << F(9,3,simtime) <<
				" interface_energy: " << F(9,3,interface_energy) <<
				" rmsd: " << F(9,3,rmsd) <<
				" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
				" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

			string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

			basic::prof_show_oneliner();

			if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
				append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
			}
			check_simtime();

		} // nstruct

		signal_that_simfile_is_done( simfile );
	}
	signal_that_job_is_done();

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
refold_and_dock_test()
{

	Real const trans_mag( 1.5 ), rot_mag( 15 );

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( core::scoring::get_score_function() ),
		cen_scorefxn( new ScoreFunction() );

	// # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
	// env     1.0
	// pair    1.0
	// cbeta   1.0
	// vdw     1.0
	// rg      3.0
	// cenpack 1.0
	// hs_pair 1.0
	// ss_pair 1.0
	// rsigma  1.0
	// sheet   1.0

	cen_scorefxn->set_weight( env, 1.0 );
	cen_scorefxn->set_weight( scoring::pair, 1.0 );
	cen_scorefxn->set_weight( cbeta, 1.0 );
	cen_scorefxn->set_weight( vdw, 1.0 );
	cen_scorefxn->set_weight( cenpack, 1.0 ); // not sure about this guy...



	strings files( start_files() );
	numeric::random::random_permutation( files, numeric::random::rg() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );
		string const simfile( shared_output_tag()+"_redock_"+filebase(filename)+".work" );

		if ( simfile_is_done( simfile ) ) continue;


		Pose pdb_pose;
		pose_from_pdb( pdb_pose, filename );

		runtime_assert( num_chains( pdb_pose ) == 2 );

		Size const nres_monomer1( chain_end( 1, pdb_pose ) );
		Size const nres_monomer2( pdb_pose.total_residue() - nres_monomer1 );

		Pose pdb_pose1( pdb_pose ), pdb_pose2( pdb_pose );
		pdb_pose1.conformation().delete_residue_range_slow( nres_monomer1+1, pdb_pose.total_residue() );
		pdb_pose2.conformation().delete_residue_range_slow( 1, nres_monomer1 );
		string const seq1( pdb_pose.sequence().substr(0,nres_monomer1) );
		string const seq2( pdb_pose.sequence().substr(nres_monomer1) );
		runtime_assert( seq1.size() == nres_monomer1 );
		runtime_assert( seq2.size() == nres_monomer2 );
		runtime_assert( seq1 == pdb_pose1.sequence() );
		runtime_assert( seq2 == pdb_pose2.sequence() );
		// {
		//  FoldTree f( pose.total_residue() );
		//  f.new_jump( 5, nres_monomer+5, nres_monomer );
		//  f.reorder(5);
		//  pose.fold_tree( f );
		// }

		//Pose const start_pose( pose );

		while ( true ) {
			string const worktag("tmp");
			Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( n > nstruct() ) break;

			// refold pose1 and pose2
			Pose pose1, pose2;
			if ( dry_run() ) pose1 = pdb_pose1;
			else abrelax_sequence( seq1, *fa_scorefxn, pose1 );
			if ( dry_run() ) pose2 = pdb_pose2;
			else abrelax_sequence( seq2, *fa_scorefxn, pose2 );

			Pose pose(  pose1 );
			pose.append_residue_by_jump( pose2.residue(1), 1 );
			pose.conformation().insert_chain_ending( pose.total_residue()-1 );
			for ( Size i=2; i<= pose2.total_residue(); ++i ) pose.append_residue_by_bond( pose2.residue(i) );

			runtime_assert( pose.total_residue() == pdb_pose.total_residue() );
			runtime_assert( pose.sequence() == pdb_pose.sequence() );

			{ // randomly re-orient the two poses
				using namespace kinematics;
				using namespace numeric;
				xyzMatrix_double R( protocols::geometry::random_reorientation_matrix( 360.0, 360.0 ) );

				/// first pre-rotate the second monomer by our random rotation
				pose.set_jump( 1, Jump( RT( R, Vector(0,0,0) ) ) );

				/// figure out where the COM of mon2 is
				Vector com1(0,0,0), com2(0,0,0);
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					if ( i<= nres_monomer1 ) com1 += pose.residue(i).xyz("CA");
					else com2 += pose.residue(i).xyz("CA");
				}
				com1/= nres_monomer1;
				com2/= nres_monomer2;

				/// so we get a center of mass of "com" with a translation of 0,0,0
				Vector const com_offset_with_zero_translation( com2-com1 );


				/// now translate

				// slide into contact
				Real const contact_dis2_threshold( 7.0 * 7.0 );
				Real distance( 100.0 );
				Vector const translation_unit(  random_unit_vector() );
				while ( true ) {
					distance -= 1.;
					Vector const desired_com_offset( distance * translation_unit );
					Stub const upstub( pose.conformation().upstream_jump_stub(1) );

					pose.set_jump( 1, Jump( RT( R, upstub.global2local( desired_com_offset-com_offset_with_zero_translation+upstub.v))));

					{ // debugging
						Vector com1(0,0,0), com2(0,0,0);
						for ( Size i=1; i<= pose.total_residue(); ++i ) {
							if ( i<= nres_monomer1 ) com1 += pose.residue(i).xyz("CA");
							else com2 += pose.residue(i).xyz("CA");
						}
						com1/= nres_monomer1;
						com2/= nres_monomer2;

						/// so we get a center of mass of "com" with a translation of 0,0,0
						Vector const actual_com_offset( com2-com1 );
						Real const com_error( actual_com_offset.distance( desired_com_offset ) );
						TR.Trace << "com_error " << F(9,3,com_error ) << endl;
					}

					Real mindis2( 1e6 );
					for ( Size i=1; i<= nres_monomer1; ++i ) {
						for ( Size j=nres_monomer1+1; j<= pose.total_residue(); ++j ) {
							Real const dis2( pose.residue(i).nbr_atom_xyz().distance_squared( pose.residue(j).nbr_atom_xyz() ) );
							mindis2 = std::min( dis2, mindis2 );
						}
					}
					TR.Trace << "mindis: " << F(9,3,distance) << F(9,3,sqrt( mindis2 ) ) << endl;
					if ( mindis2 < contact_dis2_threshold ) break;
				}

				//pose.dump_pdb("test_"+string_of(n)+".pdb");
			}


			{ /// centroid perturbations
				Size const n_outer( 5 ), n_inner( dry_run() ? 10 : 500 );

				devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

				(*cen_scorefxn)( pose );
				set_jump_rb_centers( pose );


				Real const mc_hitemp( 4 ), mc_lotemp( 0.5 );
				Real const gamma = std::pow( (mc_lotemp/mc_hitemp), 1.0/Real(n_inner*n_outer));
				Real mc_temp( mc_hitemp / gamma );

				protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *cen_scorefxn, mc_temp ) );

				for ( Size n=1; n<= n_outer; ++n ) {
					for ( Size m=1; m<= n_inner; ++m ) {
						mc_temp *= gamma;
						mc->set_temperature( mc_temp );

						protocols::rigid::gaussian_jump_move( pose, 1, trans_mag, rot_mag, 0 ); // dir=0 --> choose randomly
						(*cen_scorefxn)( pose );
						bool const mc_accept( mc->boltzmann( pose ) );
						TR.Trace << "mc_accept " << I(4,n) << I(4,m) << I(3,mc_accept) << F(9,3,mc->last_accepted_score()) <<
							F(9,3,mc->lowest_score() ) << F(9,3,mc->temperature() ) << endl;
						set_jump_rb_centers( pose );
					}
					mc->recover_low( pose );
				}



				devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );

				// recover the sidechains
				for ( Size i=1; i<= pose1.total_residue(); ++i ) {
					pose.replace_residue( i, pose1.residue(i), true ); // orient_backbone
				}
				for ( Size i=1; i<= pose2.total_residue(); ++i ) {
					pose.replace_residue( i+nres_monomer1, pose2.residue(i), true ); // orient_backbone
				}
			}


			{
				/// fullatom "refinement" -- just chi-angles/rots at interface and rb-offset

				protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

				MoveMapOP mm( new MoveMap );
				setup_interface_movemap( pose, *mm );

				fastrelax.set_movemap( mm );
				if ( !dry_run() ) fastrelax.apply( pose );
			}


			Real const final_score( (*fa_scorefxn )( pose ) );

			// compute interface energies
			EnergyMap interface_emap;
			{
				EnergyGraph const & energy_graph( pose.energies().energy_graph() );
				for ( Size pos1 = 1; pos1 <= nres_monomer1; ++pos1 ) {
					runtime_assert( pose.chain(pos1) == 1 );
					for ( utility::graph::Graph::EdgeListConstIter
							ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
							ire = energy_graph.get_node( pos1 )->const_edge_list_end();
							ir != ire; ++ir ) {
						EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
						Size const pos2( edge->get_other_ind( pos1 ) );
						if ( pose.chain(pos2) == 2 ) {
							edge->add_to_energy_map( interface_emap );
						}
					}
				}
			}
			Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );

			// figure out the rmsd between the starting model and the final model
			// probably also worth calculating interface rmsd as well
			Real const rmsd( scoring::nbr_atom_rmsd( pose, pdb_pose ) );
			Real const rmsd1( scoring::nbr_atom_rmsd( pose1, pdb_pose1 ) );
			Real const rmsd2( scoring::nbr_atom_rmsd( pose2, pdb_pose2 ) );


			bool const passed_score_filter
				( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );
			string const outfilename( output_tag() + "_dock_N"+ lead_zero_string_of( n,4 )+".pdb");

			ostringstream out;
			out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
				//" simtime: " << F(9,3,simtime) <<
				" interface_energy: " << F(9,3,interface_energy) <<
				" rmsd: " << F(9,3,rmsd) <<
				" rmsd1: " << F(9,3,rmsd1) <<
				" rmsd2: " << F(9,3,rmsd2) <<
				" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
				" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

			string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

			basic::prof_show_oneliner();

			if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
				append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
			}
			check_simtime();

		} // nstruct

		signal_that_simfile_is_done( simfile );
	}
	signal_that_job_is_done();

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
template_redock_test()
{

	Real const trans_mag( 1.5 ), rot_mag( 15 );
	Size const repeatlen( 35 );

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP fa_scorefxn( core::scoring::get_score_function() ),
		cen_scorefxn( new ScoreFunction() );

	// # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
	// env     1.0
	// pair    1.0
	// cbeta   1.0
	// vdw     1.0
	// rg      3.0
	// cenpack 1.0
	// hs_pair 1.0
	// ss_pair 1.0
	// rsigma  1.0
	// sheet   1.0

	cen_scorefxn->set_weight( env, 1.0 );
	cen_scorefxn->set_weight( scoring::pair, 1.0 );
	cen_scorefxn->set_weight( cbeta, 1.0 );
	cen_scorefxn->set_weight( vdw, 1.0 );
	cen_scorefxn->set_weight( cenpack, 1.0 ); // not sure about this guy...



	strings template_files( option[ my_options::template_pdbs ]() );
	numeric::random::random_permutation( template_files, numeric::random::rg() );

	strings files( start_files() );
	numeric::random::random_permutation( files, numeric::random::rg() );


	for ( Size fi=1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );
		string const simfile( shared_output_tag()+"_redock_"+filebase(filename)+".work" );

		if ( simfile_is_done( simfile ) ) continue;

		Pose model_pose;
		pose_from_pdb( model_pose, filename );


		runtime_assert( num_chains( model_pose ) == 2 );
		runtime_assert( model_pose.total_residue() == 2 * chain_end( 1, model_pose ) );
		Size const nres_monomer( chain_end( 1, model_pose ) );
		runtime_assert( nres_monomer%repeatlen == 0 );

		for ( Size ti=1; ti<= template_files.size(); ++ti ) {
			string const template_filename( template_files[ti] );
			Pose template_pose;
			pose_from_pdb( template_pose, template_filename );

			if ( template_pose.total_residue() > nres_monomer ) {
				template_pose.conformation().delete_residue_range_slow( nres_monomer+1, template_pose.total_residue() );
				add_upper_terminus_type_to_pose_residue( template_pose, nres_monomer );
			}


			Sizes spins( make_vector1( 0,1,2,3,4,5 ) );
			numeric::random::random_permutation( spins, numeric::random::rg() );

			for ( Size ss=1; ss<= spins.size(); ++ss ) {

				Size const spin( spins[ss] );
				// mutate template_pose to match the sequence of the design

				Pose pose( template_pose );
				Size n_mutations(0);
				for ( Size i=1; i<= nres_monomer; ++i ) {
					Size const j( (i+spin*repeatlen-1)%nres_monomer+1 );
					if ( pose.residue(i).aa() != model_pose.residue(j).aa() ) {
						++n_mutations;
						make_sequence_change( i, model_pose.residue(j).aa(), pose );
					}
				}

				// pre-relax the template pose
				{
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

					MoveMapOP mm( new MoveMap );
					mm->set_chi( true );
					fastrelax.set_movemap( mm );
					fastrelax.apply( pose );
				}

				for ( Size i=1; i<= nres_monomer; ++i ) {
					if ( i==1 ) {
						pose.append_residue_by_jump( pose.residue(i), 1 );
					} else {
						pose.append_residue_by_bond( pose.residue(i) );
					}
				}

				pose.conformation().insert_chain_ending( nres_monomer );
				runtime_assert( pose.total_residue() == 2*nres_monomer );
				runtime_assert( num_chains( pose ) == 2 );

				{
					FoldTree f( pose.total_residue() );
					f.new_jump( 5, nres_monomer+5, nres_monomer );
					f.reorder(5);
					pose.fold_tree( f );
				}

				Pose const start_pose( pose );

				while ( true ) {
					string const worktag(filebase( template_filename)+"_S"+string_of(spin));
					Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
					if ( n > nstruct() ) break;

					pose = start_pose;


					{ // randomly re-orient the two poses
						using namespace kinematics;
						using namespace numeric;
						xyzMatrix_double R( protocols::geometry::random_reorientation_matrix( 360.0, 360.0 ) );

						/// first pre-rotate the second monomer by our random rotation
						pose.set_jump( 1, Jump( RT( R, Vector(0,0,0) ) ) );

						/// figure out where the COM of mon2 is
						Vector com1(0,0,0), com2(0,0,0);
						for ( Size i=1; i<= pose.total_residue(); ++i ) {
							if ( i<= nres_monomer ) com1 += pose.residue(i).xyz("CA");
							else com2 += pose.residue(i).xyz("CA");
						}
						com1/= nres_monomer;
						com2/= nres_monomer;

						/// so we get a center of mass of "com" with a translation of 0,0,0
						Vector const com_offset_with_zero_translation( com2-com1 );


						/// now translate

						// slide into contact
						Real const contact_dis2_threshold( 7.0 * 7.0 );
						Real distance( 100.0 );
						Vector const translation_unit(  random_unit_vector() );
						while ( true ) {
							distance -= 1.;
							Vector const desired_com_offset( distance * translation_unit );
							Stub const upstub( pose.conformation().upstream_jump_stub(1) );

							pose.set_jump( 1, Jump( RT( R, upstub.global2local( desired_com_offset -
								com_offset_with_zero_translation + upstub.v ))));

							{ // debugging
								Vector com1(0,0,0), com2(0,0,0);
								for ( Size i=1; i<= pose.total_residue(); ++i ) {
									if ( i<= nres_monomer ) com1 += pose.residue(i).xyz("CA");
									else com2 += pose.residue(i).xyz("CA");
								}
								com1/= nres_monomer;
								com2/= nres_monomer;

								/// so we get a center of mass of "com" with a translation of 0,0,0
								Vector const actual_com_offset( com2-com1 );
								Real const com_error( actual_com_offset.distance( desired_com_offset ) );
								TR.Trace << "com_error " << F(9,3,com_error ) << endl;
							}

							Real mindis2( 1e6 );
							for ( Size i=1; i<= nres_monomer; ++i ) {
								for ( Size j=nres_monomer+1; j<= pose.total_residue(); ++j ) {
									Real const dis2( pose.residue(i).nbr_atom_xyz().distance_squared( pose.residue(j).nbr_atom_xyz() ) );
									mindis2 = std::min( dis2, mindis2 );
								}
							}
							TR.Trace << "mindis: " << F(9,3,distance) << F(9,3,sqrt( mindis2 ) ) << endl;
							if ( mindis2 < contact_dis2_threshold ) break;
						}

						//pose.dump_pdb("test_"+string_of(n)+".pdb");
					}


					{ /// centroid perturbations
						Size const n_outer( 5 ), n_inner( 500 );

						devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID );

						(*cen_scorefxn)( pose );
						set_jump_rb_centers( pose );


						Real const mc_hitemp( 4 ), mc_lotemp( 0.5 );
						Real const gamma = std::pow( (mc_lotemp/mc_hitemp), 1.0/Real(n_inner*n_outer));
						Real mc_temp( mc_hitemp / gamma );

						protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *cen_scorefxn, mc_temp ) );

						for ( Size n=1; n<= n_outer; ++n ) {
							for ( Size m=1; m<= n_inner; ++m ) {
								mc_temp *= gamma;
								mc->set_temperature( mc_temp );

								protocols::rigid::gaussian_jump_move( pose, 1, trans_mag, rot_mag, 0 ); // dir=0 --> choose randomly
								(*cen_scorefxn)( pose );
								bool const mc_accept( mc->boltzmann( pose ) );
								TR.Trace << "mc_accept " << I(4,n) << I(4,m) << I(3,mc_accept) << F(9,3,mc->last_accepted_score()) <<
									F(9,3,mc->lowest_score() ) << F(9,3,mc->temperature() ) << endl;
								set_jump_rb_centers( pose );
							}
							mc->recover_low( pose );
						}



						devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );

						// recover the sidechains
						for ( Size i=1; i<= pose.total_residue(); ++i ) {
							pose.replace_residue( i, start_pose.residue(i), true ); // orient_backbone
						}
					}


					{
						/// fullatom "refinement" -- just chi-angles/rots at interface and rb-offset

						protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

						MoveMapOP mm( new MoveMap );
						setup_interface_movemap( pose, *mm );

						fastrelax.set_movemap( mm );
						fastrelax.apply( pose );
					}


					Real const final_score( (*fa_scorefxn )( pose ) );

					// compute interface energies
					EnergyMap interface_emap;
					{
						EnergyGraph const & energy_graph( pose.energies().energy_graph() );
						for ( Size pos1 = 1; pos1 <= nres_monomer; ++pos1 ) {
							runtime_assert( pose.chain(pos1) == 1 );
							for ( utility::graph::Graph::EdgeListConstIter
									ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
									ire = energy_graph.get_node( pos1 )->const_edge_list_end();
									ir != ire; ++ir ) {
								EnergyEdge const * edge( static_cast< EnergyEdge const *>( *ir ) );
								Size const pos2( edge->get_other_ind( pos1 ) );
								if ( pose.chain(pos2) == 2 ) {
									edge->add_to_energy_map( interface_emap );
								}
							}
						}
					}
					Real const interface_energy( fa_scorefxn->weights().dot( interface_emap ) );

					// figure out the rmsd between the starting model and the final model
					// probably also worth calculating interface rmsd as well
					//Real const rmsd( scoring::nbr_atom_rmsd( pose, start_pose ) );
					Real rmsd; //
					{
						Vectors new_coords, old_coords;
						for ( Size ch=0; ch< 2; ++ch ) {
							Size const chain_offset( nres_monomer*ch );
							for ( Size i=1; i<= nres_monomer; ++i ) {
								Size const newpos( chain_offset + i );
								Size const oldpos( chain_offset + (i+spin*repeatlen-1)%nres_monomer+1 );
								runtime_assert( pose.residue( newpos ).aa() == model_pose.residue( oldpos ).aa() );
								new_coords.push_back( pose.residue( newpos ).nbr_atom_xyz() );
								old_coords.push_back( model_pose.residue( oldpos ).nbr_atom_xyz() );
							}
						}
						rmsd = numeric::model_quality::calc_rms( new_coords, old_coords );
					}

					bool const passed_score_filter
						( append_score_to_scorefile_and_filter( worktag, final_score, score_filter_acceptance_rate,
						score_filter_pass_early, simfile ) );
					string const outfilename( output_tag() + "_redock_T"+worktag+"_N"+ lead_zero_string_of( n,4 )+".pdb");

					ostringstream out;
					out << "final_scores " << F(9,3,final_score) << ' ' << outfilename << ' ' <<
						//" simtime: " << F(9,3,simtime) <<
						" interface_energy: " << F(9,3,interface_energy) <<
						" rmsd: " << F(9,3,rmsd) <<
						" n_mutations: " << n_mutations <<
						" spin: " << spin <<
						" model_pdb: " << filebase( filename ) <<
						" template_pdb: " << filebase(template_filename) <<
						" interface_energies: " << interface_emap.weighted_string_of( fa_scorefxn->weights()) <<
						" total_energies: " << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights()) << '\n';

					string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

					basic::prof_show_oneliner();

					if ( passed_score_filter && option[ my_options::dump_hbonds ] && pdbfilename.size()>0 ) {
						append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
					}
					check_simtime();

				} // nstruct
			} // spins
		} // template files
		signal_that_simfile_is_done( simfile );
	}
	signal_that_job_is_done();

}


///////////////////////////////////////////////////////////////////////////////

void
aa_test()
{
	strings const files( start_files() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		Pose pose;
		pose_from_pdb( pose, files[fi] );

		set_ss_from_dssp( files[fi], pose );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( pose.secstruct(i) != 'H' ) continue;
			for ( Size j=i+1; j<= i+5 && pose.total_residue(); ++j ) {
				if ( pose.secstruct(j) != 'H' ) break;
				cout << "Hsep " << j-i << ' ' << pose.residue(i).name1() << ' ' << pose.residue(j).name1() << endl;
			}
		}
	}



}



///////////////////////////////////////////////////////////////////////////////

void
phenix_rebuild_again_test()
{

	strings files( start_files() );

	string const simfile( shared_output_tag()+"_rebuild.work" );

	if ( simfile_is_done( simfile ) ) {
		signal_that_job_is_done();
		return;
	}

	numeric::random::random_permutation( files, numeric::random::rg() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );
		string const worktag( filebase( filename ) );

		while ( true ) {

			Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );

			if ( n > nstruct() ) break;

			string const outfile( output_tag() + "_" + filebase( filename ) +"_N"+lead_zero_string_of(n,4)+".pdb" );
			run_command( "cp "+filename+" "+outfile );

			run_command( "python /home/pbradley/python/phenix/run_toroid_phaser_and_phenix_refinement.py "+\
				outfile+" 0.0 -topfiles 6");

			cout << "finished_rebuild " << outfile << endl;
			fflush( stdout );
		}
	}

	signal_that_job_is_done();

}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
semet_test()
{
	Size const nrepeat( 6 ), repeatlen( 35 );

	ScoreFunctionOP fa_scorefxn( get_score_function() );

	strings files( start_files() );

	for ( Size fi=1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );

		Pose pose;
		pose_from_pdb( pose, filename );

		runtime_assert( pose.total_residue() == nrepeat * repeatlen );

		Pose const start_pose( pose );
		for ( Size seqpos=1; seqpos<= pose.total_residue(); ++seqpos ) {
			if ( pose.residue(seqpos).aa() != aa_leu ) continue;

			pose = start_pose;

			Reals scores;
			scores.push_back( (*fa_scorefxn)( pose ) );

			for ( Size i=1; i<= 2; ++i ) {

				// try mutating to met
				if ( i==2 ) make_sequence_change( seqpos, aa_met, pose );

				// now repack just the sidechain
				pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
				task->set_bump_check( false );
				task->initialize_from_command_line();
				bools is_flexible( pose.total_residue(), false ); is_flexible[ seqpos ] = true;
				task->restrict_to_repacking();
				task->restrict_to_residues( is_flexible );
				//pack::pack_rotamers( pose, *fa_scorefxn, task );
				pack::RTMin rtmin;
				rtmin.rtmin( pose, *fa_scorefxn, task );
				scores.push_back( (*fa_scorefxn)( pose ) );
			}

			string const outfilename( "semet_" + filebase( files[fi] ) + "_" + lead_zero_string_of( seqpos, 3 )+".pdb" );

			Size const rpos( (seqpos-1)%repeatlen+1 ), repeat( (seqpos-1)/repeatlen+1 );

			cout << "final_score " << F(9,3,scores[3]-scores[2])<< ' ' << outfilename <<
				" rpos: " << I(4,rpos) <<
				" repeat: " << I(3,repeat) <<
				" start_score: " << F(9,3,scores[1]) <<
				" leu_score: " << F(9,3,scores[2]) <<
				" met_score: " << F(9,3,scores[3]) <<
				endl;
			pose.dump_pdb( outfilename );
		}
	}// files
}

////////////////////////////
void
count_similar_helix_helix_frags(
	Size const helix1_begin,
	Size const helix1_end,
	Vectors const & coords1,
	Size const helix2_begin,
	Size const helix2_end,
	Vectors const & coords2,
	HelixHelixFrags const & hh_frags,
	Size & nsim,
	bool & found_contact
)
{
	static HelixHelixFrag hh_frag;
	Real const max_hh_frag_dis( 3.0 );

	found_contact = get_helix_contact_coords( helix1_begin, helix1_end, coords1, helix2_begin, helix2_end, coords2,
		hh_frag.h1_coords, hh_frag.h2_coords );
	nsim = 0;
	if ( found_contact ) {
		// TR.Trace << "fcsims: " <<
		//  " h1: " << F(9,3,hh_frag.h1_coords[1].x()) <<F(9,3,hh_frag.h1_coords[1].y()) <<F(9,3,hh_frag.h1_coords[1].z()) <<
		//  " h2: " << F(9,3,hh_frag.h2_coords[1].x()) <<F(9,3,hh_frag.h2_coords[1].y()) <<F(9,3,hh_frag.h2_coords[1].z());

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			Real const dis( helix_helix_frag_distance( hh_frag, hh_frags[i] ) );
			if ( dis <= max_hh_frag_dis ) {
				//TR.Trace << ' ' << i << F(5,2,dis);
				//runtime_assert(  fabs( dis - helix_helix_frag_distance( hh_frags[i], hh_frag ) )<1e-3 ); // HACKING
				++nsim;
			}
		}
		//TR.Trace << endl;
	}
}
///////////////////////////////////////////////////////////////////////////////

void
helix_frag_test()
{
	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++nsim;
			}
			TR.Trace << "nsim " << nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	map< string, StubFrags > all_turn_frags;
	map< Size, StubFrags > all_helix_frags;

	Sizes const helixlens( option[ my_options::helix_lens ]() );
	strings const turns( option[ my_options::turns ] );
	//strings const turns( make_vector1( string("GBB") ) );
	Size const max_frags( max( Size(1000), nstruct() ) ); // per helixlen

	strings turn1s, turn2s;
	if ( option[ my_options::turn1s ].user() ) turn1s = option[ my_options::turn1s ]();
	if ( option[ my_options::turn2s ].user() ) turn2s = option[ my_options::turn2s ]();


	build_helix_library( helixlens, max_frags, all_helix_frags );
	build_turn_library( turns, all_turn_frags );


	Size const nrepeat( option[ my_options::nrepeat ] );
	Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	Real const max_shift( 2 );

	Size const max_clash( 1 ), min_close_contacts( 5 ), min_long_contacts( 10 ), min_nsims( 10 );

	// for ( Sizes::const_iterator h1len= helixlens.begin(); h1len != helixlens.end(); ++h1len ) {
	//  StubFrags const & helix1_frags( all_helix_frags.find( *h1len )->second );
	//  for ( strings::const_iterator turn1= turns.begin(); turn1!= turns.end(); ++turn1 ) {
	//   StubFrags const & turn1_frags( all_turn_frags.find( *turn1 )->second );
	//   for ( Sizes::const_iterator h2len= helixlens.begin(); h2len != helixlens.end(); ++h2len ) {
	//    StubFrags const & helix2_frags( all_helix_frags.find(*h2len)->second );
	//    for ( strings::const_iterator turn2= turns.begin(); turn2!= turns.end(); ++turn2 ) {
	//     StubFrags const & turn2_frags( all_turn_frags.find( *turn2 )->second );

	//     Size pdbcounter(0);
	//     for ( Size n=1; n<= nstruct(); ++n ) {

	//      Stub const start_stub;

	//      StubFrag const h1( random_element( helix1_frags ) ), h2( random_element( helix2_frags ) ),
	//       t1( random_element( turn1_frags ) ), t2( random_element( turn2_frags ) );
	HelixHelixFrag hh_frag; // re-used

	for ( Size nn=1; nn<= nstruct(); ++nn ) {
		Size const helix1_len( random_element( helixlens ) ), helix2_len( random_element( helixlens ) );
		string const turn1( random_element( turns )), turn2( random_element( turns ) );

		if ( ( !turn1s.empty() && !has_element( turn1s, turn1 ) ) ||
				( !turn2s.empty() && !has_element( turn2s, turn2 ) ) ) {
			--nn;
			continue;
		}

		StubFrags const & helix1_frags( all_helix_frags.find( helix1_len )->second );
		StubFrags const & turn1_frags( all_turn_frags.find( turn1 )->second );
		StubFrags const & helix2_frags( all_helix_frags.find(helix2_len)->second );
		StubFrags const & turn2_frags( all_turn_frags.find( turn2 )->second );

		Stub const start_stub;

		StubFrag const h1( random_element( helix1_frags ) ), h2( random_element( helix2_frags ) ),
			t1( random_element( turn1_frags ) ), t2( random_element( turn2_frags ) );

		Stub const stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

		Vector center, n, t;
		Real theta;
		get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

		Real const shift( n.dot(t) );

		TR.Trace << "transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
			F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) << endl;

		Real const theta_dev( theta / target_theta );

		if ( theta_dev > 0.8 && theta_dev < 1.25 && fabs( shift ) < max_shift ) {
			// reconstruct coords
			Vectors coords;
			Stub stub( start_stub );
			Sizes seg_begins, seg_ends;
			for ( Size n=1; n<= nrepeat; ++n ) {
				for ( Size r=1; r<= 4; ++r ) {
					seg_begins.push_back( coords.size()+1 );
					StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
					for ( Size i=1; i<= f.coords.size(); ++i ) {
						coords.push_back( stub.local2global( f.coords[i] ) );
					}
					stub = f.rt.make_jump( stub );
					seg_ends.push_back( coords.size() );
				}
			}
			Stub const final_stub( stub );
			Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
			Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
			Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
			Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
			runtime_assert( coords.size() == nrepeat * repeatlen );
			runtime_assert( seg_begins.size() == nrepeat*4 );
			runtime_assert( seg_ends.size() == nrepeat*4 );


			/// look for clashes between c-alphas in first two repeats
			// h1-1 to
			Real const clash_dis2( 3.75*3.75 ), close_contact_dis2( 6.5*6.5 ), long_contact_dis2( 8*8 );
			Size n_clash(0), n_close_contacts(0), n_long_contacts(0);
			Real dis2;
			for ( Size ii=1; ii<= 8; ++ii ) { // 8 segments in the first 2 repeats
				for ( Size jj=ii+3; jj<= 8; ++jj ) {
					for ( Size i=seg_begins[ii]; i<= seg_ends[ii]; ++i ) {
						for ( Size j=seg_begins[jj]; j<= seg_ends[jj]; ++j ) {
							dis2 = coords[i].distance_squared( coords[j] );
							if ( dis2 < clash_dis2 ) ++n_clash;
							else { // only consider non-clashes as contacts
								if ( ii%2==1 && jj%2==1 ) { // both segments are helices
									if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
									if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
								}
							}
						}
					}
				}
			}

			TR.Trace << "clash_transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
				F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
				" n_clash: " << I(4,n_clash) << endl;

			if ( n_clash <= max_clash &&
					n_close_contacts >= min_close_contacts &&
					n_long_contacts >= min_long_contacts ) {

				// compute helix-helix similarities
				bools found_contacts;
				Sizes nsims;
				for ( Size r=1; r<= 2; ++r ) {
					Vectors h1_coords, h2_coords;
					Size const h1_seg_index( r == 1 ? 1 : 3 ), h2_seg_index( r == 1 ? 5 : 7 );
					for ( Size i= seg_begins[ h1_seg_index ]; i<= seg_ends[ h1_seg_index ]; ++i ) {
						h1_coords.push_back( coords[i] );
					}
					for ( Size i= seg_begins[ h2_seg_index ]; i<= seg_ends[ h2_seg_index ]; ++i ) {
						h2_coords.push_back( coords[i] );
					}

					bool const found_contact( get_helix_contact_coords( h1_coords, h2_coords,
						hh_frag.h1_coords, hh_frag.h2_coords ) );
					Size nsim( 0 );
					if ( found_contact ) {
						for ( Size i=1; i<= hh_frags.size(); ++i ) {
							Real const dis( helix_helix_frag_distance( hh_frag, hh_frags[i] ) );
							TR.Trace << "hh_frag_dis: " << F(9,3,dis) << endl;
							if ( dis <= max_hh_frag_dis ) ++nsim;
						}
					}
					nsims.push_back( nsim );
					found_contacts.push_back( found_contact );
				} // both helices




				Real helix1_dist, helix1_twist, helix2_dist, helix2_twist;
				Size const base_repeat(3);
				compute_helix_axis_angles( coords, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
					helix1_dist, helix1_twist, helix2_dist, helix2_twist );

				// write out a C-alpha PDB file...
				string const outfilename( "ca_coords_"+string_of(helix1_len)+(turn1)+string_of(helix2_len)+(turn2)+ \
					"_N"+string_of(nn)+".pdb" );
				Real const handedness( get_chirality( repeatlen, coords ) );

				Real helix1_strain, helix2_strain, turn1_strain, turn2_strain, inner_helix_strain, outer_helix_strain,
					inner_turn_strain, outer_turn_strain, inner_turn_norm_strain, outer_turn_norm_strain;

				Real tmp;
				compute_helix_strain( 1, helix1_len, coords, helix1_strain, tmp );

				compute_helix_strain( helix1_len + turn1.size() + 1,
					helix1_len + turn1.size() + helix2_len, coords, helix2_strain, tmp );

				compute_turn_strain( turn1, helix1_len+1, coords, turn1_strain, tmp );
				compute_turn_strain( turn2, helix1_len+turn1.size()+helix2_len+1, coords, turn2_strain, tmp );

				Size inner_helix_len, outer_helix_len, inner_nsim, outer_nsim;
				string inner_turn, outer_turn, inner_turn_tag, outer_turn_tag, inner_helix_tag, outer_helix_tag;
				Real inner_helix_twist, outer_helix_twist, inner_helix_dist, outer_helix_dist;
				bool inner_fc, outer_fc;

				if ( helix1_dist < helix2_dist ) {
					inner_helix_len = helix1_len;
					outer_helix_len = helix2_len;
					inner_helix_dist = helix1_dist;
					outer_helix_dist = helix2_dist;
					inner_helix_twist = helix1_twist;
					outer_helix_twist = helix2_twist;
					inner_turn = turn1;
					outer_turn = turn2;
					inner_turn_tag = t1.id;
					outer_turn_tag = t2.id;
					inner_helix_tag = h1.id;
					outer_helix_tag = h2.id;
					inner_helix_strain = helix1_strain;
					outer_helix_strain = helix2_strain;
					inner_turn_strain = turn1_strain;
					outer_turn_strain = turn2_strain;
					inner_nsim = nsims[1];
					outer_nsim = nsims[2];
					inner_fc = found_contacts[1];
					outer_fc = found_contacts[2];
					inner_turn_norm_strain = t1.norm_strain;
					outer_turn_norm_strain = t2.norm_strain;
				} else {
					inner_helix_len = helix2_len;
					outer_helix_len = helix1_len;
					inner_helix_dist = helix2_dist;
					outer_helix_dist = helix1_dist;
					inner_helix_twist = helix2_twist;
					outer_helix_twist = helix1_twist;
					inner_turn = turn2;
					outer_turn = turn1;
					inner_turn_tag = t2.id;
					outer_turn_tag = t1.id;
					inner_helix_tag = h2.id;
					outer_helix_tag = h1.id;
					inner_helix_strain = helix2_strain;
					outer_helix_strain = helix1_strain;
					inner_turn_strain = turn2_strain;
					outer_turn_strain = turn1_strain;
					inner_nsim = nsims[2];
					outer_nsim = nsims[1];
					inner_fc = found_contacts[2];
					outer_fc = found_contacts[1];
					inner_turn_norm_strain = t2.norm_strain;
					outer_turn_norm_strain = t1.norm_strain;
				}

				cout << "pdb_transform: " << inner_helix_len << ' ' << outer_helix_len <<
					' ' << inner_turn << ' ' << outer_turn <<
					F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
					" n2c_distance: "<< F(9,3,n2c_distance ) <<
					" n_clash: " << I(4,n_clash) <<
					" n_close_contacts: " << I(4,n_close_contacts) <<
					" n_long_contacts: " << I(4,n_long_contacts) <<
					" nsims: " << inner_nsim << ' ' << outer_nsim <<
					" found_contacts: " << inner_fc << ' ' << outer_fc <<
					" handedness: " << F(9,3,handedness ) <<
					" inner_twist: " << F(9,3,inner_helix_twist) <<
					" outer_twist: " << F(9,3,outer_helix_twist) <<
					" inner_dist: " << F(9,3,inner_helix_dist) <<
					" outer_dist: " << F(9,3,outer_helix_dist) <<
					" inner_turn_tag: " << inner_turn_tag <<
					" outer_turn_tag: " << outer_turn_tag <<
					" inner_helix_tag: " << inner_helix_tag <<
					" outer_helix_tag: " << outer_helix_tag <<
					" inner_turn_strain: " << F(9,3,inner_turn_strain) <<
					" outer_turn_strain: " << F(9,3,outer_turn_strain) <<
					" inner_helix_strain: " << F(9,3,inner_helix_strain) <<
					" outer_helix_strain: " << F(9,3,outer_helix_strain) <<
					" inner_turn_norm_strain: " << F(9,3,inner_turn_norm_strain) <<
					" outer_turn_norm_strain: " << F(9,3,outer_turn_norm_strain) <<
					' ' << outfilename << endl;
				if ( option[ my_options::output_pdb_files ] && ( inner_nsim > min_nsims || outer_nsim > min_nsims ) ) {
					write_ca_pdbfile( coords, outfilename );
				}
			} // clash/contacts filter
		} // theta/shift filter
	} // nstruct
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

void
helix_junction_test()
{
	// { // first go through the pdb and get some geometries of antiparallel coiled coils
	//  strings const files( start_files() );

	//  foreach_( string filename, files ) {
	//   Pose pose;
	//   pose_from_pdb( pose, filename );

	//  }
	//  exit(0);
	// }




	if ( false ) {
		Real const max_hh_frag_dis( 3.0 );
		HelixHelixFrags hh_frags;
		{ // hacking
			build_helix_transform_library( hh_frags );

			for ( Size i=1; i<= hh_frags.size(); ++i ) {
				Size nsim( 0 );
				for ( Size j=1; j<= hh_frags.size(); ++j ) {
					if ( j==i ) continue;
					Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
					if ( dis <= max_hh_frag_dis ) ++nsim;
				}
				TR.Trace << "nsim " << nsim << ' ' << i << ' ' << hh_frags.size() << endl;
			}
		}
	}

	// note: HelixStubFrag is just a StubFrag
	map< string, StubFrags > all_turn_frags;
	map< Size, HelixStubFrags > all_helix_frags;

	Size const max_frags( 5000 );

	Sizes const long_helixlens( make_vector1( 12, 13, 14, 15, 16, 17, 18 ) );
	Sizes const short_helixlens( make_vector1( 6, 7, 8, 9, 10, 11, 12 ) );
	Sizes all_helixlens;
	foreach_ ( Size hl, long_helixlens ) all_helixlens.push_back(hl);
	foreach_ ( Size hl, short_helixlens ) all_helixlens.push_back(hl);

	strings const all_turns( make_vector1( string("GB"), string("GBB"), string("BAB") ) );

	build_helix_library_v2( all_helixlens, max_frags, all_helix_frags );
	build_turn_library_v2( all_turns, all_turn_frags );


	Size const max_clash( 1 ),
		min_close_contacts_13( 4 ),
		min_close_contacts_57( 4 ),
		min_close_contacts_17( 4 ),
		min_long_contacts_13( 10 ),
		min_long_contacts_57( 10 ),
		min_long_contacts_17( 10 );
	//min_nsims( 10 );

	HelixHelixFrag hh_frag; // re-used

	HelixBarrelMultifunc barrel_func;

	//clock_t starttime( clock() ), min_clocks(0);

	for ( Size nn=1; nn<= nstruct(); ++nn ) {
		strings const turns( make_vector1( string("GB"),
			random_element( make_vector1( string("GBB"), string("BAB") ) ),
			string("GB") ) );
		Sizes const helixlens( make_vector1( random_element( long_helixlens ),
			random_element( short_helixlens ),
			random_element( short_helixlens ),
			random_element( long_helixlens ) ));

		Stub const start_stub;

		HelixStubFrag
			h1( random_element( all_helix_frags.find( helixlens[1] )->second ) ),
			h2( random_element( all_helix_frags.find( helixlens[2] )->second ) ),
			h3( random_element( all_helix_frags.find( helixlens[3] )->second ) ),
			h4( random_element( all_helix_frags.find( helixlens[4] )->second ) );
		StubFrag const
			t1( random_element( all_turn_frags.find( turns[1] )->second ) ),
			t2( random_element( all_turn_frags.find( turns[2] )->second ) ),
			t3( random_element( all_turn_frags.find( turns[3] )->second ) );

		StubFrags frags( make_vector1( h1, t1, h2, t2, h3, t3, h4 ) );
		Size const nfrags( frags.size() );

		//Stub stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));
		{
			// reconstruct coords
			Vectors coords;
			Stub stub( start_stub );
			Sizes seg_begins, seg_ends;
			foreach_ ( StubFrag const & f, frags ) {
				seg_begins.push_back( coords.size()+1 );
				for ( Size i=1; i<= f.coords.size(); ++i ) {
					coords.push_back( stub.local2global( f.coords[i] ) );
				}
				stub = f.rt.make_jump( stub );
				seg_ends.push_back( coords.size() );
			}

			Stub final_stub( stub );


			/// look for clashes between c-alphas in first two repeats
			// h1-1 to
			// reduce clash_dis2 from 4.5^2 to 3.75^2
			//
			Real const clash_dis2( 3.75*3.75 ), close_contact_dis2( 6.5*6.5 ), long_contact_dis2( 8*8 );
			Size n_clash(0),
				n_close_contacts_13(0),
				n_close_contacts_17(0),
				n_close_contacts_57(0),
				n_long_contacts_13(0),
				n_long_contacts_17(0),
				n_long_contacts_57(0);
			Real dis2;
			for ( Size ii=1; ii<= nfrags; ++ii ) {
				for ( Size jj=ii+2; jj<= nfrags; ++jj ) {
					for ( Size i=seg_begins[ii]; i<= seg_ends[ii]; ++i ) {
						for ( Size j=seg_begins[jj]; j<= seg_ends[jj]; ++j ) {
							dis2 = coords[i].distance_squared( coords[j] );
							if ( dis2 < clash_dis2 ) ++n_clash;
							else { // only consider non-clashes as contacts
								if ( ii%2==1 && jj%2==1 ) { // both segments are helices
									if ( dis2 < close_contact_dis2 ) {
										if      ( ii==1 && jj==3 ) ++n_close_contacts_13;
										else if ( ii==5 && jj==7 ) ++n_close_contacts_57;
										else if ( ii==1 && jj==7 ) ++n_close_contacts_17;
									}
									if ( dis2 < long_contact_dis2 ) {
										if      ( ii==1 && jj==3 ) ++n_long_contacts_13;
										else if ( ii==5 && jj==7 ) ++n_long_contacts_57;
										else if ( ii==1 && jj==7 ) ++n_long_contacts_17;
									}
								}
							}
						}
					}
				}
			}


			if ( n_clash <= max_clash &&
					n_close_contacts_13 >= min_close_contacts_13 &&
					n_close_contacts_57 >= min_close_contacts_57 &&
					n_close_contacts_17 >= min_close_contacts_17 &&
					n_long_contacts_13 >= min_long_contacts_13 &&
					n_long_contacts_57 >= min_long_contacts_57 &&
					n_long_contacts_17 >= min_long_contacts_17 ) {

				string const outfilename( "hj_N"+lead_zero_string_of(nn,4)+".pdb");
				write_ca_pdbfile( coords, outfilename );
				cout << "n_clash: " << n_clash <<
					" n_close_contacts: " << n_close_contacts_13 << ' ' << n_close_contacts_57 << ' ' << n_close_contacts_17 <<
					" n_long_contacts: " << n_long_contacts_13 << ' ' << n_long_contacts_57 << ' ' << n_long_contacts_17 <<
					' ' << outfilename << endl;
			}
		} // scope
	} // nstruct
}


///////////////////////////////////////////////////////////////////////////////

void
helix_frag_v2_test()
{
	bool const verbose( false );
	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++nsim;
			}
			//TR.Trace << "nsim " << nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}

	bool const depth_filter( option[ my_options::min_depth ].user() );
	Real const min_depth( depth_filter ? option[ my_options::min_depth ]() : 0.0 );

	bool const helix_angle_filter( option[ my_options::max_helix_angle ].user() );
	Real const max_helix_angle( helix_angle_filter ? option[ my_options::max_helix_angle ]() : 0.0 );

	map< string, StubFrags > all_turn_frags;
	map< Size, HelixStubFrags > all_helix_frags;


	strings bbtags;
	Sizes helixlens;
	strings turns;

	if ( option[ my_options::resample_bbs ].user() ) {
		bbtags = option[ my_options::resample_bbs ];
		for ( Size i=1; i<= bbtags.size(); ++i ) {
			strings const l( split_to_vector1( bbtags[i], ".") );
			runtime_assert( l.size() == 5 || l.size() == 6 ); // might not include nrepeat
			Size const h1len( int_of( l[1] ) ), h2len( int_of( l[3] ) );
			if ( !has_element( helixlens, h1len ) ) helixlens.push_back( h1len );
			if ( !has_element( helixlens, h2len ) ) helixlens.push_back( h2len );
			if ( !has_element( turns, l[2] ) ) turns.push_back( l[2] );
			if ( !has_element( turns, l[4] ) ) turns.push_back( l[4] );
		}
	} else {

		helixlens = option[ my_options::helix_lens ]();
		turns = option[ my_options::turns ];

		if ( option[ my_options::interpolate_helixlens ] ) {
			runtime_assert( helixlens.size() == 2 );
			Size const mn( min( helixlens ) ), mx( max( helixlens ) );
			helixlens.clear();
			for ( Size i=mn; i<= mx; ++i ) helixlens.push_back( i );
		}
	}

	bool const force_helixlen_delta( option[ my_options::helixlen_delta ].user() );
	int const helixlen_delta( option[ my_options::helixlen_delta ] );
	bool const force_max_helixlen_delta( option[ my_options::max_helixlen_delta ].user() );
	int const max_helixlen_delta( option[ my_options::max_helixlen_delta ] );

	Size const max_frags( max( Size(1000), nstruct() ) ); // per helixlen

	strings turn1s, turn2s;
	if ( option[ my_options::turn1s ].user() ) turn1s = option[ my_options::turn1s ]();
	if ( option[ my_options::turn2s ].user() ) turn2s = option[ my_options::turn2s ]();


	Sizes src_helixlens;
	if ( option[ my_options::src_helixlens ].user() ) {
		src_helixlens = option[ my_options::src_helixlens ]();
	}


	build_helix_library_v2( helixlens, max_frags, all_helix_frags, src_helixlens );
	build_turn_library_v2( turns, all_turn_frags );


	Sizes nrepeats( option[ my_options::nrepeats ]() );

	if ( nrepeats.size() == 2 ) {
		Size const mn( min( nrepeats ) ), mx( max( nrepeats ) );
		nrepeats.clear();
		for ( Size i=mn; i<= mx; ++i ) nrepeats.push_back( i );
	}

	Sizes const cmdline_nrepeats( nrepeats );

	char cmdline_force_hand('-');
	if ( option[ my_options::target_hand ].user() ) {
		runtime_assert( option[ my_options::target_hand ].size() == 1 );
		cmdline_force_hand = option[ my_options::target_hand ]()[0];
	}

	// Size const nrepeat( option[ my_options::nrepeat ] );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	Real const max_shift( 2 );

	Size const max_clash( 1 ), min_close_contacts( 5 ), min_long_contacts( 10 ), min_nsims( 10 );

	HelixHelixFrag hh_frag; // re-used

	HelixBarrelMultifunc barrel_func;

	clock_t starttime( clock() ), min_clocks(0);

	for ( Size nn=1; nn<= nstruct(); ++nn ) {
		Size helix1_len, helix2_len;
		string turn1, turn2;
		char force_hand( cmdline_force_hand );
		Sizes nrepeats;

		if ( bbtags.size() ) {
			strings const l( split_to_vector1( random_element( bbtags ), "." ) );
			helix1_len = int_of( l[1] );
			helix2_len = int_of( l[3] );
			turn1 = l[2];
			turn2 = l[4];
			force_hand = l[5][0];
			runtime_assert( force_hand == 'R' || force_hand == 'L' );
			if ( l.size() == 5 ) {
				nrepeats = cmdline_nrepeats;
			} else {
				runtime_assert( l.size() == 6 );
				nrepeats = make_vector1( Size( int_of( l[6] ) ) );
			}
		} else {
			helix1_len = random_element( helixlens );
			helix2_len = random_element( helixlens );
			turn1 = random_element( turns );
			turn2 = random_element( turns );

			if ( ( !turn1s.empty() && !has_element( turn1s, turn1 ) ) ||
					( !turn2s.empty() && !has_element( turn2s, turn2 ) ) ) {
				--nn;
				continue;
			}
			nrepeats = cmdline_nrepeats;
		}

		if ( force_helixlen_delta && int( helix1_len ) - int( helix2_len ) != helixlen_delta ) {
			--nn;
			continue;
		}

		if ( force_max_helixlen_delta && abs( int( helix1_len ) - int( helix2_len ) ) > max_helixlen_delta ) {
			--nn;
			continue;
		}

		HelixStubFrags const & helix1_frags( all_helix_frags.find( helix1_len )->second );
		StubFrags const & turn1_frags( all_turn_frags.find( turn1 )->second );
		HelixStubFrags const & helix2_frags( all_helix_frags.find(helix2_len)->second );
		StubFrags const & turn2_frags( all_turn_frags.find( turn2 )->second );

		Stub const start_stub;

		HelixStubFrag h1( random_element( helix1_frags ) ), h2( random_element( helix2_frags ) ); // may be optimized
		StubFrag const t1( random_element( turn1_frags ) ), t2( random_element( turn2_frags ) );

		Stub stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

		Real theta, shift;
		Vector axis_center, axis_vector;
		{
			Vector t;
			get_stub_transform_data( start_stub, stop_stub, axis_center, axis_vector, t, theta );
			shift = axis_vector.dot(t);
		}

		if ( verbose ) {
			TR.Trace << "transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
				F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) << endl;
		}

		Size nrepeat(0);
		Real min_dev( 1e6 );
		for ( Sizes::const_iterator nr = nrepeats.begin(); nr != nrepeats.end(); ++nr ) {
			Real const target_theta( (2*numeric::constants::d::pi)/(*nr) );
			Real const dev( fabs( log( theta / target_theta ) ) );
			if ( dev < min_dev ) {
				min_dev = dev;
				nrepeat = *nr;
			}
		}

		Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
		Real theta_dev( theta / target_theta );

		if ( theta_dev > 0.8 && theta_dev < 1.25 && fabs( shift ) < max_shift ) {
			// reconstruct coords
			Vectors coords, stub_vs;
			Stub stub( start_stub );
			Sizes seg_begins, seg_ends;
			for ( Size n=1; n<= nrepeat; ++n ) {
				for ( Size r=1; r<= 4; ++r ) {
					stub_vs.push_back( stub.v );
					seg_begins.push_back( coords.size()+1 );
					StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
					for ( Size i=1; i<= f.coords.size(); ++i ) {
						coords.push_back( stub.local2global( f.coords[i] ) );
					}
					stub = f.rt.make_jump( stub );
					seg_ends.push_back( coords.size() );
				}
			}

			stub_vs.push_back( stub.v ); // include final_stub.v

			Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );

			// compute approx helix angles
			Real h1_axis_angle(0), h2_axis_angle(0);
			{
				Vector const h1_axis( ( stub_vs[2] - stub_vs[1] ).normalized()),
					h2_axis( ( stub_vs[4]-stub_vs[3] ).normalized() );
				h1_axis_angle = degrees( std::acos( numeric::sin_cos_range( fabs( h1_axis.dot( axis_vector ) ) ) ) );
				h2_axis_angle = degrees( std::acos( numeric::sin_cos_range( fabs( h2_axis.dot( axis_vector ) ) ) ) );
			}

			if ( helix_angle_filter && ( h1_axis_angle > max_helix_angle || h2_axis_angle > max_helix_angle ) ) {
				continue;
			}


			// compute depth along symmetry axis
			Real depth(0);
			{
				Real min_proj(1e6), max_proj(-1e6);
				foreach_ ( Vector const & v, coords ) {
					min_proj = min( min_proj, axis_vector.dot( v - axis_center ) );
					max_proj = max( max_proj, axis_vector.dot( v - axis_center ) );
				}
				depth = max_proj - min_proj;
			}

			if ( depth_filter && depth < min_depth ) {
				TR.Trace << "short: " << depth << ' ' << min_depth << endl;
				continue;
			}

			// compute centroid, get twist around centroid
			runtime_assert( stub_vs.size() == nrepeat*4 + 1 );
			Real const initial_total_rotation( get_total_rotation_degrees( stub_vs, axis_center, axis_vector ) );
			{
				Vectors tmpcoords( coords );
				tmpcoords.push_back( coords.front() );
				runtime_assert( tmpcoords.size() == nrepeat*repeatlen + 1 );
				Real const initial_total_rotation_redo( get_total_rotation_degrees( tmpcoords, axis_center, axis_vector ) );
				TR.Trace << "initial_total_rotation: " << F(9,3,initial_total_rotation) <<
					" initial_total_rotation_redo: " << F(9,3,initial_total_rotation_redo) <<
					" total_theta: " << degrees( theta * nrepeat ) << endl;
			}


			if ( force_hand != '-' ) {
				Real const handedness( get_chirality( repeatlen, coords ) );
				if ( ( force_hand == 'R' && handedness<0 ) || ( force_hand == 'L' && handedness>0 ) ) continue;
			}

			Stub final_stub( stub );
			// Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
			// Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
			// Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
			Real n2c_distance( final_stub.v.distance( start_stub.v ) );
			runtime_assert( coords.size() == nrepeat * repeatlen );
			runtime_assert( seg_begins.size() == nrepeat*4 );
			runtime_assert( seg_ends.size() == nrepeat*4 );


			/// look for clashes between c-alphas in first two repeats
			// h1-1 to
			// reduce clash_dis2 from 4.5^2 to 3.75^2
			//
			Real const clash_dis2( 3.75*3.75 ), close_contact_dis2( 6.5*6.5 ), long_contact_dis2( 8*8 );
			Size n_clash(0), n_close_contacts(0), n_long_contacts(0);
			Real dis2;
			for ( Size ii=1; ii<= 8; ++ii ) { // 8 segments in the first 2 repeats
				for ( Size jj=ii+2; jj<= 8; ++jj ) { // reduce from ii+3 to ii+2
					for ( Size i=seg_begins[ii]; i<= seg_ends[ii]; ++i ) {
						for ( Size j=seg_begins[jj]; j<= seg_ends[jj]; ++j ) {
							dis2 = coords[i].distance_squared( coords[j] );
							if ( dis2 < clash_dis2 ) ++n_clash;
							else { // only consider non-clashes as contacts
								if ( ii%2==1 && jj%2==1 ) { // both segments are helices
									if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
									if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
								}
							}
						}
					}
				}
			}

			if ( verbose ) {
				TR.Trace << "clash_transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
					F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
					" n_clash: " << I(4,n_clash) << endl;
			}

			Real final_total_rotation( initial_total_rotation );
			if ( n_clash <= max_clash &&
					n_close_contacts >= min_close_contacts &&
					n_long_contacts >= min_long_contacts ) {

				{ // try optimizing the helices to get small n2c distance
					//write_ca_pdbfile( coords, "beforemin.pdb" );

					barrel_func.set_nrepeat_and_frags_and_helix_lengths( nrepeat, h1, h2, t1, t2, helix1_len, helix2_len );
					optimization::Multivec params;
					barrel_func.params_from_helix_frags( h1, h2, params );

					// this should be the square of n2c_distance, well no not really
					//Real const start_func( barrel_func( params ) );
					barrel_func( params ); // need this?

					//TR.Trace << "same_dis? " << F(9,3,n2c_distance ) << F(9,3,sqrt( start_func ) ) << endl;

					Real const tolerance( 1e-3 );
					Size iterations;
					Real final_func_value;
					clock_t min_start( clock() );
					optimization::powell( params, barrel_func, tolerance, iterations, final_func_value );
					min_clocks += ( clock() - min_start );

					barrel_func.helix_frags_from_params( params, h1, h2 ); // modifies h1,h2

					/// update everything
					stop_stub = ( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

					{
						Vector t;
						//Real theta;
						get_stub_transform_data( start_stub, stop_stub, axis_center, axis_vector, t, theta );

						shift = axis_vector.dot(t);
					}

					TR.Trace << "transform_afteropt: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
						F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) << endl;

					theta_dev = theta / target_theta;

					// reconstruct coords
					//Vectors coords;
					stub = start_stub;
					stub_vs.clear();
					Size resnum(0);
					for ( Size n=1; n<= nrepeat; ++n ) {
						for ( Size r=1; r<= 4; ++r ) {
							stub_vs.push_back( stub.v );
							StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
							for ( Size i=1; i<= f.coords.size(); ++i ) {
								++resnum;
								coords[ resnum ] = stub.local2global( f.coords[i] );
							}
							stub = f.rt.make_jump( stub );
						}
					}
					stub_vs.push_back( stub.v );
					final_stub = stub;
					n2c_distance = final_stub.v.distance( start_stub.v );
					Real const total_rotation( get_total_rotation_degrees( stub_vs, axis_center, axis_vector ) );
					final_total_rotation = total_rotation;
					{
						Real min_proj(1e6), max_proj(-1e6);
						foreach_ ( Vector const & v, coords ) {
							min_proj = min( min_proj, axis_vector.dot( v - axis_center ) );
							max_proj = max( max_proj, axis_vector.dot( v - axis_center ) );
						}
						depth = max_proj - min_proj;
					}

					{
						Vector const h1_axis( ( stub_vs[2] - stub_vs[1] ).normalized()),
							h2_axis( ( stub_vs[4]-stub_vs[3] ).normalized() );
						h1_axis_angle = degrees( std::acos( numeric::sin_cos_range( fabs( h1_axis.dot( axis_vector ) ) ) ) );
						h2_axis_angle = degrees( std::acos( numeric::sin_cos_range( fabs( h2_axis.dot( axis_vector ) ) ) ) );
					}


					// TR.Trace << "new_n2c_distance: " << F(9,3,n2c_distance) << F(9,3,new_n2c_distance) <<
					//  " final_func_value: " << F(9,3,sqrt( final_func_value ) ) << endl;

					//write_ca_pdbfile( coords, "aftermin.pdb" );

					Size const old_n_clash( n_clash ), old_n_close_contacts( n_close_contacts ),
						old_n_long_contacts( n_long_contacts);
					n_clash = 0; n_close_contacts = 0; n_long_contacts = 0;
					Real dis2;
					for ( Size ii=1; ii<= 8; ++ii ) { // 8 segments in the first 2 repeats
						for ( Size jj=ii+2; jj<= 8; ++jj ) {
							for ( Size i=seg_begins[ii]; i<= seg_ends[ii]; ++i ) {
								for ( Size j=seg_begins[jj]; j<= seg_ends[jj]; ++j ) {
									dis2 = coords[i].distance_squared( coords[j] );
									if ( dis2 < clash_dis2 ) ++n_clash;
									else { // only consider non-clashes as contacts
										if ( ii%2==1 && jj%2==1 ) { // both segments are helices
											if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
											if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
										}
									}
								}
							}
						}
					}

					Real const tot_time = ((double) clock() - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
					Real const min_time = ((double) min_clocks )/( CLOCKS_PER_SEC*60 ); // in minutes


					cout << "old2new tot_time: " << F(9,3,tot_time) << " min_time: " << F(9,3,min_time) <<
						" n_clash: " << I(3,old_n_clash) << I(3,n_clash) <<
						" n_close_contacts: " << I(3,old_n_close_contacts) << I(3,n_close_contacts) <<
						" n_long_contacts: " << I(3,old_n_long_contacts) << I(3,n_long_contacts) <<
						" total_rot: " << F(9,3,initial_total_rotation) << F(9,3,total_rotation) << endl;

					if ( !( n_clash <= max_clash &&
							n_close_contacts >= min_close_contacts &&
							n_long_contacts >= min_long_contacts ) ) continue;

				}


				// compute helix-helix similarities
				bools found_contacts;
				Sizes nsims;
				for ( Size r=1; r<= 2; ++r ) {
					Vectors h1_coords, h2_coords;
					Size const h1_seg_index( r == 1 ? 1 : 3 ), h2_seg_index( r == 1 ? 5 : 7 );
					for ( Size i= seg_begins[ h1_seg_index ]; i<= seg_ends[ h1_seg_index ]; ++i ) {
						h1_coords.push_back( coords[i] );
					}
					for ( Size i= seg_begins[ h2_seg_index ]; i<= seg_ends[ h2_seg_index ]; ++i ) {
						h2_coords.push_back( coords[i] );
					}

					bool const found_contact( get_helix_contact_coords( h1_coords, h2_coords,
						hh_frag.h1_coords, hh_frag.h2_coords ) );
					Size nsim( 0 );
					if ( found_contact ) {
						for ( Size i=1; i<= hh_frags.size(); ++i ) {
							Real const dis( helix_helix_frag_distance( hh_frag, hh_frags[i] ) );
							//TR.Trace << "hh_frag_dis: " << F(9,3,dis) << endl;
							if ( dis <= max_hh_frag_dis ) ++nsim;
						}
					}
					nsims.push_back( nsim );
					found_contacts.push_back( found_contact );
				} // both helices




				Real helix1_dist, helix1_twist, helix2_dist, helix2_twist;
				Size const base_repeat(nrepeat==3 ? 2 : 3);
				compute_helix_axis_angles( coords, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
					helix1_dist, helix1_twist, helix2_dist, helix2_twist );

				// write out a C-alpha PDB file...
				string const outfilename( "ca_coords_"+string_of(helix1_len)+(turn1)+string_of(helix2_len)+(turn2)+ \
					"_N"+string_of(nn)+".pdb" );
				Real const handedness( get_chirality( repeatlen, coords ) );

				Real helix1_strain, helix2_strain, turn1_strain, turn2_strain, inner_helix_strain, outer_helix_strain,
					inner_turn_strain, outer_turn_strain, inner_turn_norm_strain, outer_turn_norm_strain;

				Real tmp;
				compute_helix_strain( 1, helix1_len, coords, helix1_strain, tmp );

				compute_helix_strain( helix1_len + turn1.size() + 1,
					helix1_len + turn1.size() + helix2_len, coords, helix2_strain, tmp );

				compute_turn_strain( turn1, helix1_len+1, coords, turn1_strain, tmp );
				compute_turn_strain( turn2, helix1_len+turn1.size()+helix2_len+1, coords, turn2_strain, tmp );

				Size inner_helix_len, outer_helix_len, inner_nsim, outer_nsim;
				string inner_turn, outer_turn, inner_turn_tag, outer_turn_tag, inner_helix_tag, outer_helix_tag;
				Real inner_helix_twist, outer_helix_twist, inner_helix_dist, outer_helix_dist;
				bool inner_fc, outer_fc;
				HelixParams inner_helix_params, outer_helix_params;

				if ( helix1_dist < helix2_dist ) {
					inner_helix_params = h1.hparams;
					outer_helix_params = h2.hparams;
					inner_helix_len = helix1_len;
					outer_helix_len = helix2_len;
					inner_helix_dist = helix1_dist;
					outer_helix_dist = helix2_dist;
					inner_helix_twist = helix1_twist;
					outer_helix_twist = helix2_twist;
					inner_turn = turn1;
					outer_turn = turn2;
					inner_turn_tag = t1.id;
					outer_turn_tag = t2.id;
					inner_helix_tag = h1.id;
					outer_helix_tag = h2.id;
					inner_helix_strain = helix1_strain;
					outer_helix_strain = helix2_strain;
					inner_turn_strain = turn1_strain;
					outer_turn_strain = turn2_strain;
					inner_nsim = nsims[1];
					outer_nsim = nsims[2];
					inner_fc = found_contacts[1];
					outer_fc = found_contacts[2];
					inner_turn_norm_strain = t1.norm_strain;
					outer_turn_norm_strain = t2.norm_strain;
				} else {
					inner_helix_params = h2.hparams;
					outer_helix_params = h1.hparams;
					inner_helix_len = helix2_len;
					outer_helix_len = helix1_len;
					inner_helix_dist = helix2_dist;
					outer_helix_dist = helix1_dist;
					inner_helix_twist = helix2_twist;
					outer_helix_twist = helix1_twist;
					inner_turn = turn2;
					outer_turn = turn1;
					inner_turn_tag = t2.id;
					outer_turn_tag = t1.id;
					inner_helix_tag = h2.id;
					outer_helix_tag = h1.id;
					inner_helix_strain = helix2_strain;
					outer_helix_strain = helix1_strain;
					inner_turn_strain = turn2_strain;
					outer_turn_strain = turn1_strain;
					inner_nsim = nsims[2];
					outer_nsim = nsims[1];
					inner_fc = found_contacts[2];
					outer_fc = found_contacts[1];
					inner_turn_norm_strain = t2.norm_strain;
					outer_turn_norm_strain = t1.norm_strain;
				}

				using numeric::conversions::degrees;
				cout << "pdb_transform: " << nrepeat << ' ' <<
					inner_helix_len << ' ' << outer_helix_len <<
					' ' << inner_turn << ' ' << outer_turn <<
					F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
					" n2c_distance: "<< F(9,3,n2c_distance ) <<
					" n_clash: " << I(4,n_clash) <<
					" n_close_contacts: " << I(4,n_close_contacts) <<
					" n_long_contacts: " << I(4,n_long_contacts) <<
					" nsims: " << inner_nsim << ' ' << outer_nsim <<
					" found_contacts: " << inner_fc << ' ' << outer_fc <<
					" handedness: " << F(9,3,handedness ) <<
					" inner_twist: " << F(9,3,inner_helix_twist) <<
					" outer_twist: " << F(9,3,outer_helix_twist) <<
					" inner_dist: " << F(9,3,inner_helix_dist) <<
					" outer_dist: " << F(9,3,outer_helix_dist) <<
					" inner_turn_tag: " << inner_turn_tag <<
					" outer_turn_tag: " << outer_turn_tag <<
					" inner_helix_tag: " << inner_helix_tag <<
					" outer_helix_tag: " << outer_helix_tag <<
					" inner_turn_strain: " << F(9,3,inner_turn_strain) <<
					" outer_turn_strain: " << F(9,3,outer_turn_strain) <<
					" inner_helix_strain: " << F(9,3,inner_helix_strain) <<
					" outer_helix_strain: " << F(9,3,outer_helix_strain) <<
					" inner_turn_norm_strain: " << F(9,3,inner_turn_norm_strain) <<
					" outer_turn_norm_strain: " << F(9,3,outer_turn_norm_strain) <<
					" inner_helix_params: " <<
					F(9,6,inner_helix_params.rise) <<
					F(9,3,degrees(inner_helix_params.twist)) <<
					F(9,3,degrees(inner_helix_params.tilt)) <<
					F(9,3,degrees(inner_helix_params.tilt_direction)) <<
					F(9,6,inner_helix_params.ca_distance) <<
					" outer_helix_params: " <<
					F(9,6,outer_helix_params.rise) <<
					F(9,3,degrees(outer_helix_params.twist)) <<
					F(9,3,degrees(outer_helix_params.tilt)) <<
					F(9,3,degrees(outer_helix_params.tilt_direction)) <<
					F(9,6,outer_helix_params.ca_distance) <<
					" total_rotation: " << F(9,3,final_total_rotation) <<
					" depth: " << F(9,3,depth) <<
					" h1_axis_angle: " << F(9,3,h1_axis_angle) <<
					" h2_axis_angle: " << F(9,3,h2_axis_angle) <<
					" cdens_close: " << F(9,3,Real(n_close_contacts)/(inner_helix_len+outer_helix_len)) <<
					" cdens_long: " << F(9,3,Real(n_long_contacts)/(inner_helix_len+outer_helix_len)) <<
					" cdens_tot: " << F(9,3,Real(2*n_close_contacts+n_long_contacts)/(inner_helix_len+outer_helix_len)) <<
					' ' << outfilename << endl;
				if ( option[ my_options::output_pdb_files ] && ( inner_nsim > min_nsims || outer_nsim > min_nsims ) ) {
					write_ca_pdbfile( coords, outfilename );
				}
			} // clash/contacts filter
		} // theta/shift filter
	} // nstruct
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
helix_frag_v3_test()
{
	bool const verbose( true );
	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++nsim;
			}
			//TR.Trace << "nsim " << nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	map< string, StubFrags > all_turn_frags;
	map< Size, HelixStubFrags > all_helix_frags;

	Sizes const all_helixlens( option[ my_options::helix_lens ]() );
	strings const all_turns( option[ my_options::turns ] );
	//strings const turns( make_vector1( string("GBB") ) );
	Size const max_frags( max( Size(1000), nstruct() ) ); // per helixlen

	build_helix_library_v2( all_helixlens, max_frags, all_helix_frags );
	build_turn_library_v2( all_turns, all_turn_frags );

	Size const nrepeat( option[ my_options::nrepeat ] );
	Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	Real const max_shift( 2 );

	Size const max_clash( 1 ), min_close_contacts( 5 ), min_long_contacts( 10 );//, min_nsims( 10 );

	HelixHelixFrag hh_frag; // re-used

	HelixBarrelMultifunc barrel_func;

	//clock_t starttime( clock() ), min_clocks(0);

	Size const nhelices(3 ), nsegs_per_repeat( 2* nhelices );

	for ( Size nn=1; nn<= nstruct(); ++nn ) {

		Sizes helixlens;
		strings turns;

		StubFrags combined_frags, turn_frags, helix_frags;

		Size repeatlen(0);
		string bbtag;
		for ( Size i=1; i<= nhelices; ++i ) {
			Size const helixlen( random_element( all_helixlens ) );
			helixlens.push_back( helixlen );
			repeatlen += helixlen;
			bbtag += lead_zero_string_of( helixlen, 2 )+".";

			StubFrag const h( random_element( all_helix_frags.find( helixlen )->second ) );
			helix_frags.push_back( h );
			combined_frags.push_back( h );

			string const turn( random_element( all_turns ) );
			turns.push_back( turn );
			repeatlen += turn.size();
			bbtag += turn+".";

			StubFrag const t( random_element( all_turn_frags.find( turn )->second ) );
			turn_frags.push_back( t );
			combined_frags.push_back( t );
		}



		Stub const start_stub;
		Stub stop_stub( start_stub );
		for ( Size i=1; i<= nsegs_per_repeat; ++i ) stop_stub = combined_frags[i].rt.make_jump( stop_stub );

		Real theta, shift;
		Vector toroid_axis_center, toroid_axis_vector;
		{
			Vector center, n, t;
			get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );
			shift = n.dot(t);
			toroid_axis_vector = n;
			toroid_axis_center = center;
		}

		if ( verbose ) {
			TR.Trace << "transform: " << bbtag <<
				F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) << endl;
		}

		Real theta_dev( theta / target_theta );

		if ( theta_dev > 0.8 && theta_dev < 1.25 && fabs( shift ) < max_shift ) {
			// reconstruct coords
			ostringstream helix_geometry_status; // hold info on angles, distances...
			Vectors coords;
			Stub stub( start_stub );
			Sizes seg_begins, seg_ends;
			vector1< Vectors > seg_vertices;
			for ( Size n=1; n<= nrepeat; ++n ) {
				for ( Size r=1; r<= combined_frags.size(); ++r ) {
					seg_begins.push_back( coords.size()+1 );
					StubFrag const & f( combined_frags[ r ] );
					for ( Size i=1; i<= f.coords.size(); ++i ) {
						coords.push_back( stub.local2global( f.coords[i] ) );
					}
					seg_vertices.push_back( make_vector1( stub.v ) );
					stub = f.rt.make_jump( stub );
					seg_vertices.back().push_back( stub.v );
					seg_ends.push_back( coords.size() );

					if ( n==1 && r%2==1 ) {
						Vector const midpoint( 0.5 * ( seg_vertices.back()[1] + seg_vertices.back()[2] ) ),
							axis( ( seg_vertices.back()[2] - seg_vertices.back()[1] ).normalized() );

						/// distance from midpoint to the symmetry axis
						Vector radiusv( midpoint - toroid_axis_center );
						radiusv -= toroid_axis_vector.dot( radiusv ) * toroid_axis_vector; // make normal to rotation axis vector
						Real const helix_distance( radiusv.length() );
						Real const helix_angle( degrees( acos( radiusv.normalized().dot( axis.normalized() ) ) ) );

						/// angle of helix axis twist
						Vector const & n( toroid_axis_vector ), &center( toroid_axis_center );
						Vector const p1( midpoint + axis ), p2( midpoint ), p3( midpoint - radiusv ),
							p4( ( axis.dot(n)>0 ) ? p3 + n : p3 - n );
						//TR.Trace << n.length() << ' ' << (center-p3).length() << ' ' << n.dot( (center - p3).normalized() ) << endl;
						runtime_assert( (center-p3).length() < 1e-3 ||
							fabs( n.dot( (center - p3).normalized() ) )>0.99 ); // confirm p3 is on the symmetry axis
						Real const helix_twist( numeric::dihedral_degrees( p1, p2, p3, p4 ) );
						helix_geometry_status << " helix" << (r+1)/2 <<"_geometry: " <<
							F(9,3,helix_distance) << F(9,3,helix_angle) << F(9,3,helix_twist);
					}
				}
			}
			Stub final_stub( stub );
			// Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
			// Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
			// Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
			Real n2c_distance( final_stub.v.distance( start_stub.v ) );
			//Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
			//Size const nsegs_per_repeat( combined_frags.size() );
			runtime_assert( nsegs_per_repeat == combined_frags.size() );
			runtime_assert( coords.size() == nrepeat * repeatlen );
			runtime_assert( seg_begins.size() == nrepeat*nsegs_per_repeat );
			runtime_assert( seg_ends.size() == nrepeat*nsegs_per_repeat );


			/// look for clashes between c-alphas in first two repeats
			// h1-1 to
			Real const clash_dis2( 3.75*3.75 ), close_contact_dis2( 6.5*6.5 ), long_contact_dis2( 8*8 );
			Size n_clash(0), n_close_contacts(0), n_long_contacts(0);
			Real dis2;
			for ( Size ii=1; ii<= 2*nsegs_per_repeat; ++ii ) { // 8 segments in the first 2 repeats
				for ( Size jj=ii+3; jj<= 2*nsegs_per_repeat; ++jj ) {
					for ( Size i=seg_begins[ii]; i<= seg_ends[ii]; ++i ) {
						for ( Size j=seg_begins[jj]; j<= seg_ends[jj]; ++j ) {
							dis2 = coords[i].distance_squared( coords[j] );
							if ( dis2 < clash_dis2 ) ++n_clash;
							else { // only consider non-clashes as contacts
								if ( ii%2==1 && jj%2==1 ) { // both segments are helices
									if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
									if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
								}
							}
						}
					}
				}
			}

			if ( verbose ) {
				TR.Trace << "clash_transform: " << bbtag <<
					F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
					" n_clash: " << I(4,n_clash) << endl;
			}

			if ( n_clash <= max_clash &&
					n_close_contacts >= min_close_contacts &&
					n_long_contacts >= min_long_contacts ) {

				{ // try optimizing the helices to get small n2c distance
					// get this part from helix_frag_v2_test
				}


				// compute helix-helix similarities
				bools found_contacts;
				Sizes nsims;
				Size min_nsims(100000);
				Size helix_contacts(0);

				for ( Size h1=1; h1<= 2*nhelices; ++h1 ) {
					for ( Size h2=h1+2; h2<= 2*nhelices; ++h2 ) {
						Vectors h1_coords, h2_coords;
						Size const h1_seg_index( 2*h1-1 ), h2_seg_index( 2*h2-1 );
						for ( Size i= seg_begins[ h1_seg_index ]; i<= seg_ends[ h1_seg_index ]; ++i ) {
							h1_coords.push_back( coords[i] );
						}
						for ( Size i= seg_begins[ h2_seg_index ]; i<= seg_ends[ h2_seg_index ]; ++i ) {
							h2_coords.push_back( coords[i] );
						}

						bool const found_contact( get_helix_contact_coords( seg_begins[ h1_seg_index ], seg_ends[ h1_seg_index ], coords,
							seg_begins[ h2_seg_index ], seg_ends[ h2_seg_index ], coords,
							hh_frag.h1_coords, hh_frag.h2_coords ) );
						Size nsim( 0 );
						if ( found_contact ) {
							for ( Size i=1; i<= hh_frags.size(); ++i ) {
								Real const dis( helix_helix_frag_distance( hh_frag, hh_frags[i] ) );
								if ( dis <= max_hh_frag_dis ) {
									//TR.Trace << "hh_frag_dis: " << F(9,3,dis) << endl;
									++nsim;
								}
							}
							++helix_contacts;
							min_nsims = min( min_nsims, nsim );
						}
						nsims.push_back( nsim );
						found_contacts.push_back( found_contact );
						TR.Trace << "nsims " << nsim << " fc: " << found_contact << endl;
					}
				} // both helices

				// write out a C-alpha PDB file...
				string outfilename( "ca_coords_");
				for ( Size i=1; i<= nhelices; ++i ) {
					outfilename += string_of(helixlens[i])+(turns[i]);
				}
				outfilename += "_N"+string_of(nn)+".pdb";
				Real const handedness( get_chirality( repeatlen, coords ) );

				Reals helix_strains, turn_strains;
				for ( Size i=1; i<= nhelices; ++i ) {
					Size const iseg( 2*i-1 );
					runtime_assert( seg_ends[iseg]-seg_begins[iseg]+1 == helixlens[i] );
					Real helix_strain, turn_strain, tmp;
					//cout << "compute_helix_strain " << i << ' ' << seg_begins[iseg] << ' '<< helixlens[i] << endl; fflush(stdout );
					compute_helix_strain( seg_begins[iseg], seg_ends[iseg], coords, helix_strain, tmp );
					//cout << "compute_turn_strain " << i << ' ' << turns[i] << ' ' << seg_begins[2*i] << endl; fflush(stdout );
					compute_turn_strain( turns[i], seg_begins[2*i], coords, turn_strain, tmp );
					helix_strains.push_back( helix_strain );
					turn_strains.push_back( turn_strain );
					//cout << "compute strains " << i << F(9,3,helix_strain) << F(9,3,turn_strain) << endl; fflush( stdout );
				}


				using numeric::conversions::degrees;
				cout << "pdb_transform: ";
				for ( Size i=1; i<= nhelices; ++i ) cout << ' ' << helixlens[i];
				for ( Size i=1; i<= nhelices; ++i ) cout << ' ' << turns[i];
				cout << F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
					" n2c_distance: "<< F(9,3,n2c_distance ) <<
					" n_clash: " << I(4,n_clash) <<
					" n_close_contacts: " << I(4,n_close_contacts) <<
					" n_long_contacts: " << I(4,n_long_contacts) <<
					" helix_contacts: " << helix_contacts <<
					" min_nsims: " << min_nsims <<
					//" nsims: " <<
					//" found_contacts: " << inner_fc << ' ' << outer_fc <<
					" handedness: " << F(9,3,handedness ) <<
					helix_geometry_status.str();
				for ( Size i=1; i<= nhelices; ++i ) {
					HelixParams const & hparams( helix_frags[i].hparams );
					cout << " helix" << i << "_params: " <<
						F(9,6,hparams.rise) <<
						F(9,3,degrees(hparams.twist)) <<
						F(9,3,degrees(hparams.tilt)) <<
						F(9,3,degrees(hparams.tilt_direction)) <<
						F(9,6,hparams.ca_distance);
					cout << " helix" << i << "_strain: " << F(9,3,helix_strains[i]);
				}
				for ( Size i=1; i<= nhelices; ++i ) {
					cout << " turn" << i << "_tag: " << turn_frags[i].id << " turn" << i << "_strain: " << F(9,3,turn_strains[i]);
				}
				cout << ' ' << outfilename << endl;
				fflush( stdout );
				if ( min_nsims>5 ) {
					//option[ my_options::output_pdb_files ] && ( inner_nsim > min_nsims || outer_nsim > min_nsims ) ) {
					write_ca_pdbfile( coords, outfilename );
				}
			} // clash/contacts filter
		} // theta/shift filter
	} // nstruct
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct HelixFragLine {
	Size nrepeat, helix1_len, helix2_len, index;
	string turn1, turn2, helix1_tag, helix2_tag, turn1_tag, turn2_tag, line;
	HelixParams h1params, h2params;
};

typedef vector1< HelixFragLine > HelixFragLines;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StubFrag
get_stub_frag_for_id(
	string const id,
	StubFrags const & frags
)
{
	for ( Size i=1; i<= frags.size(); ++i ) {
		if ( frags[i].id == id ) return frags[i];
	}
	utility_exit_with_message( "cant find frag for id "+id );
	return frags.front();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
reconstruct_helix_frag_test()
{

	Size const nrepeat( option[ my_options::nrepeat ] );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	string const logfile( start_file() );

	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++nsim;
			}
			TR.Trace << "nsim " << nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	Sizes helixlens;
	strings turns;

	HelixFragLines hf_lines;

	{
		ifstream data( logfile.c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings const l( split_to_vector1( line ) );
			HelixFragLine hf;
			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
			}
			hf_lines.push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
			if ( !has_element( helixlens, hf.helix1_len ) ) helixlens.push_back( hf.helix1_len );
			if ( !has_element( helixlens, hf.helix2_len ) ) helixlens.push_back( hf.helix2_len );
		}
	}


	map< string, StubFrags > all_turn_frags;
	map< Size, StubFrags > all_helix_frags;

	{
		Size const max_frags( 1000000 ); // want all the helix frags!
		build_helix_library( helixlens, max_frags, all_helix_frags );
		build_turn_library( turns, all_turn_frags );
	}

	for ( Size li=1; li<= hf_lines.size(); ++li ) {
		HelixFragLine const & hfline( hf_lines[li] );

		StubFrag h1( get_stub_frag_for_id( hfline.helix1_tag, all_helix_frags.find( hfline.helix1_len )->second ) );
		StubFrag h2( get_stub_frag_for_id( hfline.helix2_tag, all_helix_frags.find( hfline.helix2_len )->second ) );
		StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
		StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

		Stub const start_stub,
			stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

		Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );
		string const turn1( hfline.turn1 ), turn2( hfline.turn2 );

		Vector center, n, t;
		Real theta;
		get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

		Real const shift( n.dot(t) );

		TR.Trace << "transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
			F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) << endl;

		//Real const theta_dev( theta / target_theta );

		// reconstruct coords
		Vectors coords;
		Stub stub( start_stub );
		Sizes seg_begins, seg_ends;
		for ( Size n=1; n<= nrepeat; ++n ) {
			for ( Size r=1; r<= 4; ++r ) {
				seg_begins.push_back( coords.size()+1 );
				StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
				for ( Size i=1; i<= f.coords.size(); ++i ) {
					coords.push_back( stub.local2global( f.coords[i] ) );
				}
				stub = f.rt.make_jump( stub );
				seg_ends.push_back( coords.size() );
			}
		}
		Stub const final_stub( stub );
		Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
		Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
		Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
		Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
		runtime_assert( coords.size() == nrepeat * repeatlen );
		runtime_assert( seg_begins.size() == nrepeat*4 );
		runtime_assert( seg_ends.size() == nrepeat*4 );


		/// look for clashes between c-alphas in first two repeats
		// h1-1 to
		Real const clash_dis2( 3.75*3.75 ), close_contact_dis2( 6.5*6.5 ), long_contact_dis2( 8*8 );
		Size n_clash(0), n_close_contacts(0), n_long_contacts(0);
		Real dis2;
		for ( Size ii=1; ii<= 8; ++ii ) { // 8 segments in the first 2 repeats
			for ( Size jj=ii+3; jj<= 8; ++jj ) {
				for ( Size i=seg_begins[ii]; i<= seg_ends[ii]; ++i ) {
					for ( Size j=seg_begins[jj]; j<= seg_ends[jj]; ++j ) {
						dis2 = coords[i].distance_squared( coords[j] );
						if ( dis2 < clash_dis2 ) ++n_clash;
						else { // only consider non-clashes as contacts
							if ( ii%2==1 && jj%2==1 ) { // both segments are helices
								if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
								if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
							}
						}
					}
				}
			}
		}

		TR.Trace << "clash_transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
			F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
			" n_clash: " << I(4,n_clash) << endl;


		// compute helix-helix similarities
		bools found_contacts;
		Sizes nsims;
		HelixHelixFrag hh_frag;
		for ( Size r=1; r<= 2; ++r ) {
			Vectors h1_coords, h2_coords;
			Size const h1_seg_index( r == 1 ? 1 : 3 ), h2_seg_index( r == 1 ? 5 : 7 );
			for ( Size i= seg_begins[ h1_seg_index ]; i<= seg_ends[ h1_seg_index ]; ++i ) {
				h1_coords.push_back( coords[i] );
			}
			for ( Size i= seg_begins[ h2_seg_index ]; i<= seg_ends[ h2_seg_index ]; ++i ) {
				h2_coords.push_back( coords[i] );
			}

			bool const found_contact( get_helix_contact_coords( h1_coords, h2_coords,
				hh_frag.h1_coords, hh_frag.h2_coords ) );
			Size nsim( 0 );
			if ( found_contact ) {
				for ( Size i=1; i<= hh_frags.size(); ++i ) {
					Real const dis( helix_helix_frag_distance( hh_frag, hh_frags[i] ) );
					TR.Trace << "hh_frag_dis: " << F(9,3,dis) << endl;
					if ( dis <= max_hh_frag_dis ) ++nsim;
				}
			}
			nsims.push_back( nsim );
			found_contacts.push_back( found_contact );
		} // both helices




		Real helix1_dist, helix1_twist, helix2_dist, helix2_twist;
		Size const base_repeat(3);
		compute_helix_axis_angles( coords, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
			helix1_dist, helix1_twist, helix2_dist, helix2_twist );

		// write out a C-alpha PDB file...
		Real const handedness( get_chirality( repeatlen, coords ) );

		Real helix1_strain, helix2_strain, turn1_strain, turn2_strain, inner_helix_strain, outer_helix_strain,
			inner_turn_strain, outer_turn_strain;

		Real tmp;
		compute_helix_strain( 1, helix1_len, coords, helix1_strain, tmp );

		compute_helix_strain( helix1_len + turn1.size() + 1,
			helix1_len + turn1.size() + helix2_len, coords, helix2_strain, tmp );

		compute_turn_strain( turn1, helix1_len+1, coords, turn1_strain, tmp );
		compute_turn_strain( turn2, helix1_len+turn1.size()+helix2_len+1, coords, turn2_strain, tmp );

		Size inner_helix_len, outer_helix_len, inner_nsim, outer_nsim;
		string inner_turn, outer_turn, inner_turn_tag, outer_turn_tag, inner_helix_tag, outer_helix_tag;
		Real inner_helix_twist, outer_helix_twist, inner_helix_dist, outer_helix_dist;
		bool inner_fc, outer_fc;

		if ( helix1_dist < helix2_dist ) {
			inner_helix_len = helix1_len;
			outer_helix_len = helix2_len;
			inner_helix_dist = helix1_dist;
			outer_helix_dist = helix2_dist;
			inner_helix_twist = helix1_twist;
			outer_helix_twist = helix2_twist;
			inner_turn = turn1;
			outer_turn = turn2;
			inner_turn_tag = t1.id;
			outer_turn_tag = t2.id;
			inner_helix_tag = h1.id;
			outer_helix_tag = h2.id;
			inner_helix_strain = helix1_strain;
			outer_helix_strain = helix2_strain;
			inner_turn_strain = turn1_strain;
			outer_turn_strain = turn2_strain;
			inner_nsim = nsims[1];
			outer_nsim = nsims[2];
			inner_fc = found_contacts[1];
			outer_fc = found_contacts[2];
		} else {
			inner_helix_len = helix2_len;
			outer_helix_len = helix1_len;
			inner_helix_dist = helix2_dist;
			outer_helix_dist = helix1_dist;
			inner_helix_twist = helix2_twist;
			outer_helix_twist = helix1_twist;
			inner_turn = turn2;
			outer_turn = turn1;
			inner_turn_tag = t2.id;
			outer_turn_tag = t1.id;
			inner_helix_tag = h2.id;
			outer_helix_tag = h1.id;
			inner_helix_strain = helix2_strain;
			outer_helix_strain = helix1_strain;
			inner_turn_strain = turn2_strain;
			outer_turn_strain = turn1_strain;
			inner_nsim = nsims[2];
			outer_nsim = nsims[1];
			inner_fc = found_contacts[2];
			outer_fc = found_contacts[1];
		}

		string const outfilename( output_tag() + "ca_coords_"+string_of(helix1_len)+(turn1)+string_of(helix2_len)+(turn2)+ \
			"_reconstruct_L"+string_of(li)+".pdb" );

		cout << "pdb_transform: " << inner_helix_len << ' ' << outer_helix_len <<
			' ' << inner_turn << ' ' << outer_turn <<
			F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
			" n2c_distance: "<< F(9,3,n2c_distance ) <<
			" n_clash: " << I(4,n_clash) <<
			" n_close_contacts: " << I(4,n_close_contacts) <<
			" n_long_contacts: " << I(4,n_long_contacts) <<
			" nsims: " << inner_nsim << ' ' << outer_nsim <<
			" found_contacts: " << inner_fc << ' ' << outer_fc <<
			" handedness: " << F(9,3,handedness ) <<
			" inner_twist: " << F(9,3,inner_helix_twist) <<
			" outer_twist: " << F(9,3,outer_helix_twist) <<
			" inner_dist: " << F(9,3,inner_helix_dist) <<
			" outer_dist: " << F(9,3,outer_helix_dist) <<
			" inner_turn_tag: " << inner_turn_tag <<
			" outer_turn_tag: " << outer_turn_tag <<
			" inner_helix_tag: " << inner_helix_tag <<
			" outer_helix_tag: " << outer_helix_tag <<
			" inner_turn_strain: " << F(9,3,inner_turn_strain) <<
			" outer_turn_strain: " << F(9,3,outer_turn_strain) <<
			" inner_helix_strain: " << F(9,3,inner_helix_strain) <<
			" outer_helix_strain: " << F(9,3,outer_helix_strain) <<
			' ' << outfilename << endl;

		write_ca_pdbfile( coords, outfilename );

	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
reconstruct_helix_frag_v2_test()
{
	using numeric::conversions::degrees;

	string const logfile( start_file() );

	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++nsim;
			}
			//TR.Trace << "nsim " << nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	//Sizes helixlens;
	strings turns;

	HelixFragLines hf_lines;

	{
		using numeric::conversions::radians;
		ifstream data( logfile.c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings l( split_to_vector1( line ) );

			HelixFragLine hf;
			if ( l[9] == "n2c_distance:" ) { // new format
				hf.nrepeat = int_of( l[2] );
				l.erase( l.begin()+1 );
			} else {
				hf.nrepeat = option[ my_options::nrepeat ];
			}

			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
				else if ( l[i] == "inner_helix_params:" ) {
					hf.h1params.rise           = float_of( l[i+1] );
					hf.h1params.twist          = radians( float_of( l[i+2] ) );
					hf.h1params.tilt           = radians( float_of( l[i+3] ) );
					hf.h1params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h1params.ca_distance    = float_of( l[i+5] );
				} else if ( l[i] == "outer_helix_params:" ) {
					hf.h2params.rise           = float_of( l[i+1] );
					hf.h2params.twist          = radians( float_of( l[i+2] ) );
					hf.h2params.tilt           = radians( float_of( l[i+3] ) );
					hf.h2params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h2params.ca_distance    = float_of( l[i+5] );
				}
			}
			hf_lines.push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
			// if ( !has_element( helixlens, hf.helix1_len ) ) helixlens.push_back( hf.helix1_len );
			// if ( !has_element( helixlens, hf.helix2_len ) ) helixlens.push_back( hf.helix2_len );
		}
	}


	map< string, StubFrags > all_turn_frags;
	//map< Size, StubFrags > all_helix_frags;

	{
		//Size const max_frags( 1000000 ); // want all the helix frags!
		//build_helix_library( helixlens, max_frags, all_helix_frags );
		build_turn_library_v2( turns, all_turn_frags );
	}

	HelixBarrelMultifunc barrel_func;
	for ( Size li=1; li<= hf_lines.size(); ++li ) {
		HelixFragLine const & hfline( hf_lines[li] );
		Size const nrepeat( hfline.nrepeat );
		//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );

		Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );

		StubFrag h1,h2;
		helix_frag_from_helix_params_and_length( helix1_len, hfline.h1params, h1 );
		helix_frag_from_helix_params_and_length( helix2_len, hfline.h2params, h2 );
		StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
		StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

		Stub const start_stub,
			stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

		string const turn1( hfline.turn1 ), turn2( hfline.turn2 );

		Vector center, n, t;
		Real theta;
		get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

		Real const shift( n.dot(t) );

		TR.Trace << "transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
			F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) << endl;

		//Real const theta_dev( theta / target_theta );

		// reconstruct coords
		Vectors coords;
		Stub stub( start_stub );
		Sizes seg_begins, seg_ends;
		for ( Size n=1; n<= nrepeat; ++n ) {
			for ( Size r=1; r<= 4; ++r ) {
				seg_begins.push_back( coords.size()+1 );
				StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
				for ( Size i=1; i<= f.coords.size(); ++i ) {
					coords.push_back( stub.local2global( f.coords[i] ) );
				}
				stub = f.rt.make_jump( stub );
				seg_ends.push_back( coords.size() );
			}
		}
		Stub const final_stub( stub );
		Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
		Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
		Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
		Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
		runtime_assert( coords.size() == nrepeat * repeatlen );
		runtime_assert( seg_begins.size() == nrepeat*4 );
		runtime_assert( seg_ends.size() == nrepeat*4 );


		/// look for clashes between c-alphas in first two repeats
		// h1-1 to
		Real const clash_dis2( 3.75*3.75 ), close_contact_dis2( 6.5*6.5 ), long_contact_dis2( 8*8 );
		Size n_clash(0), n_close_contacts(0), n_long_contacts(0);
		Real dis2;
		for ( Size ii=1; ii<= 8; ++ii ) { // 8 segments in the first 2 repeats
			for ( Size jj=ii+3; jj<= 8; ++jj ) {
				for ( Size i=seg_begins[ii]; i<= seg_ends[ii]; ++i ) {
					for ( Size j=seg_begins[jj]; j<= seg_ends[jj]; ++j ) {
						dis2 = coords[i].distance_squared( coords[j] );
						if ( dis2 < clash_dis2 ) ++n_clash;
						else { // only consider non-clashes as contacts
							if ( ii%2==1 && jj%2==1 ) { // both segments are helices
								if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
								if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
							}
						}
					}
				}
			}
		}

		TR.Trace << "clash_transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
			F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
			" n_clash: " << I(4,n_clash) << endl;


		// compute helix-helix similarities
		bools found_contacts;
		Sizes nsims;
		HelixHelixFrag hh_frag;
		for ( Size r=1; r<= 2; ++r ) {
			Vectors h1_coords, h2_coords;
			Size const h1_seg_index( r == 1 ? 1 : 3 ), h2_seg_index( r == 1 ? 5 : 7 );
			for ( Size i= seg_begins[ h1_seg_index ]; i<= seg_ends[ h1_seg_index ]; ++i ) {
				h1_coords.push_back( coords[i] );
			}
			for ( Size i= seg_begins[ h2_seg_index ]; i<= seg_ends[ h2_seg_index ]; ++i ) {
				h2_coords.push_back( coords[i] );
			}

			bool const found_contact( get_helix_contact_coords( h1_coords, h2_coords,
				hh_frag.h1_coords, hh_frag.h2_coords ) );
			Size nsim( 0 );
			if ( found_contact ) {
				for ( Size i=1; i<= hh_frags.size(); ++i ) {
					Real const dis( helix_helix_frag_distance( hh_frag, hh_frags[i] ) );
					TR.Trace << "hh_frag_dis: " << F(9,3,dis) << endl;
					if ( dis <= max_hh_frag_dis ) ++nsim;
				}
			}
			nsims.push_back( nsim );
			found_contacts.push_back( found_contact );
		} // both helices




		Real helix1_dist, helix1_twist, helix2_dist, helix2_twist;
		Size const base_repeat(nrepeat==3?2:3);
		compute_helix_axis_angles( coords, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
			helix1_dist, helix1_twist, helix2_dist, helix2_twist );

		// write out a C-alpha PDB file...
		Real const handedness( get_chirality( repeatlen, coords ) );

		Real helix1_strain, helix2_strain, turn1_strain, turn2_strain, inner_helix_strain, outer_helix_strain,
			inner_turn_strain, outer_turn_strain;

		Real tmp;
		compute_helix_strain( 1, helix1_len, coords, helix1_strain, tmp );

		compute_helix_strain( helix1_len + turn1.size() + 1,
			helix1_len + turn1.size() + helix2_len, coords, helix2_strain, tmp );

		compute_turn_strain( turn1, helix1_len+1, coords, turn1_strain, tmp );
		compute_turn_strain( turn2, helix1_len+turn1.size()+helix2_len+1, coords, turn2_strain, tmp );

		Size inner_helix_len, outer_helix_len, inner_nsim, outer_nsim;
		string inner_turn, outer_turn, inner_turn_tag, outer_turn_tag, inner_helix_tag, outer_helix_tag;
		Real inner_helix_twist, outer_helix_twist, inner_helix_dist, outer_helix_dist;
		bool inner_fc, outer_fc;
		HelixParams inner_helix_params, outer_helix_params;

		if ( helix1_dist < helix2_dist ) {
			inner_helix_params = h1.hparams;
			outer_helix_params = h2.hparams;
			inner_helix_len = helix1_len;
			outer_helix_len = helix2_len;
			inner_helix_dist = helix1_dist;
			outer_helix_dist = helix2_dist;
			inner_helix_twist = helix1_twist;
			outer_helix_twist = helix2_twist;
			inner_turn = turn1;
			outer_turn = turn2;
			inner_turn_tag = t1.id;
			outer_turn_tag = t2.id;
			inner_helix_tag = h1.id;
			outer_helix_tag = h2.id;
			inner_helix_strain = helix1_strain;
			outer_helix_strain = helix2_strain;
			inner_turn_strain = turn1_strain;
			outer_turn_strain = turn2_strain;
			inner_nsim = nsims[1];
			outer_nsim = nsims[2];
			inner_fc = found_contacts[1];
			outer_fc = found_contacts[2];
		} else {
			inner_helix_params = h2.hparams;
			outer_helix_params = h1.hparams;
			inner_helix_len = helix2_len;
			outer_helix_len = helix1_len;
			inner_helix_dist = helix2_dist;
			outer_helix_dist = helix1_dist;
			inner_helix_twist = helix2_twist;
			outer_helix_twist = helix1_twist;
			inner_turn = turn2;
			outer_turn = turn1;
			inner_turn_tag = t2.id;
			outer_turn_tag = t1.id;
			inner_helix_tag = h2.id;
			outer_helix_tag = h1.id;
			inner_helix_strain = helix2_strain;
			outer_helix_strain = helix1_strain;
			inner_turn_strain = turn2_strain;
			outer_turn_strain = turn1_strain;
			inner_nsim = nsims[2];
			outer_nsim = nsims[1];
			inner_fc = found_contacts[2];
			outer_fc = found_contacts[1];
		}

		string const outfilename( output_tag() + "ca_coords_"+string_of(helix1_len)+(turn1)+string_of(helix2_len)+(turn2)+ \
			"_reconstruct_L"+string_of(li)+".pdb" );


		{ /// try making an all-bb-atom model
			Size const nres_protein( repeatlen + helix1_len ), nsegs( 5 );
			Pose pose;
			for ( Size i=1; i<= nres_protein; ++i ) {
				pose.append_residue_by_bond( *get_vanilla_protein_residue('A'), true );
			}
			for ( Size j=1; j<= nsegs; ++j ) {
				for ( Size i=seg_begins[j]; i<= seg_ends[j]; ++i ) {
					if ( j%2==1 ) { // helix
						// set generic helical torsions...
						pose.set_phi  ( i, -62.0 );
						pose.set_psi  ( i, -43.0 );
						pose.set_omega( i, 180.0 );
					} else {
						StubFrag const & t( j==2 ? t1 : t2 );
						Size const turnpos( i-seg_begins[j]+1 );
						pose.set_phi  ( i, t.torsions[ turnpos ][1] );
						pose.set_psi  ( i, t.torsions[ turnpos ][2] );
						pose.set_omega( i, t.torsions[ turnpos ][3] );
						make_sequence_change( i, aa_from_oneletter_code( t.sequence[ turnpos-1 ] ), pose );
					}
				}
			}
			add_lower_terminus_type_to_pose_residue( pose, 1 );
			add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

			// now minimize to reduce the coord error
			ResidueOP vrtrsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRT" ) ) );
			pose.append_residue_by_jump( *vrtrsd, repeatlen/2 );
			Size const vrtpos( pose.total_residue() );

			{
				FoldTree f( pose.total_residue() );
				f.new_jump( nres_protein/2, vrtpos, vrtpos-1 );
				f.reorder( vrtpos );
				pose.fold_tree( f );
			}

			{
				using namespace id;
				using namespace scoring::constraints; using namespace scoring::func;

				ConstraintSetOP cst_set( new ConstraintSet() );

				/// add cartesian constraints on coords
				Real const slope( 1.0 ), exponent( 2.0 ), distance_tolerance(0);
				FuncOP min_max_func( new MinMaxFunc( 0, distance_tolerance, slope, exponent ) );

				for ( Size i=1; i<= nres_protein; ++i ) {
					cst_set->add_constraint( ConstraintOP( new CoordinateConstraint( AtomID( pose.residue(i).atom_index("CA"), i ),
						AtomID( 1, vrtpos ),
						coords[i],
						min_max_func->clone() ) ) );
				}

				// if ( false ) { /// helix hbond constraints
				//  FuncOP hbond_min_max_func( new MinMaxFunc( 1.4, 2.3, 1.0, 2.0 /*power*/ ) );
				//  for ( Size j=1; j<= 5; j+=2 ) {
				//   for ( Size i=seg_begins[j]; i<= seg_ends[j]-4; ++i ) {
				//    /// i-i+4 hbonds
				//    AtomID
				//     id1( pose.residue(   i ).atom_index("O"), i ),
				//     id2( pose.residue( i+4 ).atom_index("H"), i+4 );
				//    Real const hbdis( pose.xyz(id1).distance( pose.xyz(id2) ) );
				//    TR.Trace << "hbdis: " << F(9,3,hbdis) << I(4,i) << I(4,i+4) << I(4,j) << endl;
				//    cst_set->add_constraint( new AtomPairConstraint( id1, id2, hbond_min_max_func->clone() ) );
				//   }
				//  }
				// }

				/// constraint the dihedrals in the 5th segment to be the same as those in the first
				HarmonicFuncOP dihedral_tether_func( new HarmonicFunc( 0.0, 3.0 /* sdev in degrees*/ ) );
				for ( Size i=1; i<= seg_ends[1]; ++i ) {
					Size const j( i+repeatlen );
					Residue const & rsd1( pose.residue(i) ), &rsd2( pose.residue( j ) );
					AtomID const
						n1( rsd1.atom_index("N"), i ), ca1( rsd1.atom_index("CA"), i ), c1( rsd1.atom_index("C"), i ),
						n2( rsd2.atom_index("N"), j ), ca2( rsd2.atom_index("CA"), j ), c2( rsd2.atom_index("C"), j );

					if ( i>1 ) { // PHI
						AtomID const
							prevc1( pose.residue(i-1).atom_index("C"), i-1 ), prevc2( pose.residue(j-1).atom_index("C"), j-1 );
						cst_set->add_constraint
							( ConstraintOP( new DihedralPairConstraint( prevc1, n1, ca1, c1,
							prevc2, n2, ca2, c2, dihedral_tether_func ) ) );
					}
					if ( i < seg_ends[1] ) { // PSI and OMEGA
						AtomID const
							nextn1( pose.residue(i+1).atom_index("N"), i+1 ), nextca1( pose.residue(i+1).atom_index("CA"), i+1 ),
							nextn2( pose.residue(j+1).atom_index("N"), j+1 ), nextca2( pose.residue(j+1).atom_index("CA"), j+1 );
						cst_set->add_constraint
							( ConstraintOP( new DihedralPairConstraint( n1, ca1, c1, nextn1,
							n2, ca2, c2, nextn2, dihedral_tether_func ) ) );
						cst_set->add_constraint
							( ConstraintOP( new DihedralPairConstraint( ca1, c1, nextn1, nextca1,
							ca2, c2, nextn2, nextca2, dihedral_tether_func ) ) );
					}
				}

				/// add torsion tethers on the loop coords
				Real const torsion_constraint_sdev_degrees_helix( 20.0 ),torsion_constraint_sdev_degrees_turn( 5.0 );
				//torsion_constraint_weight( 1.0 );
				for ( Size j=1; j<= nsegs; ++j ) {
					for ( Size i=seg_begins[j]; i<= seg_ends[j]; ++i ) {
						for ( Size k=1; k<= 3; ++k ) {
							TorsionID const id( i, BB, k );
							DOF_ID const dof( pose.conformation().dof_id_from_torsion_id( id ) );
							if ( !dof.valid() ) {
								TR.Trace << "invalid dof: " << id << ' '<< dof << endl;
								continue;
							}
							Real const target( pose.dof( dof ) );
							Real const sdev( j%2==1 ?
								numeric::conversions::radians( torsion_constraint_sdev_degrees_helix ) :
								numeric::conversions::radians( torsion_constraint_sdev_degrees_turn ) );
							Real const weight( 1.0 );
							Size const torsion_constraint_power( 2 );
							cst_set->add_dof_constraint( dof, FuncOP( new CircularPowerFunc( target, sdev, torsion_constraint_power,
								weight ) ) );
						}
					}
				}
				pose.constraint_set( cst_set );

				ScoreFunctionOP scorefxn( new ScoreFunction() );
				scorefxn->initialize_from_file( option[ my_options::reconstruct_bb_weights_file ]() );
				// scorefxn->set_weight( omega, 0.25 );
				// //scorefxn->set_weight( atom_pair_constraint, 0.5 );
				// scorefxn->set_weight( coordinate_constraint, 10.0 );
				// scorefxn->set_weight( dof_constraint, 1.0 );
				// scorefxn->set_weight( rama, 0.2 );

				MoveMapOP mm( new MoveMap() );
				mm->set_bb( false );
				mm->set_jump( true );

				protocols::minimization_packing::MinMoverOP min_mover
					( new protocols::minimization_packing::MinMover( mm, scorefxn, "dfpmin", 0.00001, true ) );
				//true, true )); // no deriv-check right now...


				{ // RB SUP
					//pose.dump_pdb(outfilename+".start.nosup.pdb" );
					//Real const start_score( (*scorefxn)( pose ) );
					(*scorefxn)( pose );
					// cout << "beforemin rb " << F(9,3,start_score) << ' ' <<
					//  pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;
					min_mover->apply( pose ); // just rigid body
					// cout << "aftermin rb " << F(9,3,start_score) << ' ' <<
					//  pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

					//pose.dump_pdb(outfilename+".start.aftersup.pdb" );
				}

				mm->set_bb( true );

				Real const final_coordinate_constraint_weight( scorefxn->get_weight( coordinate_constraint ) );
				Real coordinate_constraint_weight( final_coordinate_constraint_weight / ( 2*2*2*2*2 ) ); // 5 times
				for ( Size i=1; i<= 5; ++i ) {
					coordinate_constraint_weight *= 2;
					scorefxn->set_weight( coordinate_constraint, coordinate_constraint_weight );

					(*scorefxn)( pose );
					// cout << "beforemin " << i << F(9,3,start_score) << ' ' <<
					//  pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

					min_mover->apply( pose );

					(*scorefxn)( pose );
					// cout << "aftermin " << i << F(9,3,final_score) << ' ' <<
					// pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;
				}


				//pose.dump_pdb(outfilename+".bb.pdb" );

				/// make a full pose
				Pose fullpose;
				for ( Size i=1; i<= nrepeat*repeatlen; ++i ) {
					fullpose.append_residue_by_bond( *get_vanilla_protein_residue('A'), true );
				}
				for ( Size i=1; i<= nrepeat*repeatlen; ++i ) {
					// what is the base position?
					// base repeat starts at 3 and ends at repeatlen+2
					Size const basepos( (i+repeatlen-3)%repeatlen+3 );
					fullpose.set_phi  ( i, pose.phi  ( basepos ) );
					fullpose.set_psi  ( i, pose.psi  ( basepos ) );
					fullpose.set_omega( i, pose.omega( basepos ) );
				}

				fullpose.dump_pdb( outfilename+".fullbb.pdb");

			}
		}


		cout << "pdb_transform: " << inner_helix_len << ' ' << outer_helix_len <<
			' ' << inner_turn << ' ' << outer_turn <<
			F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) <<
			" n2c_distance: "<< F(9,3,n2c_distance ) <<
			" n_clash: " << I(4,n_clash) <<
			" n_close_contacts: " << I(4,n_close_contacts) <<
			" n_long_contacts: " << I(4,n_long_contacts) <<
			" nsims: " << inner_nsim << ' ' << outer_nsim <<
			" found_contacts: " << inner_fc << ' ' << outer_fc <<
			" handedness: " << F(9,3,handedness ) <<
			" inner_twist: " << F(9,3,inner_helix_twist) <<
			" outer_twist: " << F(9,3,outer_helix_twist) <<
			" inner_dist: " << F(9,3,inner_helix_dist) <<
			" outer_dist: " << F(9,3,outer_helix_dist) <<
			" inner_turn_tag: " << inner_turn_tag <<
			" outer_turn_tag: " << outer_turn_tag <<
			" inner_helix_tag: " << inner_helix_tag <<
			" outer_helix_tag: " << outer_helix_tag <<
			" inner_turn_strain: " << F(9,3,inner_turn_strain) <<
			" outer_turn_strain: " << F(9,3,outer_turn_strain) <<
			" inner_helix_strain: " << F(9,3,inner_helix_strain) <<
			" outer_helix_strain: " << F(9,3,outer_helix_strain) <<
			" inner_helix_params: " <<
			F(9,6,inner_helix_params.rise) <<
			F(9,3,degrees(inner_helix_params.twist)) <<
			F(9,3,degrees(inner_helix_params.tilt)) <<
			F(9,3,degrees(inner_helix_params.tilt_direction)) <<
			F(9,6,inner_helix_params.ca_distance) <<
			" outer_helix_params: " <<
			F(9,6,outer_helix_params.rise) <<
			F(9,3,degrees(outer_helix_params.twist)) <<
			F(9,3,degrees(outer_helix_params.tilt)) <<
			F(9,3,degrees(outer_helix_params.tilt_direction)) <<
			F(9,6,outer_helix_params.ca_distance) <<
			' ' << outfilename << endl;

		write_ca_pdbfile( coords, outfilename );

	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
reconstruct_and_design_test()
{
	runtime_assert( option[ my_options::unfolded_sasas ].user() ); // make sure have datafile for buried sasa

	/// for fullatom
	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	bool const polar_pore( option[ my_options::polar_pore ].user() );

	//ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
	ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) );
	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	Real const donut_energy_weight( option[ my_options::donut_energy_weight ] );
	DonutWholeEnergy donut_whole_energy;
	Donut1B_Energy donut_1b_energy;
	if ( donut_energy_weight > 1e-3 ) {
		//runtime_assert( barrel_mode );
		fa_scorefxn->add_extra_method( donut_whole, donut_energy_weight, donut_whole_energy );
		fa_scorefxn->add_extra_method( donut_1b   , donut_energy_weight, donut_1b_energy    );
	} else {
		runtime_assert( false ); // want donut weight on...
	}


	Size const design_cycles( option[ my_options::design_cycles ] );


	Size const base_repeat( option[ my_options::base_repeat ] );
	using numeric::conversions::degrees;

	Size const nrepeat( option[ my_options::nrepeat ] );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	string const logfile( start_file() );

	// Real const max_hh_frag_dis( 3.0 );
	// HelixHelixFrags hh_frags;
	// { // hacking
	//  build_helix_transform_library( hh_frags );

	//  for ( Size i=1; i<= hh_frags.size(); ++i ) {
	//   Size nsim( 0 );
	//   for ( Size j=1; j<= hh_frags.size(); ++j ) {
	//    if ( j==i ) continue;
	//    Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
	//    if ( dis <= max_hh_frag_dis ) ++nsim;
	//   }
	//   //TR.Trace << "nsim " << nsim << ' ' << i << ' ' << hh_frags.size() << endl;
	//  }
	// }


	//Sizes helixlens;

	strings turns;
	HelixFragLines hf_lines;
	{
		using numeric::conversions::radians;
		ifstream data( logfile.c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings const l( split_to_vector1( line ) );
			HelixFragLine hf;
			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
				else if ( l[i] == "inner_helix_params:" ) {
					hf.h1params.rise           = float_of( l[i+1] );
					hf.h1params.twist          = radians( float_of( l[i+2] ) );
					hf.h1params.tilt           = radians( float_of( l[i+3] ) );
					hf.h1params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h1params.ca_distance    = float_of( l[i+5] );
				} else if ( l[i] == "outer_helix_params:" ) {
					hf.h2params.rise           = float_of( l[i+1] );
					hf.h2params.twist          = radians( float_of( l[i+2] ) );
					hf.h2params.tilt           = radians( float_of( l[i+3] ) );
					hf.h2params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h2params.ca_distance    = float_of( l[i+5] );
				}
			}
			hf.index = hf_lines.size()+1;
			hf_lines.push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
		}
	}

	map< string, StubFrags > all_turn_frags;
	build_turn_library_v2( turns, all_turn_frags );

	random_permutation( hf_lines, numeric::random::rg() );

	//HelixBarrelMultifunc barrel_func;
	for ( Size li=1; li<= hf_lines.size(); ++li ) {

		HelixFragLine const & hfline( hf_lines[li] );

		Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );

		StubFrag h1,h2;
		helix_frag_from_helix_params_and_length( helix1_len, hfline.h1params, h1 );
		helix_frag_from_helix_params_and_length( helix2_len, hfline.h2params, h2 );
		StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
		StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

		Stub const start_stub,
			stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

		string const turn1( hfline.turn1 ), turn2( hfline.turn2 );


		string const bbtag( string_of( helix1_len ) + "." + turn1 + "." + string_of( helix2_len ) + "." + turn2 );

		string const simfile( shared_output_tag() +"recdesign_"+bbtag+".work" );

		if ( simfile_is_done( simfile ) ) continue;


		// Vector center, n, t;
		// Real theta;
		// get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

		// Real const shift( n.dot(t) );

		// TR.Trace << "transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
		//  F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) << endl;

		// Real const theta_dev( theta / target_theta );

		// reconstruct coords
		Vectors coords;
		Stub stub( start_stub );
		Sizes seg_begins, seg_ends;
		for ( Size n=1; n<= nrepeat; ++n ) {
			for ( Size r=1; r<= 4; ++r ) {
				seg_begins.push_back( coords.size()+1 );
				StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
				for ( Size i=1; i<= f.coords.size(); ++i ) {
					coords.push_back( stub.local2global( f.coords[i] ) );
				}
				stub = f.rt.make_jump( stub );
				seg_ends.push_back( coords.size() );
			}
		}
		Stub const final_stub( stub );
		Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
		Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
		//Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
		Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
		runtime_assert( coords.size() == nrepeat * repeatlen );
		runtime_assert( seg_begins.size() == nrepeat*4 );
		runtime_assert( seg_ends.size() == nrepeat*4 );


		Pose start_pose;
		Real reconstruct_rmsd(0);
		{ /// try making an all-bb-atom model
			Size const nres_protein( repeatlen + helix1_len ), nsegs( 5 );
			Pose pose;
			for ( Size i=1; i<= nres_protein; ++i ) {
				pose.append_residue_by_bond( *get_vanilla_protein_residue('A'), true );
			}
			for ( Size j=1; j<= nsegs; ++j ) {
				for ( Size i=seg_begins[j]; i<= seg_ends[j]; ++i ) {
					if ( j%2==1 ) { // helix
						// set generic helical torsions...
						pose.set_phi  ( i, -62.0 );
						pose.set_psi  ( i, -43.0 );
						pose.set_omega( i, 180.0 );
					} else {
						StubFrag const & t( j==2 ? t1 : t2 );
						Size const turnpos( i-seg_begins[j]+1 );
						pose.set_phi  ( i, t.torsions[ turnpos ][1] );
						pose.set_psi  ( i, t.torsions[ turnpos ][2] );
						pose.set_omega( i, t.torsions[ turnpos ][3] );
						make_sequence_change( i, aa_from_oneletter_code( t.sequence[ turnpos-1 ] ), pose );
					}
				}
			}
			add_lower_terminus_type_to_pose_residue( pose, 1 );
			add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );

			// now minimize to reduce the coord error
			Size vrtpos(0);
			{
				ResidueOP vrtrsd
					( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRT" ) ) );
				pose.append_residue_by_jump( *vrtrsd, repeatlen/2 );

				vrtpos = pose.total_residue();

				FoldTree f( pose.total_residue() );
				f.new_jump( nres_protein/2, vrtpos, vrtpos-1 );
				f.reorder( vrtpos );
				pose.fold_tree( f );
			} // scope

			using namespace id;
			using namespace scoring::constraints; using namespace scoring::func;

			ConstraintSetOP cst_set( new ConstraintSet() );

			/// add cartesian constraints on coords
			Real const slope( 1.0 ), exponent( 2.0 ), distance_tolerance(0);
			FuncOP min_max_func( new MinMaxFunc( 0, distance_tolerance, slope, exponent ) );

			for ( Size i=1; i<= nres_protein; ++i ) {
				cst_set->add_constraint( ConstraintOP( new CoordinateConstraint( AtomID( pose.residue(i).atom_index("CA"), i ),
					AtomID( 1, vrtpos ),
					coords[i],
					min_max_func->clone() ) ) );
			}

			// if ( false ) { /// helix hbond constraints
			//  FuncOP hbond_min_max_func( new MinMaxFunc( 1.4, 2.3, 1.0, 2.0 /*power*/ ) );
			//  for ( Size j=1; j<= 5; j+=2 ) {
			//   for ( Size i=seg_begins[j]; i<= seg_ends[j]-4; ++i ) {
			//    /// i-i+4 hbonds
			//    AtomID
			//     id1( pose.residue(   i ).atom_index("O"), i ),
			//     id2( pose.residue( i+4 ).atom_index("H"), i+4 );
			//    Real const hbdis( pose.xyz(id1).distance( pose.xyz(id2) ) );
			//    TR.Trace << "hbdis: " << F(9,3,hbdis) << I(4,i) << I(4,i+4) << I(4,j) << endl;
			//    cst_set->add_constraint( new AtomPairConstraint( id1, id2, hbond_min_max_func->clone() ) );
			//   }
			//  }
			// }

			/// constraint the dihedrals in the 5th segment to be the same as those in the first
			HarmonicFuncOP dihedral_tether_func( new HarmonicFunc( 0.0, 3.0 /* sdev in degrees*/ ) );
			for ( Size i=1; i<= seg_ends[1]; ++i ) {
				Size const j( i+repeatlen );
				Residue const & rsd1( pose.residue(i) ), &rsd2( pose.residue( j ) );
				AtomID const
					n1( rsd1.atom_index("N"), i ), ca1( rsd1.atom_index("CA"), i ), c1( rsd1.atom_index("C"), i ),
					n2( rsd2.atom_index("N"), j ), ca2( rsd2.atom_index("CA"), j ), c2( rsd2.atom_index("C"), j );

				if ( i>1 ) { // PHI
					AtomID const
						prevc1( pose.residue(i-1).atom_index("C"), i-1 ), prevc2( pose.residue(j-1).atom_index("C"), j-1 );
					cst_set->add_constraint
						( ConstraintOP( new DihedralPairConstraint( prevc1, n1, ca1, c1,
						prevc2, n2, ca2, c2, dihedral_tether_func ) ) );
				}
				if ( i < seg_ends[1] ) { // PSI and OMEGA
					AtomID const
						nextn1( pose.residue(i+1).atom_index("N"), i+1 ), nextca1( pose.residue(i+1).atom_index("CA"), i+1 ),
						nextn2( pose.residue(j+1).atom_index("N"), j+1 ), nextca2( pose.residue(j+1).atom_index("CA"), j+1 );
					cst_set->add_constraint
						( ConstraintOP( new DihedralPairConstraint( n1, ca1, c1, nextn1,
						n2, ca2, c2, nextn2, dihedral_tether_func ) ) );
					cst_set->add_constraint
						( ConstraintOP( new DihedralPairConstraint( ca1, c1, nextn1, nextca1,
						ca2, c2, nextn2, nextca2, dihedral_tether_func ) ) );
				}
			}

			/// add torsion tethers on the loop coords
			Real const torsion_constraint_sdev_degrees_helix( 20.0 ),torsion_constraint_sdev_degrees_turn( 5.0 );
			//torsion_constraint_weight( 1.0 );
			for ( Size j=1; j<= nsegs; ++j ) {
				for ( Size i=seg_begins[j]; i<= seg_ends[j]; ++i ) {
					for ( Size k=1; k<= 3; ++k ) {
						TorsionID const id( i, BB, k );
						DOF_ID const dof( pose.conformation().dof_id_from_torsion_id( id ) );
						if ( !dof.valid() ) {
							TR.Trace << "invalid dof: " << id << ' '<< dof << endl;
							continue;
						}
						Real const target( pose.dof( dof ) );
						Real const sdev( j%2==1 ?
							numeric::conversions::radians( torsion_constraint_sdev_degrees_helix ) :
							numeric::conversions::radians( torsion_constraint_sdev_degrees_turn ) );
						Real const weight( 1.0 );
						Size const torsion_constraint_power( 2 );
						cst_set->add_dof_constraint( dof, FuncOP( new CircularPowerFunc( target, sdev, torsion_constraint_power,
							weight )) );
					}
				}
			}
			pose.constraint_set( cst_set );

			ScoreFunctionOP scorefxn( new ScoreFunction() );
			scorefxn->initialize_from_file( option[ my_options::reconstruct_bb_weights_file ]() );
			// scorefxn->set_weight( omega, 0.25 );
			// //scorefxn->set_weight( atom_pair_constraint, 0.5 );
			// scorefxn->set_weight( coordinate_constraint, 10.0 );
			// scorefxn->set_weight( dof_constraint, 1.0 );
			// scorefxn->set_weight( rama, 0.2 );

			MoveMapOP mm( new MoveMap() );
			mm->set_bb( false );
			mm->set_jump( true );

			protocols::minimization_packing::MinMoverOP min_mover
				( new protocols::minimization_packing::MinMover( mm, scorefxn, "dfpmin", 0.00001, true ) );
			//true, true )); // no deriv-check right now...


			// RB SUP
			min_mover->apply( pose ); // just rigid body

			mm->set_bb( true );

			Real const final_coordinate_constraint_weight( scorefxn->get_weight( coordinate_constraint ) );
			Real coordinate_constraint_weight( final_coordinate_constraint_weight / ( 2*2*2*2*2 ) ); // 5 times
			for ( Size i=1; i<= 5; ++i ) {
				coordinate_constraint_weight *= 2;
				scorefxn->set_weight( coordinate_constraint, coordinate_constraint_weight );
				min_mover->apply( pose );
			}

			Real const final_score( (*scorefxn)( pose ) );
			cout << "reconstruct_scores: " << F(9,3,final_score) <<
				" line= " << hfline.index <<
				pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;


			/// make a full pose
			Pose fullpose;
			for ( Size i=1; i<= nrepeat*repeatlen; ++i ) {
				fullpose.append_residue_by_bond( *get_vanilla_protein_residue('A'), true );
			}
			for ( Size i=1; i<= nrepeat*repeatlen; ++i ) {
				// what is the base position?
				// base repeat starts at 3 and ends at repeatlen+2
				Size const basepos( (i+repeatlen-3)%repeatlen+3 );
				fullpose.set_phi  ( i, pose.phi  ( basepos ) );
				fullpose.set_psi  ( i, pose.psi  ( basepos ) );
				fullpose.set_omega( i, pose.omega( basepos ) );
			}

			runtime_assert( seg_begins.size() == 4 * nrepeat );
			for ( Size j=1; j<= 4*nrepeat; ++j ) {
				char const ss( j%2==1 ? 'H' : 'L' );
				for ( Size i=seg_begins[j]; i<= seg_ends[j]; ++i ) {
					fullpose.set_secstruct( i, ss );
				}
			}

			Vectors newcoords;
			for ( Size i=1; i<= repeatlen*nrepeat; ++i ) {
				newcoords.push_back( fullpose.residue(i).xyz("CA") );
			}
			runtime_assert( newcoords.size() == coords.size() );
			reconstruct_rmsd = numeric::model_quality::calc_rms( coords, newcoords );

			cout << "reconstruct_rmsd: " << F(9,3,reconstruct_rmsd) << endl;

			//fullpose.dump_pdb( outfilename+".fullbb.pdb");
			start_pose = fullpose;

			/// need to append a vrtbb rsd
			/// bit of a hack: add a virtual residue at the end, fold starting from there...
			//Residue const & lastrsd( pose.residue( nres_protein ) );
			{
				ResidueOP vrtrsd
					( conformation::ResidueFactory::create_residue( start_pose.residue_type_set_for_pose()->name_map("VRTBB")));
				start_pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...
				kinematics::FoldTree f( start_pose.total_residue() );
				f.reorder( start_pose.total_residue() );
				start_pose.fold_tree( f );
			}
			add_lower_terminus_type_to_pose_residue( start_pose, 1 );
			//add_upper_terminus_type_to_pose_residue( pose, nres_protein );
			start_pose.conformation().insert_chain_ending( start_pose.total_residue()-1 );

			conformation::symmetry::SymmetryInfo symminfo;

			setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo );

			pose::symmetry::make_symmetric_pose( start_pose, symminfo );

		} // scope for reconstructing an all-bb-atom model


		// ok, now have an all-bb model for starting designs with

		// try to identify some inward pointing positions
		Sizes inward_poslist;
		if ( polar_pore ) {
			Sizes const pp_params( option[ my_options::polar_pore ]() );
			Size inward_poslist_size( max( pp_params ) );
			runtime_assert( pp_params.size() == 2 );
			Pose pose( start_pose );
			Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);
			if ( helix1_len ) {
				runtime_assert( helix1_len + turn1.size() + helix2_len + turn2.size() == repeatlen );

				compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
					helix1_dist, helix1_twist, helix2_dist, helix2_twist );

				// look from some positions in the inner helix
				Size startpos(0), endpos(0);
				if ( helix1_dist<helix2_dist ) {
					startpos = 1;
					endpos = helix1_len;
				} else {
					startpos = helix1_len + turn1.size()+1;
					endpos = startpos+helix2_len-1;
				}

				Vector axis_vector, axis_center;
				{
					Size const repeatpos(3);
					Size const pos1( (base_repeat-1)*repeatlen + repeatpos ), pos2( base_repeat * repeatlen + repeatpos );
					Residue const & rsd1( pose.residue( pos1 ) ), &rsd2( pose.residue(pos2) );
					Stub const
						stub1( rsd1.xyz("N"), rsd1.xyz("CA"), rsd1.xyz("C") ),
						stub2( rsd2.xyz("N"), rsd2.xyz("CA"), rsd2.xyz("C") );

					Real theta;
					Vector t;
					get_stub_transform_data( stub1, stub2, axis_center, axis_vector, t, theta );
				}
				vector1< std::pair< Real, Size > > ll;
				for ( Size i=startpos; i<= endpos; ++i ) {
					Residue const & rsd( pose.residue(i) );
					Vector v( rsd.aa() == aa_gly ? rsd.xyz("2HA") : rsd.xyz("CB") );
					v -= rsd.xyz("CA");
					v.normalize(); // ca-cb vector, normalized
					Vector r( axis_center - rsd.xyz("CA") );
					r = r - axis_vector.dot(r) * axis_vector;
					r.normalize(); // vector from ca to toroid axis, normalized, perpendicular to the toroid_axis_vector
					Real const inward_angle( numeric::conversions::degrees( acos( r.dot(v) ) ) ); // 0 means pointing exactly in
					ll.push_back( make_pair( inward_angle, i ) );
				}
				std::sort( ll.begin(), ll.end() );
				for ( Size i=1; i<= ll.size() && i<= inward_poslist_size; ++i ) {
					inward_poslist.push_back( ll[i].second );
					TR.Trace << "inward_angle: " << F(9,3,ll[i].first) << I(4,ll[i].second) << endl;
				}
			}
			//pose.dump_pdb("test.pdb");
			//exit(0);
		}

		while ( true ) {
			string const worktag( string_of( hfline.index ) );
			Size const n( get_next_decoy_number_and_reserve_if_not_done( worktag, nstruct(), simfile ) );
			if ( n > nstruct() ) break;

			Pose pose( start_pose );
			Size const nres_protein( pose.total_residue()-1 );
			runtime_assert( nres_protein == nrepeat * repeatlen );

			/// setup symminfo?
			Sizes force_polar_positions;
			if ( polar_pore ) {
				Size const num_to_force( min( Sizes( option[ my_options::polar_pore ]() ) ) );
				force_polar_positions = inward_poslist;
				numeric::random::random_permutation( force_polar_positions, numeric::random::rg() );
				while ( force_polar_positions.size() > num_to_force ) force_polar_positions.pop_back();
			}

			clock_t starttime = clock();
			for ( Size m=1; m<= design_cycles+1; ++m ) {
				bool const skip_relax( m == 1 && !option[ my_options::fastrelax_before_design ] );

				if ( !skip_relax ) { // relax
					protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

					MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
					movemap->set_bb (true);
					movemap->set_chi(true);
					movemap->set_bb ( pose.total_residue(), false );
					movemap->set_chi( pose.total_residue(), false );
					fastrelax.set_movemap( movemap );
					cout << "checkpoint3" << endl; fflush( stdout );
					if ( !dry_run() ) fastrelax.apply( pose );
					cout << "checkpoint4" << endl; fflush( stdout );
				}

				if ( m>design_cycles ) break; // relax is cheap, so do one extra one at the end...

				{ // design
					pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
					task->initialize_from_command_line();
					task->or_include_current( true );
					if ( option[ my_options::skip_bump_check ] ) task->set_bump_check( false );

					bools is_protein( pose.total_residue(), false );
					for ( Size i=1; i<= nres_protein; ++i ) is_protein[i] = true;
					task->restrict_to_residues( is_protein );
					//if ( nodesign ) task->restrict_to_repacking();


					if ( force_polar_positions.size() ) {
						for ( Size i=1; i<= force_polar_positions.size(); ++i ) {
							Size const repeatpos( force_polar_positions[i] );
							string const name1s( "DEKRSTNQ" );
							bools allowed_aas( num_canonical_aas, false );
							for ( Size k=0; k<name1s.size(); ++k ) {
								AA const aa( aa_from_oneletter_code( name1s[k] ) );
								allowed_aas[aa] = true;
								TR.Trace << "force_sequence_positions: allowed_aa: " << I(4,repeatpos) << ' ' << aa << endl;
							}
							for ( Size k=0; k<nrepeat; ++k ) {
								Size const seqpos( k*repeatlen + repeatpos );
								TR.Trace << "force_sequence_positions: " << seqpos << ' ' << name1s << endl;
								task->nonconst_residue_task( seqpos ).restrict_absent_canonical_aas( allowed_aas );
							}
						}
					}

					if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
						protocols::task_operations::LimitAromaChi2Operation lp_op;
						lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negatives)
						lp_op.apply( pose, *task );
					}
					Size const nloop( 25 );
					ScoreFunctionOP design_scorefxn(0);
					if ( ( option[ my_options::use_softrep_for_early_design ] && m<= design_cycles/2 ) ||
							option[ my_options::use_softrep_for_design ] ) {
						design_scorefxn = ScoreFunctionFactory::create_score_function( "soft_rep_design.wts" );
						adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
					} else {
						design_scorefxn = fa_scorefxn; // already adjusted refwts
					}
					protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
					cout << "checkpoint5" << endl; fflush( stdout );
					if ( !dry_run() ) packmover.apply( pose );
					cout << "checkpoint6" << endl; fflush( stdout );
				}

			} // cycles
			clock_t stoptime = clock();
			Real const relax_simtime = ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

			Real const final_score = (*fa_scorefxn)( pose );

			///
			bools subset( pose.total_residue(), false );
			Size const base_repeat_offset( repeatlen*(base_repeat-1));
			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;

			Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa, buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score;
			compute_sasa_scores_for_subset( subset, pose, sasapack_score, norme_score, normsasa_score,
				exposed_polar_sasa, exposed_nonpolar_sasa, buried_polar_sasa, buried_nonpolar_sasa,
				buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score );



			Real const n2cdist = pose.residue(1).xyz("N").distance( pose.residue(chain_end(1,pose)).xyz("C") );

			/// compute psipred ss using single sequence
			string repeatsspred( "-" ), fullsspred;
			Real sspred_match(0.0), sspred_turn1(0.0), sspred_turn2( 0.0 ), sspred_helix1(0.0), sspred_helix2(0.0);

			{ // run psipred to re-predict the secondary structure
				vector1< Reals > pred_HEL;
				string const sequence( pose.sequence().substr(0,chain_end(1,pose) ) );
				//string const secstruct( pose.secstruct().substr(0,chain_end(1,pose) ) );
				string sspred;
				string const repeatss( pose.secstruct().substr(0,repeatlen) ); // current, not the target

				run_psipred( sequence, sspred, pred_HEL );

				if ( sspred.size() == sequence.size() ) { // success
					fullsspred = sspred;
					runtime_assert( sspred.size() == nrepeat * repeatlen );
					repeatsspred = sspred.substr( (base_repeat-1)*repeatlen, repeatlen );
					Size const turn1_begin( helix1_len+1 ), turn1_end( helix1_len + turn1.size() ),
						turn2_begin( helix1_len+turn1.size()+helix2_len+1 );
					for ( Size i=1; i<= sspred.size(); ++i ) {
						Size const rpos( (i-1)%repeatlen+1 );
						char const ss( repeatss[ rpos-1 ] );
						Real const prediction( ss == 'H' ? pred_HEL[i][1] : ( ss == 'E' ? pred_HEL[i][2] : pred_HEL[i][3] ) );
						sspred_match += prediction;
						if ( rpos < turn1_begin                     ) sspred_helix1 += prediction;
						if ( rpos > turn1_end && rpos < turn2_begin ) sspred_helix2 += prediction;
						if ( rpos >= turn1_begin-1 && rpos <= turn1_end+1 ) sspred_turn1 += prediction;
						if ( rpos == 1 || rpos >= turn2_begin-1           ) sspred_turn2 += prediction;
					}
					sspred_match  /= (  repeatlen * nrepeat );
					sspred_helix1 /= ( helix1_len * nrepeat );
					sspred_helix2 /= ( helix2_len * nrepeat );
					sspred_turn1  /= ( ( turn1.size()+2) * nrepeat );
					sspred_turn2  /= ( ( turn2.size()+2) * nrepeat );
				}
			}

			Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);
			if ( helix1_len ) {
				runtime_assert( helix1_len + turn1.size() + helix2_len + turn2.size() == repeatlen );

				compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1.size(), helix2_len,
					helix1_dist, helix1_twist, helix2_dist, helix2_twist );
			}

			Real const
				inner_helix_dist( min( helix1_dist, helix2_dist ) ),
				outer_helix_dist( max( helix1_dist, helix2_dist ) );
			Real const inner_helix_twist( helix1_dist < helix2_dist ? helix1_twist : helix2_twist );
			Real const outer_helix_twist( helix1_dist < helix2_dist ? helix2_twist : helix1_twist );

			Size const inner_helix_len( helix1_dist < helix2_dist ? helix1_len : helix2_len );
			Size const outer_helix_len( helix1_dist < helix2_dist ? helix2_len : helix1_len );
			string const inner_helix_turn( helix1_dist < helix2_dist ? turn1 : turn2 );
			string const outer_helix_turn( helix1_dist < helix2_dist ? turn2 : turn1 );

			string const buried_unsatisfied_string( get_buried_unsatisfied_string( pose ) );

			ostringstream relax_status; // historical reasons
			relax_status <<
				" sasapack_score: " << F(9,3,sasapack_score) <<
				" norme_score: " << F(9,3,norme_score) <<
				" normsasa_score: " << F(9,3,normsasa_score) <<
				" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
				" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
				" exposed_polar_fraction: " << F(9,3,exposed_polar_sasa/(exposed_polar_sasa+exposed_nonpolar_sasa))<<
				" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
				" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
				" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc)<<
				" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc)<<
				" packstat_score: " << F(9,3,packstat_score) <<
				" repeatsspred: " << repeatsspred <<
				" sspred_match: " << F(9,3,sspred_match) <<
				" sspred_helix1: " << F(9,3,sspred_helix1) <<
				" sspred_helix2: " << F(9,3,sspred_helix2) <<
				" sspred_turn1: " << F(9,3,sspred_turn1) <<
				" sspred_turn2: " << F(9,3,sspred_turn2) <<
				" helix1_dist: " << F(9,3,helix1_dist) <<
				" helix1_twist: " << F(9,3,helix1_twist) <<
				" helix2_dist: " << F(9,3,helix2_dist) <<
				" helix2_twist: " << F(9,3,helix2_twist) <<
				" inner_helix_dist: " << F(9,3,inner_helix_dist) <<
				" inner_helix_twist: " << F(9,3,inner_helix_twist) <<
				" outer_helix_dist: " << F(9,3,outer_helix_dist) <<
				" outer_helix_twist: " << F(9,3,outer_helix_twist) <<
				" inner_helix_len: " << I(4,inner_helix_len) <<
				" inner_helix_turn: " << inner_helix_turn <<
				" outer_helix_len: " << I(4,outer_helix_len) <<
				" outer_helix_turn: " << outer_helix_turn <<
				' ' << buried_unsatisfied_string <<
				' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() );


			/// RE-compute transform: twist, axis, after relax
			Real twist, rise, min_radius, max_radius, com_radius, handedness;

			compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius );
			{
				Vectors coords;
				for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
				handedness = get_chirality( repeatlen, coords );
			}

			string scorefiltertag( string_of( nrepeat )+ "." +
				string_of( inner_helix_len ) + "." + inner_helix_turn + "." +
				string_of( outer_helix_len ) + "." + outer_helix_turn + "." );
			if ( handedness < 0 ) scorefiltertag += "L";
			else scorefiltertag += "R";

			// string worktag( resample_bbs ? "resample_bb_"+resample_bb_tag : ( barrel_mode ? "barrel" : "tal" ) ),
			//  simfile( shared_output_tag()+"unbound_frag.work" );
			// worktag += "_" + string_of( repeatlen ) + "_" + string_of( nrepeat );
			// if ( force_sequence_positions && force_seq_relaxes.size() ) {
			//  Size count(0);
			//  for ( Size i=1; i<= repeatlen; ++i ) if ( has_element( force_seq_relaxes, i ) ) ++count;
			//  worktag += "_rlxseq"+string_of(count);
			// }

			//   bool const passed_centroid_score_filter
			//    ( append_score_to_scorefile_and_filter( worktag+"_cen", final_centroid_score,
			//                        centroid_score_filter_acceptance_rate,
			//                        centroid_score_filter_pass_early, simfile ) );

			bool const passed_score_filter( append_score_to_scorefile_and_filter( scorefiltertag, final_score,
				score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );


			string const repeatseq( pose.sequence().substr(0,repeatlen) );
			//Size const base_repeat_offset( repeatlen * ( base_repeat-1 ) );
			string const repeatbb = torsion2big_bin_string( base_repeat_offset+1, base_repeat_offset+repeatlen, pose );
			string const repeatss = pose.secstruct().substr(0,repeatlen);
			ostringstream refolding_status;

			/// try refolding de novo
			Size const n_refold( option[ my_options::n_refold ] ), n_refold_sspred( option[ my_options::n_refold_sspred ] );
			if ( n_refold > 0 && ( passed_score_filter || dry_run() ) ) {
				starttime = clock();
				for ( Size RR=1; RR<= 2; ++RR ) {
					if ( RR == 2 && fullsspred.size() != nrepeat * repeatlen ) continue;
					if ( RR==1 ) refolding_status << " n_refold: " << I(3,n_refold);
					else         refolding_status << " n_refold_sspred: " << I(3,n_refold_sspred);
					Size const rr_end( RR==1 ? n_refold : n_refold_sspred );
					for ( Size rr=1; rr<= rr_end; ++rr ) {
						Pose protpose;
						string const repss( RR == 1 ? repeatss : fullsspred );
						refold_repeat_pose( repeatseq, repss, nrepeat+2, 3, protpose );
						/// calc rmsd
						Real rmsd( 0.0 ), rmsd_single_repeat( 0.0 );
						{
							using namespace core::id;
							AtomID_Map< AtomID > atom_map;
							initialize_atomid_map( atom_map, protpose, id::GLOBAL_BOGUS_ATOM_ID );
							for ( Size i=1; i<= repeatlen; ++i ) {
								atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
									AtomID( pose.residue(i).atom_index("CA"),i);
							}
							rmsd_single_repeat = rmsd_by_mapping( protpose, pose, atom_map );
							Size const nres_for_rmsd( min( chain_end( 1, protpose), nrepeat * repeatlen ) );
							for ( Size i=1; i<= nres_for_rmsd; ++i ) {
								atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
									AtomID( pose.residue(i).atom_index("CA"),i);
							}
							rmsd = rmsd_by_mapping( protpose, pose, atom_map );
						}

						Real refold_handedness( 0.0 );
						{
							Vectors coords;
							for ( Size i=1; i<= chain_end( 1, protpose); ++i ) coords.push_back( protpose.residue(i).xyz("CA") );
							refold_handedness = get_chirality( repeatlen, coords );
						}
						Real tmptwist, tmprise, tmpmin_radius, tmpmax_radius, tmpcom_radius;
						compute_repeat_params( protpose, tmptwist, tmprise, tmpmin_radius, tmpmax_radius, tmpcom_radius );
						refolding_status << " R " << rr << F(9,3,rmsd) << F(9,3,rmsd_single_repeat) << F(9,3,refold_handedness) <<
							F(9,3,tmprise) << F(9,3,tmptwist);
					}

					stoptime = clock();
				}
				Real const refolding_simtime = ((double) stoptime - starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
				refolding_status << " refolding_simtime: " << F(9,3,refolding_simtime);
			} // refold scope


			string const outfilename( output_tag() + "recdesign_L"+string_of( hfline.index )+"_"+
				string_of( helix1_len ) + turn1 + string_of( helix2_len ) + turn2 + "_"+
				lead_zero_string_of( n,4 )+".pdb" );

			// Real const final_centroid_score(0), centroid_handedness(0), centroid_rise(0), centroid_radius(0),
			//  centroid_n2cdist(0),
			//  centroid_twist(0);
			bool const relax_this_decoy( true );

			cout << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
				" handedness: " << F(9,3,handedness) <<
				" min_radius: " << F(9,3,min_radius ) <<
				" rise: " << F(9,3,rise ) <<
				" twist: " << F(9,3,twist) <<
				" final_centroid_score: "<< F(9,3,0.0) <<
				" passed_score_filter: " << passed_score_filter <<
				" relaxed: " << relax_this_decoy <<
				" min_radius: " << F(9,3,min_radius )<<
				" max_radius: " << F(9,3,max_radius )<<
				" com_radius: " << F(9,3,com_radius )<<
				" reconstruct_rmsd: " << F(9,3,reconstruct_rmsd) <<
				// " centroid_handedness: " << F(9,3,0.0) <<
				" centroid_radius: " << F(9,3,0.0) <<
				" centroid_rise: " << F(9,3,0.0) <<
				" centroid_twist: " << F(9,3,0.0) <<
				" repeatseq: " << repeatseq <<
				" repeatbb: " << repeatbb <<
				" repeatss: " << repeatss <<
				" turn1: " << turn1 <<
				" turn2: " << turn2 <<
				" helix1_len: " << helix1_len <<
				" helix2_len: " << helix2_len <<
				" repeatlen: " << repeatlen <<
				" nrepeat: " << nrepeat <<
				" n2cdist: " << F(9,3,n2cdist) <<
				" centroid_n2cdist: " << F(9,3,0.0) <<
				" target_repeatbb: " << repeatbb <<
				" target_repeatss: " << repeatss <<
				" " << relax_status.str() <<
				" " << refolding_status.str() <<
				" centroid_simtime: " << F(9,3,0.0) <<
				" relax_simtime: " << F(9,3,relax_simtime) <<
				endl;

			if ( passed_score_filter ) {
				pose.dump_pdb( outfilename );
			}
			fflush( stdout );
			check_simtime();
		} // nstruct loop
		check_if_job_is_done();
	} // lines in the input file
	signal_that_job_is_done();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Window {
	Stub stub;
	Vector p1;
	Vector p2;
};

typedef vector1< Window > Windows;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
autodock_helix_frag_test()
{

	Size const nrepeat( option[ my_options::nrepeat ] );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	string const logfile( start_file() );

	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			hh_frags[i].nsim=0;
			//Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++hh_frags[i].nsim;
			}
			TR.Trace << "nsim " << hh_frags[i].nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	Sizes helixlens;
	strings turns;

	HelixFragLines hf_lines;

	{
		ifstream data( logfile.c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings const l( split_to_vector1( line ) );
			HelixFragLine hf;
			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
			}
			hf_lines.push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
			if ( !has_element( helixlens, hf.helix1_len ) ) helixlens.push_back( hf.helix1_len );
			if ( !has_element( helixlens, hf.helix2_len ) ) helixlens.push_back( hf.helix2_len );
		}
	}


	map< string, StubFrags > all_turn_frags;
	map< Size, StubFrags > all_helix_frags;

	{
		Size const max_frags( 1000000 ); // want all the helix frags!
		build_helix_library( helixlens, max_frags, all_helix_frags );
		build_turn_library( turns, all_turn_frags );
	}

	for ( Size li=1; li<= hf_lines.size(); ++li ) {
		HelixFragLine const & hfline( hf_lines[li] );

		StubFrag h1( get_stub_frag_for_id( hfline.helix1_tag, all_helix_frags.find( hfline.helix1_len )->second ) );
		StubFrag h2( get_stub_frag_for_id( hfline.helix2_tag, all_helix_frags.find( hfline.helix2_len )->second ) );
		StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
		StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

		Stub const start_stub,
			stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

		Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );
		string const turn1( hfline.turn1 ), turn2( hfline.turn2 );

		Vector center, n, t;
		Real theta;
		get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

		Vector const toroid_axis_center( center ), toroid_axis_vector( n );

		Real const shift( n.dot(t) );

		TR.Trace << "transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
			F(9,3,numeric::conversions::degrees( theta) ) << F(9,3,shift) << endl;

		//Real const theta_dev( theta / target_theta );

		// reconstruct coords
		Vectors coords;
		Stub stub( start_stub );
		Sizes seg_begins, seg_ends;
		//vector1< Stub > seg_stubs;
		for ( Size n=1; n<= nrepeat; ++n ) {
			for ( Size r=1; r<= 4; ++r ) {
				seg_begins.push_back( coords.size()+1 );
				//seg_stubs.push_back( stub );
				StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
				for ( Size i=1; i<= f.coords.size(); ++i ) {
					coords.push_back( stub.local2global( f.coords[i] ) );
				}
				stub = f.rt.make_jump( stub );
				seg_ends.push_back( coords.size() );
			}
		}
		Stub const final_stub( stub );
		Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
		Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
		//Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
		Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
		runtime_assert( coords.size() == nrepeat * repeatlen );
		runtime_assert( seg_begins.size() == nrepeat*4 );
		runtime_assert( seg_ends.size() == nrepeat*4 );

		/// compute the helix dats
		Size const nhelices( 2 * nrepeat );
		HelixDats helixdats( nhelices );
		for ( Size i=1; i<= nhelices; ++i ) {
			Size const seg_index( 2*i-1 );
			HelixDat hd;
			setup_helix_dat( seg_begins[ seg_index ], seg_ends[ seg_index ], coords, helixdats[ i ] );
		}

		// for ( Size i=1; i<= nhelices; ++i ) {
		//  for ( Size j=i+1; j<= nhelices; ++j ) {
		//   cout << "intra_monomer_helix_dist: " << i << ' ' << j <<
		//    F(9,3,sqrt( get_helix_helix_distance_squared( helixdats[i], helixdats[j] ) ) ) << endl;
		//  }
		// }

		/// for each window of outer helix (helix2) get a helix stub and compute local coords of axis points
		///

		Size const base_repeat( 3 ), base_repeat_offset( repeatlen * ( base_repeat-1 ) );
		Size const outer_helix_begin( base_repeat_offset + helix1_len + turn1.size() + 1 ),
			outer_helix_end( outer_helix_begin+helix2_len-1 );

		Size const helix_window(7);
		runtime_assert( helix2_len >= helix_window );
		Window w;
		runtime_assert( fabs( toroid_axis_vector.length()-1)<1e-3 );
		Windows windows;
		for ( Size i=outer_helix_begin; i<= outer_helix_end-helix_window+1; ++i ) { // windows of the outer helix
			w.stub = get_helix_stub( coords, i );
			w.p1 = w.stub.global2local( toroid_axis_center );
			w.p2 = w.stub.global2local( toroid_axis_center + toroid_axis_vector );
			windows.push_back( w );
		}

		Stub newstub;
		Vector new_center, new_axis;
		Vectors tmpcoords( coords ); // reused
		tmpcoords.resize( coords.size()*2 );
		vector1< std::pair< Size, Size > > helix_contact_pairs( nhelices*nhelices ); // reused

		Size pdbcounter(0);
		for ( Size iw=1; iw<= windows.size(); ++iw ) {
			Window const & w1( windows[iw] );
			for ( Size jw=1; jw<= windows.size(); ++jw ) {
				Window const & w2( windows[jw] );
				// try out some helix docking frags

				for ( Size k=1; k<= hh_frags.size(); ++k ) {
					if ( hh_frags[k].nsim < 25 ) continue; // HACKING
					for ( Size r=1; r<= 2; ++r ) {
						RT const & rt( r == 1 ? hh_frags[k].rt_12 : hh_frags[k].rt_21 );
						// what if we use this helix-frag to dock w1 onto w2 ? where will the new axis fall?
						rt.make_jump( w1.stub, newstub );
						new_center = newstub.local2global( w2.p1 );
						new_axis = newstub.local2global( w2.p2 ) - new_center;

						// let's look for planar solutions right now
						Real const axis_angle( numeric::conversions::degrees( std::acos( new_axis.dot( toroid_axis_vector ) ) ) );

						Real const max_axis_angle( 5 );
						if ( axis_angle < max_axis_angle || axis_angle > 180 - max_axis_angle ) {
							// how do we check for horrible clashes very quickly?


							// transform coords of partner
							// whats the transform?
							// xyz0 --> ( w2.stub.global2local ) --> xyz1 --> ( newstub.local2global ) -- > xyz2
							Vector zero(0,0,0), x( 1,0,0), y( 0,1,0), z(0,0,1);
							Vector const
								translation( newstub.local2global( w2.stub.global2local( zero ) ) ),
								R_x        ( newstub.local2global( w2.stub.global2local(    x ) ) - translation ),
								R_y        ( newstub.local2global( w2.stub.global2local(    y ) ) - translation ),
								R_z        ( newstub.local2global( w2.stub.global2local(    z ) ) - translation );
							Matrix const R( Matrix::cols( R_x, R_y, R_z ) );

							// compute helix-helix distances
							Size const nres_monomer( coords.size() );
							Real mindis2( 1e6 ), dis2;
							Real const clash_dis2( 4*4 ), contact_dis2( 12*12 ); //
							Size n_helix_contact_pairs(0);
							for ( Size i=1; i<= nhelices; ++i ) {
								HelixDat const & hd1( helixdats[i] );
								for ( Size j=1; j<= nhelices; ++j ) {
									HelixDat const & hd2( helixdats[j] );
									dis2 = get_helix_helix_distance_squared( hd1.axis, hd1.begin, hd1.end,
										R * hd2.axis,
										R * hd2.begin + translation,
										R * hd2.end   + translation );
									if ( false ) {
										Real const tmpdis2
											( get_helix_helix_distance_squared_numeric( hd1.begin, hd1.end,
											R * hd2.begin + translation,
											R * hd2.end   + translation ) );
										TR.Trace << "hh_dis2_err: " << F(9,3,fabs( sqrt(dis2)-sqrt(tmpdis2))) <<
											F(9,3,sqrt(dis2)) << F(9,3,sqrt(tmpdis2)) << std::endl;
										// write some status info
										TR.Trace << "hhdis2: " << F(9,3,dis2) << " li: " << li <<
											" h1: " << seg_begins[2*i-1] << ' ' << seg_ends[2*i-1] <<
											" h2: " << seg_begins[2*j-1]+nres_monomer << ' ' << seg_ends[2*j-1]+nres_monomer << endl;
									} /// DEBUGGING
									mindis2 = min( mindis2, dis2 );
									if ( mindis2 < clash_dis2 ) break;
									if ( dis2 < contact_dis2 ) {
										//
										++n_helix_contact_pairs;
										helix_contact_pairs[ n_helix_contact_pairs ].first = i;
										helix_contact_pairs[ n_helix_contact_pairs ].second = j;
									}
								}
								if ( mindis2 < clash_dis2 ) break;
							}

							if ( mindis2 < clash_dis2 ) continue; ////////////////////////////////////// helix-helix clash check


							// compute an angle between the vectors from the dock site to the toroid axes
							Vector
								proj1( toroid_axis_center + ( w1.stub.v - toroid_axis_center ).dot(toroid_axis_vector) * toroid_axis_vector),
								proj2( new_center + ( newstub.v - new_center ).dot( new_axis ) * new_axis ),
								v1( ( proj1 - w1.stub.v ).normalized() ),
								v2( ( proj2 - newstub.v ).normalized() );

							Real const proj_angle( numeric::conversions::degrees( std::acos( v1.dot( v2 ) ) ) );


							for ( Size jj=1; jj<= nres_monomer; ++jj ) {
								tmpcoords[ nres_monomer+jj ] = ( R*coords[jj]+translation );
							}


							// get sims for helix contact pairs
							ostringstream nsims_status;
							Size nsims_count(0);
							Size min_nsim( 10000 );
							//vector1< bool > found_contacts;
							for ( Size ii=1; ii<= n_helix_contact_pairs; ++ii ) {
								Size const i( helix_contact_pairs[ii].first ), j( helix_contact_pairs[ii].second ), iseg( 2*i-1 ),
									jseg( 2*j-1 );
								Size nsim;
								bool found_contact;
								count_similar_helix_helix_frags( seg_begins[ iseg ], seg_ends[ iseg ], tmpcoords,
									seg_begins[ jseg ]+nres_monomer, seg_ends[ jseg ]+nres_monomer, tmpcoords,
									hh_frags, nsim, found_contact );
								if ( found_contact ) {
									++nsims_count;
									nsims_status << ' ' << nsim;
									min_nsim = min( min_nsim, nsim );
								}
								//nsims.push_back( nsim );
								//found_contacts.push_back( found_contact );
							}


							++pdbcounter;
							string const outfilename( output_tag()+"docktmp_L"+string_of(li)+
								"_N"+lead_zero_string_of(pdbcounter,4)+".pdb");
							string leadtag( "axis_angle:" );
							if ( nsims_count > 2 && min_nsim > 10 ) {
								write_ca_pdbfile( tmpcoords, outfilename );
								leadtag = string("pdb_")+leadtag;
							}
							TR.Trace << leadtag << ' ' << F(9,3,axis_angle) << " proj_angle: " << F(9,3,proj_angle) <<
								" min_hh_dis: " << F(9,3,sqrt(mindis2) ) <<
								" hh_frag: " << k <<
								" w1: " << outer_helix_begin+iw-1 << ' ' << outer_helix_begin+iw-1+helix_window-1 <<
								" w2: " << nres_monomer+outer_helix_begin+jw-1 << ' ' << nres_monomer+outer_helix_begin+jw-1+helix_window-1<<
								" li: " << li << // which helixline are we on?
								" nsim_count: " << nsims_count <<
								" min_nsim: " << min_nsim <<
								" nsims: " << nsims_status.str() <<
								' ' << outfilename << endl; // 'whoah!'
						} // good axis_angle
					}
				}
			}
		}

	} // loop over lines


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
autodock_helix_frag_v2_test()
{
	Size const min_min_nsim_for_pdb_output( option[ my_options::min_min_nsim ] );
	Size const min_nsims_for_pdb_output( option[ my_options::min_nsims ] );
	// Size const min_min_nsim_for_pdb_output( 11 );
	// Size const min_nsims_for_pdb_output( 3 );
	Size const max_clashes_for_pdb_output( 1 );
	bool const force_symmetry( option[ my_options::force_symmetry ] );

	Size const nrepeat( option[ my_options::nrepeat ] );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	string const logfile( start_file() );


	Reals target_axis_angles( make_vector1( 0.0, pi/2 ) );
	strings target_axis_angle_types( make_vector1( string("P"), string("P90") ) );
	{
		strings const symm_types( make_vector1( string("O"), string("T"), string("I") ) );
		for ( strings::const_iterator st = symm_types.begin(); st != symm_types.end(); ++st ) {
			Vectors vertices, nbr_vertices;
			get_symm_type_vertices( (*st)[0], vertices, nbr_vertices );
			Real axis_angle( std::acos( vertices[1].normalized().dot( nbr_vertices[1].normalized() ) ) );
			if ( force_symmetry || (*st)=="O" ) {
				target_axis_angles.push_back( axis_angle );
				target_axis_angle_types.push_back( *st );
			}
			axis_angle *= 0.5;
			target_axis_angles.push_back( axis_angle );
			target_axis_angle_types.push_back( (*st)+"2" );
		}
	}
	// Reals target_axis_angles( make_vector1( 0.0, pi/2 ) );
	// strings target_axis_angle_types( make_vector1( string("P"), string("P90") ) );
	// {
	//  strings const symm_types( make_vector1( string("O"), string("T"), string("I") ) );
	//  for ( strings::const_iterator st = symm_types.begin(); st != symm_types.end(); ++st ) {
	//   Vectors vertices, nbr_vertices;
	//   get_symm_type_vertices( (*st)[0], vertices, nbr_vertices );
	//   Real const axis_angle( 0.5 * std::acos( vertices[1].normalized().dot( nbr_vertices[1].normalized() ) ) );
	//   target_axis_angles.push_back( axis_angle );
	//   target_axis_angle_types.push_back( (*st)+"2" );
	//   if ( *st == "O" ) {
	//    target_axis_angles.push_back( 2 * axis_angle );
	//    target_axis_angle_types.push_back( "O" );
	//   }
	//  }
	// }

	for ( Size i=1; i<= target_axis_angles.size(); ++i ) {
		cout << "target_axis_angles: " << i << ' ' << F(9,3,degrees( target_axis_angles[i] ) ) << ' ' <<
			target_axis_angle_types[i] << endl;
	}


	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			hh_frags[i].nsim=0;
			//Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++hh_frags[i].nsim;
			}
			TR.Trace << "nsim " << hh_frags[i].nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	//Sizes helixlens;
	strings turns;

	HelixFragLines hf_lines;

	{
		ifstream data( logfile.c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings const l( split_to_vector1( line ) );
			HelixFragLine hf;
			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
				else if ( l[i] == "inner_helix_params:" ) {
					hf.h1params.rise           = float_of( l[i+1] );
					hf.h1params.twist          = radians( float_of( l[i+2] ) );
					hf.h1params.tilt           = radians( float_of( l[i+3] ) );
					hf.h1params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h1params.ca_distance    = float_of( l[i+5] );
				} else if ( l[i] == "outer_helix_params:" ) {
					hf.h2params.rise           = float_of( l[i+1] );
					hf.h2params.twist          = radians( float_of( l[i+2] ) );
					hf.h2params.tilt           = radians( float_of( l[i+3] ) );
					hf.h2params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h2params.ca_distance    = float_of( l[i+5] );
				}
			}
			hf_lines.push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
		}
	}


	map< string, StubFrags > all_turn_frags;
	//map< Size, StubFrags > all_helix_frags;

	{
		build_turn_library_v2( turns, all_turn_frags );
	}

	for ( Size li=1; li<= hf_lines.size(); ++li ) {
		HelixFragLine const & hfline( hf_lines[li] );

		Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );
		string const turn1( hfline.turn1 ), turn2( hfline.turn2 );

		StubFrag h1, h2;
		helix_frag_from_helix_params_and_length( helix1_len, hfline.h1params, h1 );
		helix_frag_from_helix_params_and_length( helix2_len, hfline.h2params, h2 );
		StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
		StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

		Stub const start_stub,
			stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

		Vector center, n, t;
		Real toroid_axis_theta_radians;
		get_stub_transform_data( start_stub, stop_stub, center, n, t, toroid_axis_theta_radians );

		Vector toroid_axis_center( center ), toroid_axis_vector( n );

		Real const shift( n.dot(t) );

		TR.Trace << "transform: " << I(3,helix1_len) << I(3,helix2_len) << A(5,turn1) << A(5,turn2) <<
			F(9,3,numeric::conversions::degrees( toroid_axis_theta_radians) ) << F(9,3,shift) << endl;

		//Real const theta_dev( toroid_axis_theta_radians / target_theta );

		// reconstruct coords
		Vectors coords;
		Stub stub( start_stub );
		Sizes seg_begins, seg_ends;
		//vector1< Stub > seg_stubs;
		Vector centroid(0,0,0);
		for ( Size n=1; n<= nrepeat; ++n ) {
			for ( Size r=1; r<= 4; ++r ) {
				seg_begins.push_back( coords.size()+1 );
				//seg_stubs.push_back( stub );
				StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
				for ( Size i=1; i<= f.coords.size(); ++i ) {
					coords.push_back( stub.local2global( f.coords[i] ) );
					centroid += coords.back();
				}
				stub = f.rt.make_jump( stub );
				seg_ends.push_back( coords.size() );
			}
		}
		centroid /= coords.size();
		// shift the toroid_axis_center to center on the center of mass
		toroid_axis_center += toroid_axis_vector * ( centroid - toroid_axis_center ).dot( toroid_axis_vector );
		Stub const final_stub( stub );
		Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
		Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
		//Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
		Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
		runtime_assert( coords.size() == nrepeat * repeatlen );
		runtime_assert( seg_begins.size() == nrepeat*4 );
		runtime_assert( seg_ends.size() == nrepeat*4 );

		/// compute the helix dats
		Size const nhelices( 2 * nrepeat );
		HelixDats helixdats( nhelices );
		for ( Size i=1; i<= nhelices; ++i ) {
			Size const seg_index( 2*i-1 );
			HelixDat hd;
			setup_helix_dat( seg_begins[ seg_index ], seg_ends[ seg_index ], coords, helixdats[ i ] );
		}

		// for ( Size i=1; i<= nhelices; ++i ) {
		//  for ( Size j=i+1; j<= nhelices; ++j ) {
		//   cout << "intra_monomer_helix_dist: " << i << ' ' << j <<
		//    F(9,3,sqrt( get_helix_helix_distance_squared( helixdats[i], helixdats[j] ) ) ) << endl;
		//  }
		// }

		/// for each window of outer helix (helix2) get a helix stub and compute local coords of axis points
		///

		Size const base_repeat( 3 ), base_repeat_offset( repeatlen * ( base_repeat-1 ) );
		Size const outer_helix_begin( base_repeat_offset + helix1_len + turn1.size() + 1 ),
			outer_helix_end( outer_helix_begin+helix2_len-1 );

		Size const helix_window(7);
		runtime_assert( helix2_len >= helix_window );
		Window w;
		runtime_assert( fabs( toroid_axis_vector.length()-1)<1e-3 );
		Windows windows;
		for ( Size i=outer_helix_begin; i<= outer_helix_end-helix_window+1; ++i ) { // windows of the outer helix
			w.stub = get_helix_stub( coords, i );
			w.p1 = w.stub.global2local( toroid_axis_center );
			w.p2 = w.stub.global2local( toroid_axis_center + toroid_axis_vector );
			windows.push_back( w );
		}

		Stub newstub;
		Vector new_center, new_axis;
		Vectors tmpcoords( coords ); // reused
		tmpcoords.resize( coords.size()*2 );
		vector1< std::pair< Size, Size > > helix_contact_pairs( nhelices*nhelices ); // reused

		Size pdbcounter(0);
		for ( Size iw=1; iw<= windows.size(); ++iw ) {
			Window const & w1( windows[iw] );
			for ( Size jw=1; jw<= windows.size(); ++jw ) {
				Window const & w2( windows[jw] );
				// try out some helix docking frags

				for ( Size k=1; k<= hh_frags.size(); ++k ) {
					if ( hh_frags[k].nsim < 25 ) continue; // HACKING
					for ( Size r=1; r<= 1; ++r ) { /// NOT DOING THIS REPEAT=2 ANYMORE SINCE WE'RE AUTODOCKING
						RT const & rt( r == 1 ? hh_frags[k].rt_12 : hh_frags[k].rt_21 );
						// what if we use this helix-frag to dock w1 onto w2 ? where will the new axis fall?
						rt.make_jump( w1.stub, newstub );
						new_center = newstub.local2global( w2.p1 );
						new_axis = newstub.local2global( w2.p2 ) - new_center;

						// what would the optimal axis vectors be for the various symmetry types?
						// planar is easy: the original vector
						Vector const spin_axis( ( toroid_axis_vector.cross( new_center - toroid_axis_center ) ).normalized() );
						Real const axis_angle( numeric::conversions::degrees( std::acos( new_axis.dot( toroid_axis_vector ) ) ) );
						//if ( axis_angle > 90 ) axis_angle = 180-axis_angle; // hack with this...

						Real min_angle_dev( 1e6 );
						Size best_axis_angle_match(0);

						for ( Size i=1; i<= target_axis_angles.size(); ++i ) {
							for ( Size r=1; r<= 2; ++r ) {
								//if ( force_symmetry && r == 2  ) continue;
								Real const spin_angle( r == 1 ? target_axis_angles[i] : -1 * target_axis_angles[i] );
								Vector const target_axis_vector
									( numeric::rotation_matrix( spin_axis, spin_angle ) * toroid_axis_vector );

								Real angle_dev( numeric::conversions::degrees( std::acos( new_axis.dot( target_axis_vector ) ) ) );
								if ( !force_symmetry && angle_dev > 90 ) angle_dev = 180-angle_dev;

								if ( angle_dev < min_angle_dev ) {
									min_angle_dev = angle_dev;
									best_axis_angle_match = i;
								}
							}
						}

						Real const max_axis_angle_dev( 5 );
						if ( min_angle_dev < max_axis_angle_dev ) {

							// transform coords of partner
							// whats the transform?
							// xyz0 --> ( w2.stub.global2local ) --> xyz1 --> ( newstub.local2global ) -- > xyz2
							Vector zero(0,0,0), x( 1,0,0), y( 0,1,0), z(0,0,1);
							Vector const
								translation( newstub.local2global( w2.stub.global2local( zero ) ) ),
								R_x        ( newstub.local2global( w2.stub.global2local(    x ) ) - translation ),
								R_y        ( newstub.local2global( w2.stub.global2local(    y ) ) - translation ),
								R_z        ( newstub.local2global( w2.stub.global2local(    z ) ) - translation );
							Matrix const R( Matrix::cols( R_x, R_y, R_z ) );


							// figure out whether we are a symmetric homodimer:
							Real symdev_R_angle(1e6), symdev_translation_angle(1e6);
							int best_isym( 0 );
							{
								Vector const & old_center( toroid_axis_center ); // rename

								for ( int i=-1; i<= 1; ++i ) {
									// try pre-rotating about the toroid axis
									Real theta;
									Vector R_axis;
									if ( i == 0 ) {
										R_axis = numeric::rotation_axis( R, theta ).normalized();
									} else {
										R_axis = numeric::rotation_axis( R*rotation_matrix( toroid_axis_vector,
											i * toroid_axis_theta_radians ),
											theta ).normalized();
									}

									// for a symmetric homodimer, theta should be 180.0 degrees
									Real const symdev_R_angle_i( fabs( subtract_degree_angles( degrees( theta ), 180.0 ) ) );
									TR.Trace << "symdev_R_angle_i: " << i << F(9,3,symdev_R_angle_i) << endl;
									if ( symdev_R_angle_i < symdev_R_angle || symdev_R_angle_i < 5.0 ) {

										// check angle between R_axis and the toroid-axis-center translation vector
										// R_axis and the translation vector should be perpendicular
										Real const symdev_translation_angle_i =
											fabs( subtract_degree_angles( degrees( acos( R_axis.dot( ( new_center-old_center ).normalized()))),
											90.0 ) );

										if ( symdev_R_angle_i + symdev_translation_angle_i <
												symdev_R_angle + symdev_translation_angle ) {
											symdev_R_angle = symdev_R_angle_i;
											symdev_translation_angle = symdev_translation_angle_i;
											best_isym = i;
										}
									}
								}
							}

							if ( force_symmetry && ( symdev_R_angle > 3 || symdev_translation_angle > 5 ) ) continue;

							// how do we check for horrible clashes very quickly?




							// compute helix-helix distances
							Size const nres_monomer( coords.size() );
							Real mindis2( 1e6 ), dis2;
							Real const clash_dis2( 4*4 ), contact_dis2( 12*12 ); //
							Size n_helix_contact_pairs(0);
							for ( Size i=1; i<= nhelices; ++i ) {
								HelixDat const & hd1( helixdats[i] );
								for ( Size j=1; j<= nhelices; ++j ) {
									HelixDat const & hd2( helixdats[j] );
									dis2 = get_helix_helix_distance_squared( hd1.axis, hd1.begin, hd1.end,
										R * hd2.axis,
										R * hd2.begin + translation,
										R * hd2.end   + translation );
									if ( false ) {
										Real const tmpdis2
											( get_helix_helix_distance_squared_numeric( hd1.begin, hd1.end,
											R * hd2.begin + translation,
											R * hd2.end   + translation ) );
										TR.Trace << "hh_dis2_err: " << F(9,3,fabs( sqrt(dis2)-sqrt(tmpdis2))) <<
											F(9,3,sqrt(dis2)) << F(9,3,sqrt(tmpdis2)) << std::endl;
										// write some status info
										TR.Trace << "hhdis2: " << F(9,3,dis2) << " li: " << li <<
											" h1: " << seg_begins[2*i-1] << ' ' << seg_ends[2*i-1] <<
											" h2: " << seg_begins[2*j-1]+nres_monomer << ' ' << seg_ends[2*j-1]+nres_monomer << endl;
									} /// DEBUGGING
									mindis2 = min( mindis2, dis2 );
									if ( mindis2 < clash_dis2 ) break;
									if ( dis2 < contact_dis2 ) {
										//
										++n_helix_contact_pairs;
										helix_contact_pairs[ n_helix_contact_pairs ].first = i;
										helix_contact_pairs[ n_helix_contact_pairs ].second = j;
									}
								}
								if ( mindis2 < clash_dis2 ) break;
							}

							if ( mindis2 < clash_dis2 ) continue; ////////////////////////////////////// helix-helix clash check


							// compute an angle between the vectors from the dock site to the toroid axes
							// Vector
							//  proj1( toroid_axis_center + ( w1.stub.v - toroid_axis_center ).dot(toroid_axis_vector) * toroid_axis_vector),
							//  proj2( new_center + ( newstub.v - new_center ).dot( new_axis ) * new_axis ),
							//  v1( ( proj1 - w1.stub.v ).normalized() ),
							//  v2( ( proj2 - newstub.v ).normalized() );

							// Real const proj_angle( numeric::conversions::degrees( std::acos( v1.dot( v2 ) ) ) );


							for ( Size jj=1; jj<= nres_monomer; ++jj ) {
								tmpcoords[ nres_monomer+jj ] = ( R*coords[jj]+translation );
							}


							// get sims for helix contact pairs
							ostringstream nsims_status;
							Size nsims_count(0);
							Size min_nsim( 10000 );
							//vector1< bool > found_contacts;
							for ( Size ii=1; ii<= n_helix_contact_pairs; ++ii ) {
								Size const i( helix_contact_pairs[ii].first ), j( helix_contact_pairs[ii].second ), iseg( 2*i-1 ),
									jseg( 2*j-1 );
								Size nsim;
								bool found_contact;
								count_similar_helix_helix_frags( seg_begins[ iseg ], seg_ends[ iseg ], tmpcoords,
									seg_begins[ jseg ]+nres_monomer,
									seg_ends[ jseg ]+nres_monomer, tmpcoords,
									hh_frags, nsim, found_contact );
								TR.Trace << "found_contact: " << ii << ' ' << i << ' ' << j << ' ' << found_contact << endl;
								if ( found_contact ) {
									TR.Trace << "found_nsims: " << ii << ' ' << i << ' ' << j << ' ' << nsim << endl;
									++nsims_count;
									nsims_status << ' ' << nsim;
									min_nsim = min( min_nsim, nsim );
									if ( min_nsim < min_min_nsim_for_pdb_output ) break;
								}
								//nsims.push_back( nsim );
								//found_contacts.push_back( found_contact );
							}
							if ( min_nsim < min_min_nsim_for_pdb_output ) continue;

							if ( nsims_count < min_nsims_for_pdb_output ) continue;

							/// check for clashes based on C-alpha coords
							Size n_clash(0), n_close_contacts(0), n_long_contacts( 0 );

							{ // count Ca-Ca contacts between partners
								Size const base_repeat( 3 );
								Real const kill_dis2( 2*2 ), clash_dis2( 4*4 ), close_contact_dis2( 6.5*6.5 ),
									long_contact_dis2( 8*8 );
								Size const repeatlen( nres_monomer/nrepeat);
								Size const ibegin( (base_repeat-2)*repeatlen + 1 );
								Size const jbegin( (base_repeat-2)*repeatlen + 1 + nres_monomer );
								Size const iend( (base_repeat+1)*repeatlen );
								Size const jend( (base_repeat+1)*repeatlen + nres_monomer );
								for ( Size i=ibegin; i<= iend && n_clash<= max_clashes_for_pdb_output; ++i ) {
									for ( Size j=jbegin; j<= jend && n_clash<= max_clashes_for_pdb_output; ++j ) {
										dis2 = tmpcoords[i].distance_squared( tmpcoords[j] );
										if ( dis2 < kill_dis2 ) n_clash+=100;
										else if ( dis2 < clash_dis2 ) ++n_clash;
										else {
											if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
											if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
										}
									}
								}
							}

							if ( n_clash > max_clashes_for_pdb_output ) continue;


							++pdbcounter;
							string const outfilename( output_tag()+"docktmp_L"+string_of(li)+
								"_N"+lead_zero_string_of(pdbcounter,4)+".pdb");
							string leadtag( "axis_angle:" ); // SILLY AND STUPID SINCE WE ARE ALREADY FILTERING BY NSIMS ETC ABOVE
							if ( nsims_count >= min_nsims_for_pdb_output && min_nsim >= min_min_nsim_for_pdb_output ) {
								leadtag = string("pdb_")+leadtag;
								if ( option[ my_options::output_pdb_files ] ) write_ca_pdbfile( tmpcoords, outfilename );
							}
							cout << leadtag << ' ' << F(9,3,axis_angle) <<
								" min_angle_dev: " << F(9,3,min_angle_dev) <<
								" symm_type: " << target_axis_angle_types[ best_axis_angle_match ] <<
								" min_hh_dis: " << F(9,3,sqrt(mindis2) ) <<
								" hh_frag: " << k <<
								" w1: " << outer_helix_begin+iw-1 << ' ' << outer_helix_begin+iw-1+helix_window-1 <<
								" w2: " << nres_monomer+outer_helix_begin+jw-1 << ' ' << nres_monomer+outer_helix_begin+jw-1+helix_window-1<<
								" li: " << li << // which helixline are we on?
								" n_clash: " << n_clash <<
								" n_close_contacts: " << n_close_contacts <<
								" n_long_contacts: " << n_long_contacts <<
								" nsim_count: " << nsims_count <<
								" min_nsim: " << min_nsim <<
								" symdev_R_angle: " << F(9,3,symdev_R_angle) <<
								" symdev_translation_angle: " << F(9,3,symdev_translation_angle) <<
								" best_isym: " << best_isym <<
								" hh_frag_r: " << r <<
								" nsims: " << nsims_status.str() <<
								' ' << outfilename << endl; // 'whoah!'
						} // good axis_angle
					}
				}
			}
		}

	} // loop over lines


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
toroid1_rmsd_test()
{

	Size const nrepeat( 6 );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	string const logfile( start_file() );
	string const template_pdb( option[ my_options::template_pdb ] );

	Vectors template_coords;
	Size const nres_template( nrepeat * 35 - 2 );
	{
		Pose template_pose;
		pose_from_pdb( template_pose, template_pdb );
		runtime_assert( template_pose.sequence().substr(0,2) == "VS" );
		for ( Size i=1; i<= nres_template ; ++i ) {
			template_coords.push_back( template_pose.residue(i+2).xyz("CA") );
		}
	}


	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			hh_frags[i].nsim=0;
			//Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++hh_frags[i].nsim;
			}
			TR.Trace << "nsim " << hh_frags[i].nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	Sizes helixlens;
	strings turns;

	HelixFragLines hf_lines;

	{
		ifstream data( logfile.c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings const l( split_to_vector1( line ) );
			HelixFragLine hf;
			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];
			hf.line = line;

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
			}
			hf_lines.push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
			if ( !has_element( helixlens, hf.helix1_len ) ) helixlens.push_back( hf.helix1_len );
			if ( !has_element( helixlens, hf.helix2_len ) ) helixlens.push_back( hf.helix2_len );
		}
	}


	map< string, StubFrags > all_turn_frags;
	map< Size, StubFrags > all_helix_frags;

	{
		Size const max_frags( 1000000 ); // want all the helix frags!
		build_helix_library( helixlens, max_frags, all_helix_frags );
		build_turn_library( turns, all_turn_frags );
	}

	for ( Size li=1; li<= hf_lines.size(); ++li ) {
		HelixFragLine const & hfline( hf_lines[li] );

		StubFrag h1( get_stub_frag_for_id( hfline.helix1_tag, all_helix_frags.find( hfline.helix1_len )->second ) );
		StubFrag h2( get_stub_frag_for_id( hfline.helix2_tag, all_helix_frags.find( hfline.helix2_len )->second ) );
		StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
		StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

		Stub const start_stub,
			stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

		Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );
		string const turn1( hfline.turn1 ), turn2( hfline.turn2 );

		Vector center, n, t;
		Real theta;
		get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

		Vector const toroid_axis_center( center ), toroid_axis_vector( n );

		//Real const shift( n.dot(t) );

		//Real const theta_dev( theta / target_theta );

		// reconstruct coords
		Vectors coords;
		Stub stub( start_stub );
		Sizes seg_begins, seg_ends;
		//vector1< Stub > seg_stubs;
		for ( Size n=1; n<= nrepeat; ++n ) {
			for ( Size r=1; r<= 4; ++r ) {
				seg_begins.push_back( coords.size()+1 );
				//seg_stubs.push_back( stub );
				StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
				for ( Size i=1; i<= f.coords.size(); ++i ) {
					coords.push_back( stub.local2global( f.coords[i] ) );
				}
				stub = f.rt.make_jump( stub );
				seg_ends.push_back( coords.size() );
			}
		}
		Stub const final_stub( stub );
		Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
		Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
		//Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
		Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
		runtime_assert( coords.size() == nrepeat * repeatlen );
		runtime_assert( seg_begins.size() == nrepeat*4 );
		runtime_assert( seg_ends.size() == nrepeat*4 );


		/// compute rmsd to template pose
		coords.resize( nres_template );
		Real const rmsd( numeric::model_quality::calc_rms( coords, template_coords ) );

		cout << F(9,3,rmsd) << ' ' << hfline.line << endl;
	} // loop over lines


}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
toroid1_rmsd_v2_test()
{

	Size const nrepeat( 6 );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	string const logfile( start_file() );
	string const template_pdb( option[ my_options::template_pdb ] );

	Vectors template_coords;
	Size const nres_template( nrepeat * 35 - 2 );
	{
		Pose template_pose;
		pose_from_pdb( template_pose, template_pdb );
		runtime_assert( template_pose.sequence().substr(0,2) == "VS" );
		for ( Size i=1; i<= nres_template ; ++i ) {
			template_coords.push_back( template_pose.residue(i+2).xyz("CA") );
		}
	}


	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			hh_frags[i].nsim=0;
			//Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++hh_frags[i].nsim;
			}
			TR.Trace << "nsim " << hh_frags[i].nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	Sizes helixlens;
	strings turns;

	HelixFragLines hf_lines;

	{
		ifstream data( logfile.c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings const l( split_to_vector1( line ) );
			HelixFragLine hf;
			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];
			hf.line = line;

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
				else if ( l[i] == "inner_helix_params:" ) {
					hf.h1params.rise           = float_of( l[i+1] );
					hf.h1params.twist          = radians( float_of( l[i+2] ) );
					hf.h1params.tilt           = radians( float_of( l[i+3] ) );
					hf.h1params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h1params.ca_distance    = float_of( l[i+5] );
				} else if ( l[i] == "outer_helix_params:" ) {
					hf.h2params.rise           = float_of( l[i+1] );
					hf.h2params.twist          = radians( float_of( l[i+2] ) );
					hf.h2params.tilt           = radians( float_of( l[i+3] ) );
					hf.h2params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h2params.ca_distance    = float_of( l[i+5] );
				}
			}
			hf_lines.push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
			if ( !has_element( helixlens, hf.helix1_len ) ) helixlens.push_back( hf.helix1_len );
			if ( !has_element( helixlens, hf.helix2_len ) ) helixlens.push_back( hf.helix2_len );
		}
	}


	map< string, StubFrags > all_turn_frags;
	//map< Size, StubFrags > all_helix_frags;

	{
		//Size const max_frags( 1000000 ); // want all the helix frags!
		//build_helix_library( helixlens, max_frags, all_helix_frags );
		build_turn_library_v2( turns, all_turn_frags );
	}

	for ( Size li=1; li<= hf_lines.size(); ++li ) {
		HelixFragLine const & hfline( hf_lines[li] );

		Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );
		string const turn1( hfline.turn1 ), turn2( hfline.turn2 );
		StubFrag h1, h2;
		helix_frag_from_helix_params_and_length( helix1_len, hfline.h1params, h1 );
		helix_frag_from_helix_params_and_length( helix2_len, hfline.h2params, h2 );
		StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
		StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

		Stub const start_stub,
			stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));


		Vector center, n, t;
		Real theta;
		get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

		Vector const toroid_axis_center( center ), toroid_axis_vector( n );

		//Real const shift( n.dot(t) );

		//Real const theta_dev( theta / target_theta );

		// reconstruct coords
		Vectors coords;
		Stub stub( start_stub );
		Sizes seg_begins, seg_ends;
		//vector1< Stub > seg_stubs;
		for ( Size n=1; n<= nrepeat; ++n ) {
			for ( Size r=1; r<= 4; ++r ) {
				seg_begins.push_back( coords.size()+1 );
				//seg_stubs.push_back( stub );
				StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
				for ( Size i=1; i<= f.coords.size(); ++i ) {
					coords.push_back( stub.local2global( f.coords[i] ) );
				}
				stub = f.rt.make_jump( stub );
				seg_ends.push_back( coords.size() );
			}
		}
		Stub const final_stub( stub );
		Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
		Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
		//Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
		Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
		runtime_assert( coords.size() == nrepeat * repeatlen );
		runtime_assert( seg_begins.size() == nrepeat*4 );
		runtime_assert( seg_ends.size() == nrepeat*4 );


		/// compute rmsd to template pose
		coords.resize( nres_template );
		Real const rmsd( numeric::model_quality::calc_rms( coords, template_coords ) );

		cout << F(9,3,rmsd) << ' ' << hfline.line << endl;
	} // loop over lines


}

///////////////////////////////////////////////////////////////////////////////
void
setup_from_hfline(
	HelixFragLine const & hfline,
	map< Size, StubFrags > const & all_helix_frags,
	map< string, StubFrags > const & all_turn_frags,
	Vectors & coords,
	Windows & windows,
	Sizes & seg_begins,
	Sizes & seg_ends,
	HelixDats & helixdats,
	Vector & toroid_axis_vector,
	Vector & toroid_axis_center
)
{
	Size const nrepeat( option[ my_options::nrepeat ] );
	coords.clear(); windows.clear(); seg_begins.clear(); seg_ends.clear(); helixdats.clear();


	StubFrag h1( get_stub_frag_for_id( hfline.helix1_tag, all_helix_frags.find( hfline.helix1_len )->second ) );
	StubFrag h2( get_stub_frag_for_id( hfline.helix2_tag, all_helix_frags.find( hfline.helix2_len )->second ) );
	StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
	StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

	Stub const start_stub,
		stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

	Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );
	string const turn1( hfline.turn1 ), turn2( hfline.turn2 );

	Vector center, n, t;
	Real theta;
	get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

	toroid_axis_center =center;
	toroid_axis_vector = n;

	// reconstruct coords
	Stub stub( start_stub );
	for ( Size n=1; n<= nrepeat; ++n ) {
		for ( Size r=1; r<= 4; ++r ) {
			seg_begins.push_back( coords.size()+1 );
			StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
			for ( Size i=1; i<= f.coords.size(); ++i ) {
				coords.push_back( stub.local2global( f.coords[i] ) );
			}
			stub = f.rt.make_jump( stub );
			seg_ends.push_back( coords.size() );
		}
	}
	Stub const final_stub( stub );
	// Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
	// Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
	// Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
	Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
	runtime_assert( coords.size() == nrepeat * repeatlen );
	runtime_assert( seg_begins.size() == nrepeat*4 );
	runtime_assert( seg_ends.size() == nrepeat*4 );

	/// compute the helix dats
	Size const nhelices( 2 * nrepeat );
	helixdats.resize( nhelices );
	for ( Size i=1; i<= nhelices; ++i ) {
		Size const seg_index( 2*i-1 );
		HelixDat hd;
		setup_helix_dat( seg_begins[ seg_index ], seg_ends[ seg_index ], coords, helixdats[ i ] );
	}


	/// for each window of outer helix (helix2) get a helix stub and compute local coords of axis points
	///

	Size const base_repeat( 3 ), base_repeat_offset( repeatlen * ( base_repeat-1 ) );
	Size const outer_helix_begin( base_repeat_offset + helix1_len + turn1.size() + 1 ),
		outer_helix_end( outer_helix_begin+helix2_len-1 );

	Size const helix_window(7);
	runtime_assert( helix2_len >= helix_window );
	Window w;
	runtime_assert( fabs( toroid_axis_vector.length()-1)<1e-3 );
	//Windows windows;
	for ( Size i=outer_helix_begin; i<= outer_helix_end-helix_window+1; ++i ) { // windows of the outer helix
		w.stub = get_helix_stub( coords, i );
		w.p1 = w.stub.global2local( toroid_axis_center );
		w.p2 = w.stub.global2local( toroid_axis_center + toroid_axis_vector );
		windows.push_back( w );
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
setup_from_hfline_v2(
	HelixFragLine const & hfline,
	//map< Size, StubFrags > const & all_helix_frags,
	map< string, StubFrags > const & all_turn_frags,
	Vectors & coords,
	Windows & windows,
	Sizes & seg_begins,
	Sizes & seg_ends,
	HelixDats & helixdats,
	Vector & toroid_axis_vector,
	Vector & toroid_axis_center
)
{
	Size const nrepeat( option[ my_options::nrepeat ] );
	coords.clear(); windows.clear(); seg_begins.clear(); seg_ends.clear(); helixdats.clear();

	Size const helix1_len( hfline.helix1_len ), helix2_len( hfline.helix2_len );
	string const turn1( hfline.turn1 ), turn2( hfline.turn2 );

	StubFrag h1,h2;
	helix_frag_from_helix_params_and_length( helix1_len, hfline.h1params, h1 );
	helix_frag_from_helix_params_and_length( helix2_len, hfline.h2params, h2 );
	StubFrag t1( get_stub_frag_for_id( hfline.turn1_tag, all_turn_frags.find( hfline.turn1 )->second ) );
	StubFrag t2( get_stub_frag_for_id( hfline.turn2_tag, all_turn_frags.find( hfline.turn2 )->second ) );

	Stub const start_stub,
		stop_stub( t2.rt.make_jump( h2.rt.make_jump( t1.rt.make_jump( h1.rt.make_jump( start_stub )))));

	Vector center, n, t;
	Real theta;
	get_stub_transform_data( start_stub, stop_stub, center, n, t, theta );

	toroid_axis_center =center;
	toroid_axis_vector = n;

	// reconstruct coords
	Stub stub( start_stub );
	for ( Size n=1; n<= nrepeat; ++n ) {
		for ( Size r=1; r<= 4; ++r ) {
			seg_begins.push_back( coords.size()+1 );
			StubFrag const & f( r == 1 ? h1 : ( r == 2 ? t1 : ( r == 3 ? h2 : t2 ) ) );
			for ( Size i=1; i<= f.coords.size(); ++i ) {
				coords.push_back( stub.local2global( f.coords[i] ) );
			}
			stub = f.rt.make_jump( stub );
			seg_ends.push_back( coords.size() );
		}
	}
	Stub const final_stub( stub );
	// Vector const overlap_coord1( start_stub.local2global( h1.coords[1] ) );
	// Vector const overlap_coord2( final_stub.local2global( h1.coords[1] ) );
	// Real const n2c_distance( overlap_coord1.distance( overlap_coord2 ) );
	Size const repeatlen( (helix1_len) + turn1.size() + (helix2_len) + turn2.size() );
	runtime_assert( coords.size() == nrepeat * repeatlen );
	runtime_assert( seg_begins.size() == nrepeat*4 );
	runtime_assert( seg_ends.size() == nrepeat*4 );

	/// compute the helix dats
	Size const nhelices( 2 * nrepeat );
	helixdats.resize( nhelices );
	for ( Size i=1; i<= nhelices; ++i ) {
		Size const seg_index( 2*i-1 );
		HelixDat hd;
		setup_helix_dat( seg_begins[ seg_index ], seg_ends[ seg_index ], coords, helixdats[ i ] );
	}


	/// for each window of outer helix (helix2) get a helix stub and compute local coords of axis points
	///

	Size const base_repeat( 3 ), base_repeat_offset( repeatlen * ( base_repeat-1 ) );
	Size const outer_helix_begin( base_repeat_offset + helix1_len + turn1.size() + 1 ),
		outer_helix_end( outer_helix_begin+helix2_len-1 );

	Size const helix_window(7);
	runtime_assert( helix2_len >= helix_window );
	Window w;
	runtime_assert( fabs( toroid_axis_vector.length()-1)<1e-3 );
	//Windows windows;
	for ( Size i=outer_helix_begin; i<= outer_helix_end-helix_window+1; ++i ) { // windows of the outer helix
		w.stub = get_helix_stub( coords, i );
		w.p1 = w.stub.global2local( toroid_axis_center );
		w.p2 = w.stub.global2local( toroid_axis_center + toroid_axis_vector );
		windows.push_back( w );
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
dock_helix_frag_test()
{
	Size const min_min_nsim_for_pdb_output( 11 );
	Size const min_nsims_for_pdb_output( 3 );
	Size const max_clashes_for_pdb_output( 1 );


	strings const files( start_files() );
	runtime_assert( files.size() == 2);

	Size const nrepeat( option[ my_options::nrepeat ] );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	//string const logfile( start_file() );

	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			hh_frags[i].nsim=0;
			//Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++hh_frags[i].nsim;
			}
			//TR.Trace << "nsim " << hh_frags[i].nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	Sizes helixlens;
	strings turns;

	vector1< HelixFragLines > hf_lines(2);

	for ( Size fi=1; fi<= 2; ++fi ) {
		ifstream data( (files[fi]).c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings const l( split_to_vector1( line ) );
			HelixFragLine hf;
			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
			}
			hf_lines[fi].push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
			if ( !has_element( helixlens, hf.helix1_len ) ) helixlens.push_back( hf.helix1_len );
			if ( !has_element( helixlens, hf.helix2_len ) ) helixlens.push_back( hf.helix2_len );
		}
	}

	map< string, StubFrags > all_turn_frags;
	map< Size, StubFrags > all_helix_frags;

	{
		Size const max_frags( 1000000 ); // want all the helix frags!
		build_helix_library( helixlens, max_frags, all_helix_frags );
		build_turn_library( turns, all_turn_frags );
	}

	for ( Size li1=1; li1<= hf_lines[1].size(); ++li1 ) {
		Vectors coords1;
		HelixDats helixdats1;
		Windows windows1;
		Sizes seg_begins1, seg_ends1;
		Vector toroid_axis_vector1, toroid_axis_center1;
		setup_from_hfline( hf_lines[1][ li1 ], all_helix_frags, all_turn_frags,
			coords1, windows1, seg_begins1, seg_ends1, helixdats1,
			toroid_axis_vector1, toroid_axis_center1 );
		Size const nres_monomer1( coords1.size() );

		Vectors tmpcoords( coords1 ); // reused
		Size const nhelices( 2*nrepeat );
		vector1< std::pair< Size, Size > > helix_contact_pairs( nhelices*nhelices ); // reused

		for ( Size li2=1; li2<= hf_lines[2].size(); ++li2 ) {
			Vectors coords2;
			HelixDats helixdats2;
			Windows windows2;
			Sizes seg_begins2, seg_ends2;
			Vector toroid_axis_vector2, toroid_axis_center2;
			setup_from_hfline( hf_lines[2][ li2 ], all_helix_frags, all_turn_frags,
				coords2, windows2, seg_begins2, seg_ends2, helixdats2,
				toroid_axis_vector2, toroid_axis_center2 );
			Size const nres_monomer2( coords2.size() );

			tmpcoords.resize( nres_monomer1 + nres_monomer2 );

			Size pdbcounter(0);
			Stub newstub;
			Vector new_center, new_axis;
			for ( Size iw=1; iw<= windows1.size(); ++iw ) {
				Window const & w1( windows1[iw] );
				for ( Size jw=1; jw<= windows2.size(); ++jw ) {
					Window const & w2( windows2[jw] );
					// try out some helix docking frags

					for ( Size k=1; k<= hh_frags.size(); ++k ) {
						if ( hh_frags[k].nsim < 25 ) continue; // HACKING

						for ( Size r=1; r<= 2; ++r ) {
							RT const & rt( r == 1 ? hh_frags[k].rt_12 : hh_frags[k].rt_21 );
							// what if we use this helix-frag to dock w1 onto w2 ? where will the new axis fall?
							rt.make_jump( w1.stub, newstub );
							new_center = newstub.local2global( w2.p1 );
							new_axis = newstub.local2global( w2.p2 ) - new_center;

							// let's look for planar solutions right now
							Real const axis_angle( numeric::conversions::degrees( std::acos( new_axis.dot( toroid_axis_vector1 ) ) ) );

							Real const max_axis_angle( 5 );
							if ( axis_angle < max_axis_angle || axis_angle > 180 - max_axis_angle ) {
								// how do we check for horrible clashes very quickly?


								// transform coords of partner
								// whats the transform?
								// xyz0 --> ( w2.stub.global2local ) --> xyz1 --> ( newstub.local2global ) -- > xyz2
								Vector zero(0,0,0), x( 1,0,0), y( 0,1,0), z(0,0,1);
								Vector const
									translation( newstub.local2global( w2.stub.global2local( zero ) ) ),
									R_x        ( newstub.local2global( w2.stub.global2local(    x ) ) - translation ),
									R_y        ( newstub.local2global( w2.stub.global2local(    y ) ) - translation ),
									R_z        ( newstub.local2global( w2.stub.global2local(    z ) ) - translation );
								Matrix const R( Matrix::cols( R_x, R_y, R_z ) );

								// compute helix-helix distances
								//Size const nres_monomer1( coords1.size() );
								Real mindis2( 1e6 ), dis2;
								Real const clash_dis2( 4*4 ), contact_dis2( 12*12 ); //
								Size n_helix_contact_pairs(0);
								for ( Size i=1; i<= nhelices; ++i ) {
									HelixDat const & hd1( helixdats1[i] );
									for ( Size j=1; j<= nhelices; ++j ) {
										HelixDat const & hd2( helixdats2[j] );
										dis2 = get_helix_helix_distance_squared( hd1.axis, hd1.begin, hd1.end,
											R * hd2.axis,
											R * hd2.begin + translation,
											R * hd2.end   + translation );
										mindis2 = min( mindis2, dis2 );
										if ( mindis2 < clash_dis2 ) break;
										if ( dis2 < contact_dis2 ) {
											//
											++n_helix_contact_pairs;
											helix_contact_pairs[ n_helix_contact_pairs ].first = i;
											helix_contact_pairs[ n_helix_contact_pairs ].second = j;
										}
									}
									if ( mindis2 < clash_dis2 ) break;
								}

								if ( mindis2 < clash_dis2 ) continue; ////////////////////////////////////// helix-helix clash check


								// reconstruct full coords
								for ( Size jj=1; jj<= nres_monomer2; ++jj ) {
									tmpcoords[ nres_monomer1+jj ] = ( R*coords2[jj]+translation );
								}


								// get sims for helix contact pairs
								ostringstream nsims_status;
								Size nsims_count(0);
								Size min_nsim( 10000 );
								//vector1< bool > found_contacts;
								for ( Size ii=1; ii<= n_helix_contact_pairs; ++ii ) {
									Size const i( helix_contact_pairs[ii].first ), j( helix_contact_pairs[ii].second ), iseg( 2*i-1 ),
										jseg( 2*j-1 );
									Size nsim;
									bool found_contact;
									count_similar_helix_helix_frags( seg_begins1[ iseg ], seg_ends1[ iseg ], tmpcoords,
										seg_begins2[ jseg ]+nres_monomer1, seg_ends2[ jseg ]+nres_monomer1,
										tmpcoords,
										hh_frags, nsim, found_contact );
									if ( found_contact ) {
										++nsims_count;
										nsims_status << ' ' << nsim;
										min_nsim = min( min_nsim, nsim );
										if ( min_nsim < min_min_nsim_for_pdb_output ) break;
									}
								}

								if ( min_nsim < min_min_nsim_for_pdb_output ) continue;

								if ( nsims_count < min_nsims_for_pdb_output ) continue;


								Size n_clash(0), n_close_contacts(0), n_long_contacts( 0 );

								{ // count Ca-Ca contacts between partners
									Size const base_repeat( 3 );
									Real const kill_dis2( 2*2 ), clash_dis2( 4*4 ), close_contact_dis2( 6.5*6.5 ),
										long_contact_dis2( 8*8 );
									Size const repeatlen1( nres_monomer1/nrepeat), repeatlen2( nres_monomer2/nrepeat);
									Size const ibegin( (base_repeat-2)*repeatlen1 + 1 );
									Size const jbegin( (base_repeat-2)*repeatlen2 + 1 + nres_monomer1 );
									Size const iend( (base_repeat+1)*repeatlen1 );
									Size const jend( (base_repeat+1)*repeatlen2 + nres_monomer1 );
									for ( Size i=ibegin; i<= iend && n_clash<= max_clashes_for_pdb_output; ++i ) {
										for ( Size j=jbegin; j<= jend && n_clash<= max_clashes_for_pdb_output; ++j ) {
											dis2 = tmpcoords[i].distance_squared( tmpcoords[j] );
											if ( dis2 < kill_dis2 ) n_clash+=100;
											else if ( dis2 < clash_dis2 ) ++n_clash;
											else {
												if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
												if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
											}
										}
									}
								}

								if ( n_clash > max_clashes_for_pdb_output ) continue;

								++pdbcounter;
								string const outfilename( output_tag()+"docktmp_L"+string_of(li1)+"_"+string_of(li2)+
									"_N"+lead_zero_string_of(pdbcounter,4)+".pdb");
								string leadtag( "axis_angle:" );
								if ( nsims_count >= min_nsims_for_pdb_output && min_nsim >= min_min_nsim_for_pdb_output ) {
									write_ca_pdbfile( tmpcoords, outfilename );
									leadtag = string("pdb_")+leadtag;
								}
								TR.Trace << leadtag << ' ' << F(9,3,axis_angle) <<
									" min_hh_dis: " << F(9,3,sqrt(mindis2) ) <<
									" hh_frag: " << k <<
									" li1: " << li1 <<
									" li2: " << li2 <<
									" n_clash: " << n_clash <<
									" n_close_contacts: " << n_close_contacts <<
									" n_long_contacts: " << n_long_contacts <<
									" nsim_count: " << nsims_count <<
									" min_nsim: " << min_nsim <<
									" nsims: " << nsims_status.str() <<
									' ' << outfilename << endl; // 'whoah!'
								fflush( stdout );
							} // good axis_angle
						} // use hh-frag stubs for both directions
					} // loop over hh-frags
				} // jw=1,windows2.size()
			} // iw=1,windows1.size()
		} // li2
	} // li1

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
dock_helix_frag_v2_test()
{
	Size const min_min_nsim_for_pdb_output( option[ my_options::min_min_nsim ] );
	Size const min_nsims_for_pdb_output( option[ my_options::min_nsims ] );
	Size const max_clashes_for_pdb_output( 1 );

	Reals target_axis_angles( make_vector1( 0.0, pi/2 ) );
	strings target_axis_angle_types( make_vector1( string("P"), string("P90") ) );
	{
		strings const symm_types( make_vector1( string("O"), string("T"), string("I") ) );
		for ( strings::const_iterator st = symm_types.begin(); st != symm_types.end(); ++st ) {
			Vectors vertices, nbr_vertices;
			get_symm_type_vertices( (*st)[0], vertices, nbr_vertices );
			Real const axis_angle( 0.5 * std::acos( vertices[1].normalized().dot( nbr_vertices[1].normalized() ) ) );
			target_axis_angles.push_back( axis_angle );
			target_axis_angle_types.push_back( (*st)+"2" );
			if ( *st == "O" ) {
				target_axis_angles.push_back( 2 * axis_angle );
				target_axis_angle_types.push_back( "O" );
			}
		}
	}

	for ( Size i=1; i<= target_axis_angles.size(); ++i ) {
		cout << "target_axis_angles: " << i << ' ' << F(9,3,degrees( target_axis_angles[i] ) ) << ' ' <<
			target_axis_angle_types[i] << endl;
	}


	strings const files( start_files() );
	runtime_assert( files.size() == 2);

	Size const nrepeat( option[ my_options::nrepeat ] );
	//Real const target_theta( (2*numeric::constants::d::pi)/nrepeat );
	//string const logfile( start_file() );

	Real const max_hh_frag_dis( 3.0 );
	HelixHelixFrags hh_frags;
	{ // hacking
		build_helix_transform_library( hh_frags );

		for ( Size i=1; i<= hh_frags.size(); ++i ) {
			hh_frags[i].nsim=0;
			//Size nsim( 0 );
			for ( Size j=1; j<= hh_frags.size(); ++j ) {
				if ( j==i ) continue;
				Real const dis( helix_helix_frag_distance( hh_frags[i], hh_frags[j] ) );
				if ( dis <= max_hh_frag_dis ) ++hh_frags[i].nsim;
			}
			//TR.Trace << "nsim " << hh_frags[i].nsim << ' ' << i << ' ' << hh_frags.size() << endl;
		}
	}


	//Sizes helixlens;
	strings turns;

	vector1< HelixFragLines > hf_lines(2);

	for ( Size fi=1; fi<= 2; ++fi ) {
		ifstream data( (files[fi]).c_str() );
		string line;
		while ( getline( data, line ) ) {
			if ( line.find( "pdb_transform:" ) == string::npos ) continue;
			strings const l( split_to_vector1( line ) );
			HelixFragLine hf;
			hf.helix1_len = int_of( l[2] );
			hf.helix2_len = int_of( l[3] );
			hf.turn1 = l[4];
			hf.turn2 = l[5];

			for ( Size i=1; i< l.size(); ++i ) {
				if      ( l[i] == "inner_helix_tag:" ) hf.helix1_tag = l[i+1];
				else if ( l[i] == "outer_helix_tag:" ) hf.helix2_tag = l[i+1];
				else if ( l[i] == "inner_turn_tag:" ) hf.turn1_tag = l[i+1];
				else if ( l[i] == "outer_turn_tag:" ) hf.turn2_tag = l[i+1];
				else if ( l[i] == "inner_helix_params:" ) {
					hf.h1params.rise           = float_of( l[i+1] );
					hf.h1params.twist          = radians( float_of( l[i+2] ) );
					hf.h1params.tilt           = radians( float_of( l[i+3] ) );
					hf.h1params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h1params.ca_distance    = float_of( l[i+5] );
				} else if ( l[i] == "outer_helix_params:" ) {
					hf.h2params.rise           = float_of( l[i+1] );
					hf.h2params.twist          = radians( float_of( l[i+2] ) );
					hf.h2params.tilt           = radians( float_of( l[i+3] ) );
					hf.h2params.tilt_direction = radians( float_of( l[i+4] ) );
					hf.h2params.ca_distance    = float_of( l[i+5] );
				}
			}
			hf_lines[fi].push_back( hf );
			if ( !has_element( turns, hf.turn1 ) ) turns.push_back( hf.turn1 );
			if ( !has_element( turns, hf.turn2 ) ) turns.push_back( hf.turn2 );
		}
	}

	map< string, StubFrags > all_turn_frags;
	//map< Size, StubFrags > all_helix_frags;

	{
		//Size const max_frags( 1000000 ); // want all the helix frags!
		//build_helix_library( helixlens, max_frags, all_helix_frags );
		build_turn_library_v2( turns, all_turn_frags );
	}

	for ( Size li1=1; li1<= hf_lines[1].size(); ++li1 ) {
		Vectors coords1;
		HelixDats helixdats1;
		Windows windows1;
		Sizes seg_begins1, seg_ends1;
		Vector toroid_axis_vector1, toroid_axis_center1;
		setup_from_hfline_v2( hf_lines[1][ li1 ], all_turn_frags, coords1, windows1, seg_begins1, seg_ends1, helixdats1,
			toroid_axis_vector1, toroid_axis_center1 );
		Size const nres_monomer1( coords1.size() );

		Vectors tmpcoords( coords1 ); // reused
		Size const nhelices( 2*nrepeat );
		vector1< std::pair< Size, Size > > helix_contact_pairs( nhelices*nhelices ); // reused

		for ( Size li2=1; li2<= hf_lines[2].size(); ++li2 ) {
			Vectors coords2;
			HelixDats helixdats2;
			Windows windows2;
			Sizes seg_begins2, seg_ends2;
			Vector toroid_axis_vector2, toroid_axis_center2;
			setup_from_hfline_v2( hf_lines[2][ li2 ], all_turn_frags, coords2, windows2, seg_begins2, seg_ends2, helixdats2,
				toroid_axis_vector2, toroid_axis_center2 );
			Size const nres_monomer2( coords2.size() );

			tmpcoords.resize( nres_monomer1 + nres_monomer2 );

			Size pdbcounter(0);
			Stub newstub;
			Vector new_center, new_axis;
			for ( Size iw=1; iw<= windows1.size(); ++iw ) {
				Window const & w1( windows1[iw] );
				for ( Size jw=1; jw<= windows2.size(); ++jw ) {
					Window const & w2( windows2[jw] );
					// try out some helix docking frags

					for ( Size k=1; k<= hh_frags.size(); ++k ) {
						if ( hh_frags[k].nsim < 25 ) continue; // HACKING

						for ( Size r=1; r<= 2; ++r ) {
							RT const & rt( r == 1 ? hh_frags[k].rt_12 : hh_frags[k].rt_21 );
							// what if we use this helix-frag to dock w1 onto w2 ? where will the new axis fall?
							rt.make_jump( w1.stub, newstub );
							new_center = newstub.local2global( w2.p1 );
							new_axis = newstub.local2global( w2.p2 ) - new_center;


							// what would the optimal axis vectors be for the various symmetry types?
							// planar is easy: the original vector
							Vector const spin_axis( ( toroid_axis_vector1.cross( new_center - toroid_axis_center1 ) ).normalized() );
							Real axis_angle( numeric::conversions::degrees( std::acos( new_axis.dot( toroid_axis_vector1 ) ) ) );
							if ( axis_angle > 90 ) axis_angle = 180-axis_angle;

							Real min_angle_dev( 1e6 );
							Size best_axis_angle_match(0);

							for ( Size i=1; i<= target_axis_angles.size(); ++i ) {
								for ( Size r=1; r<= 2; ++r ) {
									Real const spin_angle( r == 1 ? target_axis_angles[i] : -1 * target_axis_angles[i] );
									Vector const target_axis_vector
										( numeric::rotation_matrix( spin_axis, spin_angle ) * toroid_axis_vector1 );

									Real angle_dev( numeric::conversions::degrees( std::acos( new_axis.dot( target_axis_vector ) ) ) );
									if ( angle_dev > 90 ) angle_dev = 180-angle_dev;

									if ( angle_dev < min_angle_dev ) {
										min_angle_dev = angle_dev;
										best_axis_angle_match = i;
									}
								}
							}

							Real const max_axis_angle_dev( 5 );
							if ( min_angle_dev < max_axis_angle_dev ) {



								// transform coords of partner
								// whats the transform?
								// xyz0 --> ( w2.stub.global2local ) --> xyz1 --> ( newstub.local2global ) -- > xyz2
								Vector zero(0,0,0), x( 1,0,0), y( 0,1,0), z(0,0,1);
								Vector const
									translation( newstub.local2global( w2.stub.global2local( zero ) ) ),
									R_x        ( newstub.local2global( w2.stub.global2local(    x ) ) - translation ),
									R_y        ( newstub.local2global( w2.stub.global2local(    y ) ) - translation ),
									R_z        ( newstub.local2global( w2.stub.global2local(    z ) ) - translation );
								Matrix const R( Matrix::cols( R_x, R_y, R_z ) );

								// compute helix-helix distances
								//Size const nres_monomer1( coords1.size() );
								Real mindis2( 1e6 ), dis2;
								Real const clash_dis2( 4*4 ), contact_dis2( 12*12 ); //
								Size n_helix_contact_pairs(0);
								for ( Size i=1; i<= nhelices; ++i ) {
									HelixDat const & hd1( helixdats1[i] );
									for ( Size j=1; j<= nhelices; ++j ) {
										HelixDat const & hd2( helixdats2[j] );
										dis2 = get_helix_helix_distance_squared( hd1.axis, hd1.begin, hd1.end,
											R * hd2.axis,
											R * hd2.begin + translation,
											R * hd2.end   + translation );
										mindis2 = min( mindis2, dis2 );
										if ( mindis2 < clash_dis2 ) break;
										if ( dis2 < contact_dis2 ) {
											//
											++n_helix_contact_pairs;
											helix_contact_pairs[ n_helix_contact_pairs ].first = i;
											helix_contact_pairs[ n_helix_contact_pairs ].second = j;
										}
									}
									if ( mindis2 < clash_dis2 ) break;
								}

								if ( mindis2 < clash_dis2 ) continue; ////////////////////////////////////// helix-helix clash check


								// reconstruct full coords
								for ( Size jj=1; jj<= nres_monomer2; ++jj ) {
									tmpcoords[ nres_monomer1+jj ] = ( R*coords2[jj]+translation );
								}


								// get sims for helix contact pairs
								ostringstream nsims_status;
								Size nsims_count(0);
								Size min_nsim( 10000 );
								//vector1< bool > found_contacts;
								for ( Size ii=1; ii<= n_helix_contact_pairs; ++ii ) {
									Size const i( helix_contact_pairs[ii].first ), j( helix_contact_pairs[ii].second ), iseg( 2*i-1 ),
										jseg( 2*j-1 );
									Size nsim;
									bool found_contact;
									count_similar_helix_helix_frags( seg_begins1[ iseg ], seg_ends1[ iseg ], tmpcoords,
										seg_begins2[ jseg ]+nres_monomer1, seg_ends2[ jseg ]+nres_monomer1,
										tmpcoords,
										hh_frags, nsim, found_contact );
									if ( found_contact ) {
										++nsims_count;
										nsims_status << ' ' << nsim;
										min_nsim = min( min_nsim, nsim );
										if ( min_nsim < min_min_nsim_for_pdb_output ) break;
									}
								}

								if ( min_nsim < min_min_nsim_for_pdb_output ) continue;

								if ( nsims_count < min_nsims_for_pdb_output ) continue;


								Size n_clash(0), n_close_contacts(0), n_long_contacts( 0 );

								{ // count Ca-Ca contacts between partners
									Size const base_repeat( 3 );
									Real const kill_dis2( 2*2 ), clash_dis2( 4*4 ), close_contact_dis2( 6.5*6.5 ),
										long_contact_dis2( 8*8 );
									Size const repeatlen1( nres_monomer1/nrepeat), repeatlen2( nres_monomer2/nrepeat);
									Size const ibegin( (base_repeat-2)*repeatlen1 + 1 );
									Size const jbegin( (base_repeat-2)*repeatlen2 + 1 + nres_monomer1 );
									Size const iend( (base_repeat+1)*repeatlen1 );
									Size const jend( (base_repeat+1)*repeatlen2 + nres_monomer1 );
									for ( Size i=ibegin; i<= iend && n_clash<= max_clashes_for_pdb_output; ++i ) {
										for ( Size j=jbegin; j<= jend && n_clash<= max_clashes_for_pdb_output; ++j ) {
											dis2 = tmpcoords[i].distance_squared( tmpcoords[j] );
											if ( dis2 < kill_dis2 ) n_clash+=100;
											else if ( dis2 < clash_dis2 ) ++n_clash;
											else {
												if ( dis2 < close_contact_dis2 ) ++n_close_contacts;
												if ( dis2 < long_contact_dis2 ) ++n_long_contacts;
											}
										}
									}
								}

								if ( n_clash > max_clashes_for_pdb_output ) continue;

								++pdbcounter;
								string const outfilename( output_tag()+"docktmp_L"+string_of(li1)+"_"+string_of(li2)+
									"_N"+lead_zero_string_of(pdbcounter,4)+".pdb");
								string leadtag( "axis_angle:" );
								if ( nsims_count >= min_nsims_for_pdb_output && min_nsim >= min_min_nsim_for_pdb_output ) {
									leadtag = string("pdb_")+leadtag;
									if ( option[ my_options::output_pdb_files ] ) write_ca_pdbfile( tmpcoords, outfilename );
								}
								cout << leadtag << ' ' << F(9,3,axis_angle) <<
									" min_angle_dev: " << F(9,3,min_angle_dev) <<
									" symm_type: " << target_axis_angle_types[ best_axis_angle_match ] <<
									" min_hh_dis: " << F(9,3,sqrt(mindis2) ) <<
									" hh_frag: " << k <<
									" li1: " << li1 <<
									" li2: " << li2 <<
									" n_clash: " << n_clash <<
									" n_close_contacts: " << n_close_contacts <<
									" n_long_contacts: " << n_long_contacts <<
									" nsim_count: " << nsims_count <<
									" min_nsim: " << min_nsim <<
									" nsims: " << nsims_status.str() <<
									' ' << outfilename << endl; // 'whoah!'
								fflush( stdout );
							} // good axis_angle
						} // use hh-frag stubs for both directions
					} // loop over hh-frags
				} // jw=1,windows2.size()
			} // iw=1,windows1.size()
		} // li2
	} // li1

}

///////////////////////////////////////////////////////////////////////////////

void
helical_params_test()
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;
	using namespace core::optimization;

	HelixParams hparams;

	hparams.rise = 1.5;
	hparams.twist = radians( 99.0 );
	hparams.tilt = 0;
	hparams.tilt_direction = 0;

	hparams.ca_distance = 1.5; // guess

	Vectors coords;
	Stub start_stub;
	generate_helix_coords( start_stub, hparams, 20, coords );

	write_ca_pdbfile( coords, "helix20.pdb" );

	Multivec params( 3 );
	params[1] = 1.5;
	params[2] = radians( 99.0 );
	params[3] = 2.0; // radius, guess


	Vectors target_coords;
	{
		Pose pose;
		pose_from_pdb( pose, "./input/perfect_helix_1elwA_A108_A115.pdb");
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			target_coords.push_back( pose.residue(i).xyz("CA") );
		}
	}

	HelicalParamsFitMultifunc func( target_coords );

	Size iterations;
	Real const tolerance( option[ basic::options::OptionKeys::run::min_tolerance ] );
	Real final_func_value;
	Real start_func_value( func( params ) );
	optimization::powell( params, func, tolerance, iterations, final_func_value );

	cout << "final_params " << F(9,3,start_func_value) << F(9,3,final_func_value) <<
		F(9,3,params[1]) << F(9,3,degrees(params[2])) <<F(9,3,params[3]) << endl;

	Multivec const idl_params( params );

	vector1< Vectors > all_helix_coords;
	strings all_helix_tags;
	{
		string const pdb_coords_file( option[ my_options::pdb_coords_file ] );
		ifstream data( pdb_coords_file.c_str() );
		string line;
		while ( getline( data, line ) ) {
			strings const l( split_to_vector1( line ) );
			// if ( l[1] == "helix_coords:" ) {
			if ( l[1] == "clean_helix_coords:" ) {
				Size const helixlen( int_of( l[2] ) );
				Vectors coords;
				for ( Size i=1; i<= helixlen; ++i ) {
					coords.push_back( Vector( float_of( l[ 3*i+1 ] ), float_of( l[ 3*i+2 ] ), float_of( l[3*i+3 ] ) ) );
				}
				string const tag( filebase( l.back() )+"_"+ l[ l.size()-1 ]+"_"+l[ l.size()-3 ] );
				all_helix_coords.push_back( coords );
				all_helix_tags.push_back( tag );
			}
		}
		data.close();
		TR.Trace << "Read " << all_helix_coords.size() << " helix coords from file " << pdb_coords_file << endl;
	}

	Size const helix_window( 7 );
	target_coords.resize( helix_window );
	for ( Size ii=1; ii<= all_helix_coords.size(); ++ii ) {
		Vectors const & pdb_coords( all_helix_coords[ii] );
		if ( pdb_coords.size() < helix_window ) continue;

		for ( Size i=1; i<= pdb_coords.size()-helix_window+1; ++i ) {
			func.set_target_coords( i, i+helix_window-1, pdb_coords );

			params = idl_params;

			Real start_func_value( func( params ) );
			optimization::powell( params, func, tolerance, iterations, final_func_value );

			cout << "final_params_helix_window " << I(4,ii) << I(4,i) << F(9,3,start_func_value) << F(9,3,final_func_value) <<
				" rise: " << F(9,3,params[1]) <<
				" twist: " << F(9,3,degrees(params[2])) <<
				" radius: " << F(9,3,params[3]) << ' ' << all_helix_tags[ii] << endl;


			/// now try optimizing tilt params
			func.optimize_tilt( true );
			params.resize( 5 );
			params[1] = idl_params[1];
			params[2] = idl_params[2];
			params[3] = 0;
			params[4] = 0;
			params[5] = idl_params[3];

			start_func_value = func( params );
			optimization::powell( params, func, tolerance, iterations, final_func_value );

			cout << "final_params_helix_window_tilt " << I(4,ii) << I(4,i) <<
				F(9,3,start_func_value) << F(9,3,final_func_value) <<
				" rise: " << F(9,3,params[1]) <<
				" twist: " << F(9,3,degrees(params[2])) <<
				" tilt: " << F(9,3,degrees(params[3])) <<
				" tilt_direction: " << F(9,3,degrees(params[4])) <<
				" radius: " << F(9,3,params[5]) << ' ' << all_helix_tags[ii] << endl;

			func.optimize_tilt( false );
		} // window loop

		/// now the whole deal
		func.set_target_coords( 1, pdb_coords.size(), pdb_coords );

		params = idl_params;

		Real start_func_value( func( params ) );
		optimization::powell( params, func, tolerance, iterations, final_func_value );

		cout << "final_params_helix_full " << I(4,ii) << I(4,pdb_coords.size()) <<
			F(9,3,start_func_value) << F(9,3,final_func_value) <<
			" rise: " << F(9,3,params[1]) <<
			" twist: " << F(9,3,degrees(params[2])) <<
			" radius: " << F(9,3,params[3]) << ' ' << all_helix_tags[ii] << endl;

		/// now try optimizing tilt params
		func.optimize_tilt( true );
		params.resize( 5 );
		params[1] = idl_params[1];
		params[2] = idl_params[2];
		params[3] = 0;
		params[4] = 0;
		params[5] = idl_params[3];

		start_func_value = func( params );
		optimization::powell( params, func, tolerance, iterations, final_func_value );

		cout << "final_params_helix_full_tilt " << I(4,ii) << I(4,pdb_coords.size()) <<
			F(9,3,start_func_value) << F(9,3,final_func_value) <<
			" rise: " << F(9,3,params[1]) <<
			" twist: " << F(9,3,degrees(params[2])) <<
			" tilt: " << F(9,3,degrees(params[3])) <<
			" tilt_direction: "  << F(9,3,periodic_range( degrees(params[4]), 360.0 ) ) <<
			" radius: " << F(9,3,params[5]) << ' ' << all_helix_tags[ii] << endl;

		func.optimize_tilt( false );

		/// now try optimizing tilt params plus precession
		func.optimize_tilt( true );
		func.optimize_tilt_precession( true );
		params.resize( 6 );
		params[1] = idl_params[1];
		params[2] = idl_params[2];
		params[3] = 0;
		params[4] = 0;
		params[5] = 0;
		params[6] = idl_params[3];

		start_func_value = func( params );
		optimization::powell( params, func, tolerance, iterations, final_func_value );

		cout << "final_params_helix_full_tilt_precession " << I(4,ii) << I(4,pdb_coords.size()) <<
			F(9,3,start_func_value) << F(9,3,final_func_value) <<
			" rise: " << F(9,3,params[1]) <<
			" twist: " << F(9,3,degrees(params[2])) <<
			" tilt: " << F(9,3,degrees(params[3])) <<
			" tilt_direction: "  << F(9,3,periodic_range( degrees(params[4]), 360.0 ) ) <<
			" tilt_precession: " << F(9,3,periodic_range( degrees(params[5]), 360.0 ) ) <<
			" radius: " << F(9,3,params[6]) << ' ' << all_helix_tags[ii] << endl;

		func.optimize_tilt( false );
		func.optimize_tilt_precession( false );


	} // all-coords loop


}
///////////////////////////////////////////////////////////////////////////////
void
tal_repeat_params_test()
{
	string const filename( start_file() );
	Pose pose;
	pose_from_pdb( pose, filename );

	Size const startpos( pose.pdb_info()->pdb2pose( 'A', 334 ) ), repeatlen( 34 ), nrepeat( 8 );

	SizePairs mapping;
	for ( Size i=1; i<= nrepeat; ++i ) {
		Size const repeatbegin( startpos+repeatlen*(i-1));
		cout << "rep " << I(2,i) << pose.sequence().substr( repeatbegin-1, repeatlen ) <<
			I(4,pose.pdb_info()->number(repeatbegin)) <<endl;
		if ( i<nrepeat ) {
			for ( Size j=0; j< repeatlen; ++j ) mapping.push_back( make_pair( repeatbegin+j, repeatbegin+repeatlen+j ) );
		}
	}

	// these looked good, so let's superimpose
	Pose movpose( pose );
	Real const rmsd( superimpose_pose_using_CA_mapping( movpose, pose, mapping ) );

	Stub const istub( pose.residue(1).xyz("CA"), pose.residue(2).xyz("CA"), pose.residue(3).xyz("CA") ),
		jstub( movpose.residue(1).xyz("CA"), movpose.residue(2).xyz("CA"), movpose.residue(3).xyz("CA") );
	Vector center, t, n;
	Real theta;
	get_stub_transform_data( istub, jstub, center, n, t, theta );

	Real const twist = numeric::conversions::degrees( theta );
	Real const rise = t.dot( n );

	cout << "rise: " << F(9,3,rise) << " twist: " << F(9,3,twist) << " rmsd: "<< F(9,3,rmsd) << ' ' <<
		filename << endl;


}



///////////////////////////////////////////////////////////////////////////////
void
annotate_pdb_test()
{
	//bool const skip_termini_bb_

	using namespace id;
	Size const nrepeat( 10 );
	Real const probe_radius_for_unsat_calcs( 1.0 );
	strings const files( start_files() );


	for ( Size fi=1; fi<= files.size(); ++fi ) {
		string const filename( files[fi] );
		Pose pose;
		pose_from_pdb( pose, filename );

		if ( !pose.residue( pose.total_residue() ).is_protein() ) {
			pose.conformation().delete_residue_slow( pose.total_residue());
		}

		AtomID_Map< Real > unsat_frequency;
		initialize_atomid_map( unsat_frequency, pose, 0.0 );

		Reals avg_rsd_sasapack_scores( pose.total_residue(), 0.0 );

		vector1< Reals > all_rsd_sasapack_scores;
		Real sasapack_score(0.0), donbb_per_35aa(0), donsc_per_35aa(0), accbb_per_35aa(0), accsc_per_35aa(0),
			buried_nonpolar_sasa(0), buried_nonpolar_sasa_sc(0);
		for ( Size n=1; n<= nrepeat; ++n ) {

			vector1< id::AtomID > buried_unsatisfied_donors, buried_unsatisfied_donors_all,
				buried_unsatisfied_acceptors, buried_unsatisfied_acceptors_all;

			Pose posetmp( pose );
			randomly_shift_and_tilt_pose( posetmp );
			bools const subset( pose.total_residue(), true );

			find_buried_unsatisfied_polars( subset, posetmp, buried_unsatisfied_donors_all, buried_unsatisfied_acceptors_all,
				filename, probe_radius_for_unsat_calcs );
			// ignore specific backbone atoms at the termini
			foreach_ ( AtomID atm, buried_unsatisfied_donors_all ) {
				if ( atm.rsd() == 1 && pose.residue( atm.rsd() ).atom_is_backbone( atm.atomno() ) ) continue;
				buried_unsatisfied_donors.push_back( atm );
			}
			foreach_ ( AtomID atm, buried_unsatisfied_acceptors_all ) {
				if ( atm.rsd() == pose.total_residue() && pose.residue( atm.rsd() ).atom_is_backbone( atm.atomno() ) ) continue;
				buried_unsatisfied_acceptors.push_back( atm );
			}
			Real const increment( 100.0 / Real(nrepeat) );
			foreach_ ( AtomID atm, buried_unsatisfied_donors    ) unsat_frequency[atm] += increment;
			foreach_ ( AtomID atm, buried_unsatisfied_acceptors ) unsat_frequency[atm] += increment;

			Size donbb, donsc, accbb, accsc;
			count_bb_and_sc_atomids( buried_unsatisfied_donors   , pose, donbb, donsc );
			count_bb_and_sc_atomids( buried_unsatisfied_acceptors, pose, accbb, accsc );

			donbb_per_35aa += ( 35.0 * donbb ) / ( pose.total_residue() * nrepeat );
			donsc_per_35aa += ( 35.0 * donsc ) / ( pose.total_residue() * nrepeat );
			accbb_per_35aa += ( 35.0 * accbb ) / ( pose.total_residue() * nrepeat );
			accsc_per_35aa += ( 35.0 * accsc ) / ( pose.total_residue() * nrepeat );

			Real exposed_polar_sasa, exposed_nonpolar_sasa, buried_polar_sasa, tmp_buried_nonpolar_sasa,
				tmp_buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc;
			get_exposed_and_buried_sasas_by_atomtype( pose, subset, exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, tmp_buried_nonpolar_sasa, tmp_buried_nonpolar_sasa_sc,
				buried_nonpolar_rsd_sasa_sc );
			buried_nonpolar_sasa    += tmp_buried_nonpolar_sasa    / ( pose.total_residue() * nrepeat );
			buried_nonpolar_sasa_sc += tmp_buried_nonpolar_sasa_sc / ( pose.total_residue() * nrepeat );

			Reals rsd_sasapack_scores, rsd_sasa14_normalized, rsd_norme_scores;

			Real average_sasapack, average_normsasa;
			protocols::sasa_scores::compute_sasapack_scores( posetmp, rsd_sasapack_scores, rsd_sasa14_normalized,
				average_sasapack, average_normsasa );
			TR.Trace << "average_sasapack: " << F(9,3,average_sasapack) << I(4,n) << ' ' << filename << endl;
			for ( Size i=1; i<= pose.total_residue(); ++i ) avg_rsd_sasapack_scores[i] += rsd_sasapack_scores[i]/nrepeat;
			all_rsd_sasapack_scores.push_back( rsd_sasapack_scores );
			sasapack_score += average_sasapack/nrepeat;
		}

		cout << "final_scores " <<
			" bsasa: "    << F(5,1,buried_nonpolar_sasa ) <<
			" bsasa_sc: " << F(5,1,buried_nonpolar_sasa_sc ) <<
			" sasapack: " << F(5,2,sasapack_score) <<
			" donbb35: "  << F(5,2,donbb_per_35aa) <<
			" donsc35: "  << F(5,2,donsc_per_35aa) <<
			" accbb35: "  << F(5,2,accbb_per_35aa) <<
			" accsc35: "  << F(5,2,accsc_per_35aa) <<
			' ' << filename << endl;

		{ // now try writing a bfactor pdb for unsat hbonds
			string const outfilename( output_tag() + "anno_hb_"+filebase( filename ) );
			ofstream out( outfilename.c_str() );
			utility_exit();
			//io::pdb::dump_bfactor_pdb( pose, unsat_frequency, out, "NO_MODEL_LINE_IN_OUTPUT" ); // stupid tag
			out.close();
		}

		{ // now try writing a bfactor pdb for sasapack
			AtomID_Map< Real > pack_bfactor;
			initialize_atomid_map( pack_bfactor, pose, 0.0 );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				for ( Size j=1; j<= pose.residue(i).natoms(); ++j ) {
					pack_bfactor[ AtomID( j,i) ] = 50.0 + avg_rsd_sasapack_scores[i];
				}
			}

			string const outfilename( output_tag() + "anno_pack_"+filebase( filename ) );
			ofstream out( outfilename.c_str() );
			utility_exit();
			//io::pdb::dump_bfactor_pdb( pose, pack_bfactor, out, "NO_MODEL_LINE_IN_OUTPUT" ); // stupid tag
			out.close();

			// what are the stddevs of the scores?
			Reals devs( pose.total_residue(), 0.0 );
			for ( Size n=1; n<= nrepeat; ++n ) {
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					devs[i] += numeric::square( all_rsd_sasapack_scores[n][i] - avg_rsd_sasapack_scores[i] );
				}
			}
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				TR.Trace << "sasapack_stdev: " << I(4,i) << ' ' << pose.residue(i).name1() <<
					" score: " << F(9,3,avg_rsd_sasapack_scores[i]) <<
					" stddev: " << F(9,3,sqrt( devs[i]/nrepeat )) << ' ' << filename << endl;
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
db_turn_test()
{
	// {
	//  Pose pose;
	//  pose.append_residue_by_bond( *get_vanilla_dna_residue('a') );
	//  pose.dump_pdb("ADE.pdb");
	// }
	// {
	//  Pose pose;
	//  pose.append_residue_by_bond( *get_vanilla_dna_residue('c') );
	//  pose.dump_pdb("CYT.pdb");
	// }
	// exit(0);
	strings const turns( option[ my_options::turns ] );//split_to_vector1( string("B BAAB BAB BB E G GABB GB GBB GBBBB") ) );
	map< string, StubFrags > frags;
	build_turn_library_v2( turns, frags );

}


///////////////////////////////////////////////////////////////////////////////

void
native_hand_test()
{
	strings const files( start_files() );


	foreach_ ( string const filename, files ) {
		Pose pose;
		pose_from_pdb( pose, filename );

		// parse filename
		Size nrepeats(0), median_repeatlen(0);
		{
			strings const l( split_to_vector1( filebase(filename), ".") );
			nrepeats = int_of( l[2].substr(1) );
			median_repeatlen = int_of( l[3].substr(1) );
		}

		Vectors coords;
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			coords.push_back( pose.residue(i).xyz("CA") );
		}

		if ( coords.size() <= median_repeatlen ) {
			cout << "SKIP too short " << filename << ' ' << coords.size() << ' ' << median_repeatlen << endl;
			cerr << "SKIP too short " << filename << ' ' << coords.size() << ' ' << median_repeatlen << endl;
			continue;
		}

		Real const hand( get_chirality( median_repeatlen, coords ) );

		cout << "hand: " << F(9,3,hand) <<
			" nrepeats: " << nrepeats <<
			" median_repeatlen: " << median_repeatlen <<
			" nres: " << coords.size() <<
			' ' << filename << endl;

	}




}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
class VectorsDistanceMetricWithShiftingVariableRepeatlen { // but constant number of repeats
public:
	VectorsDistanceMetricWithShiftingVariableRepeatlen(
																										 Size const nrepeat,
																										 Size const ncluster
																										 ):
		nrepeat_( nrepeat ),
		ncluster_( ncluster ),
		best_shift_( 0 )
	{
		coords1_.resize( (nrepeat-1) * ncluster );
		coords2_.resize( (nrepeat-1) * ncluster );
	}

public:
	void
	setup_coords( Vectors const & coordsfull, Size const offset, Vectors & coords ) const
	{
		pbassert( coordsfull.size()%nrepeat_ == 0 );
		Size const repeatlenfull( coordsfull.size() / nrepeat_ );
		pbassert( repeatlenfull >= ncluster_ );
		pbassert( coords.size() == ( nrepeat_-1)*ncluster_ );
		for ( Size i=0; i< nrepeat_-1; ++i ) {
			for ( Size j=1; j<= ncluster_; ++j ) {
				Size const pos( i*ncluster_ + j ), posfull( i*repeatlenfull + offset + j );
				pbassert( pos<= coords.size() && posfull <= coordsfull.size() );
				coords[ i*ncluster_ + j ] = coordsfull[ i*repeatlenfull + offset + j ];
			}
		}
	}


public:

	Real
	operator()( Vectors const & coords1full_in, Vectors const & coords2full_in ) const
	{
		Vectors const & coords1full( coords1full_in.size() <= coords2full_in.size() ? // shorter
																 coords1full_in : coords2full_in );
		Vectors const & coords2full( coords1full_in.size() <= coords2full_in.size() ? // longer
																 coords2full_in : coords1full_in );

		Size const repeatlen2( coords2full.size()/nrepeat_ );

		Real min_rmsd( 1e6 );
		for ( Size shift=0; shift< repeatlen2; ++shift ) {
			setup_coords( coords1full, 0    , coords1_ );
			setup_coords( coords2full, shift, coords2_ );
			Real const rmsd( numeric::model_quality::calc_rms( coords1_, coords2_ ) );
			if ( rmsd < min_rmsd ) {
				min_rmsd = rmsd;
				best_shift_ = shift; // save the optimal shift value
			}
		}

		return min_rmsd;
	}

	Size
	best_shift() const { return best_shift_; } // get best alignment from previous operator() call

private:
	Size const nrepeat_;
	Size const ncluster_;
	mutable Size best_shift_;
	mutable Vectors coords1_;
	mutable Vectors coords2_;
};


///////////////////////////////////////////////////////////////////////////////////////
/// protein only clustering
///
void
unbound_cluster_test()
{
	if ( !option[ my_options::clustering_params ].user() || option[ my_options::clustering_params].size() != 8 ) {
		cout << "-clustering_params min_threshold max_threshold " <<
			" min_cluster_size min_top_cluster_size try_top_cluster_size max_top_cluster_size " <<
			" max_clusters max_decoys_per_cluster_pdbfile" << endl;
		utility_exit();
	}

	Real min_threshold, max_threshold;
	Size min_cluster_size, min_top_cluster_size, try_top_cluster_size, max_top_cluster_size,
		max_clusters, max_decoys_per_cluster_pdbfile;

	parse_clustering_params( min_threshold, max_threshold, min_cluster_size, min_top_cluster_size, try_top_cluster_size,
		max_top_cluster_size, max_clusters, max_decoys_per_cluster_pdbfile );

	//ScoreFunctionOP fa_scorefxn( get_score_function_from_command_line() );

	bool const force_nrepeat( option[ my_options::nrepeat ].user() );

	// Reals const clustering_params( option[ my_options::clustering_params ] );
	// pbassert( clustering_params.size() == 8 );
	// Real const min_threshold( clustering_params[1] ), max_threshold( clustering_params[2] );
	// Size const min_cluster_size( clustering_params[3] ),
	// 	min_top_cluster_size( clustering_params[4] ),
	// 	try_top_cluster_size( clustering_params[5] ),
	// 	max_top_cluster_size( clustering_params[6] ),
	// 	max_clusters( clustering_params[7] ),
	// 	max_decoys_per_cluster_pdbfile( clustering_params[8] );

	// how many coords to cluster over?
	// will use window consisting of anchorpos-clusterwindow --> anchorpos+clusterwindow
	//Size const clusterwindow( option[ my_options::clusterwindow ] );

// 	Size const cluster_coords_begin( option[ my_options::cluster_coords ]()[1] ),
// 		cluster_coords_end( option[ my_options::cluster_coords ]()[2] );

	//Size const min_overlap( 30 );
	// coords for each decoy and for the natives will be shifted so that anchorpos is position "centerpos", ie
	// the middle of the clustering window
	Size nrepeat( 0 ); //, repeatlen( 0 );//, base_repeat( 3 );
	if ( force_nrepeat ) nrepeat = option[ my_options::nrepeat ];
	Size ncluster( 10000 );


	Size counter(0);

	//VectorsDistanceMetricNoSuper distance_metric;

	strings all_filenames, all_scorelines;
	vector1< Vectors > all_coords;
// 	strings all_filenames, all_seq, all_ss, all_startfiles, all_bb, all_sspred, all_scorelines;
// 	vector1< Vectors > all_coords;
// 	Reals all_rmsds, all_scores, all_interface_scores, all_anchorpos_distances, all_sasapack_scores, all_norme_scores,
// 		all_normsasa_scores, all_hands, all_sspred_matches, all_refold_rmsds;
// 	Sizes all_repeatlens;


	/// read the logfile

	string const logfile( start_file() );
	//bool rescored( false ), new_format( false ), have_sspred( option[ my_options::pred_ss ] ), have_refold_rmsds( false );

	ifstream data( logfile.c_str() );


// rhino1 symdes$ showfields tmp9 | more
// 1 ./output/job9_1_000_hyraxB50_br2norm/job9_1_000_0.log:final_scores
// 2 -70.994
// 3 job9_1_000_0symfragrebuild_symmpose_02.pdb_Srevcendesign_N0001.pdb
// 4 relax_rmsd:
// 5 13.652
// 6 passed_score_filter:
// 7 0
// 8 centroid_score:
// 9 44.430
// 10 centroid_rmsd:
// 11 13.076
// 12 centroid_simtime:
// 13 1.242
// 14 relax_simtime:
// 15 26.753
// 16 anchorpos_loop_begin:
// 17 14
// 18 anchorpos_loop_end:
// 19 16
// 20 cutpoint_loop_begin:
// 21 34
// 22 cutpoint_loop_end:
// 23 2
// 24 anchorseq:
// 25 D
// 26 repeatseq:
// 27 GIAYPQLVAGGRASDSTGASQAERLDNWLDKWGA
// 28 repeatss:
// 29 LLLHHHHHHHHHHLLLHHHHHHHHHHHHHHHHLL
// 30 na:
// 31 c
// 32 unbound_score:
// 33 -69.305
// 34 interface_score:
// 35 -1.551
// 36 weighted_energies:

	string line;
	while ( getline( data, line ) ) {
		if ( line.find("final_scores") == string::npos ) continue;
		++counter;
		if ( counter%10000 == 0 ) cout << counter << ' ' << all_coords.size() << endl;

		//istringstream l( line );
		strings const l( split_to_vector1( line ) );
		//Size const NF( l.size() );
		string dirtag( l[1] ), filename( dirtag.find(".out") == string::npos ?
			dirtag.substr( 0, dirtag.find_last_of( "/" ) ) +"/" + l[3] :
			dirtag.substr( 0, dirtag.find( ".out" ) ) +".pdbs/" + l[3] );

		//Real const score( float_of( l[2] ) );
		bool passed_score_filter( false ), relaxed( false );

		Size nrepeat_this_decoy( 0 ), repeatlen_this_decoy(0);
		if ( l[14] == "passed_score_filter:" ) {
			passed_score_filter = ( l[15] == "1" );
			relaxed = ( l[17] == "1" );

			if ( l[40] == "repeatlen:" ) {
				pbassert( l[42] == "nrepeat:" );
				repeatlen_this_decoy = int_of( l[41] );
				nrepeat_this_decoy   = int_of( l[43] );
			} else {
				pbassert( l[46] == "repeatlen:" );
				pbassert( l[48] == "nrepeat:" );
				repeatlen_this_decoy = int_of( l[47] );
				nrepeat_this_decoy   = int_of( l[49] );
			}
		} else if ( l[14] == "passed_fullatom_score_filter:" ) {
			passed_score_filter = ( l[15] == "1" );
			relaxed = ( l[19] == "1" );
			pbassert( l[42] == "repeatlen:" );
			pbassert( l[44] == "nrepeat:" );
			repeatlen_this_decoy = int_of( l[43] );
			nrepeat_this_decoy   = int_of( l[45] );
		} else utility_exit_with_message("bad line "+line );

		if ( force_nrepeat && nrepeat_this_decoy != nrepeat ) continue;

		if ( !nrepeat ) nrepeat = nrepeat_this_decoy;
		else { pbassert( nrepeat == nrepeat_this_decoy ); }

// 		if ( !repeatlen ) repeatlen = repeatlen_this_decoy;
// 		else { pbassert( repeatlen == repeatlen_this_decoy ); }

		ncluster = min( ncluster, repeatlen_this_decoy );

		if ( ! ( passed_score_filter && relaxed ) ) continue;

		if ( !utility::file::file_exists( filename ) && utility::file::file_exists( filename+".gz" ) ) {
			filename += ".gz";
		}

		if ( !utility::file::file_exists( filename ) ) {
			cout << "MISSING file: " << filename << endl;
			continue;
		}

		/// read the clustering coords
		Vectors coords( read_CA_coords_from_file( 1, nrepeat*repeatlen_this_decoy, filename ) );

		if ( coords.size() != nrepeat*repeatlen_this_decoy ) {
			TR.Trace << "bad coords size: "<< coords.size() << ' ' << filename << endl;
		}
		//Real const handedness( get_chirality( repeatlen_this_decoy, coords ) );


		all_coords.push_back( coords );
		all_filenames.push_back( filename );
		all_scorelines.push_back( line );

		if ( all_coords.size()%25 == 0 ) cout << "read " << I(4,all_coords.size()) << " files, nlines= " << I(6,counter) <<endl;
	}


	// now do the clustering
	//// now cluster
	VectorsDistanceMetricWithShiftingVariableRepeatlen distance_metric( nrepeat, ncluster );
	//VectorsDistanceMetricWithShifting distance_metric( repeatlen );
	cout << "clustering... ndecoys= " << all_coords.size() << endl;
	Real threshold;
	vector1< Sizes > cluster_members;
	Sizes nbr_counts;
	devel::blab::cluster::simple_cluster( all_coords, distance_metric, min_cluster_size, min_top_cluster_size, try_top_cluster_size,
		max_top_cluster_size, min_threshold, max_threshold,
		threshold, nbr_counts, cluster_members );

	cout << "simple_cluster_status: ndecoys: " << I(6,all_coords.size()) <<
		" n_clusters: " << I(4,cluster_members.size()) <<
		" cluster_sizes(big,small): " << I(4,cluster_members.front().size()) << I(4,cluster_members.back().size()) <<
		" threshold: " << F(9,3,threshold) << endl;

	string const prefix( output_tag() + filebase( logfile ) +
											 "_nrep"+string_of(nrepeat)+
											 //"_CW" + string_of( clusterwindow ) +
											 "_N" + string_of( all_coords.size() ) +
											 "_C" + string_of( cluster_members.front().size() ) +
											 "_T" + string_of( int( 100 * threshold ) ) + "_" );

	string const logfilename( prefix+".log");
	std::ofstream out( logfilename.c_str() );
	out << "simple_cluster_status: n_files: " << I(8,counter) <<
		" ndecoys: " << I(6,all_coords.size()) <<
		" n_clusters: " << I(4,cluster_members.size()) <<
		" cluster_sizes(big,small): " << I(4,cluster_members.front().size()) << I(4,cluster_members.back().size()) <<
		" threshold: " << F(9,3,threshold) << endl;

	for ( Size i=1; i<= all_coords.size(); ++i ) {
		out << "nbr_count: " << nbr_counts[i] << ' ' << all_filenames[i] << '\n';
	}

	vector1< Vectors > cluster_center_coords;

	for ( Size clusterno=1; clusterno<= cluster_members.size() && clusterno <= max_clusters; ++clusterno ) {
		Sizes const & members( cluster_members[ clusterno ] );


		string const clusterfilename( prefix + lead_zero_string_of( clusterno, 3 ) +"_" +
																	lead_zero_string_of( members.size(),3 ) +".pdb" );
		ofstream pdbout( clusterfilename.c_str());

		// write all the scorelines to the top of the cluster file
		for ( Size i=1; i<= members.size(); ++i ) {
			pdbout << "REMARK SCORELINE " << i << ' ' << members.size() << ' ' << all_scorelines[ members[i] ] << '\n';
		}

		Sizes const pdb_members( random_subset( max_decoys_per_cluster_pdbfile, members ) );
		for ( Size ii=1; ii<= pdb_members.size(); ++ii ) {
			Size const m( ii == 1 ? members[ii] : pdb_members[ii] ); // always include the cluster center first
			//Pose pose;
			//pose_from_pdb( pose, all_filenames[ members[ii] ] );
			//Stub const & T( symm_pose_transforms[ all_startfiles[ m ] ] );
			Stub T;

			if ( ii>1 ) {
				// superimpose to fit cluster center
				Size const center( members.front() );
				Vectors const & coords1full( all_coords[ center ] ), &coords2full( all_coords[m] );
				Real const rms( distance_metric( coords1full, coords2full ) );
				if ( rms > threshold + 1e-3 ) cerr << "threshold error? " << F(9,3,threshold) << F(9,3,rms) << endl;
				Size const best_shift( distance_metric.best_shift() );
				Size const shift1( coords1full.size() <= coords2full.size() ? 0 : best_shift );
				Size const shift2( coords1full.size() <= coords2full.size() ? best_shift : 0 );
				Vectors coords1( (nrepeat-1)*ncluster ), coords2( (nrepeat-1)*ncluster );
				distance_metric.setup_coords( coords1full, shift1, coords1 );
				distance_metric.setup_coords( coords2full, shift2, coords2 );
				// now rmsfit
				Vectors const coords2_before( coords2 );
				superimpose_coords( coords1, coords2 );
				Real const rms_redo2( numeric::model_quality::calc_rms( coords1, coords2 ) );
				if ( abs( rms - rms_redo2 )>1e-2 ) cerr << "big rms_redo2 error: " << F(9,3,rms) << F(9,3,rms_redo2) << endl;
				//if ( abs( rms - rms_redo )>1e-2 ) cerr << "big rms_redo error: " << F(9,3,rms) << F(9,3,rms_redo) << endl;
				// now figure out what the transformation is
				Stub const stub1( coords2_before[1], coords2_before[2], coords2_before[3] ),
					stub2( coords2[1], coords2[2], coords2[3] );
				// R * stub1.M = stub2.M
				Matrix const R( stub2.M * stub1.M.transposed() );
				// R * x_old + v = x_new
				Vector const v( coords2[1] - R * coords2_before[1] );
				T.M = R;
				T.v = v;
			}
// 			if ( ii>1 ) {
// 				// superimpose to fit cluster center
// 				Size const center( members.front() );
// 				Real const rms( distance_metric( all_coords[ center ], all_coords[ m ] ) );
// 				if ( rms > threshold + 1e-3 ) cerr << "threshold error? " << F(9,3,threshold) << F(9,3,rms) << endl;
// 				Size const shift( distance_metric.best_shift() );
// 				Vectors coords1( all_coords[center] ), coords2( all_coords[ m ] );
// 				Vectors const & coords2_full( all_coords[m] );
// 				// trim off last repeat, shift coords2 by shift
// 				coords1.erase( coords1.end() - repeatlen , coords1.end() );
// 				coords2.erase( coords2.end() - repeatlen , coords2.end() );
// 				for ( Size i=1; i<= shift; ++i ) {
// 					coords2.erase( coords2.begin() );
// 					Size const newpos( coords2_full.size() - repeatlen + i );
// 					pbassert( coords2_full[ newpos-1 ].distance_squared( coords2.back() ) < 1e-3 );
// 					coords2.push_back( coords2_full[ newpos ] );
// 				}
// 				// now rmsfit
// 				Vectors const coords2_before( coords2 );
// 				Real const rms_redo( numeric::model_quality::rmsfit_wrapper( coords1, coords2 ) );
// 				Real const rms_redo2( numeric::model_quality::rms_wrapper( coords1, coords2 ) );
// 				if ( abs( rms - rms_redo2 )>1e-2 ) cerr << "big rms_redo2 error: " << F(9,3,rms) << F(9,3,rms_redo2) << endl;
// 				if ( abs( rms - rms_redo )>1e-2 ) cerr << "big rms_redo error: " << F(9,3,rms) << F(9,3,rms_redo) << endl;
// 				// now figure out what the transformation is
// 				Stub const stub1( coords2_before[1], coords2_before[2], coords2_before[3] ),
// 					stub2( coords2[1], coords2[2], coords2[3] );
// 				// R * stub1.M = stub2.M
// 				Matrix const R( stub2.M * stub1.M.transposed() );
// 				// R * x_old + v = x_new
// 				Vector const v( coords2[1] - R * coords2_before[1] );
// 				T.M = R;
// 				T.v = v;
// 			}

			//pose.apply_transform_Rx_plus_v( T.M, T.v );
			pdbout << "REMARK MODELNO " << I(4,ii) << ' ' << all_filenames[m] << '\n';
			pdbout << "MODEL" << I(9,ii) << '\n';
			//io::pdb::FileData::dump_pdb( pose, pdbout );
			// should be faster than pose i/o
			Size const nres(all_coords[m].size());
			transform_pdbfile_coords_and_append_to_stream( all_filenames[ m ], T.M, T.v, 1, nres, pdbout );
			pdbout << "ENDMDL\n";
		}
		pdbout.close();
		cout << "cluster " << I(4,clusterno) << I(4,members.size()) << ' ' << clusterfilename << endl;


// 		// write the seq/ss for this cluster
// 		out << "CLUSTER_SEQ: " << I(4,clusterno) << I(6,members.size());
// 		for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << all_seq[ members[ii] ];
// 		out << '\n';

// 		out << "CLUSTER_SS: " << I(4,clusterno) << I(6,members.size());
// 		for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << all_ss[ members[ii] ];
// 		out << '\n';

		out << "CLUSTER_FILENAMES: " << I(4,clusterno) << I(6,members.size());
		for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << all_filenames[ members[ii] ];
		out << '\n';

// 		if ( rescored || new_format ) {
// 			out << "CLUSTER_BB: " << I(4,clusterno) << I(6,members.size());
// 			for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << all_bb[ members[ii] ];
// 			out << '\n';
// 			out << "CLUSTER_SASAPACK_SCORES: " << I(4,clusterno) << I(6,members.size());
// 			for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_sasapack_scores[ members[ii] ] );
// 			out << '\n';
// 			out << "CLUSTER_NORME_SCORES: " << I(4,clusterno) << I(6,members.size());
// 			for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_norme_scores[ members[ii] ] );
// 			out << '\n';
// 			if ( new_format ) {
// 				out << "CLUSTER_NORMSASA_SCORES: " << I(4,clusterno) << I(6,members.size());
// 				for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_normsasa_scores[ members[ii] ] );
// 				out << '\n';
// 				out << "CLUSTER_REPEATLEN_SCORES: " << I(4,clusterno) << I(6,members.size());
// 				for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << I(4,all_repeatlens[ members[ii] ] );
// 				out << '\n';
// 			}
// 		}


// 		out << "CLUSTER_DIRECTIONS: " << I(4,clusterno) << I(6,members.size());
// 		for ( Size ii=1; ii<= members.size(); ++ii ) {
// 			out << ' ' << ( all_filenames[ members[ii] ].find("rev") == string::npos ? 1 : -1 );
// 		}
// 		out << '\n';

// 		out << "CLUSTER_HANDS: " << I(4,clusterno) << I(6,members.size());
// 		for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_hands[ members[ii] ] );
// 		out << '\n';

// 		out << "CLUSTER_NAT_RMSDS: " << I(4,clusterno) << I(6,members.size());
// 		for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_rmsds[ members[ii] ] );
// 		out << '\n';

// 		out << "CLUSTER_SCORES: " << I(4,clusterno) << I(6,members.size());
// 		for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_scores[ members[ii] ] );
// 		out << '\n';

// 		if ( have_sspred ) {
// 			out << "CLUSTER_SSPREDMATCHES: " << I(4,clusterno) << I(6,members.size());
// 			for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_sspred_matches[ members[ii] ] );
// 			out << '\n';
// 			out << "CLUSTER_SSPRED: " << I(4,clusterno) << I(6,members.size());
// 			for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << all_sspred[ members[ii] ];
// 			out << '\n';
// 		}

// 		if ( have_refold_rmsds ) {
// 			out << "CLUSTER_REFOLD_RMSDS: " << I(4,clusterno) << I(6,members.size());
// 			for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_refold_rmsds[ members[ii] ] );
// 			out << '\n';
// 		}

// 		out << "CLUSTER_INTERFACE_SCORES: " << I(4,clusterno) << I(6,members.size());
// 		for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_interface_scores[ members[ii] ] );
// 		out << '\n';

// 		out << "CLUSTER_ANCHORPOS_DISTANCES: " << I(4,clusterno) << I(6,members.size());
// 		for ( Size ii=1; ii<= members.size(); ++ii ) out << ' ' << F(9,3,all_anchorpos_distances[ members[ii] ] );
// 		out << '\n';

		Vectors const & center_coords( all_coords[ members[1] ] );
		cluster_center_coords.push_back( center_coords );

		out << "CLUSTER_RMSDS: " << I(4,clusterno) << I(6,members.size());
		for ( Size other_cluster=1; other_cluster<= clusterno; ++other_cluster ) {
			out << F(9,3,distance_metric( center_coords, cluster_center_coords[ other_cluster] ));
		}
		out << '\n';

	}
	out.close();
}

///////////////////////////////////////////////////////////////////////////////
void
unbound_frag_redesign_test()
{
	strings const files( start_files() );


	foreach_( string const & filename, files ) {
		Pose pose;
		pose_from_pdb( pose, filename );


		Reals rsd_sasa;
		get_residue_sasa_for_layers( pose, rsd_sasa );

		// Real const core_threshold_sse    ( option[ my_options::core_threshold_sse ] );
		// Real const core_threshold_loop   ( option[ my_options::core_threshold_loop ] );
		// Real const surface_threshold_sse ( option[ my_options::surface_threshold_sse ] );
		// Real const surface_threshold_loop( option[ my_options::surface_threshold_loop ] );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			cout << "layers_sasa: " << F(9,3,rsd_sasa[i]) << I(4,i) << ' ' << pose.residue(i).name1() << ' ' <<
				filebase( filename ) << endl;
		}

	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for biangle triangle type sims
//

string
get_single_helix_stats(
	Size const nrepeat,
	string repeatbb,
	string repeatseq,
	Vectors coords
)
{
	Size const repeatlen( repeatbb.size() );
	runtime_assert( coords.size() == repeatlen * nrepeat );
	runtime_assert( repeatbb.size() == repeatseq.size() );

	// figure out the helix len...
	while ( repeatbb[0] != 'A' ) { // move from beginning to end
		repeatbb.push_back( repeatbb[0] );
		repeatbb.erase( repeatbb.begin());
		repeatseq.push_back( repeatseq[0] );
		repeatseq.erase( repeatseq.begin());
		coords.push_back( coords.front() );
		coords.erase( coords.begin() );
	}
	while ( repeatbb[ repeatbb.size()-1 ] == 'A' ) { // move from end to beginning
		repeatbb.insert( repeatbb.begin(), repeatbb[ repeatbb.size()-1 ] );
		repeatbb.pop_back();
		repeatseq.insert( repeatseq.begin(), repeatseq[ repeatseq.size()-1 ] );
		repeatseq.pop_back();
		coords.insert( coords.begin(), coords.back() );
		coords.pop_back();
	}
	Size const hlen( repeatbb.find_first_not_of("A") );
	string const turn( repeatbb.substr(hlen));

	if ( hlen < 7 ) {
		cout << "SKIP short helix: " << repeatbb << ' ' << hlen << endl;
		return string("");
	}



	// these are not normed helix strains
	Real helix_strain, turn_strain, turn_strain_normed, tmp;
	compute_helix_strain( repeatlen+1, repeatlen+hlen, coords, true, helix_strain, tmp );

	// these are normed values, running from 0-1
	compute_turn_strain( turn, repeatlen+hlen+1, coords, turn_strain, turn_strain_normed, tmp );

	string const helix_stats( analyze_helix_coords( "helix1", repeatlen+1, repeatlen+hlen, coords ) );

	string const new_bbtag( string_of(nrepeat)+".4.-."+string_of(hlen-4)+"."+turn );

	ostringstream out;
	out << " new_repeatbb: " << repeatbb <<
		" new_repeatseq: " << repeatseq <<
		" new_bbtag: " << new_bbtag <<
		" new_helix1_len: " << hlen <<
		" new_turn1: " << turn <<
		" helix1_strain: " << helix_strain <<
		" turn1_strain: " << turn_strain <<
		" turn1_strain_normed: " << turn_strain_normed <<
		' ' << helix_stats;

	return out.str();

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for biangle triangle type sims
//

string
get_double_helix_stats(
	Size const nrepeat,
	string repeatbb,
	string repeatseq,
	Vectors coords
)
{
	Size const min_helix_len(7);

	Size const repeatlen( repeatbb.size() );
	runtime_assert( coords.size() == repeatlen * nrepeat );
	runtime_assert( repeatbb.size() == repeatseq.size() );

	if ( repeatbb.find('A') == string::npos ) {
		cout << "SKIP no A: " << repeatbb << endl;
		return string("");
	}

	// figure out the helix len...
	while ( repeatbb[0] != 'A' ) { // move from beginning to end
		repeatbb.push_back( repeatbb[0] );
		repeatbb.erase( repeatbb.begin());
		repeatseq.push_back( repeatseq[0] );
		repeatseq.erase( repeatseq.begin());
		coords.push_back( coords.front() );
		coords.erase( coords.begin() );
	}
	while ( repeatbb[ repeatbb.size()-1 ] == 'A' ) { // move from end to beginning
		repeatbb.insert( repeatbb.begin(), repeatbb[ repeatbb.size()-1 ] );
		repeatbb.pop_back();
		repeatseq.insert( repeatseq.begin(), repeatseq[ repeatseq.size()-1 ] );
		repeatseq.pop_back();
		coords.insert( coords.begin(), coords.back() );
		coords.pop_back();
	}


	Size const h1_len( repeatbb.find_first_not_of("A") );

	// where does the second helix start? There could be 'A's in the turn...
	string const hbb( min_helix_len,'A' );

	if ( h1_len < min_helix_len || repeatbb.substr(h1_len).find(hbb) == string::npos ) {
		cout << "SKIP short helices: " << repeatbb << ' ' << h1_len << endl;
		return string("");
	}

	Size const h1_begin(1), t1_begin( h1_len+1 ), h2_begin( repeatbb.substr(h1_len).find(hbb)+h1_len+1 ), // 1-indexed
		t1_len( h2_begin-t1_begin ), h2_len( repeatbb.substr(h2_begin-1).find_first_not_of("A") ),
		t2_begin( h2_begin+h2_len ), t2_len( repeatlen-t2_begin+1 );
	runtime_assert( h2_len >= min_helix_len );
	runtime_assert( repeatbb[ h2_begin-1 ] == 'A' && repeatbb[ h2_begin-2] != 'A' );

	string const turn1( repeatbb.substr(t1_begin-1,t1_len));
	string const turn2( repeatbb.substr(t2_begin-1,t2_len));


	// these are not normed helix strains
	Real h1_strain, h2_strain, t1_strain, t1_strain_normed, t2_strain, t2_strain_normed, tmp;
	Size const offset( repeatlen ); // coords are wonky near the terms due to shifting
	compute_helix_strain( offset+h1_begin, offset+h1_begin+h1_len-1, coords, true, h1_strain, tmp );
	compute_helix_strain( offset+h2_begin, offset+h2_begin+h2_len-1, coords, true, h2_strain, tmp );

	// these are normed values, running from 0-1
	compute_turn_strain( turn1, offset+t1_begin, coords, t1_strain, t1_strain_normed, tmp );
	compute_turn_strain( turn2, offset+t2_begin, coords, t2_strain, t2_strain_normed, tmp );

	string const h1_stats( analyze_helix_coords( "helix1", offset+h1_begin, offset+h1_begin+h1_len-1, coords ) );
	string const h2_stats( analyze_helix_coords( "helix2", offset+h2_begin, offset+h2_begin+h2_len-1, coords ) );

	string const new_bbtag( string_of(nrepeat)+"."+string_of(h1_len)+"."+turn1+"."+string_of(h2_len)+"."+turn2 );

	ostringstream out;
	out << " new_repeatbb: " << repeatbb <<
		" new_repeatseq: " << repeatseq <<
		" new_bbtag: " << new_bbtag <<
		" new_helix1_len: " << h1_len <<
		" new_turn1: " << turn1 <<
		" new_helix2_len: " << h2_len <<
		" new_turn2: " << turn2 <<
		" helix1_strain: " << h1_strain <<
		" helix2_strain: " << h2_strain <<
		" turn1_strain: " << t1_strain <<
		" turn1_strain_normed: " << t1_strain_normed <<
		" turn2_strain: " << t2_strain <<
		" turn2_strain_normed: " << t2_strain_normed <<
		' ' << h1_stats << ' ' << h2_stats;

	return out.str();

}



///////////////////////////////////////////////////////////////////////////////
//
// don't simulate the full number of repeats
//
void
new_unbound_frag_test()
{
	runtime_assert( !option[ my_options::local_simfile ] );
	runtime_assert( option[ my_options::unfolded_sasas ].user() ); // make sure have datafile for buried sasa
	runtime_assert( option[ my_options::resample_centroid_seqs ].user() ||
		option[ my_options::centroid_helixseq ].user() );

	bool const force_biangles( option[ my_options::force_biangles ] );
	bool const force_triangles( option[ my_options::force_triangles ] );
	bool const force_twistrise( option[ my_options::target_twist ].user() );

	//runtime_assert( option[ my_options::local_simfile ] );
	bool const resample_centroid_seqs( option[ my_options::resample_centroid_seqs ].user() );
	map< string, strings > bbtag2centroid_seq;
	if ( resample_centroid_seqs ) {
		strings l( option[ my_options::resample_centroid_seqs ] );
		string bbtag;
		while ( !l.empty() ) {
			if ( is_int( l.front().substr(0,1) ) ) {
				bbtag = l.front();
			} else {
				runtime_assert( bbtag.size() );
				TR.Trace << "bbtag2centroid_seq: " << bbtag << ' ' << l.front() << endl;
				bbtag2centroid_seq[ bbtag ].push_back( l.front() );
			}
			l.erase( l.begin() );
		}
	}


	//Real const min_helix2_len( 5 );

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
 	ScoreFunctionOP fa_scorefxn(0);
	if ( option[ OptionKeys::dna::specificity::score_function ].user() ) {
		fa_scorefxn = get_score_function_from_command_line();
	} else {
		fa_scorefxn = ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag );
	}
	//ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) ); // 09/16/15
	runtime_assert( pose::symmetry::is_symmetric( *score3 ) );

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	Real const donut_energy_weight( option[ my_options::donut_energy_weight ] );
	bool const use_donut_energy( donut_energy_weight>1e-3 );
	bool const use_fullatom_twistrise_energy( force_biangles || force_triangles || force_twistrise );
	if ( use_fullatom_twistrise_energy ) {
		runtime_assert( !use_donut_energy );
	} else {
		runtime_assert( use_donut_energy );
	}
	DonutWholeEnergy donut_whole_energy;
	Donut1B_Energy donut_1b_energy;
	donut_whole_energy.pose_is_subset( true );
	donut_1b_energy.pose_is_subset( true );
	if ( use_donut_energy ) {
		fa_scorefxn->add_extra_method( donut_whole, donut_energy_weight, donut_whole_energy );
		fa_scorefxn->add_extra_method( donut_1b   , donut_energy_weight, donut_1b_energy    );
	}


	clock_t search_starttime( clock() );

	Size const nstruct( option[ OptionKeys::out::nstruct ].user() ? option[ OptionKeys::out::nstruct ]() : 100000 );

	string const outdir( shared_output_dir_create_if_necessary() ); // has trailing /

	for ( Size n=1; n<= nstruct; ++n ) {
		// undo any fiddling
		symminfo_hack::nrepeat_ = 0;
		symminfo_hack::repeatlen_ = 0;

		clock_t sim_starttime = clock();

		Size const base_repeat( option[ my_options::base_repeat ].user() ? option[ my_options::base_repeat ]: 3 ),
			nrepeat( option[ my_options::nrepeat_sim ].user() ? option[ my_options::nrepeat_sim ] : 5 );

		Pose pose;

		// create a symmetric pose
		Size helix1_len(0), helix2_len(0);
		string turn1, turn2;
		Size repeatlen(0 );
		char expected_hand( '-' );
		string resample_bb_tag;
		//string tal_twist_type;

		Size nrepeat_full(0);


		if ( option[ my_options::resample_bbs ].user() ) { // NOTE: NO convention on which is the inner helix!!!
			// EXAMPLE: " -resample_bbs  8-9,13,GB,13,GB,L  9,14-17,GBB,11,GBB,L "
			string const bb( random_element( option[ my_options::resample_bbs ]() ) );
			strings l;
			if ( bb.find(',') != string::npos ) l = split_to_vector1( bb, "," );
			else l = split_to_vector1( bb, "." );
			runtime_assert( l.size() == 5 || l.size() == 6 );
			nrepeat_full = get_int_from_resample_bb_tag( l[1] );
			helix1_len   = get_int_from_resample_bb_tag( l[2] );
			helix2_len   = get_int_from_resample_bb_tag( l[4] );
			turn1 = l[3];
			turn2 = l[5];
			if ( l.size() == 6 ) {
				runtime_assert( l[6].size() == 1 );
				expected_hand = l[6][0];
				runtime_assert( expected_hand == 'R' || expected_hand == 'L' || expected_hand == 'U' );
			}
			repeatlen = helix1_len + get_turnlen( turn1 ) + helix2_len + get_turnlen( turn2 );
			// resample_bb_tag = l[1];
			// for ( Size j=2; j<= l.size(); ++j ) resample_bb_tag += "." + l[j];
		} else { //if ( option[ my_options::turn_pairs ].user() ) {
			nrepeat_full = random_element( option[ my_options::nrepeats ]() );
			//runtime_assert( nrepeat_full==24); // tmp hack
			strings const tp( split_to_vector1( random_element( option[ my_options::turn_pairs ]() ), "." ) );
			runtime_assert( tp.size() == 2 );
			turn1 = tp[1];
			turn2 = tp[2];
			helix1_len = random_element( option[ my_options::helix1_lens ]() );
			helix2_len = helix1_len + random_element( option[ my_options::helixlen_deltas ]() );
			repeatlen = helix1_len + helix2_len + get_turnlen( turn1 ) + get_turnlen( turn2 );
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// end setup
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		TR.Trace << "nrepeat: " << nrepeat << " repeatlen: " << repeatlen << ' ' << turn1 << ' ' << turn2 << ' ' <<
			helix1_len << ' ' << helix2_len << endl;


		Size const turn1_len( get_turnlen( turn1 ) ), turn2_len( get_turnlen( turn2 ) );
		string repeatbb(
			string( helix1_len, 'A' ) + ( turn1_len==0 ? string() : turn1 ) +
			string( helix2_len, 'A' ) + ( turn2_len==0 ? string() : turn2 ) );
		string repeatss(
			string( helix1_len, 'H' ) + string( turn1_len, 'L' ) +
			string( helix2_len, 'H' ) + string( turn2_len, 'L' ) );

		string const target_repeatbb( repeatbb ), target_repeatss( repeatss );

		runtime_assert( repeatbb.size() == repeatss.size() && repeatbb.size() == repeatlen );

		Size nres_protein( nrepeat * repeatlen );
		for ( Size i=1; i<= nres_protein; ++i ) {
			pose.append_residue_by_bond( *get_vanilla_protein_residue( 'L' ), true ); // build_ideal_geometry
		}


		/// setup the centroid repeat sequence
		string const bbtag
			( string_of( nrepeat_full ) + "." +
				string_of( helix1_len ) + "." + turn1 + "." +
				string_of( helix2_len ) + "." + turn2 );

		string centroid_repeatseq;

		if ( resample_centroid_seqs ) {
			runtime_assert( bbtag2centroid_seq.count( bbtag ) );
			centroid_repeatseq = random_element( bbtag2centroid_seq[ bbtag ] );
			runtime_assert( centroid_repeatseq.size() == repeatlen );
		} else {
			for ( Size i=1; i<= repeatlen; ++i ) {
				if ( repeatss[i-1] != 'L' ) { // stupid backwards compatibility...
					centroid_repeatseq.push_back( oneletter_code_from_aa( aa_from_name( random_element(
									option[ my_options::centroid_helixseq ]() ) ) ) );
				} else centroid_repeatseq.push_back( 'G' );
			}
		}

		// now mutate to centroid_repeatseq
		for ( Size i=1; i<= repeatlen; ++i ) {
			AA const aa( aa_from_oneletter_code( centroid_repeatseq[i-1] ) );
			for ( Size j=0; j< nrepeat; ++j ) {
				Size const seqpos( j*repeatlen + i );
				make_sequence_change( seqpos, aa, pose );
			}
		}

		/// bit of a hack: add a virtual residue at the end, fold starting from there...
		{
			//Residue const & lastrsd( pose.residue( nres_protein ) );
			ResidueOP vrtrsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
			pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...

			kinematics::FoldTree f( pose.total_residue() );
			f.reorder( pose.total_residue() );
			pose.fold_tree( f );
		}

		/// I think that these terminus types may be lost during replace residue calls?
		add_lower_terminus_type_to_pose_residue( pose, 1 );
		//add_upper_terminus_type_to_pose_residue( pose, nres_protein );
		pose.conformation().insert_chain_ending( nres_protein );

		{
			// Size const turn1_begin( helix1_len+1 ), turn1_end( helix1_len+turn1_len ),
			// 	turn2_begin( helix1_len+turn1_len+helix2_len+1 );
			for ( Size i=1; i<= nres_protein; ++i ) {
				pose.set_phi  ( i, init_phi   );
				pose.set_psi  ( i, init_psi   );
				pose.set_omega( i, init_omega );
				pose.set_secstruct( i, 'L' );
			}
		}

		/// pick some fragments
		FragLib fraglib;
		string fullss, fullbb;
		for ( Size i=1; i<= nrepeat; ++i ) {
			fullss += repeatss;
			fullbb += repeatbb;
		}
		TR.Trace << "pick_frags: " << nrepeat << ' ' << repeatlen << ' ' << repeatss << ' ' << repeatbb << ' ' <<
			fullss << ' ' << fullbb << endl;
		// default for resample_centroid_seqs_fragweight is 1.0
		//
		Real const seq_weight( resample_centroid_seqs ? 1.0 : 0.0 ), ss_weight( 10.0 ), bb_weight( 100.0 );
		Sizes const fragsizes( make_vector1( 3, 9 ) );
		Size const nfrags( 200 );
		kinematics::MoveMap mm;
		for ( Size i=1; i<= chain_end( 1, pose ); ++i ) mm.set_bb( i, true );

		Sizes const homs_to_exclude;
		devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, pose, mm, fullss, seq_weight, ss_weight,
			fraglib, homs_to_exclude, bb_weight, fullbb );


		{ // get rid of frags that violate the bb constraints
			Sizes const fragsizes( make_vector1( 3, 9 ) );
			Sizes const min_nns( make_vector1( 50, 25 ) );

			for ( Size si=1; si<= 2; ++si ) {
				Size const fragsize( fragsizes[si] );
				Size const min_nn( min_nns[si] );
				devel::blab::classic_frags::TorsionFragmentLibrary & lib( fraglib.library( fragsize ) );

				for ( Size fragpos=1; fragpos<= lib.size(); ++fragpos ) {
					for ( Size nn=lib[ fragpos ].size(); nn> min_nn; --nn ) {
						bool badfrag( false );

						devel::blab::classic_frags::TorsionFragment const frag( lib[ fragpos ][ nn ] );

						for ( Size k=1; k<= fragsize; ++k ) {
							//char const ss( frag.get_secstruct(k) );
							Real const phi( frag.get_torsion(k,1 ) ), psi( frag.get_torsion( k, 2 ) ), omega( frag.get_torsion(k,3 ) );
							char const bigbin( torsion2big_bin( phi, psi, omega ) );
							Size const repeatpos( ( ( fragpos +k-1 ) -1 )%repeatlen + 1 );
							char const desired_bigbin( repeatbb[ repeatpos-1 ] );
							if ( desired_bigbin != 'X' && desired_bigbin != bigbin ) {
								runtime_assert( string("ABGEO").find( desired_bigbin ) != string::npos );
								badfrag = true;
								break;
							}
						}

						if ( badfrag ) {
							TR.Trace << "delete_fragment: " << I(4,fragsize) << I(4,fragpos) << I(4,nn) << endl;
							lib[ fragpos ].erase( nn );
						}
					} // nn
					TR.Trace << "final_nn: " << I(4,fragsize) << I(4,fragpos) << I(4,(fragpos-1)%repeatlen+1) <<
						I(4,lib[ fragpos ].size()) << endl;
				} // fragpos
			} // fragsize
		} // scope


		/// setup symminfo?
		conformation::symmetry::SymmetryInfo symminfo;

		setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo );

		/// switch to centroid
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID ); // was CENTROID_DNA

		/// now make symmetric
		pose::symmetry::make_symmetric_pose( pose, symminfo );



		Real target_twist( 360.0 / nrepeat_full ), target_rise( 0.0 );
		if ( force_biangles ) {
			target_twist = 180.0;
			target_rise = 1000.0; // signal not to score rise
		} else if ( force_triangles ) {
			target_twist = 120.0;
			target_rise = 1000.0; // signal not to score rise
		} else if ( force_twistrise ) {
			target_twist = option[ my_options::target_twist ];
			target_rise  = option[ my_options::target_rise ];
		}

		{
			Real const twistrise_energy_weight( 1.0 );
			use_linear_twist_rise_energy = true; /// HACKY GLOBAL VAR
			Sizes fragseq_poslist;
			bool use_twistrise_energy( true ); // default
			simple_fold_abinitio( fraglib, fragseq_poslist, pose, use_twistrise_energy, target_twist, target_rise, 0.0,
				twistrise_energy_weight );
		}

		//basic::prof_show_oneliner();
		Real const centroid_simtime( ( (double) clock() - sim_starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

		Real const final_centroid_score( (*score3)( pose ) );
		use_linear_twist_rise_energy = false; /// HACKY GLOBAL VAR

		/// get handedness
		Real handedness(0.0);
		{
			Vectors coords;
			for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
			handedness = get_chirality( repeatlen, coords );
		}
		Real const centroid_handedness( handedness );

		/// compute transform: twist, axis,
		Real twist, rise, min_radius, max_radius, com_radius, depth;
		compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius, depth );
		Real const centroid_twist( twist ), centroid_rise( rise ), centroid_radius( min_radius );


		Real final_score( final_centroid_score ), relax_simtime( 0.0 ), res_score(0);

		bool relax_this_decoy( false );

		// Size const n2cdistpos( repeatlen*(nrepeat-2) ); // hacking
		// Real n2cdist( pose.residue(1).xyz("N").distance( pose.residue(n2cdistpos).xyz("C") ) ),
		// 	centroid_n2cdist( n2cdist );
		if ( force_biangles ) {
			Real const max_twist_dev(10.0);
			relax_this_decoy = ( fabs( twist - target_twist ) < max_twist_dev );
		} else if ( force_triangles ) {
			Real const max_twist_dev(5.0);
			relax_this_decoy = ( fabs( twist - target_twist ) < max_twist_dev );
		} else if ( force_twistrise ) {
			Real const max_twist_dev(5.0), max_rise_dev( 2.0 );
			relax_this_decoy = ( fabs( twist - target_twist ) < max_twist_dev ) &&
				( fabs( rise - target_rise ) < max_rise_dev );
		} else {
			/// old way: ( rise > -2 && rise < 2 && twist > 0.9 * target_twist && twist < 1.1 * target_twist );
			Real const max_twist_dev(10.0);
			relax_this_decoy = ( rise > -2 && rise < 2 && fabs( twist - target_twist ) < max_twist_dev );
		}

		if ( ( expected_hand == 'R' && centroid_handedness < 0 ) || ( expected_hand == 'L' && centroid_handedness > 0 ) ) {
			relax_this_decoy = false;
		}
		cout << "relax_this_decoy: final " << relax_this_decoy << endl;

		//Size nrepeat_this_decoy( nrepeat );

		ostringstream relax_status; // some scores only get written if we relaxed this decoy

		string fullsspred( "-" );
		Real searchtime(0);
		bool refold_this_decoy( true );
		if ( relax_this_decoy || dry_run() ) {
			// not sure we need this here... but make sure fullatom donut scoring works OK
			symminfo_hack::nrepeat_ = nrepeat_full;
			symminfo_hack::repeatlen_ = repeatlen;

			searchtime = ((double) clock() - search_starttime )/( CLOCKS_PER_SEC*60 ); // records time other than in relax,refold

			// convert to fullatom
			conformation::symmetry::SymmetryInfo symminfo( *pose::symmetry::symmetry_info( pose ) );
			pose::symmetry::make_asymmetric_pose( pose );
			devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
			pose::symmetry::make_symmetric_pose( pose, symminfo );

			clock_t const relax_starttime = clock();

			SequenceConstraints sequence_constraints;

			if ( option[ my_options::force_helix_capping ] ) {
				Sizes ncaps( 1, repeatlen );
				Sizes ccaps( 1, helix1_len+1 );
				ncaps.push_back( helix1_len + turn1_len );
				ccaps.push_back( helix1_len + turn1_len + helix2_len+1 );
				for ( Size repeatpos=1; repeatpos<= repeatlen; ++repeatpos ) {
					string seqcst;
					if ( has_element( repeatpos, ccaps ) && repeatbb[ repeatpos-1 ] == 'G' ) seqcst = "G";
					if ( has_element( repeatpos, ncaps ) && repeatbb[ repeatpos-1 ] == 'B' ) seqcst = "STND";
					if ( seqcst.empty() || sequence_constraints.count( repeatpos ) ) continue; // earlier constraints have prio
					TR.Trace << "force_helix_capping: " << repeatpos << ' ' << seqcst << endl;
					for ( Size k=0; k<nrepeat; ++k ) {
						Size const seqpos( k*repeatlen + repeatpos );
						sequence_constraints[ seqpos ] = seqcst;
					}
				}
			}


			if ( option[ my_options::layer_design ] ) {
				///
				push_sequence_constraints_to_base_repeat( pose, sequence_constraints );
				string fullss;
				for ( Size i=1; i<= nrepeat; ++i ) fullss += target_repeatss;
				fullss.push_back( 'X' );
				runtime_assert( fullss.size() == pose.total_residue() );
				add_layer_design_constraints( pose, fullss, sequence_constraints );
			}

			{
				bool const add_jump_flex( false ), use_atom_pair_constraints( false ), use_coordinate_constraints( false );
				//bool const use_twistrise_energy( option[ my_options::fullatom_twistrise_energy ] );
				symmetric_design_and_relax( repeatlen, nrepeat_full, option[ my_options::use_softrep_for_design ],
					use_donut_energy, option[ my_options::design_cycles ], sequence_constraints, pose,
					add_jump_flex, use_atom_pair_constraints, use_coordinate_constraints,
					use_fullatom_twistrise_energy, target_twist, target_rise );
			}

			/// restore these since they are zeroed by symmetric_design_and_relax
			symminfo_hack::nrepeat_ = nrepeat_full;
			symminfo_hack::repeatlen_ = repeatlen;


			relax_simtime = ((double) clock() - relax_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

			final_score = (*fa_scorefxn)( pose );
			res_score = final_score / repeatlen;

			///
			bools subset( pose.total_residue(), false );
			Size const base_repeat_offset( repeatlen*(base_repeat-1));
			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;

			Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa, buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score;
			compute_sasa_scores_for_subset_slow( 5, subset, pose, sasapack_score, norme_score, normsasa_score,
				exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa,
				buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score );

			Real n_buried_unsatisfied_donors_bb_per_repeat, n_buried_unsatisfied_acceptors_bb_per_repeat,
				n_buried_unsatisfied_donors_sc_per_repeat, n_buried_unsatisfied_acceptors_sc_per_repeat;
			{ // compute unsats over full pose
				Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
					n_buried_unsatisfied_acceptors_sc;
				bools subset_full( pose.total_residue(), true ); subset_full[ pose.total_residue() ] = false;
				get_buried_unsatisfied_counts_real_slow( subset_full, pose,
					n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc,
					n_buried_unsatisfied_acceptors_bb, n_buried_unsatisfied_acceptors_sc );

				n_buried_unsatisfied_donors_bb_per_repeat = Real( n_buried_unsatisfied_donors_bb )/nrepeat;
				n_buried_unsatisfied_donors_sc_per_repeat = Real( n_buried_unsatisfied_donors_sc )/nrepeat;
				n_buried_unsatisfied_acceptors_bb_per_repeat = Real( n_buried_unsatisfied_acceptors_bb )/nrepeat;
				n_buried_unsatisfied_acceptors_sc_per_repeat = Real( n_buried_unsatisfied_acceptors_sc )/nrepeat;
			}

			if ( option[ my_options::pick_decoys ] ) {
				Real const exposed_polar_fraction( exposed_polar_sasa / ( exposed_polar_sasa + exposed_nonpolar_sasa ) );
				if ( exposed_polar_fraction > 0.55 || exposed_polar_fraction < 0.35 ) refold_this_decoy = false;
				if ( sasapack_score > 0 ) refold_this_decoy = false;

				if ( refold_this_decoy ) {
					Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
						n_buried_unsatisfied_acceptors_sc;
					get_buried_unsatisfied_counts_real_slow( subset, pose,
						n_buried_unsatisfied_donors_bb,
						n_buried_unsatisfied_donors_sc,
						n_buried_unsatisfied_acceptors_bb,
						n_buried_unsatisfied_acceptors_sc );
					if ( ( n_buried_unsatisfied_donors_bb + n_buried_unsatisfied_donors_sc > 1 ) ||
							( n_buried_unsatisfied_acceptors_bb + n_buried_unsatisfied_acceptors_sc > 0 ) ) refold_this_decoy = false;
				}
			}

			//n2cdist = pose.residue(1).xyz("N").distance( pose.residue(chain_end(1,pose)).xyz("C") );

			/// compute psipred ss using single sequence
			string repeatsspred( "-" );
			Real sspred_match(0.0), sspred_turn1(0.0), sspred_turn2( 0.0 ), sspred_helix1(0.0), sspred_helix2(0.0);
			{ // run psipred to re-predict the secondary structure
				vector1< Reals > pred_HEL;
				string const sequence( pose.sequence().substr(0,chain_end(1,pose) ) );
				//string const secstruct( pose.secstruct().substr(0,chain_end(1,pose) ) );
				string sspred;
				string const repeatss( pose.secstruct().substr(0,repeatlen) ); // current, not the target

				run_psipred( sequence, sspred, pred_HEL );

				if ( sspred.size() == sequence.size() ) { // success
					fullsspred = sspred;
					runtime_assert( sspred.size() == nrepeat * repeatlen );
					repeatsspred = sspred.substr( (base_repeat-1)*repeatlen, repeatlen );
					Size const turn1_begin( helix1_len+1 ), turn1_end( helix1_len + turn1_len ),
						turn2_begin( helix1_len+turn1_len+helix2_len+1 );
					for ( Size i=1; i<= sspred.size(); ++i ) {
						Size const rpos( (i-1)%repeatlen+1 );
						char const ss( repeatss[ rpos-1 ] );
						Real const prediction( ss == 'H' ? pred_HEL[i][1] : ( ss == 'E' ? pred_HEL[i][2] : pred_HEL[i][3] ) );
						sspred_match += prediction;
						if ( rpos < turn1_begin                     ) sspred_helix1 += prediction;
						if ( rpos > turn1_end && rpos < turn2_begin ) sspred_helix2 += prediction;
						if ( rpos >= turn1_begin-1 && rpos <= turn1_end+1 ) sspred_turn1 += prediction;
						if ( rpos == 1 || rpos >= turn2_begin-1           ) sspred_turn2 += prediction;
					}
					sspred_match  /= (  repeatlen * nrepeat );
					sspred_helix1 /= ( helix1_len * nrepeat );
					sspred_helix2 /= ( helix2_len * nrepeat );
					sspred_turn1  /= ( ( turn1_len+2) * nrepeat );
					sspred_turn2  /= ( ( turn2_len+2) * nrepeat );
				}
			}

			Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);
			if ( helix1_len ) {
				runtime_assert( helix1_len + turn1_len + helix2_len + turn2_len == repeatlen );

				compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1_len, helix2_len,
					helix1_dist, helix1_twist, helix2_dist, helix2_twist );
			}

			Real const
				inner_helix_dist( min( helix1_dist, helix2_dist ) ),
				outer_helix_dist( max( helix1_dist, helix2_dist ) );
			Real const inner_helix_twist( helix1_dist < helix2_dist ? helix1_twist : helix2_twist );
			Real const outer_helix_twist( helix1_dist < helix2_dist ? helix2_twist : helix1_twist );

			Size const inner_helix_len( helix1_dist < helix2_dist ? helix1_len : helix2_len );
			Size const outer_helix_len( helix1_dist < helix2_dist ? helix2_len : helix1_len );
			string const inner_helix_turn( helix1_dist < helix2_dist ? turn1 : turn2 );
			string const outer_helix_turn( helix1_dist < helix2_dist ? turn2 : turn1 );

			string const buried_unsatisfied_string( get_buried_unsatisfied_string( pose ) );

			if ( option[ my_options::max_sasapack_score ].user() && sasapack_score > option[ my_options::max_sasapack_score]){
				refold_this_decoy = false;
			}

			relax_status <<
				" sasapack_score: " << F(9,3,sasapack_score) <<
				" norme_score: " << F(9,3,norme_score) <<
				" normsasa_score: " << F(9,3,normsasa_score) <<
				" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
				" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
				" exposed_polar_fraction: " << F(9,3,exposed_polar_sasa/(exposed_polar_sasa+exposed_nonpolar_sasa))<<
				" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
				" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
				" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc)<<
				" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc)<<
				" packstat_score: " << F(9,3,packstat_score) <<
				" res_score: " << F(9,3,res_score ) <<
				" n_buried_unsatisfied_donors_bb_per_repeat: " << F(9,3,n_buried_unsatisfied_donors_bb_per_repeat) <<
				" n_buried_unsatisfied_donors_sc_per_repeat: " << F(9,3,n_buried_unsatisfied_donors_sc_per_repeat) <<
				" n_buried_unsatisfied_acceptors_bb_per_repeat: " <<F(9,3,n_buried_unsatisfied_acceptors_bb_per_repeat)<<
				" n_buried_unsatisfied_acceptors_sc_per_repeat: " <<F(9,3,n_buried_unsatisfied_acceptors_sc_per_repeat)<<
				" repeatsspred: " << repeatsspred <<
				" sspred_match: " << F(9,3,sspred_match) <<
				" sspred_helix1: " << F(9,3,sspred_helix1) <<
				" sspred_helix2: " << F(9,3,sspred_helix2) <<
				" sspred_turn1: " << F(9,3,sspred_turn1) <<
				" sspred_turn2: " << F(9,3,sspred_turn2) <<
				" helix1_dist: " << F(9,3,helix1_dist) <<
				" helix1_twist: " << F(9,3,helix1_twist) <<
				" helix2_dist: " << F(9,3,helix2_dist) <<
				" helix2_twist: " << F(9,3,helix2_twist) <<
				" inner_helix_dist: " << F(9,3,inner_helix_dist) <<
				" inner_helix_twist: " << F(9,3,inner_helix_twist) <<
				" outer_helix_dist: " << F(9,3,outer_helix_dist) <<
				" outer_helix_twist: " << F(9,3,outer_helix_twist) <<
				" inner_helix_len: " << I(4,inner_helix_len) <<
				" inner_helix_turn: " << inner_helix_turn <<
				" outer_helix_len: " << I(4,outer_helix_len) <<
				" outer_helix_turn: " << outer_helix_turn <<
				' ' << buried_unsatisfied_string <<
				' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() );


			/// RE-compute transform: twist, axis, after relax
			compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius );
			{
				Vectors coords;
				for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
				handedness = get_chirality( repeatlen, coords );
			}
		}

		string const repeatseq( pose.sequence().substr(0,repeatlen) );




		Real ala_fraction(0);
		foreach_ ( char const s, repeatseq ) { ala_fraction += ( s == 'A' ); }
		ala_fraction /= repeatseq.size();

		if ( option[ my_options::max_ala_fraction ].user() && ala_fraction - 1e-6 > option[ my_options::max_ala_fraction]){
			refold_this_decoy = false;
		}

		if ( use_donut_energy ) { // problem with donut energy!
			Real const max_twist_dev( 5.0 );
			refold_this_decoy = refold_this_decoy && ( rise > -2 && rise < 2 && fabs( twist - target_twist ) < max_twist_dev );
		} else if ( force_biangles || force_triangles ) {
			Real const max_twist_dev(2.5);
			refold_this_decoy = refold_this_decoy && ( fabs( twist - target_twist ) < max_twist_dev );
		} else if ( force_twistrise ) {
			Real const max_twist_dev(2.0), max_rise_dev( 1.0 );
			refold_this_decoy = ( fabs( twist - target_twist ) < max_twist_dev ) &&
				( fabs( rise - target_rise ) < max_rise_dev );
		}


		Size const base_repeat_offset( repeatlen * ( base_repeat-1 ) );
		repeatbb = torsion2big_bin_string( base_repeat_offset+1, base_repeat_offset+repeatlen, pose );
		repeatss = pose.secstruct().substr(0,repeatlen);

		Vectors ca_coords;
		for ( Size i=1; i<= nres_protein; ++i ) ca_coords.push_back( pose.residue(i).xyz("CA") );
		string helix_stats;
		if ( force_biangles || force_triangles ) {
			helix_stats = get_single_helix_stats( nrepeat, repeatbb, repeatseq, ca_coords );
		} else {
			helix_stats = get_double_helix_stats( nrepeat, repeatbb, repeatseq, ca_coords );
		}

		cout << "refold_this_decoy: final " << refold_this_decoy << ' ' << relax_this_decoy <<
			" ala_fraction: " << ala_fraction << ' ' << bbtag << endl;

		string const simfile( outdir+"unbound_frag_"+string_of( repeatlen/4 )+".work" );

		/// how do we want to score filter?? by bbtags? seems too granular.
		bool const passed_score_filter( relax_this_decoy && refold_this_decoy && // short-circuit evaluation
			append_score_to_scorefile_and_filter( "tmp", res_score, // was final_score (doh!)
				score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );

		ostringstream refolding_status;

		/// try refolding de novo
		Size const n_refold( option[ my_options::n_refold ] ), n_refold_sspred( option[ my_options::n_refold_sspred ] );
		Real refolding_simtime(0);
		if ( refold_this_decoy && n_refold+n_refold_sspred > 0 && ( passed_score_filter || dry_run() ) ) {
			clock_t const refold_starttime = clock();
			for ( Size RR=1; RR<= 2; ++RR ) {
				Real min_refold_rmsd(1e6);
				Size n_under_4(0);
				if ( RR == 2 && fullsspred.size() != nrepeat * repeatlen ) continue;
				if ( RR==1 ) refolding_status << " n_refold: " << I(3,n_refold);
				else         refolding_status << " n_refold_sspred: " << I(3,n_refold_sspred);
				Size const rr_end( RR==1 ? n_refold : n_refold_sspred );
				for ( Size rr=1; rr<= rr_end; ++rr ) {
					Pose protpose;
					string const repss( RR == 1 ? repeatss : fullsspred );
					Size const nrepeat_refold( nrepeat );
					Size const base_repeat_refold( base_repeat );
					refold_repeat_pose( repeatseq, repss, nrepeat_refold, base_repeat_refold, protpose );
					/// calc rmsd
					Real rmsd( 0.0 ), rmsd_single_repeat( 0.0 );
					{
						using namespace core::id;
						AtomID_Map< AtomID > atom_map;
						initialize_atomid_map( atom_map, protpose, id::GLOBAL_BOGUS_ATOM_ID );
						for ( Size i=1; i<= repeatlen; ++i ) {
							atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
								AtomID( pose.residue(i).atom_index("CA"),i);
						}
						rmsd_single_repeat = rmsd_by_mapping( protpose, pose, atom_map );
						Size nrepeat_for_refolding_rmsd( min( nrepeat, nrepeat_refold ) );
						if ( option[ my_options::nrepeat_for_refolding_rmsd ].user() ) {
							nrepeat_for_refolding_rmsd = option[ my_options::nrepeat_for_refolding_rmsd ];
						}
						Size const nres_for_rmsd( min( chain_end( 1, protpose), nrepeat_for_refolding_rmsd * repeatlen ) );
						for ( Size i=1; i<= nres_for_rmsd; ++i ) {
							atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
								AtomID( pose.residue(i).atom_index("CA"),i);
						}
						rmsd = rmsd_by_mapping( protpose, pose, atom_map );
						min_refold_rmsd = min( min_refold_rmsd, rmsd );
						if ( rmsd < 4 ) ++n_under_4;
					}

					Real refold_handedness( 0.0 );
					{
						Vectors coords;
						for ( Size i=1; i<= chain_end( 1, protpose); ++i ) coords.push_back( protpose.residue(i).xyz("CA") );
						refold_handedness = get_chirality( repeatlen, coords );
					}
					Real tmptwist, tmprise, tmpmin_radius, tmpmax_radius, tmpcom_radius;
					compute_repeat_params( protpose, tmptwist, tmprise, tmpmin_radius, tmpmax_radius, tmpcom_radius );
					refolding_status << " R " << rr << F(9,3,rmsd) << F(9,3,rmsd_single_repeat) << F(9,3,refold_handedness) <<
						F(9,3,tmprise) << F(9,3,tmptwist);
				}
				refolding_status << " min_refold" << RR << "_rmsd: "<< F(9,3,min_refold_rmsd) <<
					" n_refold" << RR << "_under_4: " << n_under_4;

			}
			refolding_simtime = ((double) clock() - refold_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
		}


		string const outfilename( output_tag() + "unbound_frag_"+bbtag+"_"+lead_zero_string_of( n,4 )+".pdb" );

		Real const simtime = ((double) clock() - sim_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

		ostringstream out;
		out << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" handedness: " << F(9,3,handedness) <<
			" min_radius: " << F(9,3,min_radius ) <<
			" rise: " << F(9,3,rise ) <<
			" twist: " << F(9,3,twist) <<
			" final_centroid_score: "<< F(9,3,final_centroid_score) <<
			" passed_score_filter: " << passed_score_filter <<
			" relaxed: " << relax_this_decoy <<
			" min_radius: " << F(9,3,min_radius )<<
			" max_radius: " << F(9,3,max_radius )<<
			" com_radius: " << F(9,3,com_radius )<<
			" centroid_handedness: " << F(9,3,centroid_handedness) <<
			" centroid_radius: " << F(9,3,centroid_radius ) <<
			" centroid_rise: " << F(9,3,centroid_rise ) <<
			" centroid_twist: " << F(9,3,centroid_twist ) <<
			" repeatseq: " << repeatseq <<
			" repeatbb: " << repeatbb <<
			" repeatss: " << repeatss <<
			" turn1: " << turn1 <<
			" turn2: " << turn2 <<
			" helix1_len: " << helix1_len <<
			" helix2_len: " << helix2_len <<
			" repeatlen: " << repeatlen <<
			" nrepeat: " << nrepeat <<
			//" n2cdist: " << F(9,3,n2cdist) <<
			//" centroid_n2cdist: " << F(9,3,centroid_n2cdist) <<
			" target_repeatbb: " << target_repeatbb <<
			" target_repeatss: " << target_repeatss <<
			" " << relax_status.str() <<
			" ala_fraction: " << F(9,3,ala_fraction) <<
			" bbtag: " << bbtag <<
			" " << refolding_status.str() <<
			" centroid_simtime: " << F(9,3,centroid_simtime) <<
			" relax_simtime: " << F(9,3,relax_simtime) <<
			" refolding_simtime: " << F(9,3,refolding_simtime) <<
			" searchtime: " << F(9,3,searchtime) <<
			" simtime: " << F(9,3,simtime) <<
			' ' << helix_stats << '\n';

		string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

		// if ( passed_score_filter || dry_run() ) {
		// 	pose.dump_pdb( outfilename );
		// 	run_command("gzip "+outfilename ); // NEW NEW NEW
		// } // passed_score_filter
		if ( passed_score_filter && pdbfilename.size()>0 ) {
			append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
			run_command("gzip "+pdbfilename );
		}

		fflush( stdout );
		check_simtime();

		search_starttime = clock(); // restart the search timer
	}




}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
add_helix_stats_test()
{
	strings const files( start_files() );

	foreach_( string filename, files ) {
		string scoreline;
		{ //
			utility::io::izstream data( filename );
			string line;
			while ( getline( data, line ) ) {
				if ( line.substr(0,16) == "REMARK SCORELINE" ) {
					scoreline = line.substr(17); // starts "final_scores"
					if ( filename.find(".pdbs") != string::npos ) {
						/// make it look like we grepped from .out files:
						scoreline = filename.substr(0,filename.find(".pdbs"))+".out:"+scoreline;
					}
					break;
				}
			}
			data.close();
		}

		if ( scoreline.empty() ) {
			cout << "SKIP missing scoreline: " << filename << endl;
			continue;
		}


		// parse the scoreline
		Size nrepeat(0);
		string repeatbb, repeatseq;
		{
			strings const l( split_to_vector1( scoreline ) );
			for ( Size i=1; i< l.size(); ++i ) {
				string const tag( l[i] ), val( l[i+1] );
				if ( tag == "nrepeat:" ) nrepeat = int_of(val);
				else if ( tag == "repeatbb:" ) repeatbb = val;
				else if ( tag == "repeatseq:" ) repeatseq = val;
			}
		}

		if ( repeatbb.find_first_not_of("A") == string::npos ) { // all A
			cout << "SKIP bad repeatbb: " << repeatbb << ' ' << filename << endl;
			continue;
		}

		Vectors const coords( read_CA_coords_from_file( filename ) );

		string helix_stats;
		if ( option[ my_options::force_biangles ] ||option[ my_options::force_triangles ] ) {
			helix_stats = get_single_helix_stats( nrepeat, repeatbb, repeatseq, coords );
		} else {
			helix_stats = get_double_helix_stats( nrepeat, repeatbb, repeatseq, coords );
		}

		if ( helix_stats.empty() ) {
			cout << "SKIP bad stats: " << helix_stats << ' ' << filename << endl;
			continue;
		}

		cout << scoreline << " filename: " << filename << ' ' << helix_stats << '\n';
	}
}

///////////////////////////////////////////////////////////////////////////////
//
void
design_stats_test()
{

	ScoreFunctionOP fa_scorefxn( get_score_function_from_command_line() );

	strings const files( start_files() );

	foreach_( string const filename, files ) {

		Pose pose;
		cout << "reading: " << filename << endl;

		pose_from_pdb( pose, filename );

		set_ss_from_dssp( filename, pose ); // phil_io.hh

		bools const residue_changed( remove_funny_variants( pose, filename ) ); // hopefully not doing anything...

		Size const nres( pose.total_residue() );

		for ( Size i=1; i<= nres;++i ) { runtime_assert( pose.residue(i).is_protein() ); }

		Real const total_score( (*fa_scorefxn)( pose ) ), res_score( total_score / nres );

		/// sasa stuff
		///
		bools subset( pose.total_residue(), true );

		Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
			buried_polar_sasa, buried_nonpolar_sasa, buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
			packstat_score;
		compute_sasa_scores_for_subset_slow( 5, subset, pose, sasapack_score, norme_score, normsasa_score,
			exposed_polar_sasa, exposed_nonpolar_sasa,
			buried_polar_sasa, buried_nonpolar_sasa,
			buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
			packstat_score );

		Real const exposed_polar_fraction( exposed_polar_sasa / ( exposed_polar_sasa + exposed_nonpolar_sasa ) );

		Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
			n_buried_unsatisfied_acceptors_sc;
		get_buried_unsatisfied_counts_real_slow( subset, pose,
			n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc,
			n_buried_unsatisfied_acceptors_bb, n_buried_unsatisfied_acceptors_sc );

		/// normalize these to 100 rsds length
		Real const scale( 100.0 / nres );
		n_buried_unsatisfied_acceptors_bb *= scale;
		n_buried_unsatisfied_acceptors_sc *= scale;
		n_buried_unsatisfied_donors_bb *= scale;
		n_buried_unsatisfied_donors_sc *= scale;

		Real ala_fraction(0), H_fraction(0), E_fraction(0), L_fraction(0);
		foreach_( char const s, pose.sequence() ) {
			ala_fraction += (s=='A');
		}
		foreach_( char const s, pose.secstruct() ) {
			H_fraction += (s=='H');
			E_fraction += (s=='E');
			L_fraction += (s=='L');
		}

		ala_fraction /= nres;
		H_fraction /= nres;
		E_fraction /= nres;
		L_fraction /= nres;

		cout << "design_stats: " << filename <<
			" nres: " << nres <<
			" res_score: " << F(9,3,res_score ) <<
			" ala_fraction: " << F(9,3,ala_fraction) <<
			" H_fraction: " << F(9,3,H_fraction) <<
			" E_fraction: " << F(9,3,E_fraction) <<
			" L_fraction: " << F(9,3,L_fraction) <<
			" sasapack_score: " << F(9,3,sasapack_score) <<
			" norme_score: " << F(9,3,norme_score) <<
			" normsasa_score: " << F(9,3,normsasa_score) <<
			" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
			" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
			" exposed_polar_fraction: " << F(9,3,exposed_polar_fraction) <<
			" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
			" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
			" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc)<<
			" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc)<<
			" packstat_score: " << F(9,3,packstat_score) <<
			" n_buried_unsatisfied_donors_bb_per_100: " << F(9,3,n_buried_unsatisfied_donors_bb) <<
			" n_buried_unsatisfied_donors_sc_per_100: " << F(9,3,n_buried_unsatisfied_donors_sc) <<
			" n_buried_unsatisfied_acceptors_bb_per_100: " <<F(9,3,n_buried_unsatisfied_acceptors_bb)<<
			" n_buried_unsatisfied_acceptors_sc_per_100: " <<F(9,3,n_buried_unsatisfied_acceptors_sc)<<
			" secstruct: " << pose.secstruct() <<
			endl;

	} // filename

}


///////////////////////////////////////////////////////////////////////////////
void
extend_design_pdbs_test()
{
	strings const files( start_files() );

	foreach_( string const filename, files ) {

		Pose pose;
		pose_from_pdb( pose, filename );

		while ( !pose.residue( pose.total_residue() ).is_protein() ) {
			pose.conformation().delete_residue_slow( pose.total_residue() );
		}
		remove_upper_terminus_type_from_pose_residue( pose, pose.total_residue() );

		Size const nres( pose.total_residue() ), repeatlen( deduce_repeatlen( pose.sequence() ) ),
			nrepeat( nres/repeatlen );


		Vectors fixcoords, movcoords;

		runtime_assert( nrepeat>=3 );
		Vector com(0,0,0);
		Size const nsup( nres-repeatlen );
		for ( Size i=1; i<= nsup; ++i ) {
			Size const movpos( i ), fixpos( repeatlen+i );
			fixcoords.push_back( pose.residue(fixpos).xyz("CA"));
			movcoords.push_back( pose.residue(movpos).xyz("CA"));
			com += movcoords.back();
		}
		com /= nsup;

		Stub const stub1( movcoords[1], movcoords[2], movcoords[3] );
		superimpose_coords( fixcoords, movcoords );
		Stub const stub2( movcoords[1], movcoords[2], movcoords[3] );

		Real theta;
		Vector axis, center, t;

		get_stub_transform_data( stub1, stub2, center, axis, t, theta );
		runtime_assert( fabs( axis.length_squared() - 1.0 )<1e-2 ); // confirm that axis is normal vector

		// figure out full nrepeat
		Size const nrepeat_full( int( 0.5 + 360. / degrees(theta) ) );

		TR.Trace << "extend_design_pdbs_test: expected_theta: " << 360.0/nrepeat_full <<
			" actual_theta: " << degrees( theta ) << endl;

		center += ( com-center ).dot( axis ) * axis; // slide center along axis to align with com of repeat2
		Vector const new_axis( 0,0,1 ), new_center(0,0,0); // z-axis

		Real const rotangle( acos( axis.dot( new_axis ) ) );
		Vector const rotaxis( axis.cross( new_axis ) );

		numeric::xyzMatrix< Real > R( numeric::rotation_matrix( rotaxis, rotangle ) );
		runtime_assert( is_small( new_axis.distance_squared( R*axis ) ) );

		Vector const v( new_center - R*center );
		runtime_assert( is_small( new_center.distance_squared( R*center + v ) ) );

		pose.apply_transform_Rx_plus_v( R, v );


		// num extra repeats
		Size const num_extra_repeats( nrepeat_full - nrepeat );
		Pose spin_pose( pose );
		for ( Size i=1; i<= num_extra_repeats; ++i ) {
			numeric::xyzMatrix< Real > R2( numeric::rotation_matrix( new_axis, 2 * pi / nrepeat_full ) );
			spin_pose.apply_transform_Rx_plus_v( R2, Vector(0.,0.,0.) );
			for ( Size j=1; j<= repeatlen; ++j ) {
				pose.append_residue_by_bond( spin_pose.residue( nres-repeatlen+j ) );
			}
		}

		runtime_assert( pose.total_residue() == nrepeat_full * repeatlen );

		pose.dump_pdb( filename+".extended.pdb" );
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// don't simulate the full number of repeats
//
void
triple_unbound_frag_test()
{
	runtime_assert( !option[ my_options::local_simfile ] );
	runtime_assert( option[ my_options::unfolded_sasas ].user() ); // make sure have datafile for buried sasa
	runtime_assert( option[ my_options::resample_centroid_seqs ].user() ||
		option[ my_options::centroid_helixseq ].user() );

	//runtime_assert( option[ my_options::local_simfile ] );
	bool const resample_centroid_seqs( option[ my_options::resample_centroid_seqs ].user() );
	map< string, strings > bbtag2centroid_seq;
	if ( resample_centroid_seqs ) {
		strings l( option[ my_options::resample_centroid_seqs ] );
		string bbtag;
		while ( !l.empty() ) {
			if ( is_int( l.front().substr(0,1) ) ) {
				bbtag = l.front();
			} else {
				runtime_assert( bbtag.size() );
				TR.Trace << "bbtag2centroid_seq: " << bbtag << ' ' << l.front() << endl;
				bbtag2centroid_seq[ bbtag ].push_back( l.front() );
			}
			l.erase( l.begin() );
		}
	}


	//Real const min_helix2_len( 5 );

	Real const score_filter_acceptance_rate( option[ my_options::score_filter_acceptance_rate ] );
	bool const score_filter_pass_early( option[ my_options::score_filter_pass_early ] );

	ScoreFunctionOP score3( ScoreFunctionFactory::create_score_function( "score3.wts" ) );
 	ScoreFunctionOP fa_scorefxn(0);
	if ( option[ OptionKeys::dna::specificity::score_function ].user() ) {
		fa_scorefxn = get_score_function_from_command_line();
	} else {
		fa_scorefxn = ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag );
	}
	//ScoreFunctionOP fa_scorefxn( ScoreFunctionFactory::create_score_function( "score12prime" ) ); // 09/16/15
	runtime_assert( pose::symmetry::is_symmetric( *score3 ) );

	adjust_ref_weights_from_command_line( *fa_scorefxn );

	Real const donut_energy_weight( option[ my_options::donut_energy_weight ] );
	bool const use_donut_energy( donut_energy_weight>1e-3 );
	// bool const use_fullatom_twistrise_energy( force_biangles || force_triangles || force_twistrise );
	// if ( use_fullatom_twistrise_energy ) {
	// 	runtime_assert( !use_donut_energy );
	// } else {
	// 	runtime_assert( use_donut_energy );
	// }
	DonutWholeEnergy donut_whole_energy;
	Donut1B_Energy donut_1b_energy;
	donut_whole_energy.pose_is_subset( true );
	donut_1b_energy.pose_is_subset( true );
	if ( use_donut_energy ) {
		fa_scorefxn->add_extra_method( donut_whole, donut_energy_weight, donut_whole_energy );
		fa_scorefxn->add_extra_method( donut_1b   , donut_energy_weight, donut_1b_energy    );
	}


	clock_t search_starttime( clock() );

	Size const nstruct( option[ OptionKeys::out::nstruct ].user() ? option[ OptionKeys::out::nstruct ]() : 100000 );

	string const outdir( shared_output_dir_create_if_necessary() ); // has trailing /

	for ( Size n=1; n<= nstruct; ++n ) {
		// undo any fiddling
		symminfo_hack::nrepeat_ = 0;
		symminfo_hack::repeatlen_ = 0;

		clock_t sim_starttime = clock();

		Size const base_repeat( option[ my_options::base_repeat ].user() ? option[ my_options::base_repeat ]: 3 ),
			nrepeat( option[ my_options::nrepeat_sim ].user() ? option[ my_options::nrepeat_sim ] : 5 );

		Pose pose;

		// create a symmetric pose
		Size helix1_len(0), helix2_len(0), helix3_len(0);
		string turn1, turn2, turn3;
		Size repeatlen(0 );
		// char expected_hand( '-' );
		//string resample_bb_tag;
		//string tal_twist_type;

		Size nrepeat_full(0);


		if ( option[ my_options::resample_bbs ].user() ) { // NOTE: NO convention on which is the inner helix!!!
			// EXAMPLE: " -resample_bbs  8-9,13,GB,13,GB,L  9,14-17,GBB,11,GBB,L "
			string const bb( random_element( option[ my_options::resample_bbs ]() ) );
			strings l;
			if ( bb.find(',') != string::npos ) l = split_to_vector1( bb, "," );
			else l = split_to_vector1( bb, "." );
			runtime_assert( l.size() == 7 );
			nrepeat_full = get_int_from_resample_bb_tag( l[1] );
			helix1_len   = get_int_from_resample_bb_tag( l[2] );
			helix2_len   = get_int_from_resample_bb_tag( l[4] );
			helix3_len   = get_int_from_resample_bb_tag( l[6] );
			turn1 = l[3];
			turn2 = l[5];
			turn3 = l[7];
			repeatlen = helix1_len + helix2_len + helix3_len + \
				get_turnlen( turn1 ) + get_turnlen( turn2 ) + get_turnlen( turn3 );
		} else { //if ( option[ my_options::turn_pairs ].user() ) {
			nrepeat_full = random_element( option[ my_options::nrepeats ]() );
			// strings const tp( split_to_vector1( random_element( option[ my_options::turn_pairs ]() ), "." ) );
			// runtime_assert( tp.size() == 3 );
			turn1 = random_element( option[ my_options::turn1s ]() );
			turn2 = random_element( option[ my_options::turn2s ]() );
			turn3 = random_element( option[ my_options::turn3s ]() );
			helix1_len = random_element( option[ my_options::helix1_lens ]() );
			helix2_len = random_element( option[ my_options::helix2_lens ]() );
			helix3_len = random_element( option[ my_options::helix3_lens ]() );
			repeatlen = helix1_len + helix2_len + helix3_len + \
				get_turnlen( turn1 ) + get_turnlen( turn2 ) + get_turnlen( turn3 );
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// end setup
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		TR.Trace << "nrepeat: " << nrepeat << " repeatlen: " << repeatlen << ' ' << turn1 << ' ' << turn2 << ' ' <<
			turn2 << ' ' << helix1_len << ' ' << helix2_len << ' ' << helix3_len << endl;


		Size const turn1_len( get_turnlen( turn1 ) ), turn2_len( get_turnlen( turn2 ) ), turn3_len( get_turnlen( turn3 ) );
		string repeatbb(
			string( helix1_len, 'A' ) + ( turn1_len==0 ? string() : turn1 ) +
			string( helix2_len, 'A' ) + ( turn2_len==0 ? string() : turn2 ) +
			string( helix3_len, 'A' ) + ( turn3_len==0 ? string() : turn3 ) );
		string repeatss(
			string( helix1_len, 'H' ) + string( turn1_len, 'L' ) +
			string( helix2_len, 'H' ) + string( turn2_len, 'L' ) +
			string( helix3_len, 'H' ) + string( turn3_len, 'L' ) );

		string const target_repeatbb( repeatbb ), target_repeatss( repeatss );

		runtime_assert( repeatbb.size() == repeatss.size() && repeatbb.size() == repeatlen );

		Size nres_protein( nrepeat * repeatlen );
		for ( Size i=1; i<= nres_protein; ++i ) {
			pose.append_residue_by_bond( *get_vanilla_protein_residue( 'L' ), true ); // build_ideal_geometry
		}


		/// setup the centroid repeat sequence
		string const bbtag
			( string_of( nrepeat_full ) + "." +
				string_of( helix1_len ) + "." + turn1 + "." +
				string_of( helix2_len ) + "." + turn2 + "." +
				string_of( helix3_len ) + "." + turn3 );

		string centroid_repeatseq;

		if ( resample_centroid_seqs ) {
			runtime_assert( bbtag2centroid_seq.count( bbtag ) );
			centroid_repeatseq = random_element( bbtag2centroid_seq[ bbtag ] );
			runtime_assert( centroid_repeatseq.size() == repeatlen );
		} else {
			for ( Size i=1; i<= repeatlen; ++i ) {
				if ( repeatss[i-1] != 'L' ) { // stupid backwards compatibility...
					centroid_repeatseq.push_back( oneletter_code_from_aa( aa_from_name( random_element(
									option[ my_options::centroid_helixseq ]() ) ) ) );
				} else centroid_repeatseq.push_back( 'G' );
			}
		}

		// now mutate to centroid_repeatseq
		for ( Size i=1; i<= repeatlen; ++i ) {
			AA const aa( aa_from_oneletter_code( centroid_repeatseq[i-1] ) );
			for ( Size j=0; j< nrepeat; ++j ) {
				Size const seqpos( j*repeatlen + i );
				make_sequence_change( seqpos, aa, pose );
			}
		}

		/// bit of a hack: add a virtual residue at the end, fold starting from there...
		{
			//Residue const & lastrsd( pose.residue( nres_protein ) );
			ResidueOP vrtrsd
				( conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose()->name_map( "VRTBB" ) ) );
			pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...

			kinematics::FoldTree f( pose.total_residue() );
			f.reorder( pose.total_residue() );
			pose.fold_tree( f );
		}

		/// I think that these terminus types may be lost during replace residue calls?
		add_lower_terminus_type_to_pose_residue( pose, 1 );
		//add_upper_terminus_type_to_pose_residue( pose, nres_protein );
		pose.conformation().insert_chain_ending( nres_protein );

		{
			for ( Size i=1; i<= nres_protein; ++i ) {
				pose.set_phi  ( i, init_phi   );
				pose.set_psi  ( i, init_psi   );
				pose.set_omega( i, init_omega );
				pose.set_secstruct( i, 'L' );
			}
		}

		/// pick some fragments
		FragLib fraglib;
		string fullss, fullbb;
		for ( Size i=1; i<= nrepeat; ++i ) {
			fullss += repeatss;
			fullbb += repeatbb;
		}
		TR.Trace << "pick_frags: " << nrepeat << ' ' << repeatlen << ' ' << repeatss << ' ' << repeatbb << ' ' <<
			fullss << ' ' << fullbb << endl;
		// default for resample_centroid_seqs_fragweight is 1.0
		//
		Real const seq_weight( resample_centroid_seqs ? 1.0 : 0.0 ), ss_weight( 10.0 ), bb_weight( 100.0 );
		Sizes const fragsizes( make_vector1( 3, 9 ) );
		Size const nfrags( 200 );
		kinematics::MoveMap mm;
		for ( Size i=1; i<= chain_end( 1, pose ); ++i ) mm.set_bb( i, true );

		Sizes const homs_to_exclude;
		devel::blab::classic_frags::add_vall_fragments( fragsizes, nfrags, pose, mm, fullss, seq_weight, ss_weight,
			fraglib, homs_to_exclude, bb_weight, fullbb );


		{ // get rid of frags that violate the bb constraints
			Sizes const fragsizes( make_vector1( 3, 9 ) );
			Sizes const min_nns( make_vector1( 50, 25 ) );

			for ( Size si=1; si<= 2; ++si ) {
				Size const fragsize( fragsizes[si] );
				Size const min_nn( min_nns[si] );
				devel::blab::classic_frags::TorsionFragmentLibrary & lib( fraglib.library( fragsize ) );

				for ( Size fragpos=1; fragpos<= lib.size(); ++fragpos ) {
					for ( Size nn=lib[ fragpos ].size(); nn> min_nn; --nn ) {
						bool badfrag( false );

						devel::blab::classic_frags::TorsionFragment const frag( lib[ fragpos ][ nn ] );

						for ( Size k=1; k<= fragsize; ++k ) {
							//char const ss( frag.get_secstruct(k) );
							Real const phi( frag.get_torsion(k,1 ) ), psi( frag.get_torsion( k, 2 ) ), omega( frag.get_torsion(k,3 ) );
							char const bigbin( torsion2big_bin( phi, psi, omega ) );
							Size const repeatpos( ( ( fragpos +k-1 ) -1 )%repeatlen + 1 );
							char const desired_bigbin( repeatbb[ repeatpos-1 ] );
							if ( desired_bigbin != 'X' && desired_bigbin != bigbin ) {
								runtime_assert( string("ABGEO").find( desired_bigbin ) != string::npos );
								badfrag = true;
								break;
							}
						}

						if ( badfrag ) {
							TR.Trace << "delete_fragment: " << I(4,fragsize) << I(4,fragpos) << I(4,nn) << endl;
							lib[ fragpos ].erase( nn );
						}
					} // nn
					TR.Trace << "final_nn: " << I(4,fragsize) << I(4,fragpos) << I(4,(fragpos-1)%repeatlen+1) <<
						I(4,lib[ fragpos ].size()) << endl;
				} // fragpos
			} // fragsize
		} // scope


		/// setup symminfo?
		conformation::symmetry::SymmetryInfo symminfo;

		setup_unbound_frag_symminfo( repeatlen, nrepeat, base_repeat, symminfo );

		/// switch to centroid
		devel::blab::loops::switch_pose_to_residue_type_set( pose, CENTROID ); // was CENTROID_DNA

		/// now make symmetric
		pose::symmetry::make_symmetric_pose( pose, symminfo );



		Real target_twist( 360.0 / nrepeat_full ), target_rise( 0.0 );
		// if ( force_biangles ) {
		// 	target_twist = 180.0;
		// 	target_rise = 1000.0; // signal not to score rise
		// } else if ( force_triangles ) {
		// 	target_twist = 120.0;
		// 	target_rise = 1000.0; // signal not to score rise
		// } else if ( force_twistrise ) {
		// 	target_twist = option[ my_options::target_twist ];
		// 	target_rise  = option[ my_options::target_rise ];
		// }

		{
			Real const twistrise_energy_weight( 1.0 );
			use_linear_twist_rise_energy = true; /// HACKY GLOBAL VAR
			Sizes fragseq_poslist;
			bool use_twistrise_energy( true ); // default
			simple_fold_abinitio( fraglib, fragseq_poslist, pose, use_twistrise_energy, target_twist, target_rise, 0.0,
				twistrise_energy_weight );
		}

		//basic::prof_show_oneliner();
		Real const centroid_simtime( ( (double) clock() - sim_starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

		Real const final_centroid_score( (*score3)( pose ) );
		use_linear_twist_rise_energy = false; /// HACKY GLOBAL VAR

		/// get handedness
		Real handedness(0.0);
		{
			Vectors coords;
			for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
			handedness = get_chirality( repeatlen, coords );
		}
		Real const centroid_handedness( handedness );

		/// compute transform: twist, axis,
		Real twist, rise, min_radius, max_radius, com_radius, depth;
		compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius, depth );
		Real const centroid_twist( twist ), centroid_rise( rise ), centroid_radius( min_radius );


		Real final_score( final_centroid_score ), relax_simtime( 0.0 ), res_score(0);

		bool relax_this_decoy( false );

		// Size const n2cdistpos( repeatlen*(nrepeat-2) ); // hacking
		// Real n2cdist( pose.residue(1).xyz("N").distance( pose.residue(n2cdistpos).xyz("C") ) ),
		// 	centroid_n2cdist( n2cdist );
		// if ( force_biangles ) {
		// 	Real const max_twist_dev(10.0);
		// 	relax_this_decoy = ( fabs( twist - target_twist ) < max_twist_dev );
		// } else if ( force_triangles ) {
		// 	Real const max_twist_dev(5.0);
		// 	relax_this_decoy = ( fabs( twist - target_twist ) < max_twist_dev );
		// } else if ( force_twistrise ) {
		// 	Real const max_twist_dev(5.0), max_rise_dev( 2.0 );
		// 	relax_this_decoy = ( fabs( twist - target_twist ) < max_twist_dev ) &&
		// 		( fabs( rise - target_rise ) < max_rise_dev );
		// } else {
		{
			/// old way: ( rise > -2 && rise < 2 && twist > 0.9 * target_twist && twist < 1.1 * target_twist );
			Real const max_twist_dev(10.0);
			relax_this_decoy = ( rise > -2 && rise < 2 && fabs( twist - target_twist ) < max_twist_dev );
		}

		{ // check if the three helices are bundled:
			Real const dis2_threshold( 9*9 );
			Vectors coords;
			for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
			// look for contacts between helices over 4 rsd windows
			// helix3 should contact both other helices over its full length
			Size const h1_begin(1),               h1_end( h1_begin + helix1_len - 1 ),
				h2_begin( h1_end + turn1_len + 1 ), h2_end( h2_begin + helix2_len - 1 ),
				h3_begin( h2_end + turn2_len + 1 ), h3_end( h3_begin + helix3_len - 1 );
			Sizes h1_posl, h2_posl;
			for ( Size i=h1_begin; i<= h1_end; ++i ) h1_posl.push_back(i);
			for ( Size i=h2_begin; i<= h2_end; ++i ) h2_posl.push_back(i);
			bool all_contacts( true );
			Size const window_size(4);
			for ( Size i=h3_begin; i<= h3_end-window_size+1; ++i ) {
				bool h1_contact(false), h2_contact(false);
				for ( Size j=i; j<i+window_size; ++j ) {
					for ( Size k : h1_posl ) {
						if ( coords[j].distance_squared( coords[k] ) <= dis2_threshold ) {
							h1_contact = true;
						}
					}
					for ( Size k : h2_posl ) {
						if ( coords[j].distance_squared( coords[k] ) <= dis2_threshold ) {
							h2_contact = true;
						}
					}
				}
				if ( !( h1_contact && h2_contact ) ) {
					all_contacts = false;
					break;
				}
			} // loop over window_size windows in helix3
			if ( !all_contacts ) {
				cout << "unbundled" << endl;
				relax_this_decoy = false;
			}
		}

		cout << "relax_this_decoy: final " << relax_this_decoy << endl;

		//Size nrepeat_this_decoy( nrepeat );

		ostringstream relax_status; // some scores only get written if we relaxed this decoy

		string fullsspred( "-" );
		Real searchtime(0);
		bool refold_this_decoy( true );
		if ( relax_this_decoy || dry_run() ) {
			// not sure we need this here... but make sure fullatom donut scoring works OK
			symminfo_hack::nrepeat_ = nrepeat_full;
			symminfo_hack::repeatlen_ = repeatlen;

			searchtime = ((double) clock() - search_starttime )/( CLOCKS_PER_SEC*60 ); // records time other than in relax,refold

			// convert to fullatom
			conformation::symmetry::SymmetryInfo symminfo( *pose::symmetry::symmetry_info( pose ) );
			pose::symmetry::make_asymmetric_pose( pose );
			devel::blab::loops::switch_pose_to_residue_type_set( pose, FA_STANDARD );
			pose::symmetry::make_symmetric_pose( pose, symminfo );

			clock_t const relax_starttime = clock();

			SequenceConstraints sequence_constraints;

			if ( option[ my_options::force_helix_capping ] ) {
				Sizes const ncaps( make_vector1( repeatlen, helix1_len+turn1_len, helix1_len+turn1_len+helix2_len+turn2_len ));
				Sizes const ccaps( make_vector1( helix1_len+1, helix1_len+turn1_len+helix2_len+1,
						helix1_len+turn1_len+helix2_len+turn2_len+helix3_len+1 ) );
				for ( Size repeatpos=1; repeatpos<= repeatlen; ++repeatpos ) {
					string seqcst;
					if ( has_element( repeatpos, ccaps ) && repeatbb[ repeatpos-1 ] == 'G' ) seqcst = "G";
					if ( has_element( repeatpos, ncaps ) && repeatbb[ repeatpos-1 ] == 'B' ) seqcst = "STND";
					if ( seqcst.empty() || sequence_constraints.count( repeatpos ) ) continue; // earlier constraints have prio
					TR.Trace << "force_helix_capping: " << repeatpos << ' ' << seqcst << endl;
					for ( Size k=0; k<nrepeat; ++k ) {
						Size const seqpos( k*repeatlen + repeatpos );
						sequence_constraints[ seqpos ] = seqcst;
					}
				}
			}


			if ( option[ my_options::layer_design ] ) {
				///
				push_sequence_constraints_to_base_repeat( pose, sequence_constraints );
				string fullss;
				for ( Size i=1; i<= nrepeat; ++i ) fullss += target_repeatss;
				fullss.push_back( 'X' );
				runtime_assert( fullss.size() == pose.total_residue() );
				add_layer_design_constraints( pose, fullss, sequence_constraints );
			}

			{
				bool const add_jump_flex( false ), use_atom_pair_constraints( false ), use_coordinate_constraints( false );
				//bool const use_twistrise_energy( option[ my_options::fullatom_twistrise_energy ] );
				symmetric_design_and_relax( repeatlen, nrepeat_full, option[ my_options::use_softrep_for_design ],
					use_donut_energy, option[ my_options::design_cycles ], sequence_constraints, pose,
					add_jump_flex, use_atom_pair_constraints, use_coordinate_constraints,
					false );//use_fullatom_twistrise_energy, target_twist, target_rise );
			}

			/// restore these since they are zeroed by symmetric_design_and_relax
			symminfo_hack::nrepeat_ = nrepeat_full;
			symminfo_hack::repeatlen_ = repeatlen;


			relax_simtime = ((double) clock() - relax_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

			final_score = (*fa_scorefxn)( pose );
			res_score = final_score / repeatlen;

			///
			bools subset( pose.total_residue(), false );
			Size const base_repeat_offset( repeatlen*(base_repeat-1));
			for ( Size i=base_repeat_offset+1; i<= base_repeat_offset+repeatlen; ++i ) subset[i] = true;

			Real sasapack_score, norme_score, normsasa_score, exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa, buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score;
			compute_sasa_scores_for_subset_slow( 5, subset, pose, sasapack_score, norme_score, normsasa_score,
				exposed_polar_sasa, exposed_nonpolar_sasa,
				buried_polar_sasa, buried_nonpolar_sasa,
				buried_nonpolar_sasa_sc, buried_nonpolar_rsd_sasa_sc,
				packstat_score );

			Real n_buried_unsatisfied_donors_bb_per_repeat, n_buried_unsatisfied_acceptors_bb_per_repeat,
				n_buried_unsatisfied_donors_sc_per_repeat, n_buried_unsatisfied_acceptors_sc_per_repeat;
			{ // compute unsats over full pose
				Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
					n_buried_unsatisfied_acceptors_sc;
				bools subset_full( pose.total_residue(), true ); subset_full[ pose.total_residue() ] = false;
				get_buried_unsatisfied_counts_real_slow( subset_full, pose,
					n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc,
					n_buried_unsatisfied_acceptors_bb, n_buried_unsatisfied_acceptors_sc );

				n_buried_unsatisfied_donors_bb_per_repeat = Real( n_buried_unsatisfied_donors_bb )/nrepeat;
				n_buried_unsatisfied_donors_sc_per_repeat = Real( n_buried_unsatisfied_donors_sc )/nrepeat;
				n_buried_unsatisfied_acceptors_bb_per_repeat = Real( n_buried_unsatisfied_acceptors_bb )/nrepeat;
				n_buried_unsatisfied_acceptors_sc_per_repeat = Real( n_buried_unsatisfied_acceptors_sc )/nrepeat;
			}

			if ( option[ my_options::pick_decoys ] ) {
				Real const exposed_polar_fraction( exposed_polar_sasa / ( exposed_polar_sasa + exposed_nonpolar_sasa ) );
				if ( exposed_polar_fraction > 0.55 || exposed_polar_fraction < 0.35 ) refold_this_decoy = false;
				if ( sasapack_score > 0 ) refold_this_decoy = false;

				if ( refold_this_decoy ) {
					Real n_buried_unsatisfied_donors_bb, n_buried_unsatisfied_donors_sc, n_buried_unsatisfied_acceptors_bb,
						n_buried_unsatisfied_acceptors_sc;
					get_buried_unsatisfied_counts_real_slow( subset, pose,
						n_buried_unsatisfied_donors_bb,
						n_buried_unsatisfied_donors_sc,
						n_buried_unsatisfied_acceptors_bb,
						n_buried_unsatisfied_acceptors_sc );
					if ( ( n_buried_unsatisfied_donors_bb + n_buried_unsatisfied_donors_sc > 1 ) ||
							( n_buried_unsatisfied_acceptors_bb + n_buried_unsatisfied_acceptors_sc > 0 ) ) refold_this_decoy = false;
				}
			}

			//n2cdist = pose.residue(1).xyz("N").distance( pose.residue(chain_end(1,pose)).xyz("C") );

			/// compute psipred ss using single sequence
			string repeatsspred( "-" );
			//Real sspred_match(0.0), sspred_turn1(0.0), sspred_turn2( 0.0 ), sspred_helix1(0.0), sspred_helix2(0.0);
			{ // run psipred to re-predict the secondary structure
				vector1< Reals > pred_HEL;
				string const sequence( pose.sequence().substr(0,chain_end(1,pose) ) );
				//string const secstruct( pose.secstruct().substr(0,chain_end(1,pose) ) );
				string sspred;
				string const repeatss( pose.secstruct().substr(0,repeatlen) ); // current, not the target

				run_psipred( sequence, sspred, pred_HEL );

				if ( sspred.size() == sequence.size() ) { // success
					fullsspred = sspred;
					runtime_assert( sspred.size() == nrepeat * repeatlen );
					repeatsspred = sspred.substr( (base_repeat-1)*repeatlen, repeatlen );
					// Size const turn1_begin( helix1_len+1 ), turn1_end( helix1_len + turn1_len ),
					// 	turn2_begin( helix1_len+turn1_len+helix2_len+1 );
					// for ( Size i=1; i<= sspred.size(); ++i ) {
					// 	Size const rpos( (i-1)%repeatlen+1 );
					// 	char const ss( repeatss[ rpos-1 ] );
					// 	Real const prediction( ss == 'H' ? pred_HEL[i][1] : ( ss == 'E' ? pred_HEL[i][2] : pred_HEL[i][3] ) );
					// 	sspred_match += prediction;
					// 	if ( rpos < turn1_begin                     ) sspred_helix1 += prediction;
					// 	if ( rpos > turn1_end && rpos < turn2_begin ) sspred_helix2 += prediction;
					// 	if ( rpos >= turn1_begin-1 && rpos <= turn1_end+1 ) sspred_turn1 += prediction;
					// 	if ( rpos == 1 || rpos >= turn2_begin-1           ) sspred_turn2 += prediction;
					// }
					// sspred_match  /= (  repeatlen * nrepeat );
					// sspred_helix1 /= ( helix1_len * nrepeat );
					// sspred_helix2 /= ( helix2_len * nrepeat );
					// sspred_turn1  /= ( ( turn1_len+2) * nrepeat );
					// sspred_turn2  /= ( ( turn2_len+2) * nrepeat );
				}
			}

			// Real helix1_dist(0), helix1_twist(0), helix2_dist(0), helix2_twist(0);
			// if ( helix1_len ) {
			// 	runtime_assert( helix1_len + turn1_len + helix2_len + turn2_len == repeatlen );

			// 	compute_helix_axis_angles( pose, base_repeat, repeatlen, helix1_len, turn1_len, helix2_len,
			// 		helix1_dist, helix1_twist, helix2_dist, helix2_twist );
			// }

			// Real const
			// 	inner_helix_dist( min( helix1_dist, helix2_dist ) ),
			// 	outer_helix_dist( max( helix1_dist, helix2_dist ) );
			// Real const inner_helix_twist( helix1_dist < helix2_dist ? helix1_twist : helix2_twist );
			// Real const outer_helix_twist( helix1_dist < helix2_dist ? helix2_twist : helix1_twist );

			// Size const inner_helix_len( helix1_dist < helix2_dist ? helix1_len : helix2_len );
			// Size const outer_helix_len( helix1_dist < helix2_dist ? helix2_len : helix1_len );
			// string const inner_helix_turn( helix1_dist < helix2_dist ? turn1 : turn2 );
			// string const outer_helix_turn( helix1_dist < helix2_dist ? turn2 : turn1 );

			string const buried_unsatisfied_string( get_buried_unsatisfied_string( pose ) );

			if ( option[ my_options::max_sasapack_score ].user() && sasapack_score > option[ my_options::max_sasapack_score]){
				refold_this_decoy = false;
			}

			relax_status <<
				" sasapack_score: " << F(9,3,sasapack_score) <<
				" norme_score: " << F(9,3,norme_score) <<
				" normsasa_score: " << F(9,3,normsasa_score) <<
				" exposed_polar_sasa: " << F(9,3,exposed_polar_sasa)<<
				" exposed_nonpolar_sasa: " << F(9,3,exposed_nonpolar_sasa)<<
				" exposed_polar_fraction: " << F(9,3,exposed_polar_sasa/(exposed_polar_sasa+exposed_nonpolar_sasa))<<
				" buried_polar_sasa: " << F(9,3,buried_polar_sasa)<<
				" buried_nonpolar_sasa: " << F(9,3,buried_nonpolar_sasa)<<
				" buried_nonpolar_sasa_sc: " << F(9,3,buried_nonpolar_sasa_sc)<<
				" buried_nonpolar_rsd_sasa_sc: " << F(9,3,buried_nonpolar_rsd_sasa_sc)<<
				" packstat_score: " << F(9,3,packstat_score) <<
				" res_score: " << F(9,3,res_score ) <<
				" n_buried_unsatisfied_donors_bb_per_repeat: " << F(9,3,n_buried_unsatisfied_donors_bb_per_repeat) <<
				" n_buried_unsatisfied_donors_sc_per_repeat: " << F(9,3,n_buried_unsatisfied_donors_sc_per_repeat) <<
				" n_buried_unsatisfied_acceptors_bb_per_repeat: " <<F(9,3,n_buried_unsatisfied_acceptors_bb_per_repeat)<<
				" n_buried_unsatisfied_acceptors_sc_per_repeat: " <<F(9,3,n_buried_unsatisfied_acceptors_sc_per_repeat)<<
				" repeatsspred: " << repeatsspred <<
				// " sspred_match: " << F(9,3,sspred_match) <<
				// " sspred_helix1: " << F(9,3,sspred_helix1) <<
				// " sspred_helix2: " << F(9,3,sspred_helix2) <<
				// " sspred_turn1: " << F(9,3,sspred_turn1) <<
				// " sspred_turn2: " << F(9,3,sspred_turn2) <<
				// " helix1_dist: " << F(9,3,helix1_dist) <<
				// " helix1_twist: " << F(9,3,helix1_twist) <<
				// " helix2_dist: " << F(9,3,helix2_dist) <<
				// " helix2_twist: " << F(9,3,helix2_twist) <<
				// " inner_helix_dist: " << F(9,3,inner_helix_dist) <<
				// " inner_helix_twist: " << F(9,3,inner_helix_twist) <<
				// " outer_helix_dist: " << F(9,3,outer_helix_dist) <<
				// " outer_helix_twist: " << F(9,3,outer_helix_twist) <<
				// " inner_helix_len: " << I(4,inner_helix_len) <<
				// " inner_helix_turn: " << inner_helix_turn <<
				// " outer_helix_len: " << I(4,outer_helix_len) <<
				// " outer_helix_turn: " << outer_helix_turn <<
				' ' << buried_unsatisfied_string <<
				' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() );


			/// RE-compute transform: twist, axis, after relax
			compute_repeat_params( pose, twist, rise, min_radius, max_radius, com_radius );
			{
				Vectors coords;
				for ( Size i=1; i<= nres_protein; ++i ) coords.push_back( pose.residue(i).xyz("CA") );
				handedness = get_chirality( repeatlen, coords );
			}
		}

		string const repeatseq( pose.sequence().substr(0,repeatlen) );




		Real ala_fraction(0);
		foreach_ ( char const s, repeatseq ) { ala_fraction += ( s == 'A' ); }
		ala_fraction /= repeatseq.size();

		if ( option[ my_options::max_ala_fraction ].user() && ala_fraction - 1e-6 > option[ my_options::max_ala_fraction]){
			refold_this_decoy = false;
		}

		if ( use_donut_energy ) { // problem with donut energy!
			Real const max_twist_dev( 5.0 );
			refold_this_decoy = refold_this_decoy && ( rise > -2 && rise < 2 && fabs( twist - target_twist ) < max_twist_dev );
		}
		// } else if ( force_biangles || force_triangles ) {
		// 	Real const max_twist_dev(2.5);
		// 	refold_this_decoy = refold_this_decoy && ( fabs( twist - target_twist ) < max_twist_dev );
		// } else if ( force_twistrise ) {
		// 	Real const max_twist_dev(2.0), max_rise_dev( 1.0 );
		// 	refold_this_decoy = ( fabs( twist - target_twist ) < max_twist_dev ) &&
		// 		( fabs( rise - target_rise ) < max_rise_dev );
		// }


		Size const base_repeat_offset( repeatlen * ( base_repeat-1 ) );
		repeatbb = torsion2big_bin_string( base_repeat_offset+1, base_repeat_offset+repeatlen, pose );
		repeatss = pose.secstruct().substr(0,repeatlen);

		Vectors ca_coords;
		for ( Size i=1; i<= nres_protein; ++i ) ca_coords.push_back( pose.residue(i).xyz("CA") );
		string helix_stats;
		// if ( force_biangles || force_triangles ) {
		// 	helix_stats = get_single_helix_stats( nrepeat, repeatbb, repeatseq, ca_coords );
		// } else {
		{
			helix_stats = get_double_helix_stats( nrepeat, repeatbb, repeatseq, ca_coords );
		}

		cout << "refold_this_decoy: final " << refold_this_decoy << ' ' << relax_this_decoy <<
			" ala_fraction: " << ala_fraction << ' ' << bbtag << endl;

		string const simfile( outdir+"unbound_frag_"+string_of( repeatlen/4 )+".work" );

		/// how do we want to score filter?? by bbtags? seems too granular.
		bool const passed_score_filter( relax_this_decoy && refold_this_decoy && // short-circuit evaluation
			append_score_to_scorefile_and_filter( "tmp", res_score, // was final_score (doh!)
				score_filter_acceptance_rate,
				score_filter_pass_early, simfile ) );

		ostringstream refolding_status;

		/// try refolding de novo
		Size const n_refold( option[ my_options::n_refold ] ), n_refold_sspred( option[ my_options::n_refold_sspred ] );
		Real refolding_simtime(0);
		if ( refold_this_decoy && n_refold+n_refold_sspred > 0 && ( passed_score_filter || dry_run() ) ) {
			clock_t const refold_starttime = clock();
			for ( Size RR=1; RR<= 2; ++RR ) {
				Real min_refold_rmsd(1e6);
				Size n_under_4(0);
				if ( RR == 2 && fullsspred.size() != nrepeat * repeatlen ) continue;
				if ( RR==1 ) refolding_status << " n_refold: " << I(3,n_refold);
				else         refolding_status << " n_refold_sspred: " << I(3,n_refold_sspred);
				Size const rr_end( RR==1 ? n_refold : n_refold_sspred );
				for ( Size rr=1; rr<= rr_end; ++rr ) {
					Pose protpose;
					string const repss( RR == 1 ? repeatss : fullsspred );
					Size const nrepeat_refold( nrepeat );
					Size const base_repeat_refold( base_repeat );
					refold_repeat_pose( repeatseq, repss, nrepeat_refold, base_repeat_refold, protpose );
					/// calc rmsd
					Real rmsd( 0.0 ), rmsd_single_repeat( 0.0 );
					{
						using namespace core::id;
						AtomID_Map< AtomID > atom_map;
						initialize_atomid_map( atom_map, protpose, id::GLOBAL_BOGUS_ATOM_ID );
						for ( Size i=1; i<= repeatlen; ++i ) {
							atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
								AtomID( pose.residue(i).atom_index("CA"),i);
						}
						rmsd_single_repeat = rmsd_by_mapping( protpose, pose, atom_map );
						Size nrepeat_for_refolding_rmsd( min( nrepeat, nrepeat_refold ) );
						if ( option[ my_options::nrepeat_for_refolding_rmsd ].user() ) {
							nrepeat_for_refolding_rmsd = option[ my_options::nrepeat_for_refolding_rmsd ];
						}
						Size const nres_for_rmsd( min( chain_end( 1, protpose), nrepeat_for_refolding_rmsd * repeatlen ) );
						for ( Size i=1; i<= nres_for_rmsd; ++i ) {
							atom_map[ AtomID( protpose.residue(i).atom_index("CA"), i ) ] =
								AtomID( pose.residue(i).atom_index("CA"),i);
						}
						rmsd = rmsd_by_mapping( protpose, pose, atom_map );
						min_refold_rmsd = min( min_refold_rmsd, rmsd );
						if ( rmsd < 4 ) ++n_under_4;
					}

					Real refold_handedness( 0.0 );
					{
						Vectors coords;
						for ( Size i=1; i<= chain_end( 1, protpose); ++i ) coords.push_back( protpose.residue(i).xyz("CA") );
						refold_handedness = get_chirality( repeatlen, coords );
					}
					Real tmptwist, tmprise, tmpmin_radius, tmpmax_radius, tmpcom_radius;
					compute_repeat_params( protpose, tmptwist, tmprise, tmpmin_radius, tmpmax_radius, tmpcom_radius );
					refolding_status << " R " << rr << F(9,3,rmsd) << F(9,3,rmsd_single_repeat) << F(9,3,refold_handedness) <<
						F(9,3,tmprise) << F(9,3,tmptwist);
				}
				refolding_status << " min_refold" << RR << "_rmsd: "<< F(9,3,min_refold_rmsd) <<
					" n_refold" << RR << "_under_4: " << n_under_4;

			}
			refolding_simtime = ((double) clock() - refold_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes
		}


		string const outfilename( output_tag() + "unbound_frag_"+bbtag+"_"+lead_zero_string_of( n,4 )+".pdb" );

		Real const simtime = ((double) clock() - sim_starttime )/( CLOCKS_PER_SEC*60 ); // in minutes

		ostringstream out;
		out << "final_scores " << F(9,3,final_score) << ' ' << outfilename <<
			" handedness: " << F(9,3,handedness) <<
			" min_radius: " << F(9,3,min_radius ) <<
			" rise: " << F(9,3,rise ) <<
			" twist: " << F(9,3,twist) <<
			" final_centroid_score: "<< F(9,3,final_centroid_score) <<
			" passed_score_filter: " << passed_score_filter <<
			" relaxed: " << relax_this_decoy <<
			" min_radius: " << F(9,3,min_radius )<<
			" max_radius: " << F(9,3,max_radius )<<
			" com_radius: " << F(9,3,com_radius )<<
			" centroid_handedness: " << F(9,3,centroid_handedness) <<
			" centroid_radius: " << F(9,3,centroid_radius ) <<
			" centroid_rise: " << F(9,3,centroid_rise ) <<
			" centroid_twist: " << F(9,3,centroid_twist ) <<
			" repeatseq: " << repeatseq <<
			" repeatbb: " << repeatbb <<
			" repeatss: " << repeatss <<
			" turn1: " << turn1 <<
			" turn2: " << turn2 <<
			" turn3: " << turn3 <<
			" helix1_len: " << helix1_len <<
			" helix2_len: " << helix2_len <<
			" helix3_len: " << helix3_len <<
			" repeatlen: " << repeatlen <<
			" nrepeat: " << nrepeat <<
			//" n2cdist: " << F(9,3,n2cdist) <<
			//" centroid_n2cdist: " << F(9,3,centroid_n2cdist) <<
			" target_repeatbb: " << target_repeatbb <<
			" target_repeatss: " << target_repeatss <<
			" " << relax_status.str() <<
			" ala_fraction: " << F(9,3,ala_fraction) <<
			" bbtag: " << bbtag <<
			" " << refolding_status.str() <<
			" centroid_simtime: " << F(9,3,centroid_simtime) <<
			" relax_simtime: " << F(9,3,relax_simtime) <<
			" refolding_simtime: " << F(9,3,refolding_simtime) <<
			" searchtime: " << F(9,3,searchtime) <<
			" simtime: " << F(9,3,simtime) <<
			' ' << helix_stats << '\n';

		string const pdbfilename( create_output( pose, passed_score_filter, outfilename, out.str(), simfile ) );

		// if ( passed_score_filter || dry_run() ) {
		// 	pose.dump_pdb( outfilename );
		// 	run_command("gzip "+outfilename ); // NEW NEW NEW
		// } // passed_score_filter
		if ( passed_score_filter && pdbfilename.size()>0 ) {
			append_hbond_info_to_pdb_file( pose, *fa_scorefxn, pdbfilename );
			run_command("gzip "+pdbfilename );
		}

		fflush( stdout );
		check_simtime();

		search_starttime = clock(); // restart the search timer
	}




}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	using namespace basic::options;

	add_my_options();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::cout << "ARGS ";
	for ( int i=0; i< argc; ++i ) {
		std::cout << ' ' <<  argv[i];
	}
	std::cout << std::endl;

	devel::init(argc, argv);

	check_if_job_is_done();

	cout << "wow cooler1 exe= symdes.cc" << endl;

	string const mode( basic::options::option[ my_options::mode ] );

	cout << "wow cooler2 mode= " << mode << endl;

	runtime_assert( option[ basic::options::OptionKeys::out::file::output_virtual ].user() );

	modify_atom_properties_from_command_line(); // in case we are stupid and using set_atom_property

	run_command("hostname");

	if ( mode == "dna" ) {
		dna_test();
	} else if ( mode == "tal" ) {
		tal_test();
	} else if ( mode == "design_stats" ) {
		design_stats_test();
	} else if ( mode == "tal_repeat_params" ) {
		tal_repeat_params_test();
	} else if ( mode == "interface" ) {
		interface_test();
	} else if ( mode == "shift" ) {
		shift_test();
	} else if ( mode == "extend_design_pdbs" ) {
		extend_design_pdbs_test();
	} else if ( mode == "native_hand" ) {
		native_hand_test();
	} else if ( mode == "helical_params" ) {
		helical_params_test();
	} else if ( mode == "spin" ) {
		spin_test();
	} else if ( mode == "helix_frag" ) {
		helix_frag_test();
	} else if ( mode == "helix_frag_v2" ) {
		helix_frag_v2_test();
	} else if ( mode == "helix_frag_v3" ) {
		helix_frag_v3_test();
	} else if ( mode == "reconstruct_helix_frag" ) {
		reconstruct_helix_frag_test();
	} else if ( mode == "reconstruct_helix_frag_v2" ) {
		reconstruct_helix_frag_v2_test();
	} else if ( mode == "dock_helix_frag_v2" ) {
		dock_helix_frag_v2_test();
	} else if ( mode == "dock_helix_frag" ) {
		dock_helix_frag_test();
	} else if ( mode == "autodock_helix_frag" ) {
		autodock_helix_frag_test();
	} else if ( mode == "autodock_helix_frag_v2" ) {
		autodock_helix_frag_v2_test();
	} else if ( mode == "cluster_2c" ) {
		cluster_2c_test();
	} else if ( mode == "symdock_dna" ) {
		symdock_dna_test();
	} else if ( mode == "regan_tpr" ) {
		regan_tpr_test();
	} else if ( mode == "w3b" ) {
		w3b_test();
	} else if ( mode == "redock" ) {
		redock_test();
	} else if ( mode == "template_redock" ) {
		template_redock_test();
	} else if ( mode == "donut_deriv" ) {
		donut_deriv_test();
	} else if ( mode == "talspecsym" ) {
		talspecsym_test();
	} else if ( mode == "phenix_rebuild_again" ) {
		phenix_rebuild_again_test();
	} else if ( mode == "frag" ) {
		frag_test();
	} else if ( mode == "refold_toroid" ) {
		refold_toroid_test();
	} else if ( mode == "refold_logfile" ) {
		refold_logfile_test();
	} else if ( mode == "refold_logfile_asym" ) {
		refold_logfile_asym_test();
	} else if ( mode == "unbound_frag" ) {
		unbound_frag_test();
	} else if ( mode == "new_unbound_frag" ) {
		new_unbound_frag_test();
	} else if ( mode == "triple_unbound_frag" ) {
		triple_unbound_frag_test();
	} else if ( mode == "unbound_frag_extend" ) {
		unbound_frag_extend_test();
	} else if ( mode == "sasa" ) {
		sasa_test();
	} else if ( mode == "read_int_designs" ) {
		read_int_designs_test();
	} else if ( mode == "read_cage_designs" ) {
		read_cage_designs_test();
	} else if ( mode == "cage" ) {
		cage_test();
	} else if ( mode == "two_component" ) {
		two_component_test();
	} else if ( mode == "planar_two_component" ) {
		planar_two_component_test();
	} else if ( mode == "planar_one_component" ) {
		planar_one_component_test();
	} else if ( mode == "cage_design" ) {
		cage_design_test();
	} else if ( mode == "generic_cage_design" ) {
		generic_cage_design_test();
	} else if ( mode == "sasa_fractions" ) {
		sasa_fractions_test();
	} else if ( mode == "sasapack" ) {
		sasapack_test();
	// } else if ( mode == "sasapack_silent" ) {
	// 	sasapack_silent_test();
	} else if ( mode == "design" ) {
		design_test();
	} else if ( mode == "unsatisfied" ) {
		unsatisfied_test();
	} else if ( mode == "rescore" ) {
		rescore_test();
	} else if ( mode == "bb_strain" ) {
		bb_strain_test();
	} else if ( mode == "reconstruct_2c" ) {
		reconstruct_2c_test();
	} else if ( mode == "reconstruct_and_design" ) {
		reconstruct_and_design_test();
	} else if ( mode == "reconstruct_planar_2c" ) {
		reconstruct_planar_2c_test();
	} else if ( mode == "rescore_logfile" ) {
		rescore_logfile_test();
	} else if ( mode == "rescore_unbound_logfile" ) {
		rescore_unbound_logfile_test();
	} else if ( mode == "redesign_unbound_logfile" ) {
		redesign_unbound_logfile_test();
	} else if ( mode == "centroid" ) {
		centroid_test();
	} else if ( mode == "refold" ) {
		refold_test();
	} else if ( mode == "toroid1_rmsd" ) {
		toroid1_rmsd_test();
	} else if ( mode == "toroid1_rmsd_v2" ) {
		toroid1_rmsd_v2_test();
	} else if ( mode == "refold_and_dock" ) {
		refold_and_dock_test();
	} else if ( mode == "refold_9x2" ) {
		refold_9x2_test();
	} else if ( mode == "refold_3x2" ) {
		refold_3x2_test();
	} else if ( mode == "refold_tordim1" ) {
		refold_tordim1_test();
	} else if ( mode == "refold_tri" ) {
		refold_tri_test();
	} else if ( mode == "symmetric_refold_tri" ) {
		symmetric_refold_tri_test();
	} else if ( mode == "rephase_tri" ) {
		rephase_tri_test();
	} else if ( mode == "extend_tri" ) {
		extend_tri_test();
	} else if ( mode == "rsd_sasa" ) {
		rsd_sasa_test();
	} else if ( mode == "analyze_interface" ) {
		analyze_interface_test();
	} else if ( mode == "helix_loop" ) {
		helix_loop_test();
	} else if ( mode == "helix_junction" ) {
		helix_junction_test();
	} else if ( mode == "semet" ) {
		semet_test();
	} else if ( mode == "relax_designs" ) {
		relax_designs_test();
	} else if ( mode == "centroid_build" ) {
		centroid_build_test();
	} else if ( mode == "aa" ) {
		aa_test();
	} else if ( mode == "add_helix_stats" ) {
		add_helix_stats_test();
	} else if ( mode == "unbound_cluster" ) {
		unbound_cluster_test();
	} else if ( mode == "tri_rmsd" ) {
		tri_rmsd_test();
	} else if ( mode == "refold_generic" ) {
		refold_generic_test();
	} else if ( mode == "annotate_pdb" ) {
		annotate_pdb_test();
	} else if ( mode == "strand" ) {
		strand_test();
	} else if ( mode == "hairpin" ) {
		hairpin_test();
	} else if ( mode == "strand_turn" ) {
		strand_turn_test();
	} else if ( mode == "db_turn" ) {
		db_turn_test();
	} else if ( mode == "unbound_frag_redesign" ) {
		unbound_frag_redesign_test();
	} else {
		cout << "Unrecognized mode " << mode << endl;
	}


	exit(0);

}
