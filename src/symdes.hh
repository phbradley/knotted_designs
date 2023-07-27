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

// libRosetta headers
#ifndef INCLUDED_apps_pilot_phil_symdes_HH
#define INCLUDED_apps_pilot_phil_symdes_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/phil_options.hh>
#include <apps/pilot/phil/phil_io.hh>
#include <apps/pilot/phil/types.hh>
#include <apps/pilot/phil/constraints.hh>
#include <apps/pilot/phil/rotamer_correctness.hh>

#include <apps/pilot/phil/symscores.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <protocols/dna/water.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/task_operations/LimitAromaChi2Operation.hh>
// #include <protocols/task_operations/RestrictToNativelikeRotamers.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/annealer/AnnealerFactory.hh>
#include <core/pack/annealer/SimAnnealerBase.hh>
#include <core/pack/annealer/FixbbSimAnnealer.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>

static basic::Tracer TR_SYMDES_HH( "apps.pilot.phil.symdes_hh" );



// perhaps we want to make these configurable somehow
string relax_scorefxn_weights_tag( "ref2015" );
string soft_design_scorefxn_weights_tag( "ref2015_soft" );
// string relax_scorefxn_weights_tag( "talaris2013" );
// string soft_design_scorefxn_weights_tag( "talaris2013_soft" );


////////////////////////////////////////////////////////////////////////////////////////////////////////////
Pose const empty_pose_for_symdes;
void
symmetric_design_and_relax(
	Size const repeatlen,
	Size const nrepeat_full,
	bool const use_softrep_for_design,
	bool const use_donut_energy,
	Size n_cycles,
	//vector1< std::pair< Size, string > > const & sequence_constraints,
	SequenceConstraints sequence_constraints,
	Pose & pose,
	bool const add_jump_flex = false,
	bool const use_atom_pair_constraints = false,
	bool const use_coordinate_constraints = false,
	bool const use_twistrise_energy = false,
	Real const target_twist = 0.0,
	Real const target_rise = 0.0,
	bool const freeze_backbone = false,
	bool const flex_all_independent_jumps = false,
	Pose const & native_pose = empty_pose_for_symdes
)
{
	bool const debugging( false );

	bool const favor_native_residue( option[ my_options::favor_native_residue ].user() );

	/// add the residue constraints; we assume that the starting pose has the native sequence !!!
	if ( favor_native_residue ) {
		using namespace core::scoring::constraints;

		Size base_repeat(0);
		conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
		for ( Size i=1; i<= nrepeat_full; ++i ) {
			if ( i*repeatlen <= pose.total_residue() && symminfo->bb_is_independent( i*repeatlen ) ) {
				runtime_assert( !base_repeat );
				base_repeat = i;
			}
		}

		/// add some residue type constraints
		Pose const & restype_pose( native_pose.total_residue()>0 ? native_pose : pose );
		vector1< ConstraintCOP > favor_native_constraints;
		Real const bonus( option[ my_options::favor_native_residue ]() );
		Size const base_offset( (base_repeat-1)*repeatlen );
		for ( Size i = base_offset+1; i <= base_offset+repeatlen; ++i ) {
			runtime_assert( restype_pose.total_residue() >= i );
			favor_native_constraints.push_back( ConstraintOP( new ResidueTypeConstraint( restype_pose, i, bonus ) ) );
			TR_SYMDES_HH.Trace << "favor_native_residue: " << I(4,i) << ' ' << restype_pose.residue(i).name3() << ' ' <<
				bonus << endl;
		}
		//adds to pose and scorefxn
		pose.add_constraints( favor_native_constraints );
	}


	protocols::moves::PyMOLMoverOP pymol(0);
	if ( option[ my_options::show_pymol ] ) {
		pymol = protocols::moves::PyMOLMoverOP( new protocols::moves::PyMOLMover());
		pymol->print("symmetric_design_and_relax: start");
		pymol->pymol_name("symmetric_design_and_relax");
		pymol->keep_history(true );
		pymol->apply( pose );
	}

	///////////////
	/// establish new convention for SequenceConstraints:: only applied at independent positions ///////////////////
	push_sequence_constraints_to_base_repeat( pose, sequence_constraints );

	ScoreFunctionOP fa_scorefxn(0);
	if ( option[ OptionKeys::dna::specificity::score_function ].user() ) {
		fa_scorefxn = get_score_function_from_command_line();
	} else {
		fa_scorefxn = ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag );
	}
	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );
	adjust_ref_weights_from_command_line( *fa_scorefxn ); // since this is brand new

	if ( use_atom_pair_constraints ) fa_scorefxn->set_weight( atom_pair_constraint, 1.0 );
	if ( use_coordinate_constraints ) fa_scorefxn->set_weight( coordinate_constraint, 1.0 );

	DonutWholeEnergy donut_whole_energy;
	Donut1B_Energy donut_1b_energy;
	donut_whole_energy.pose_is_subset( true );
	donut_1b_energy.pose_is_subset( true );
	if ( use_donut_energy ) {
		Real const donut_energy_weight( 0.5 );//option[ my_options::donut_energy_weight ] );
		fa_scorefxn->add_extra_method( donut_whole, donut_energy_weight, donut_whole_energy );
		fa_scorefxn->add_extra_method( donut_1b   , donut_energy_weight, donut_1b_energy    );
		symminfo_hack::nrepeat_ = nrepeat_full;
		symminfo_hack::repeatlen_ = repeatlen;
	}
	TwistRiseEnergy    twist_rise_energy   ( target_twist, target_rise );
	TwistRise1B_Energy twist_rise_1b_energy( target_twist, target_rise );
	if ( use_twistrise_energy ) {
		fa_scorefxn->add_extra_method( twistrise   , 0.5, twist_rise_energy );
		fa_scorefxn->add_extra_method( twistrise_1b, 0.5, twist_rise_1b_energy );
	}
	// always safe to have SequenceBonusEnergy turned on (?)
	SequenceBonusEnergy sequence_bonus_energy;
	Real const sequence_bonus_weight( 1.0 );
	fa_scorefxn->add_extra_method( sequence_bonus, sequence_bonus_weight, sequence_bonus_energy );
	if ( favor_native_residue ) fa_scorefxn->set_weight( res_type_constraint, 1.0 );

	// { // hacking: test out closure minimization; seems to work fine
	//  MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
	//  movemap->set_bb (true);
	//  movemap->set_bb ( pose.total_residue(), false );
	//  //ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( relax_scorefxn_weights_tag ) );
	//  scoring::symmetry::SymmetricScoreFunctionOP scorefxn( new scoring::symmetry::SymmetricScoreFunction() );
	//  scorefxn->add_extra_method( donut_whole, donut_energy_weight, donut_whole_energy );
	//  scorefxn->add_extra_method( donut_1b   , donut_energy_weight, donut_1b_energy    );

	//  protocols::minimization_packing::symmetry::SymMinMoverOP min_mover
	//   ( new protocols::minimization_packing::symmetry::SymMinMover(movemap, scorefxn, "dfpmin", 0.00001, true,
	//                              true, true ));
	//  min_mover->apply( pose );
	//  exit(0);
	// }

	Real const startscore( (*fa_scorefxn)( pose ) );
	TR_SYMDES_HH.Trace << "symmetric_design_and_relax:: startscore " << F(9,3,startscore) <<
		' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << endl;
	//if ( use_atom_pair_constraints ) pose.constraint_set()->show( TR_SYMDES_HH.Trace );
	for ( Size n=1; n<= n_cycles; ++n ) {
		Real const startcyclescore( (*fa_scorefxn)( pose ) );
		static Size counter(0);
		++counter;
		string outfilename( "start_"+string_of(counter)+"_"+string_of(n)+".pdb");
		if ( debugging ) pose.dump_pdb(outfilename);
		TR_SYMDES_HH.Trace << "symmetric_design_and_relax:: startcyclescore " << F(9,3,startcyclescore) << ' ' << outfilename <<
			' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << endl;
		//if ( use_atom_pair_constraints ) pose.constraint_set()->show( TR_SYMDES_HH.Trace );

		{ // design calculation
			pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
			task->initialize_from_command_line();
			task->or_include_current( true );

			bools is_protein( pose.total_residue(), true ); is_protein[ pose.total_residue() ] = false;
			task->restrict_to_residues( is_protein );

			{ /// add sequence forcing here?
				conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
				for ( SequenceConstraints::const_iterator
						it= sequence_constraints.begin(); it!= sequence_constraints.end(); ++it ) {
					// if ( pose::symmetry::symmetry_info( pose )->contiguous_monomers() ) {
					//  Size const repeatpos( it->first );
					//  string const name1s( it->second );
					//  bools allowed_aas( num_canonical_aas, false );
					//  for ( Size k=0; k<name1s.size(); ++k ) {
					//   AA const aa( aa_from_oneletter_code( name1s[k] ) );
					//   allowed_aas[aa] = true;
					//  }
					//  runtime_assert( pose.total_residue()%repeatlen == 1 ); // includes vrt pos
					//  Size const nrepeat( ( pose.total_residue()-1 )/repeatlen );
					//  for ( Size k=0; k<nrepeat; ++k ) {
					//   Size const seqpos( k*repeatlen + repeatpos );
					//   TR_SYMDES_HH.Trace << "sequence_constraints: " << repeatpos << ' ' << seqpos << ' ' << name1s << endl;
					//   task->nonconst_residue_task( seqpos ).restrict_absent_canonical_aas( allowed_aas );
					//  }
					// } else {
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
						TR_SYMDES_HH << "sequence_constraints dont match current aa: " << base_pos << ' ' <<
							pose.residue( base_pos ).aa() << ' ' << name1s << endl;
						//make_sequence_change( base_pos, aa_from_oneletter_code( name1s[0] ), pose );
					}
					//runtime_assert( pose.total_residue()%repeatlen == 1 ); // includes vrt pos
					for ( Size seqpos=1; seqpos<= pose.total_residue(); ++seqpos ) {
						if ( seqpos == base_pos || symminfo->bb_follows( seqpos ) == base_pos ) {
							task->nonconst_residue_task( seqpos ).restrict_absent_canonical_aas( allowed_aas );
						}
						if ( seqpos == base_pos ) {
							TR_SYMDES_HH.Trace << "sequence_constraints: " << base_pos << ' ' << seqpos << ' ' << name1s << endl;
						}
					}
				}
			}

			/// for xtal design, limit the positions to potential interface positions
			if ( option[ my_options::limit_design_positions_by_chain_contacts ].user() ) {
				Size base_repeat(0);
				conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
				Size const max_design( min( Size( option[ my_options::limit_design_positions_by_chain_contacts ] ),
					repeatlen ) ); // otherwise we get in infinite loop below...
				for ( Size i=1; i<= nrepeat_full; ++i ) {
					if ( symminfo->bb_is_independent( i*repeatlen ) ) {
						runtime_assert( !base_repeat );
						base_repeat = i;
					}
				}
				Size const nres_protein( repeatlen * nrepeat_full );
				// we don't really know what the right number is for this distance threshold...
				vector1< std::pair< Real, Size > > contactlist;
				Real dis_threshold( 7.5 + 3 * uniform() );
				while ( contactlist.empty() ) {
					Real const dis2_threshold( dis_threshold * dis_threshold );// ( 7.5, 10.5 )
					// Real const dis2_threshold( 9 * 9 ); // 9 Angstroms
					for ( Size i=1; i<= repeatlen; ++i ) {
						Size const seqpos( ( base_repeat-1)*repeatlen+i );
						runtime_assert( symminfo->bb_is_independent( seqpos ) );
						Size icontacts(0);
						for ( Size j=1; j<= nres_protein; ++j ) {
							Size const jmon( (j-1)/repeatlen+1 );
							if ( jmon == base_repeat ) continue;
							if ( pose.residue(seqpos).nbr_atom_xyz().distance_squared(pose.residue(j).nbr_atom_xyz())<dis2_threshold ) {
								++icontacts;
							}
						}
						if ( icontacts ) contactlist.push_back( make_pair( Real(icontacts)+0.1*numeric::random::uniform(), i ));
					}
					if ( contactlist.size() < max_design ) { // didnt find enough contacts...
						contactlist.clear(); // go through the loop again
						dis_threshold += 1.0;
					}
				}
				runtime_assert( contactlist.size() >= max_design );
				std::sort( contactlist.begin(), contactlist.end() );
				std::reverse( contactlist.begin(), contactlist.end() );
				bools allow_design( repeatlen, false );
				Size n_allowed(0);
				for ( Size i=1; i<= contactlist.size(); ++i ) {
					Size const rpos( contactlist[i].second );
					Size const seqpos( (base_repeat-1)*repeatlen + rpos );
					string seqcst("-");
					if ( sequence_constraints.count(seqpos) ) seqcst = sequence_constraints.find( seqpos )->second;
					allow_design[ rpos ] = ( n_allowed < max_design && ( seqcst == "-" || seqcst.size() > 1 ) );
					if ( allow_design[ rpos ] ) ++n_allowed;
					TR_SYMDES_HH.Trace << "chain_contacts: " << I(3,i) << F(9,3,contactlist[i].first) <<
						I(4,rpos) << I(5,seqpos) <<
						" allow_design: " << allow_design[rpos] <<
						" seqcst: " << seqcst << endl;
				}
				for ( Size i=1; i<= repeatlen; ++i ) {
					Size const seqpos( ( base_repeat-1)*repeatlen+i );
					runtime_assert( symminfo->bb_is_independent( seqpos ) );
					if ( !allow_design[i] ) {
						task->nonconst_residue_task( seqpos ).restrict_to_repacking();
					}
				}
			} // restrict design positions by inter chain contacts

			if ( option[ basic::options::OptionKeys::relax::limit_aroma_chi2 ] ) {
				protocols::task_operations::LimitAromaChi2Operation lp_op;
				lp_op.chi2min( 55.0 ); // default is 70, seems aggressive (THIS IS STUPID, given adding 180 to negs)
				lp_op.apply( pose, *task );
			}
			Size const nloop( 25 );
			runtime_assert( !option[ OptionKeys::packing::linmem_ig ].user() );
			runtime_assert( pose::symmetry::is_symmetric( pose ) );
			runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );

			ScoreFunctionOP design_scorefxn(0);
			if ( n == n_cycles && option[ my_options::final_round_design_score_function ].user() ) {
				design_scorefxn = setup_score_function( option[ my_options::final_round_design_score_function ] );
				adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
				// always safe to have SequenceBonusEnergy turned on (?)
				design_scorefxn->add_extra_method( sequence_bonus, sequence_bonus_weight, sequence_bonus_energy );
				if ( favor_native_residue ) design_scorefxn->set_weight( res_type_constraint, 1.0 );
			} else if ( option[ my_options::design_score_function ].user() ) {
				design_scorefxn = setup_score_function( option[ my_options::design_score_function ] );
				adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
				// always safe to have SequenceBonusEnergy turned on (?)
				design_scorefxn->add_extra_method( sequence_bonus, sequence_bonus_weight, sequence_bonus_energy );
				if ( favor_native_residue ) design_scorefxn->set_weight( res_type_constraint, 1.0 );
			} else if ( use_softrep_for_design ) {
				design_scorefxn = ScoreFunctionFactory::create_score_function( soft_design_scorefxn_weights_tag );
				adjust_ref_weights_from_command_line( *design_scorefxn ); // since this is brand new
				// always safe to have SequenceBonusEnergy turned on (?)
				design_scorefxn->add_extra_method( sequence_bonus, sequence_bonus_weight, sequence_bonus_energy );
				if ( favor_native_residue ) design_scorefxn->set_weight( res_type_constraint, 1.0 );
			} else {
				runtime_assert( false );
				design_scorefxn = fa_scorefxn; // already adjusted refwts
			}
			protocols::minimization_packing::symmetry::SymPackRotamersMover packmover( design_scorefxn, task, nloop );
			if ( !dry_run() ) {
				packmover.apply( pose );
				if ( pymol ) {
					pymol->apply( pose );
					pymol->print("Finished round "+string_of(n)+" design.");
				}
			}
		}

		++counter;
		outfilename = ( "mid_"+string_of(counter)+"_"+string_of(n)+".pdb");
		if ( debugging ) pose.dump_pdb(outfilename);
		Real const midcyclescore( (*fa_scorefxn)( pose ) );
		TR_SYMDES_HH.Trace << "symmetric_design_and_relax:: midcyclescore " << F(9,3,midcyclescore) << ' ' << outfilename <<
			' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << endl;
		//if ( use_atom_pair_constraints ) pose.constraint_set()->show( TR_SYMDES_HH.Trace );

		{ // now relax
			protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

			MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
			movemap->set_bb ( !freeze_backbone );
			movemap->set_chi(true);
			movemap->set_bb ( pose.total_residue(), false );
			movemap->set_chi( pose.total_residue(), false );
			if ( flex_all_independent_jumps ) { runtime_assert( add_jump_flex ); }

			if ( add_jump_flex ) {
				if ( flex_all_independent_jumps ) {
					FoldTree const & f( pose.fold_tree() );
					conformation::symmetry::SymmetryInfo const & symminfo( *pose::symmetry::symmetry_info( pose ) );
					for ( Size i=1; i<= f.num_jump(); ++i ) {
						if ( symminfo.jump_is_independent(i) ) movemap->set_jump( i, true );
					}
				} else if ( freeze_backbone ) {
					// xtal test
					movemap->set_jump( true );
				} else if ( !pose::symmetry::symmetry_info( pose )->contiguous_monomers() ) {
					// modules setup
					FoldTree const & f( pose.fold_tree() );
					conformation::symmetry::SymmetryInfo const & symminfo( *pose::symmetry::symmetry_info( pose ) );
					for ( Size i=1; i<= f.num_jump(); ++i ) {
						if ( symminfo.jump_is_independent(i) ) movemap->set_jump( i, true );
					}
				} else if ( num_chains( pose )>2 && pose::symmetry::symmetry_info( pose )->num_virtuals() > 1 ) {
					// braid test
					Size const nbraid( num_chains(pose)-1 ), nres_protein( nrepeat_full * repeatlen ),
						nrepeat_per_braid( nrepeat_full/nbraid ), num_virtuals( pose.total_residue() - nres_protein );
					int orientation(0);
					if ( num_virtuals == nrepeat_full ) orientation = 1;
					else if ( num_virtuals == nrepeat_full + 2 ) orientation = -1;
					else utility_exit_with_message("symmetric_design_and_relax: unable to parse braid foldtree/symmetry");
					runtime_assert( chain_end( nbraid,pose) == nres_protein );
					runtime_assert( nrepeat_full%nbraid == 0 );
					runtime_assert( pose::symmetry::symmetry_info( pose )->num_virtuals() == num_virtuals );
					//FoldTree const & f( pose.fold_tree() );
					conformation::symmetry::SymmetryInfo const & symminfo( *pose::symmetry::symmetry_info( pose ) );
					if ( orientation > 0 ) { //parallel
						// flex the independent monomer jump and the independent vertical jump
						Size nflex(0);
						for ( Size i=1; i<= nrepeat_full + nrepeat_per_braid; ++i ) {
							if ( symminfo.jump_is_independent(i) ) {
								movemap->set_jump( i, true );
								TR_SYMDES_HH.Trace << "flex_jump: " << i << ' ' << pose.fold_tree().upstream_jump_residue(i) <<
									" --> " << pose.fold_tree().downstream_jump_residue(i) <<
									" nres_protein: " << nbraid * nrepeat_per_braid * repeatlen << endl;
								++nflex;
							}
						}
						runtime_assert( nflex == 2 );
					} else {
						runtime_assert( orientation == -1 && nbraid == 2 );
						Size nflex(0);
						for ( Size i=1; i<= nrepeat_full + 2*(nrepeat_per_braid-1) + 1; ++i ) {
							if ( symminfo.jump_is_independent(i) ) {
								movemap->set_jump( i, true );
								TR_SYMDES_HH.Trace << "flex_jump: " << i << ' ' << pose.fold_tree().upstream_jump_residue(i) <<
									" --> " << pose.fold_tree().downstream_jump_residue(i) <<
									" nres_protein: " << nres_protein << endl;
								++nflex;
							}
						}
						// if we are targeting a specific twist/rise then the symdofs may not allow some of these to move...
						// this happens inside pose/symmetry/util.cc:make_symmetric_movemap which is called by FastRelax
						// (I think, yeah, must be since the symdof z-angle_z type restrictions are being followed...)
						runtime_assert( nflex == 3 );
					}
				} else if ( pose::symmetry::symmetry_info( pose )->num_virtuals() > 1 ) {
					// symmetric resampling
					FoldTree const & f( pose.fold_tree() );
					conformation::symmetry::SymmetryInfo const & symminfo( *pose::symmetry::symmetry_info( pose ) );
					for ( Size i=1; i<= f.num_jump(); ++i ) {
						if ( symminfo.jump_is_independent(i) && pose.residue( f.downstream_jump_residue(i) ).is_protein() ) {
							TR_SYMDES_HH.Trace << "setting jump flexible: " << i << endl;
							movemap->set_jump( i, true );
						}
					}

				} else {
					Size upper_jumpno, lower_jumpno;
					get_base_repeat_flanking_jumps( pose, lower_jumpno, upper_jumpno, true ); // allow_single_failure
					for ( Size i=lower_jumpno+1; i<= upper_jumpno; ++i ) {
						runtime_assert( pose::symmetry::symmetry_info( pose )->jump_is_independent( i ) );
						movemap->set_jump( i, true );
					}
				}
			}
			{ // status
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					TR_SYMDES_HH.Trace << "symmetric_design_and_relax: " << i <<
						" flexbb: " << movemap->get_bb(i) <<
						" flexchi: " << movemap->get_chi(i) << endl;
				}
				for ( Size i=1; i<= pose.fold_tree().num_jump(); ++i ) {
					TR_SYMDES_HH.Trace << "symmetric_design_and_relax: " << i <<
						" flexjump: " << movemap->get_jump(i) << endl;
				}
			}
			fastrelax.set_movemap( movemap );
			if ( !dry_run() ) {
				fastrelax.apply( pose );
				if ( pymol ) {
					pymol->apply( pose );
					pymol->print("Finished round "+string_of(n)+" relax.");
				}
			}
		}
		Real const endcyclescore( (*fa_scorefxn)( pose ) );
		++counter;
		outfilename = ( "end_"+string_of(counter)+"_"+string_of(n)+".pdb");
		if ( debugging ) pose.dump_pdb(outfilename);
		TR_SYMDES_HH.Trace << "symmetric_design_and_relax:: endcyclescore " << F(9,3,endcyclescore) << ' ' << outfilename <<
			' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << endl;

		bool found_a_large_exposed_nonpolar( false );
		if ( option[ my_options::disallow_large_exposed_nonpolar ] ) { ///////////////////////  exposed nonpolar filter
			conformation::symmetry::SymmetryInfoCOP symminfo( pose::symmetry::symmetry_info( pose ) );
			// ignore edge repeats, maybe
			bools is_edge_repeat( pose.total_residue(), false );
			if ( symminfo->contiguous_monomers() && (pose.total_residue()-1)%repeatlen==0 && !use_donut_energy ) {
				for ( Size i=1; i<= repeatlen; ++i ) {
					is_edge_repeat[i] = is_edge_repeat[ pose.total_residue()-i ] = true;
				}
			}
			bools is_large_exposed_nonpolar;
			find_large_exposed_nonpolar( pose, option[ my_options::exposed_nonpolar_sasa_fraction ], is_large_exposed_nonpolar );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( !( symminfo->bb_is_independent(i) && pose.residue(i).is_protein() ) ) continue;
				Size count(0), total(0);
				for ( Size j=1; j<= pose.total_residue(); ++j ) {
					if ( is_edge_repeat[j] ) continue;
					if ( j==i || symminfo->bb_follows(j) == i ) {
						++total;
						count += ( is_large_exposed_nonpolar[j] );
					}
				}
				TR_SYMDES_HH.Trace << "disallow_large_exposed_nonpolar: " << i << ' ' << pose.residue(i).aa() << ' ' <<
					count << ' ' << total << endl;
				if ( total && Real(count)/total > 0.5 ) { // at least half the positions are large+exposed+nonpolar
					// add sequence constraints at position i
					if ( sequence_constraints.count(i) ) { // already have constraints
						string new_csts, old_csts( sequence_constraints[i] );
						for ( Size j=0; j< old_csts.size(); ++j ) {
							if ( string("FILM").find( old_csts[j] ) == string::npos ) new_csts.push_back( old_csts[j] );
						}
						if ( new_csts.empty() ) {
							sequence_constraints[i] = old_csts; // no luck, this nonpolar is forced !
						} else {
							sequence_constraints[i] = new_csts;
							found_a_large_exposed_nonpolar = true;
						}
					} else {
						sequence_constraints[i] = "DEKRHNQSTY";
						found_a_large_exposed_nonpolar = true;
					}
				}
			}
		}
		if ( option[ my_options::disallow_large_exposed_nonpolar ] && found_a_large_exposed_nonpolar && n == n_cycles ) {
			TR_SYMDES_HH << "found a large nonpolar on the last cycle! " << n << ' ' << n_cycles << endl;
			if ( !dry_run() ) ++n_cycles;
		}

	} // n=1,n_cycles

	Real const endscore( (*fa_scorefxn)( pose ) );
	TR_SYMDES_HH.Trace << "symmetric_design_and_relax:: endscore " << F(9,3,endscore) <<
		' ' << pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << endl;

	if ( use_donut_energy ) {
		// reset the hack
		// note that this means we can't use donut energy anymore without these set!!
		symminfo_hack::nrepeat_ = symminfo_hack::repeatlen_ = 0;
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// taken from xtal.cc, entropy_test

void
symmetric_entropies_calc(
	Size const nres_monomer, // tmp hack: assume indep subunit comes first, is contiguous, runs from 1-->nres_monomer
	ScoreFunctionOP fa_scorefxn,
	bools const & is_flexible,
	Pose const & pose_in,
	bool const include_current,
	Reals & rsd_entropy,
	Reals & rsd_current_rotprob,
	Reals & rsd_current_rotprob_withtol_40,
	Reals & rsd_current_rotprob_withtol_20,
	Reals & snapshot_energies,
	Size const extrachi_cutoff = 0, // to handle bound vs unbound comparisons
	Size const inneriterations_scaling = 250,
	Real const packing_temperature = 0.8,
	bool const start_with_current = true
)
{
	using namespace pack;
	using namespace pack::interaction_graph;
	using namespace pack::rotamer_set;
	using namespace pack::annealer;

	Pose pose( pose_in );

	runtime_assert( pose::symmetry::is_symmetric( pose ) );
	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	rsd_entropy.resize( nres_monomer, 0.0 );
	rsd_current_rotprob.resize( nres_monomer, 0.0 );
	rsd_current_rotprob_withtol_40.resize( nres_monomer, 0.0 );
	rsd_current_rotprob_withtol_20.resize( nres_monomer, 0.0 );
	snapshot_energies.clear();


	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ) );
	bools allow_repack( pose.total_residue(), false );
	for ( Size i=1; i<= nres_monomer; ++i ) {
		if ( is_flexible[i] && pose::symmetry::symmetry_info( pose )->bb_is_independent(i) ) {
			allow_repack[i] = true;
		}
	}
	task->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
	task->or_include_current( include_current );
	task->set_bump_check( false );

	// extrachi_cutoff -- this is important if we are doing a bound vs unbound comparison where nbr number can change!!!
	// the below only does something if the passed in value is LOWER than the cmdline/default value
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( allow_repack[i] ) task->nonconst_residue_task(i).and_extrachi_cutoff( extrachi_cutoff );
	}


	task::PackerTaskCOP non_symmetric_task( task->clone() );
	task = make_new_symmetric_PackerTask_by_requested_method( pose, non_symmetric_task ); // what does this do?

	/// would do next line in asym case
	// rotamer_set::RotamerSetsOP rotsets( new rotamer_set::RotamerSets() );
	rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets( new rotamer_set::symmetry::SymmetricRotamerSets() );

	InteractionGraphBaseOP ig = NULL;

	pack_scorefxn_pose_handshake( pose, *fa_scorefxn ); // this just scores the pose

	pose.update_residue_neighbors();

	fa_scorefxn->setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );

	{ // new
		runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );
		runtime_assert( utility::pointer::dynamic_pointer_cast<scoring::symmetry::SymmetricScoreFunction >(fa_scorefxn));
		scoring::symmetry::SymmetricScoreFunction const & symm_scfxn
			( dynamic_cast< scoring::symmetry::SymmetricScoreFunction const & >( *fa_scorefxn ) );
		symm_scfxn.correct_arrays_for_symmetry( pose );
	}


	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, *fa_scorefxn, task );

	rotsets->set_task( task );
	rotsets->initialize_pose_for_rotsets_creation(pose); // new, calculates Tsymm_, see e.g. pack_rotamers_setup
	rotsets->build_rotamers( pose, *fa_scorefxn, packer_neighbor_graph );
	rotsets->prepare_sets_for_packing( pose, *fa_scorefxn );

	ig = InteractionGraphFactory::create_interaction_graph( *task, *rotsets, pose, *fa_scorefxn, *packer_neighbor_graph );

	TR_SYMDES_HH.Trace << "built " << rotsets->nrotamers() << " rotamers at " <<
		rotsets->nmoltenres() << " positions." << std::endl;

	for ( Size i=1; i<= pose.total_residue(); ++i ) { runtime_assert( !task->design_residue(i) ); }

	rotsets->compute_energies( pose, *fa_scorefxn, packer_neighbor_graph, ig );

	TR_SYMDES_HH.Trace << "IG: " << ig->getTotalMemoryUsage() << " bytes" << std::endl;

	/// NOW the pack_rotamers_run code:
	utility::vector0< int > rot_to_pack;
	ObjexxFCL::FArray1D_int bestrotamer_at_seqpos( pose.total_residue() );
	core::PackerEnergy bestenergy;

	//bool start_with_current = true; //option[ my_options::swc ];
	ObjexxFCL::FArray1D_int current_rot_index( rotsets->nmoltenres(), -1 );

	// fill out the current_rot_index array
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Size const mres( rotsets->resid_2_moltenres(i) );
		if ( mres > 0 ) {
			Size const nrotamer_offset( rotsets->nrotamer_offset_for_moltenres( mres ) );
			Size rotset_current_id(0);
			RotamerSet const & rotset( *rotsets->rotamer_set_for_moltenresidue( mres ) );
			if ( include_current ) {
				rotset_current_id = rotset.id_for_current_rotamer();
			} else {
				// look for the closest rotamer to the designed rotamer
				Residue const & rsd( pose.residue(i) );
				Real min_angle_dev(1e6);
				for ( Size j=1; j<= rotset.num_rotamers(); ++j ) {
					Real const dev( rotamer_chi_angles_dev( rsd, *rotset.rotamer(j) ) );
					if ( dev < min_angle_dev ) {
						min_angle_dev = dev;
						rotset_current_id = j;
					}
				}
				TR_SYMDES_HH.Trace << "best_rotamer_match: " << F(9,3,min_angle_dev) << I(4,i) << ' ' << rsd.name1() << ' ' <<
					rotset_current_id << endl;
			}

			runtime_assert( rotset_current_id );
			current_rot_index(mres) = nrotamer_offset + rotset_current_id;
			TR_SYMDES_HH.Trace << "current_rot_index: " << i << ' ' << mres << ' ' << nrotamer_offset << ' ' <<
				rotset_current_id << ' ' << current_rot_index(mres) << endl;
		}
	}


	bool calc_rot_freq = true;
	ObjexxFCL::FArray1D< core::PackerEnergy > rot_freq( ig->get_num_total_states(), 0.0 );

	SimAnnealerBaseOP annealer_op = AnnealerFactory::create_annealer(
		task, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
		rotsets, current_rot_index, calc_rot_freq, rot_freq );

	/// this is really cheesy:
	runtime_assert( utility::pointer::dynamic_pointer_cast< FixbbSimAnnealer >( annealer_op ) );
	FixbbSimAnnealer & annealer( dynamic_cast< FixbbSimAnnealer & >( *annealer_op ) );

	// calc_freq_temp is hard-coded to 1.0 in SimAnnealerBase.cc
	annealer.set_lowtemp ( packing_temperature );
	annealer.set_hightemp( packing_temperature );
	annealer.set_disallow_quench( true );//task->disallow_quench() );
	annealer.n_snapshots( 100 );
	annealer.scale_inneriterations( dry_run() ? 10 : inneriterations_scaling ); // dflt is 250

	annealer.run();

	Reals const & packer_snapshot_energies( annealer.snapshot_energies() );
	vector1< vector1< int > > const & snapshot_states( annealer.snapshot_states() );
	Size const n_snapshots_actual( packer_snapshot_energies.size() );

	/// the rotamer numbers in the snapshot_states are based on the "state_on_node" array in FixbbSimAnnealer
	/// hence they are numbered internally with respect to each molten residue's rotamer set
	//Real prev_edev(0);
	for ( Size i=1; i<= n_snapshots_actual; ++i ) {
		vector1< int > const & i_state( snapshot_states[i] );
		bool all_assigned( true );
		foreach_ ( int j, i_state ) all_assigned = ( all_assigned && j != 0 );
		if ( !all_assigned ) continue;
		runtime_assert( i_state.size() == rotsets->nmoltenres() );
		// is this the same as the last state?
		ostringstream state_changes;
		if ( i>1 ) {
			vector1< int > const & prev_state( snapshot_states[i-1] );
			for ( Size j=1; j<= i_state.size(); ++j ) {
				if ( i_state[j] != prev_state[j] ) {
					state_changes << ' ' << rotsets->moltenres_2_resid( j ) << ':' << prev_state[j] << ':' <<
						i_state[j];
					// state_changes.push_back( rotsets->moltenres_2_resid( j ) );
				}
			}
			if ( state_changes.str().empty() ) continue;
		}

		// copy this state to the pose
		for ( Size mres=1; mres<= rotsets->nmoltenres(); ++mres ) {
			Size const resid( rotsets->moltenres_2_resid( mres ) );
			if ( resid > nres_monomer ) continue;
			Size const nrotamer_offset( rotsets->nrotamer_offset_for_moltenres( mres ) );
			Residue const & rot_at_resid_in_i_state( *rotsets->rotamer( nrotamer_offset + i_state[mres] ) );
			pose.replace_residue( resid, rot_at_resid_in_i_state, false ); // will symm-replace
		}
		Real const state_score( (*fa_scorefxn)( pose ) ), state_packer_energy( packer_snapshot_energies[i] ),
			edev( state_score - state_packer_energy );

		TR_SYMDES_HH.Trace << "B_snapshot: " << I(4,i) << " state_packer_energy: " << F(9,3,state_packer_energy) <<
			" state_score: " << F(9,3,state_score) <<
			" edev: " << F(9,3,edev) <<
			" state_changes: " << state_changes.str() << endl;
		snapshot_energies.push_back( state_score );
		//prev_edev = edev;
	}

	// now show the rotamer probabilities
	Size const nsteps_approximate( annealer.get_inneriterations() * annealer.get_outeriterations() );

	for ( Size mres=1; mres<= rotsets->nmoltenres(); ++mres ) {
		Size const resid( rotsets->moltenres_2_resid( mres ) );
		if ( resid > nres_monomer ) continue;
		//Residue const & rsd( pose.residue(resid) );
		Size const nrotamer_offset( rotsets->nrotamer_offset_for_moltenres( mres ) );
		Size const nrot( rotsets->nrotamers_for_moltenres( mres ) );
		Size const current_rotid( current_rot_index(mres) - nrotamer_offset );
		// Size const current_rotid( rotsets->rotamer_set_for_moltenresidue( mres )->id_for_current_rotamer() );
		runtime_assert( current_rotid );
		rsd_entropy[ resid ] = rsd_current_rotprob[ resid ] = rsd_current_rotprob_withtol_40[ resid ] =
			rsd_current_rotprob_withtol_20[ resid ] = 0.0;
		Residue const & current_rot( pose_in.residue( resid ) );
		Residue const & current_rot_approx( *rotsets->rotamer( nrotamer_offset + current_rotid ) );
		TR_SYMDES_HH.Trace << "current_rot_approx: " << F(9,3,rotamer_chi_angles_dev( current_rot, current_rot_approx ) ) <<
			I(4,resid) << endl;

		for ( Size i=1; i<= nrot; ++i ) {
			Size const rotid( nrotamer_offset + i );
			Residue const & rot( *rotsets->rotamer( rotid ) );
			Real const nsteps_with_this_rot( rot_freq(rotid)*nsteps_approximate );
			if ( nsteps_with_this_rot > 0.5 ) rsd_entropy[resid] -= rot_freq(rotid) * log( rot_freq(rotid) );
			if ( i == current_rotid ) rsd_current_rotprob[resid] = rot_freq(rotid);
			if ( rotamer_chi_angles_match( rot, current_rot, 40.0 ) ) {
				rsd_current_rotprob_withtol_40[resid] += rot_freq(rotid);
			}
			if ( rotamer_chi_angles_match( rot, current_rot, 20.0 ) ) {
				rsd_current_rotprob_withtol_20[resid] += rot_freq(rotid);
			}
		}
	}
} // scope

//Size const current_rotid( current_rot_index(mres) - nrotamer_offset );

// figure out how the rotamers for this position break down into the different
//RotamerSet const & rotset( *rotsets->rotamer_set_for_residue(resid) );
// Size const n_restypes( rotset.get_n_residue_types() );
// vector1< Sizes > rotids_for_aa( 20 );
// for ( Size i_restype=1; i_restype<= n_restypes; ++i_restype ) {
// runtime_assert( current_rotid );
// unbound_entropy[ resid ] = 0.0;
// Residue const & current_rot( input_pose.residue( resid ) );
// Residue const & current_rot_approx( *rotsets->rotamer( nrotamer_offset + current_rotid ) );
// TR_SYMDES_HH.Trace << "current_rot_approx: " << F(9,3,rotamer_chi_angles_dev( current_rot, current_rot_approx ) ) <<
//  I(4,resid) << endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
