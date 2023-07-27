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
#ifndef INCLUDED_apps_pilot_phil_fragments_HH
#define INCLUDED_apps_pilot_phil_fragments_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/types.hh>
#include <apps/pilot/phil/phil_options.hh>
#include <apps/pilot/phil/phil_io.hh>

#include <core/scoring/rms_util.hh>
#include <core/pose/util.hh>
#include <core/fragment/ConstantLengthFragSet.hh>


#include <devel/blab/classic_frags/TorsionFragment.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <fstream>

// static numeric::random::RandomGenerator RG(12321); // <- Magic number, do not change it!!!

static basic::Tracer TR_FRAGMENTS_HH( "apps.pilot.phil.fragments_hh" );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
optionally_cleanup_tmp()
{
	if ( option[ my_options::cleanup_tmp ] ) {
		string const cmd("python /home/pbradley/python/cleanup_tmp_pbradley.py");
		run_command( cmd );
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pick_nnmake_fragments_from_single_sequence(
	string const & seq,
	core::fragment::FragSetOP & small_frags,
	core::fragment::FragSetOP & large_frags,
	bool const use_nohoms = false
)
{
	optionally_cleanup_tmp();

	string const prefix( option[ my_options::tmpdir ]()+"/pbradley_"+filebase( output_tag())+"_"+string_of(uniform())),
		small_frags_filename( prefix+"_3mers.txt" ),
		large_frags_filename( prefix+"_9mers.txt" ),
		logfile( prefix+"_log.txt" ),
		errfile( prefix+"_err.txt" );

	/// NOTE: we are using -nohoms flag:
	///
	string extra_flags;
	if ( use_nohoms ) extra_flags += " -nohoms ";
	//if ( option[ my_options::ssblast ] ) extra_flags += " -ssblast ";
	string const cmd( "python /home/pbradley/nnmake/create_fragment_files_from_sequence.py "+
		small_frags_filename+" "+large_frags_filename +" "+seq+" "+extra_flags+
		" > " + logfile + " 2> " + errfile );
	run_command( cmd );

	core::fragment::ConstantLengthFragSetOP small_frags_local( new core::fragment::ConstantLengthFragSet );
	small_frags_local->read_fragment_file( small_frags_filename,
		option[ OptionKeys::abinitio::number_3mer_frags ](), 1, false );

	core::fragment::ConstantLengthFragSetOP large_frags_local( new core::fragment::ConstantLengthFragSet );
	large_frags_local->read_fragment_file( large_frags_filename,
		option[ OptionKeys::abinitio::number_9mer_frags ](), 1, false );

	small_frags = small_frags_local;
	large_frags = large_frags_local;

	run_command( "rm "+prefix+"*" );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
pick_nnmake_fragments_from_single_sequence(
	string const & seq,
	devel::blab::classic_frags::FragLib & fraglib,
	bool const use_nohoms_flag = false
)
{
	optionally_cleanup_tmp();

	string const prefix( option[ my_options::tmpdir ]()+"/pbradley_"+filebase( output_tag())+"_"+string_of(uniform())),
		small_frags_filename( prefix+"_3mers.txt" ),
		large_frags_filename( prefix+"_9mers.txt" ),
		logfile( prefix+"_log.txt" ),
		errfile( prefix+"_err.txt" );

	/// NOTE: we are using -nohoms flag:
	///
	string extra_flags;
	if ( use_nohoms_flag ) extra_flags += " -nohoms ";
	//if ( option[ my_options::ssblast ] ) extra_flags += " -ssblast ";
	string const cmd( "python /home/pbradley/nnmake/create_fragment_files_from_sequence.py "+
		small_frags_filename+" "+large_frags_filename +" "+seq+" "+extra_flags+
		" > " + logfile + " 2> " + errfile );
	run_command( cmd );

	if ( !utility::file::file_exists( small_frags_filename ) || !utility::file::file_exists( small_frags_filename ) ) {
		cout << "nnmake failed " << get_hostname() << endl;
		run_command( "date" );
		run_command( "cat "+logfile );
		run_command( "cat "+errfile );
	} else {
		fraglib.library( 3 ).read_file( small_frags_filename, 3, 3 );
		fraglib.library( 9 ).read_file( large_frags_filename, 9, 3 );
	}

	run_command( "rm "+prefix+"*" );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
run_psipred(
	string const seq,
	string & ss,
	vector1< Reals > & pred_HEL
)
{
	optionally_cleanup_tmp();

	string prefix( option[ my_options::tmpdir ]() +"/psipred_tmpfile_"+filebase( output_tag() ) + string_of( uniform() ) ),
		fastafile( prefix+".fasta" ),
		logfile  ( prefix+".log"   ),
		errfile  ( prefix+".err"   ),
		ssfile   ( prefix+".ss"    ),
		ss2file  ( prefix+".ss2"   ),
		horizfile( prefix+".horiz" );

	// silly: psipred puts them in current working directory
	ssfile = filebase( ssfile );
	ss2file = filebase( ss2file );
	horizfile = filebase( horizfile );

	string psipred_directory( option[ my_options::psipred_directory ]() );
	if ( psipred_directory[-1] != '/' ) psipred_directory += '/';
	string const psipredexe( psipred_directory+"runpsipred_single" );

	if ( !utility::file::file_exists( psipredexe ) ) {
		utility_exit_with_message("cant find psipredexe "+psipredexe );
	}

	ofstream out( fastafile.c_str() );
	out << ">tmp\n" << seq << '\n';
	out.close();

	run_command( psipredexe + " " + fastafile + " > " + logfile + " 2> " + errfile );

	if ( !utility::file::file_exists( ss2file ) ) {
		cout << "psipred failed " << logfile << endl;
		cerr << "psipred failed " << logfile << endl;
		cout << "psipred failedseq " << seq << endl;
		run_command( "hostname" );
		run_command( "cat " + logfile );
		run_command( "cat " + errfile );
		run_command( "cat " + ssfile );
		run_command( "cat " + horizfile );
		ss.clear();
		pred_HEL.clear();
		return;
	}


	string line;
	ifstream data( ss2file.c_str() );
	// george1 symdes$ head avrb3_xaneu.ss2
	// # PSIPRED VFORMAT (PSIPRED V3.3)
	//
	//    1 M C   1.000  0.000  0.000
	getline( data, line );
	getline( data, line );
	//  104 W H   0.362  0.621  0.023
	//  110 G C   0.534  0.429  0.020
	//  121 M E   0.342  0.187  0.512
	ss.clear();
	pred_HEL.clear();

	for ( Size i=1; i<= seq.size(); ++i ) {
		getline( data, line );
		istringstream l( line );
		Size pos;
		char seq1, ss1;
		Real pred_L, pred_H, pred_E;
		l >> pos >> seq1 >> ss1 >> pred_L >> pred_H >> pred_E;
		if ( l.fail() ) {
			cout << "bad psipred line: " << line << endl;
			break;
		}
		if ( ss1 == 'C' ) ss1 = 'L';

		ss.push_back( ss1 );
		pred_HEL.push_back( make_vector1( pred_H, pred_E, pred_L ) );
	}

	data.close();

	if ( ss.size() != seq.size() ) {
		cout << "psipred failed " << logfile << endl;
		cerr << "psipred failed " << logfile << endl;
		ss.clear();
		pred_HEL.clear();
		return;
	}

	utility::file::file_delete( fastafile );
	utility::file::file_delete( logfile   );
	utility::file::file_delete( errfile   );
	utility::file::file_delete( ssfile    );
	utility::file::file_delete( ss2file   );
	utility::file::file_delete( horizfile );

}
/// alt version
void
run_psipred(
	string const seq,
	string & ss,
	vector1< map< char, Real > > & pred_ss
)
{
	vector1< Reals > pred_HEL;
	run_psipred( seq, ss, pred_HEL );

	pred_ss.clear();
	for ( Size i=1; i<= seq.size(); ++i ) {
		map< char, Real > sspred;
		sspred[ 'H' ] = pred_HEL[i][1];
		sspred[ 'E' ] = pred_HEL[i][2];
		sspred[ 'L' ] = pred_HEL[i][3];
		pred_ss.push_back( sspred );
	}
}

///////////////////////////////////////////////////////////////////////////////////////


/// 200 9-mers, 3-mers, and 1-mers
///
string const empty_ss;

devel::blab::classic_frags::FragLibOP
setup_cheating_fragments(
	Pose const & pose,
	Sizes const & vall_homs,
	string const & ss_in = empty_ss,
	bools const & is_flexible = empty_bools,
	Size const big_frag_size = 9
)
{
	Real const min_torsion_dev( option[ my_options::cheating_frags_min_torsion_dev ] ),
		max_torsion_dev( option[ my_options::cheating_frags_max_torsion_dev ] ); // for fragment picking...
	string ss( ss_in );
	if ( ss.empty() ) ss = pose.secstruct();
	{ // confirm that ss has been initialized
		bool all_L( true );
		for ( Size i=0; i<ss.size(); ++i ) if ( ss[i] != 'L' ) all_L = false;
		if ( all_L ) utility_exit_with_message("Need to initialize pose's secondary structure.");
	}

	devel::blab::classic_frags::FragLibOP protein_fraglib( 0 );
	kinematics::MoveMap mm; mm.set_bb( true );
	if ( is_flexible.size() > 0 ) {
		for ( Size i=1; i<= is_flexible.size(); ++i ) {
			if ( !is_flexible[i] ) mm.set_bb( i, false );
		}
	}
	vector1< Size > frag_sizes;
	frag_sizes.push_back( big_frag_size );
	frag_sizes.push_back( 3 );
	Size const nfrags( 200 );
	Real const seq_weight( 1.0 ), ss_weight( 3.0 ), torsion_weight( 5 );

	TR_FRAGMENTS_HH.Trace << "setup_cheating_fragments: ss= " << ss << endl;

	protein_fraglib =
		devel::blab::classic_frags::setup_vall_cheating_fragments( frag_sizes, nfrags, pose, mm, ss,
		seq_weight, ss_weight, torsion_weight,
		min_torsion_dev, max_torsion_dev,
		vall_homs );
	// gets 1mers from 3mers
	protein_fraglib->library( 1 ).derive_from_src_lib( 1, 3, protein_fraglib->library( 3 ) );

	return protein_fraglib;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// each 9-mer is considered as covering the 3 central positions (4-6)
// for each position in the repeat, we compute the fraction of covering 9mers that have rmsd < 1
// thats the sub1_fraction
//
void
get_fragment_quality(
	Size const nrepeat,
	Size const repeatlen,
	devel::blab::classic_frags::FragLib const & fraglib,
	Pose const & pose,
	string & fragss,
	string & fragbb,
	Reals & bb1_recovery,
	Reals & bb3_recovery,
	Reals & sub1_fraction,
	vector1< Reals > & frag_rmsds // 9mers
)
{

	Sizes const fragsizes( fraglib.frag_sizes() );

	runtime_assert( has_element( fragsizes, Size(3) ) );
	runtime_assert( has_element( fragsizes, Size(9) ) );

	string const fullbb( torsion2big_bin_string( 1, pose.total_residue(), pose, true ) );

	Sizes good_bb1_count( repeatlen, 0 ), total_bb1_count( repeatlen, 0 );
	Sizes good_bb3_count( repeatlen, 0 ), total_bb3_count( repeatlen, 0 );
	Sizes sub1_count( repeatlen, 0 ), sub1_total( repeatlen, 0 );

	vector1< map< char, Size > > ss_counts( repeatlen );
	vector1< map< char, Size > > bb_counts( repeatlen );

	for ( Sizes::const_iterator sz= fragsizes.begin(); sz != fragsizes.end(); ++sz ) {
		Size const fragsize( *sz );

		Pose fragpose;
		for ( Size i=1; i<= fragsize+2; ++i ) {
			bool const build_ideal_geometry( true );
			fragpose.append_residue_by_bond( *(get_vanilla_protein_residue('A') ), build_ideal_geometry );
		}
		Size const insertpos( 2 );

		frag_rmsds.clear(); frag_rmsds.resize( repeatlen );

		devel::blab::classic_frags::TorsionFragmentLibrary const & lib( fraglib.library(fragsize) );

		for ( Size fragbegin=1; fragbegin<= lib.size() && fragbegin+fragsize-1 <= pose.total_residue(); ++fragbegin ) {
			Size const nn( lib[fragbegin].size() );
			for ( Size n=1; n<= nn; ++n ) {
				devel::blab::classic_frags::TorsionFragment const & frag( lib[ fragbegin ][n] );
				string frag_fullbb;
				for ( Size i=1; i<= fragsize; ++i ) {
					char const fbb( torsion2big_bin( frag.get_torsion( i, 1 ),
						frag.get_torsion( i, 2 ),
						frag.get_torsion( i, 3 ) ) );
					frag_fullbb.push_back( fbb );
					Size const seqpos( fragbegin+i-1 ), repeatpos( (seqpos-1)%repeatlen+1 );
					if ( nrepeat>1 && ( seqpos == 1 || seqpos == nrepeat*repeatlen ) ) continue; // bb weird

					++ss_counts[ repeatpos ][ frag.get_secstruct(i) ];
					++bb_counts[ repeatpos ][ fbb ];

					++total_bb1_count[ repeatpos ];
					if ( fbb == torsion2big_bin( pose.phi(seqpos), pose.psi(seqpos), pose.omega(seqpos) ) ) {
						++good_bb1_count[ repeatpos ];
					}
				}
				/// now compute rmsd
				if ( fragsize == 9 ) {
					frag.insert( fragpose, insertpos );
					using namespace core::id;
					AtomID_Map< AtomID > atom_map;
					initialize_atomid_map( atom_map, fragpose, id::GLOBAL_BOGUS_ATOM_ID );
					for ( Size i=0; i< fragsize; ++i ) {
						atom_map[ AtomID( fragpose.residue( insertpos+i ).atom_index("N"), insertpos+i ) ] =
							AtomID( pose.residue( fragbegin+i ).atom_index("N"), fragbegin+i );
						atom_map[ AtomID( fragpose.residue( insertpos+i ).atom_index("CA"), insertpos+i ) ] =
							AtomID( pose.residue( fragbegin+i ).atom_index("CA"), fragbegin+i );
						atom_map[ AtomID( fragpose.residue( insertpos+i ).atom_index("C"), insertpos+i ) ] =
							AtomID( pose.residue( fragbegin+i ).atom_index("C"), fragbegin+i );
					}
					Real const rmsd( rmsd_by_mapping( fragpose, pose, atom_map ) );
					Size const midpoint( fragbegin+4 ), repeatpos( (midpoint-1)%repeatlen+1 );
					frag_rmsds[ repeatpos ].push_back( rmsd );
					for ( Size fpos=4; fpos <= 6; ++fpos ) {
						Size const coveredpos( fragbegin+fpos-1 ), covered_repeatpos( ( coveredpos-1)%repeatlen+1 );
						++sub1_total[ covered_repeatpos ];
						if ( rmsd <= 1 ) ++sub1_count[ covered_repeatpos ];
					}
				}

				if ( fragsize == 3 ) {
					/// count fragments that have all three bb's correct
					Size const centerpos( fragbegin+1 ), repeatpos( ( centerpos-1)%repeatlen + 1 );
					if ( fragbegin == 1 || fragbegin+2 >= nrepeat*repeatlen ) continue; // includes terminal pos
					string const posebb( torsion2big_bin_string( 1, pose.total_residue(), pose, true ).substr(fragbegin-1,3));
					++total_bb3_count[ repeatpos ];
					if ( frag_fullbb == posebb ) ++good_bb3_count[ repeatpos ];
				}

			} // nn
		} // fragbegin
	} // fragsizes

	// hack for single repeat
	if ( nrepeat==1 ) {
		for ( Size i=1; i<= repeatlen; ++i ) {
			if ( !sub1_total[i] ) ++sub1_total[i];
			if ( !total_bb1_count[i] ) ++total_bb1_count[i];
			if ( !total_bb3_count[i] ) ++total_bb3_count[i];
		}
	}


	//
	sub1_fraction.clear();
	for ( Size i=1; i<= repeatlen; ++i ) {
		runtime_assert( sub1_total[i] );
		sub1_fraction.push_back( Real( sub1_count[i] )/sub1_total[i] );
	}


	bb1_recovery.clear();
	bb3_recovery.clear();
	fragss.clear();
	fragbb.clear();
	for ( Size i=1; i<= repeatlen; ++i ) {
		bb1_recovery.push_back( Real( good_bb1_count[i] )/total_bb1_count[i] );
		bb3_recovery.push_back( Real( good_bb3_count[i] )/total_bb3_count[i] );

		{
			Real H( ss_counts[i]['H'] ), E( ss_counts[i]['E'] ), L( ss_counts[i]['L'] ), total( H+E+L );
			H /= total; E/= total; L/= total;
			if ( H > E && H > L ) {
				if ( H>0.666 ) fragss += 'H';
				else fragss += 'h';
			} else if ( E > H && E > L ) {
				if ( E>0.666 ) fragss += 'E';
				else fragss += 'e';
			} else {
				if ( L>0.666 ) fragss += 'L';
				else fragss += 'l';
			}
		}

		{ /// bb
			string const abego("ABEGO");
			Size total(0),maxcount(0);
			for ( Size j=0; j< abego.size(); ++j ) {
				Size const count( bb_counts[i][abego[j]] );
				total += count;
				maxcount = max( maxcount, count );
			}
			for ( Size j=0; j< abego.size(); ++j ) {
				Size const count( bb_counts[i][abego[j]] );
				if ( count == maxcount ) { // best bb
					if ( Real(count)/total > 0.666 ) {
						fragbb += abego[j];
					} else {
						fragbb += ObjexxFCL::lowercased( abego.substr(j,1) );
					}
					break; // take first, if tied
				}
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
pick_vall_fragments_from_ss_and_bb(
	Pose const & pose, // unused I think but needed for add_vall_fragments
	string const & fullss,
	string const & fullbb,
	devel::blab::classic_frags::FragLib & fraglib
)
{
	runtime_assert( pose.total_residue() >= fullss.size() );

	Real const seq_weight( 0.0 ), ss_weight( 10.0 ), bb_weight( 100.0 );
	Sizes const fragsizes( make_vector1( 3, 9 ) );
	Size const nfrags( 200 );
	kinematics::MoveMap mm;
	for ( Size i=1; i<= fullss.size(); ++i ) mm.set_bb( i, true );

	Sizes const homs_to_exclude;
	TR_FRAGMENTS_HH.Trace << "add_vall_fragments: " << nfrags << " seq_weight: " << seq_weight <<
		" ss_weight: " << ss_weight << " bb_weight: " << bb_weight << " fullss.size(): " << fullss.size() << endl;
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
				if ( fragpos > fullss.size() ) continue;
				for ( Size nn=lib[ fragpos ].size(); nn> min_nn; --nn ) {
					bool badfrag( false );

					devel::blab::classic_frags::TorsionFragment const frag( lib[ fragpos ][ nn ] );

					for ( Size k=1; k<= fragsize; ++k ) {
						//char const ss( frag.get_secstruct(k) );
						Real const phi( frag.get_torsion(k,1 ) ), psi( frag.get_torsion( k, 2 ) ), omega( frag.get_torsion(k,3 ) );
						char const bigbin( torsion2big_bin( phi, psi, omega ) );
						Size const pos( fragpos +k-1 );
						if ( pos > fullbb.size() ) continue;
						char const desired_bigbin( fullbb[ pos-1 ] );
						// Size const repeatpos( ( ( fragpos +k-1 ) -1 )%repeatlen + 1 );
						// char const desired_bigbin( repeatbb[ repeatpos-1 ] );
						if ( desired_bigbin != 'X' && desired_bigbin != bigbin ) {
							runtime_assert( string("ABGEO").find( desired_bigbin ) != string::npos );
							badfrag = true;
							break;
						}
					}

					if ( badfrag ) {
						TR_FRAGMENTS_HH.Trace << "delete_fragment: " << I(4,fragsize) << I(4,fragpos) << I(4,nn) << endl;
						lib[ fragpos ].erase( nn );
					}
				} // nn
				TR_FRAGMENTS_HH.Trace << "final_nn: " << I(4,fragsize) << I(4,fragpos) <<
					I(4,lib[ fragpos ].size()) << endl;
			} // fragpos
		} // fragsize
	} // scope
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
