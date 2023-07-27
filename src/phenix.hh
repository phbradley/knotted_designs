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
#ifndef INCLUDED_apps_pilot_phil_phenix_HH
#define INCLUDED_apps_pilot_phil_phenix_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/phil_io.hh>
#include <apps/pilot/phil/phil_options.hh>
#include <utility/file/file_sys_util.hh>
// #include <core/scoring/rms_util.hh>
// #include <core/pose/util.hh>
// #include <numeric/model_quality/rms.hh>
// #include <numeric/model_quality/maxsub.hh>

#include <basic/Tracer.hh>
// #include <fstream>


static basic::Tracer TR_PHENIX( "apps.pilot.phil.phenix_hh" );



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
confirm_phenix_version()
{
	// if ( option[ my_options::skip_phenix_version_check ] ) return;
	//string const expected_version( phenix_version() );
	char * p;
	p = getenv("PHENIX_VERSION");
	if ( p == NULL ) utility_exit_with_message("environment variable PHENIX_VERSION is not set");
	string const actual_version( p );
	cout << "PHENIX_VERSION= " << actual_version << endl;
	// if ( expected_version != actual_version ) {
	//  utility_exit_with_message("PHENIX_VERSION mismatch: "+expected_version+" "+actual_version);
	// }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
run_phaser_and_phenix_refinement_generic(
	string const model_file,
	string const mtz_file,
	string const fasta_file,
	string const target_space_group,
	Size const target_num_models,
	Size const composition,
	Real const llg_threshold_for_refinement,
	Size & actual_num_models,
	Real & llg,
	Real & tfz,
	Real & rwork,
	Real & rfree,
	string & actual_space_group,
	string & refined_pdbfile
)
{
	runtime_assert( utility::file::file_exists( fasta_file ) );
	runtime_assert( utility::file::file_exists( model_file ) );
	runtime_assert( utility::file::file_exists( mtz_file ) );

	string cmd( string("python /home/pbradley/python/phenix/run_phaser_and_phenix_refinement_generic.py ")+
		" --model_file "+model_file+
		" --mtz_file "+mtz_file+
		" --fasta_file "+fasta_file+
		" --space_group "+target_space_group+
		" --num_models "+string_of(target_num_models)+
		" --composition "+string_of(composition)+
		" --llg_threshold "+string_of(llg_threshold_for_refinement) );

	if ( dry_run() ) cmd += " --dry_run ";
	if ( option[ my_options::skip_tncs ] ) cmd += " --skip_tncs ";
	if ( option[ my_options::sgalt_all ] ) cmd += " --sgalt_all ";
	if ( option[ my_options::sgalt_none ] ) cmd += " --sgalt_none ";
	if ( option[ my_options::full_search ] ) cmd += " --full_search ";
	if ( option[ my_options::generate_r_free_flags ] ) cmd += " --generate_r_free_flags ";
	if ( option[ my_options::autobuild_rwork_threshold ].user() ) {
		cmd += " --autobuild --autobuild_rwork_threshold "+string_of( option[ my_options::autobuild_rwork_threshold ] );
	}
	run_command( cmd );

	string const scoresfile( model_file +"_MR.scores" );
	if ( utility::file::file_exists( scoresfile ) ) {
		ifstream data( scoresfile.c_str() );
		string line, tmp;
		getline( data, line );
		istringstream l(line );
		l >> tmp >> llg >> tmp >> tfz >> tmp >> rwork >> tmp >> rfree >> tmp >> actual_space_group >> tmp >> actual_num_models;
		if ( l.fail() ) {
			rwork = rfree = 1.0;
			llg = tfz = 0.0;
			actual_space_group = "PARSE_ERROR";
		}
		data.close();
		refined_pdbfile = model_file +"_MR_refine_001.pdb";
	} else {
		rwork = rfree = 1.0;
		llg = tfz  = 0.0;
		actual_space_group = "NO_SCORESFILE";
		refined_pdbfile  = "NONE";
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
