// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
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

#ifndef INCLUDED_apps_pilot_phil_phil_io_HH
#define INCLUDED_apps_pilot_phil_phil_io_HH

#include <apps/pilot/phil/phil.hh>
#include <apps/pilot/phil/safe_unix_file_io.hh>
#include <devel/blab/loops/fragments.hh>
#include <devel/blab/motif/MotifData.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <basic/prof.hh>

basic::Tracer TR_IO( "apps.pilot.phil.phil_io" );


// bad form: global data
namespace cmdline {
string cmdline_;
}

namespace decoy_timer {
clock_t starttime;
map<string,Size> decoy_counter;
string last_filename;
}




void
set_cmdline(
	int argc,
	char * argv []
)
{
	ostringstream out;
	for ( int i=0; i< argc; ++i ) out << ' ' << argv[i];
	cmdline::cmdline_ = out.str().substr(1); // remove leading ' '
}

string
get_cmdline()
{
	return cmdline::cmdline_;
}

void
run_command( string const & cmd, bool const terse = false )
{
	if ( !terse ) {
		std::cout << "Run: " << cmd << std::endl;
		fflush( stdout );
	}
	int const retval( system( cmd.c_str() ) );
	if ( !terse ) TR_IO.Trace << "run_command: retval= "<< retval << " cmd= " << cmd << std::endl;
}



string
run_command_and_get_output( string const & cmd )
{
	FILE* cmdpopen = popen( cmd.c_str(),"r");
	char buf[5000];
	if ( fgets(buf,5000,cmdpopen) == NULL ) {
		string hostname("unknown");
		if ( cmd != "hostname" ) { // no recursion
			hostname = run_command_and_get_output( "hostname" );
		}
		std::cout << "run_command_and_get_output: NULL pointer returned by fgets. cmd= " << cmd << " host= " << hostname <<
			std::endl;
		std::cerr << "run_command_and_get_output: NULL pointer returned by fgets. cmd= " << cmd << " host= " << hostname <<
			std::endl;
		fflush( stdout );
	}
	buf[strlen(buf)-1]='\0';
	string const result( buf );
	pclose( cmdpopen );
	return result;
}

string
get_hostname()
{
	string const result( run_command_and_get_output("hostname") );
	if ( result.empty() ) return "unknown";
	return result;
}



Size
get_next_decoy_number_and_reserve_if_not_done(
	string const worktag,
	Size const nstruct,
	string const & filename
)
{
	using decoy_timer::decoy_counter;


	string const host( get_hostname() );
	decoy_timer::starttime = clock();

	TR_IO << "get_next_decoy_number_and_reserve_if_not_done: "<< worktag << ' ' << nstruct << ' ' <<
		filename << endl;

	if ( filename == decoy_timer::last_filename && decoy_counter.count( worktag ) && decoy_counter[ worktag ] >= nstruct ){
		TR_IO << "get_next_decoy_number_and_reserve_if_not_done:: use old data: " << worktag << ' ' <<
			decoy_counter[worktag] << ' ' << nstruct << endl;
		return nstruct+1;
	}

	decoy_timer::last_filename = filename;
	decoy_counter.clear();

	safely_lock_file( filename ); // LOCK the file

	std::ifstream data( filename.c_str() );
	string line, tag;
	Size last_decoy_made( 0 ), n;
	while ( getline( data, line ) ) {
		std::istringstream l( line );
		l >> tag >> n;
		map<string,Size>::iterator it( decoy_counter.find( tag ) );
		if ( it == decoy_counter.end() ) {
			decoy_counter[tag];
			it = decoy_counter.find( tag );
		}
		// store in our map
		it->second = std::max( it->second, n );
		TR_IO << "reading " << filename << " tag: " << tag << " n: " << n << endl;
		if ( l.fail() ) {
			std::cout << "get_next_decoy_number_and_reserve_if_not_done: parse error: " << filename << ' ' << line << endl;
			std::cerr << "get_next_decoy_number_and_reserve_if_not_done: parse error: " << filename << ' ' << line << endl;
			continue;
		}
		if ( tag == worktag ) last_decoy_made = std::max( last_decoy_made, n );
	}
	data.close();

	Size decoy_number( last_decoy_made + 1 );
	if ( decoy_number <= nstruct ) {
		/// write my decoy_number to the file, thereby reserving that workunit
		std::ofstream out( filename.c_str(), std::ios_base::out | std::ios_base::app );

		out << worktag << ' ' << decoy_number << ' ' << host << '\n';
		out.close();
	}


	safely_unlock_file( filename ); // UNLOCK the file

	cout << "START " << filename << ' ' << worktag << ' ' << decoy_number << ' ' << host << std::endl;
	cerr << "START " << filename << ' ' << worktag << ' ' << decoy_number << ' ' << host << std::endl;

	fflush( stdout );

	return decoy_number;
}

void
signal_that_work_is_done_and_check_simtime(
	string const worktag,
	string const & simfilename
)
{
	clock_t stoptime = clock();
	Real const simtime( ( (double) stoptime - decoy_timer::starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes

	static string const hostname( get_hostname() );

	std::string const filename( simfilename+".completed" );

	safely_lock_file( filename ); // LOCK the file

	//// how many work units have been done for this worktag?
	std::ifstream data( filename.c_str() );
	string line, tag;
	Size count( 0 ), n;
	while ( getline( data, line ) ) {
		std::istringstream l( line );
		l >> tag >> n;
		if ( l.fail() ) {
			std::cout << "signal_that_work_is_done_and_check_simtime: parse error: " << filename << ' ' << line << endl;
			std::cerr << "signal_that_work_is_done_and_check_simtime: parse error: " << filename << ' ' << line << endl;
			continue;
		}
		if ( tag == worktag ) ++count;
	}
	data.close();

	std::ofstream out( filename.c_str(), std::ios_base::out | std::ios_base::app );
	out << worktag << ' ' << count+1 << ' ' << hostname << ' ' << F(9,3,simtime) << '\n';
	out.close();
	safely_unlock_file( filename ); // UNLOCK the file

	basic::prof_show_oneliner(); // NEW

	check_simtime();
}


bool // passed score filter?
append_score_to_scorefile_and_filter(
	string const worktag,
	Real const score,
	Real const score_filter_acceptance_rate,
	bool const pass_if_threshold_index_is_below_1,
	string const & workfilename
)
{
	string const filename( workfilename + ".scores" );
	string const host( get_hostname() );

	TR_IO << "append_score_to_scorefile_and_filter: "<< worktag << ' ' << score << ' ' <<
		filename << endl;

	safely_lock_file( filename ); // LOCK the file


	std::ifstream data( filename.c_str() );
	string line, tag;
	Reals all_scores;
	while ( getline( data, line ) ) {
		Real other_score;
		std::istringstream l( line );
		l >> tag >> other_score;
		TR_IO << "reading " << filename << " tag: " << tag << " other_score: " << other_score << endl;
		if ( l.fail() ) {
			std::cout << "append_score_to_scorefile_and_filter: parse error: " << filename << ' ' << line << endl;
			std::cerr << "append_score_to_scorefile_and_filter: parse error: " << filename << ' ' << line << endl;
			continue;
		}
		if ( tag == worktag ) all_scores.push_back( other_score );
	}
	data.close();

	/// write my score to the file
	std::ofstream out( filename.c_str(), std::ios_base::out | std::ios_base::app );
	out << worktag << ' ' << ObjexxFCL::format::F(20,6,score) << ' ' << host << '\n';
	out.close();

	safely_unlock_file( filename ); // UNLOCK the file

	Real const simtime( ( (double) clock() - decoy_timer::starttime )/( CLOCKS_PER_SEC*60 ) ); // in minutes
	cout << "FINISH " << filename << ' ' << worktag << ' ' << host << ' ' << simtime << std::endl;
	cerr << "FINISH " << filename << ' ' << worktag << ' ' << host << ' ' << simtime << std::endl;

	fflush( stdout );

	/// always accept if desired:
	if ( score_filter_acceptance_rate >= 0.999 ) return true;


	/// figure out whether this decoy passed
	std::sort( all_scores.begin(), all_scores.end() );

	Size threshold_index( std::min( all_scores.size(), Size( score_filter_acceptance_rate * all_scores.size() ) ) );

	if ( threshold_index < 1 ) {
		if ( !pass_if_threshold_index_is_below_1 ) {
			TR_IO << "threshold_index<1 ==> failed score_filter" << endl;
			return false; ///// NOTE EARLY RETURN
		}
		if ( all_scores.empty() ) return true; //// NOTE EARLY RETURN
		threshold_index = 1;
	}

	Real const epsilon( 1e-6 );
	bool const passed_score_filter( score < all_scores[ threshold_index ] - epsilon ); // new 11/02/12
	//  bool const passed_score_filter( score <= all_scores[ threshold_index ] );

	if ( all_scores.empty() ) {
		TR_IO << "all_scores.empty()!" << endl;
	} else {
		TR_IO << "passed_score_filter: " << passed_score_filter << " threshold_index: " << threshold_index <<
			" all_scores.size(): " << all_scores.size() << " score: " << score <<
			" all_scores.front(): " << all_scores.front() << " all_scores.back(): " << all_scores.back() << std::endl;
	}

	return passed_score_filter;
}

void
set_ss_from_dssp( string const & filename, Pose & pose )
{
	string const dssp_file( filename.substr( 0, filename.size()-3 ) + "dssp" );
	// new logic
	if ( !utility::file::file_exists( dssp_file ) ) {
		string const dssp_exe( "/home/pbradley/dssp" );
		run_command( dssp_exe + " " + filename + " > " + dssp_file + " 2> " + dssp_file + ".err" );
	}

	devel::blab::loops::set_secstruct_from_dssp( pose, dssp_file );
}


void
dump_hbonds(
	utility::vector1< bool > const & subset,
	pose::Pose const & pose,
	scoring::hbonds::HBondSet const & hbond_set,
	std::ostream & os
)
{
	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		scoring::hbonds::HBond const & hb( hbond_set.hbond(i) );
		if ( subset[ hb.don_res() ] || subset[ hb.acc_res() ] ) {
			os << "REMARK HBOND " <<
				I(5,hb.don_res()) << ' ' << pose.residue( hb.don_res() ).atom_name( hb.don_hatm() ) << ' ' <<
				I(5,hb.acc_res()) << ' ' << pose.residue( hb.acc_res() ).atom_name( hb.acc_atm() ) << ' ' <<
				hbond_set.allow_hbond(i) << F(9,3,hb.energy() ) << F(9,3,hb.weight() ) <<
				F(9,3,hb.energy() * hb.weight() ) << std::endl;
		}
	}
}

void
dump_hbonds(
	pose::Pose const & pose,
	scoring::hbonds::HBondSet const & hbond_set,
	std::ostream & os
)
{
	utility::vector1< bool > subset( pose.total_residue(), true );
	dump_hbonds( subset, pose, hbond_set, os );
}




void
dump_hbond_info_to_stream(
	bools const & subset,
	Pose & pose,
	scoring::ScoreFunction const & scorefxn,
	std::ostream & out
)
{
	scorefxn( pose ); // ensure nbrs up-to-date

	scoring::methods::EnergyMethodOptions const & energy_method_options( scorefxn.energy_method_options() );

	/// create and fill an HBondSet
	scoring::hbonds::HBondSet hbond_set( energy_method_options.hbond_options() );
	//  hbond_set.exclude_DNA_DNA  ( energy_method_options.exclude_DNA_DNA() );
	//  hbond_set.use_hb_env_dep   ( energy_method_options.use_hb_env_dep() );
	//  hbond_set.smooth_hb_env_dep( energy_method_options.smooth_hb_env_dep() );
	//hbond_set.setup_for_residue_pair_energies( pose );
	scoring::hbonds::fill_hbond_set( pose, false, hbond_set );

	dump_hbonds( subset, pose, hbond_set, out );
	//  io::pdb::dump_hbonds( subset, pose, hbond_set, out );
}

void
append_hbond_info_to_pdb_file(
	Pose & pose,
	scoring::ScoreFunction const & scorefxn,
	string const & filename
)
{
	bools const subset( pose.total_residue(), true );
	std::ofstream out( filename.c_str(), std::ios_base::out | std::ios_base::app );
	dump_hbond_info_to_stream( subset, pose, scorefxn, out );
	out.close();
}

string // return the pdb filename
create_output(
	Pose const & pose,
	bool const passed_score_filter,
	string const & outfilename, // just basename, no slashes
	string const & scoreline,
	string const & simfile
)
{
	runtime_assert( outfilename.find("/") == string::npos );

	string const scoresfile( simfile + ".out" ), pdboutdir( simfile + ".pdbs/" );

	safely_append_lines_to_file( make_vector1( scoreline ), scoresfile );

	if ( passed_score_filter ) {
		if ( !utility::file::file_exists( pdboutdir ) ) utility::file::create_directory( pdboutdir );
		// add the scoreline
		string const real_outfilename( pdboutdir + outfilename );
		std::ofstream out( real_outfilename.c_str() );
		out << "REMARK SCORELINE " << scoreline;
		if ( devel::blab::motif::has_motif_data( pose ) ) {
			out << "REMARK " << devel::blab::motif::get_motif_data( pose ) << '\n';
		}
		out << "REMARK " << pose.fold_tree() << '\n';
		pose.dump_pdb( out );
		out.close();
		return real_outfilename;
	}
	return string();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
fill_bools_from_resfile(
	Pose const & pose,
	string const filename,
	bools & is_flexible
)
{
	is_flexible.clear();
	is_flexible.resize( pose.total_residue(), false );

	ifstream data( filename.c_str() );
	string line;
	while ( getline( data, line ) ) {
		strings const l( split_to_vector1( line ) );
		if ( l.size() == 4 ) {
			int pdbpos(0);
			char insertcode(' ');
			if ( is_int( l[1] ) ) {
				pdbpos = int_of(l[1] );
			} else {
				Size const len( l[1].size() );
				insertcode = l[1][ len-1 ];
				runtime_assert( is_int( l[1].substr(0,len-1) ) );
				pdbpos = int_of(l[1].substr(0,len-1) );
				TR_IO.Trace << "fill_bools_from_resfile: reading insertcode " << l[1] << endl;
			}
			runtime_assert( l[2].size() == 1 );
			char const chain( l[2][0] );
			Size const posepos( pose.pdb_info()->pdb2pose( chain, pdbpos, insertcode ) );
			TR_IO.Trace << filename << ": " << line << endl;
			TR_IO.Trace << posepos << ' ' << pose.total_residue() << endl;
			runtime_assert( posepos && posepos <= pose.total_residue() );
			is_flexible[ posepos ] = true;
		} else {
			TR_IO.Trace << "resfile line " << line << endl;
		}
	}
	data.close();

}


// void
// cenpose_from_pdb( Pose & pose, string const filename )
// {
//  ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );
//  pose_from_pdb( pose, *rsd_set, filename );
// }


#endif
