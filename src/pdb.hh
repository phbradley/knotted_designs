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
#ifndef INCLUDED_apps_pilot_phil_pdb_HH
#define INCLUDED_apps_pilot_phil_pdb_HH

#include <apps/pilot/phil/phil.hh>
#include <core/pose/util.hh>
#include <devel/blab/motif/MotifData.hh>
#include <devel/dna/util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
// #include <basic/Tracer.hh>
#include <fstream>


static basic::Tracer TR_PDB( "apps.pilot.phil.pdb_hh" );


///////////////////////////////////////////////////////////////////////////////////////
Real
get_resolution_from_pdb_file( string const & filename )
{
	string const restag( "REMARK   2 RES" );
	utility::io::izstream data( filename );
	string line, tag;
	while ( getline( data, line ) ) {
		if ( line.size() > restag.size() && line.substr(0,restag.size()) == restag ) {
			istringstream l(line );
			Real resolution;
			l >> tag >> tag >> tag >> resolution;
			if ( !l.fail() ) return resolution;
		}
	}
	return 9.9;
}

///////////////////////////////////////////////////////////////////////////////////////
//

void
read_backbone_torsions_from_file(
	string const & filename,
	string & sequence,
	vector1< Reals > & bb_torsions
)
{
	ifstream data( filename.c_str());

	string line;
	sequence.clear();
	Vectors bbcoords;
	while ( getline( data, line ) ) {
		if ( line.substr(0,6) == "ATOM  " ) {
			string const atomname( line.substr(12,4) );
			if ( atomname == " N  " ||
					atomname == " CA " ||
					atomname == " C  " ) {
				bbcoords.push_back( Vector( float_of( line.substr( 30,8 ) ),
					float_of( line.substr( 38,8 ) ),
					float_of( line.substr( 46,8 ) ) ) );
				if ( atomname == " CA " ) {
					string const resname( line.substr(17,3) );
					sequence.push_back( oneletter_code_from_aa( aa_from_name( resname ) ));
				}
			}
		}
	}
	data.close();

	if ( bbcoords.size() != 3 * sequence.size() ) {
		cout << "messed up bbcoords: " << bbcoords.size() << ' ' << sequence.size() << ' ' << filename << endl;
		sequence.clear();
		bb_torsions.clear();
		return;
	}

	bb_torsions.clear();
	for ( Size i=1; i<= sequence.size(); ++i ) {
		bb_torsions.push_back( Reals( 3, 0.0 ) );
	}

	for ( Size i=1; i<= bbcoords.size()-3; ++i ) {
		Real const torsion( numeric::dihedral_degrees( bbcoords[i], bbcoords[i+1], bbcoords[i+2], bbcoords[i+3] ) );
		Size const tor_index( i%3+1 ), seq_index( i/3+1 );
		bb_torsions[ seq_index ][ tor_index ] = torsion;
	}
}


///////////////////////////////////////////////////////////////////////////////
string tmpstring; // hack

vector1< Vectors >
read_CA_coords_from_complex_file(
	string const & filename,
	string & sequence = tmpstring
)
{
	vector1< Vectors > coords( 2 );

	utility::io::izstream data( filename ); //.c_str());
	// ifstream data( filename.c_str());

	string line;
	sequence.clear();
	while ( getline( data, line ) ) {
		if ( line.substr(0,6) == "ATOM  " && line.substr(12,4) == " CA " ) {
			char const chain( line[21] );
			runtime_assert( chain == 'A' || chain == 'B' );
			Size const chain_index( chain == 'A' ? 1 : 2 );
			coords[ chain_index ].push_back( Vector( float_of( line.substr( 30,8 ) ),
				float_of( line.substr( 38,8 ) ),
				float_of( line.substr( 46,8 ) ) ) );
			string const resname( line.substr(17,3) );
			sequence.push_back( oneletter_code_from_aa( aa_from_name( resname ) ));

		}
	}
	data.close();

	return coords;
}


vector1< Vectors > // RVO ?
read_CA_coords_multichain(
	string const & filename
)
{
	vector1< Vectors > coords;

	utility::io::izstream data( filename ); //.c_str());
	// ifstream data( filename.c_str());

	string line;

	char prev_chain('?');

	while ( getline( data, line ) ) {
		if ( line.substr(0,6) == "ATOM  " && line.substr(12,4) == " CA " ) {
			char const chain( line[21] );
			if ( chain != prev_chain ) {
				TR_PDB.Trace << "read_CA_coords_multichain: new_chain old: " << prev_chain << " new: " << chain << endl;
				prev_chain = chain;
				coords.push_back( Vectors() );
			}
			coords.back().push_back( Vector( float_of( line.substr( 30,8 ) ),
				float_of( line.substr( 38,8 ) ),
				float_of( line.substr( 46,8 ) ) ) );
		}
	}
	data.close();

	return coords;
}


/// this one doesnt care about chains
Vectors
read_CA_coords_from_file(
	string const & filename
)
{
	Vectors coords;
	utility::io::izstream data( filename ); //.c_str());
	// ifstream data( filename.c_str());
	runtime_assert( data.good() );
	string line;
	while ( getline( data, line ) ) {
		if ( line.substr(0,6) == "ATOM  " && line.substr(12,4) == " CA " ) {
			//char const chain( line[21] );
			coords.push_back( Vector( float_of( line.substr( 30,8 ) ),
				float_of( line.substr( 38,8 ) ),
				float_of( line.substr( 46,8 ) ) ) );
		}
	}
	data.close();
	return coords;
}


///////////////////////////////////////////////////////////////////////////////////////
Vectors
read_CA_coords_from_file(
												 Size const startpos,
												 Size const  stoppos,
												 string const & filename
												 )
{
	Vectors coords;

	utility::io::izstream data( filename );

	string line;
	while ( getline( data, line ) ) {
		if ( line.substr(0,6) == "ATOM  " && line.substr(12,4) == " CA " ) {
			Size const pos( int_of( line.substr( 22,4 ) ) );
			if ( pos < startpos )continue;
			else if ( pos > stoppos ) break;
			coords.push_back( Vector( float_of( line.substr( 30,8 ) ),
																float_of( line.substr( 38,8 ) ),
																float_of( line.substr( 46,8 ) ) ) );

		}
	}
	data.close();

	return coords;
}


///////////////////////////////////////////////////////////////////////////////
/// this one doesnt care about chains
/// returns an empty vector1 if there's a problem...
///
Vectors
read_bb_coords_from_file(
	string const & filename
)
{
	Vectors n_coords, ca_coords, c_coords, coords;
	utility::io::izstream data( filename );
	if ( !data.good() ) return coords; // empty

	string line;
	while ( getline( data, line ) ) {
		if ( line.substr(0,6) == "ATOM  " ) {
			string const atomname( line.substr(12,4) );
			if ( atomname == " N  " ) {
				//char const chain( line[21] );
				n_coords.push_back( Vector( float_of( line.substr( 30,8 ) ),
					float_of( line.substr( 38,8 ) ),
					float_of( line.substr( 46,8 ) ) ) );
			} else if ( atomname == " CA " ) {
				//char const chain( line[21] );
				ca_coords.push_back( Vector( float_of( line.substr( 30,8 ) ),
					float_of( line.substr( 38,8 ) ),
					float_of( line.substr( 46,8 ) ) ) );
			} else if ( atomname == " C  " ) {
				//char const chain( line[21] );
				c_coords.push_back( Vector( float_of( line.substr( 30,8 ) ),
					float_of( line.substr( 38,8 ) ),
					float_of( line.substr( 46,8 ) ) ) );
			}
		}
	}
	data.close();

	if ( n_coords.size() == ca_coords.size() && n_coords.size() == c_coords.size() ) {
		for ( Size i=1; i<= n_coords.size(); ++i ) {
			coords.push_back(  n_coords[i] );
			coords.push_back( ca_coords[i] );
			coords.push_back(  c_coords[i] );
		}
	}
	return coords;
}


///////////////////////////////////////////////////////////////////////////////
void
count_atom_occupancy(
	string const & filename,
	Pose const & pose,
	id::AtomID_Map< Size > & count
)
{
	core::pose::initialize_atomid_map( count, pose, Size(0) );
	utility::io::izstream data( filename );

	string line;
	while ( getline( data, line ) ) {
		if ( line.substr(0,6) == "ATOM  " ) {
			string const atomname( line.substr(12,4) );
			char const chain( line[21] );
			Size const resnum( int_of( line.substr(22,4 ) ) );
			char const insertcode( line[26] );
			Size const pos( pose.pdb_info()->pdb2pose( chain, resnum, insertcode  ) );
			if ( !pos ) {
				TR_PDB.Trace << "no pose pos: "<< chain << ' ' << resnum << ' '<< insertcode << endl;
				continue;
			}
			Residue const & rsd( pose.residue( pos  ) );
			if ( !rsd.has( atomname ) ) {
				TR_PDB.Trace << "rsd missing atom " << rsd.name() << ' ' << atomname << endl;
				continue;
			}
			++count[ id::AtomID( rsd.atom_index( atomname ), pos ) ];
		}
	}
	data.close();
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
bools
identify_residues_with_all_single_occupancy(
	Pose const & pose,
	string const filename
)
{
	bools goodrsd( pose.total_residue(), false );

	id::AtomID_Map< Size > atom_occupancy;
	count_atom_occupancy( filename, pose, atom_occupancy );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		bool good( true );
		for ( Size ii=1; ii<= rsd.nheavyatoms(); ++ii ) {
			if ( rsd.atom_type(ii).is_virtual() ) continue;
			if ( atom_occupancy[ id::AtomID( ii, i) ] != 1 ) good = false;
		}
		goodrsd[i] = good;
		TR_PDB.Trace << "goodrsd: " << I(4,i) << ' ' << rsd.name1() << ' ' << good << ' ' << filename << endl;
	}

	return goodrsd;
}


///////////////////////////////////////////////////////////////////////////////
void
transform_pdbfile_coords_and_append_to_stream(
	string const & filename,
	Matrix const & R,
	Vector const & v,
	bools const & output_mask,
	ostream & out
)
{
	utility::io::izstream data( filename );//.c_str() );
	// ifstream data( filename.c_str() );

	string line;
	while ( getline( data, line  ) ) {
		if ( line.substr( 0,6 ) == "ATOM  " || line.substr( 0,6 ) == "HETATM" ) {
			int const resnum( int_of( line.substr(22,4) ) );
			if ( resnum < 1 || resnum > int(output_mask.size()) || ( !output_mask[resnum] ) ) continue;
			Vector xyz( float_of( line.substr(30,8) ), float_of( line.substr( 38,8 ) ), float_of( line.substr( 46,8 ) ) );
			xyz = R * xyz + v;
			out << line.substr(0,30) << F(8,3,xyz.x() ) << F(8,3,xyz.y()) << F(8,3,xyz.z()) << line.substr( 54 ) << '\n';
		} else {
			out << line << '\n';
		}
	}
	data.close();

}

///////////////////////////////////////////////////////////////////////////////////////

void
transform_pdbfile_coords_and_append_to_stream(
	string const & filename,
	Matrix const & R,
	Vector const & v,
	Size const min_resnum,
	Size const max_resnum,
	ostream & out
)
{
	bools output_mask( max_resnum, false );
	for ( Size i=min_resnum; i<= max_resnum; ++i ) output_mask[i] = true;
	transform_pdbfile_coords_and_append_to_stream( filename, R, v, output_mask, out );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
read_waters_from_pdb_file(
	std::string const & filename,
	strings & names,
	Vectors & waters // appends to this vector, does not clear at the start...
)
{
	std::ifstream data( filename.c_str());
	//Vectors waters;
	string line;
	while ( getline( data,line) ) {
		if ( ( line.substr( 0,6) == "HETATM" || line.substr(0,6) == "ATOM  " ) &&
				( stripped( line.substr(12,4) ) == "O" ) &&
				( line.substr(17,3) == "HOH" || line.substr(17,3) == "WAT" || line.substr(17,3) == "WOT" ) ) {
			char const altloc( line[16] );
			string const ok_altlocs(" aA1" );
			if ( ok_altlocs.find( altloc ) == std::string::npos ) {
				TR_PDB.Warning << "SKIPPING water altloc: " << altloc << ' ' << line.substr(16,10) << endl;
				continue;
			}

			Vector const xyz( float_of( line.substr(30,8) ), float_of( line.substr(38,8) ), float_of( line.substr(46,8) ) );
			waters.push_back( xyz );
			names.push_back( line.substr(17,10) );
		}
	}
	data.close();
	TR_PDB.Trace << "Read " << waters.size() << " xtal waters from the file: " << filename << endl;
	if ( utility::file::file_exists( filename+".waters" ) ) {
		// RECURSION
		read_waters_from_pdb_file( filename+".waters", names, waters );
	}
}

void
dump_motif_pdb_w_base_partner( Pose const & pose, string const & filename )
{
	using namespace devel::blab::motif;

	utility::io::ozstream out( filename );

	TR_PDB.Trace << "dump_motif_pdb: nres= " << pose.total_residue() << ' ' << filename << std::endl;
	// mechanism to stash the BasePartner info in the MotifData
	if ( has_motif_data( pose ) ) {
		MotifData md( get_motif_data( pose ) );
		if ( scoring::dna::has_base_partner( pose ) ) {
			md.segment("BASE_PARTNER").clear();
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				md.segment("BASE_PARTNER").push_back( scoring::dna::retrieve_base_partner_from_pose( pose )[i]);
			}
		}
		out << "REMARK " << md << '\n';
	}

	pose.dump_pdb( out );
	//io::pdb::dump_complete_pdb( pose, out );
	out.close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
write_ca_pdbfile(
	Vectors const & coords,
	string const & outfilename
)
{
	ofstream out( outfilename.c_str() );

	for ( Size i=1; i<= coords.size(); ++i ) {
		char const chain( 'A' );


		out << "ATOM  " << I(5,i) << ' ' << " CA " << ' ' <<
			"ALA" << ' ' << chain << I(4,i ) << "    " <<
			F(8,3,coords[i].x()) <<
			F(8,3,coords[i].y()) <<
			F(8,3,coords[i].z()) <<
			F(6,2,1.0) << F(6,2,1.0) << '\n';
	}
	out.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
parse_defs_and_make_mutations(
	strings const & defs,
	Pose & pose,
	bools & subset
)
{
	subset.clear(); subset.resize( pose.total_residue(), false );
	for ( Size ii=1; ii<= defs.size(); ++ii ) {
		strings const l( split_to_vector1( defs[ii],"." ) );
		if ( l.size() == 2 ) {
			// chain.pos
			string const chain( l[1] );
			if ( chain.size() != 1 || !is_int( l[2] ) ) utility_exit_with_message("bad def: "+defs[ii] );
			Size const seqpos( pose.pdb_info()->pdb2pose( chain[0], int_of( l[2] ) ) );
			if ( !seqpos || seqpos > pose.total_residue() ) utility_exit_with_message("bad def: "+defs[ii] );
			subset[seqpos] = true;
		} else if ( l.size() == 3 || l.size() == 4 ) {
			runtime_assert(l[3].size() == 1 );
			// chain.pos.newname1 OR chain.pos.oldname1.newname1
			string const chain( l[1] );
			if ( chain.size() != 1 || !is_int( l[2] ) ) utility_exit_with_message("bad def: "+defs[ii] );
			Size const seqpos( pose.pdb_info()->pdb2pose( chain[0], int_of( l[2] ) ) );
			if ( !seqpos || seqpos > pose.total_residue() ) utility_exit_with_message("bad def: "+defs[ii] );
			subset[seqpos] = true;
			// possibly mutate to aa
			char const oldaa( pose.residue(seqpos).name1()), newaa( l[l.size()][0] );
			if ( l.size() == 4 ) { // sequence checking of oldaa
				runtime_assert(l[4].size() == 1 );
				runtime_assert( oldaa == l[3][0] );
			}
			AA aa;
			bool const use_methyl_cytosine        ( newaa == 'd' );
			bool const use_methyl_cytosine_partner( newaa == 'h' );
			if ( use_methyl_cytosine ) {
				aa = na_cyt;
			} else if ( use_methyl_cytosine_partner ) {
				aa = na_gua;
			} else if ( 'a'<= newaa && newaa<='z' ) { // lowercase ==> DNA
				aa = devel::dna::dna_aa_from_oneletter_code( newaa );
			} else {
				aa = aa_from_oneletter_code( newaa );
			}
			Residue const & rsd( pose.residue(seqpos) );
			if ( rsd.aa() != aa ) {
				TR_PDB.Trace << "Making mutation: " << defs[ii] << " from " << rsd.aa() << " to " << aa << ' ' << newaa << ' ' <<
					use_methyl_cytosine << ' ' << use_methyl_cytosine_partner << endl;
				if ( rsd.is_protein() ) {
					make_sequence_change( seqpos, aa, pose );
				} else if ( rsd.is_DNA() ) {
					devel::dna::make_base_pair_mutation( pose, seqpos, aa );
					if ( use_methyl_cytosine ) {
						add_variant_type_to_pose_residue( pose, chemical::METHYLATION, seqpos );
					} else if ( use_methyl_cytosine_partner ) {
						Size const seqposp( retrieve_base_partner_from_pose( pose )[ seqpos] );
						runtime_assert( seqposp );
						add_variant_type_to_pose_residue( pose, chemical::METHYLATION, seqposp );
					}
					// { // test
					// 	cout << "HACKING1: " << seqpos << ' ' << pose.residue(seqpos).name() << endl;
					// 	add_variant_type_to_pose_residue( pose, chemical::METHYLATION, seqpos );
					// 	cout << "HACKING2: " << seqpos << ' ' << pose.residue(seqpos).name() << endl;
					// }
				} else {
					utility_exit_with_message("unknown mutrsd type "+rsd.name() );
				}
			}
		}
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
parse_defs_and_make_mutations_design(
	strings const & defs,
	Pose & pose,
	bools & subset,
	bools & is_mutable // newaa == 'X' or 'x'; no muts performed here
)
{
	subset.clear(); subset.resize( pose.total_residue(), false );
	is_mutable = subset;
	for ( Size ii=1; ii<= defs.size(); ++ii ) {
		strings const l( split_to_vector1( defs[ii],"." ) );
		if ( l.size() == 2 ) {
			// chain.pos
			string const chain( l[1] );
			if ( chain.size() != 1 || !is_int( l[2] ) ) utility_exit_with_message("bad def: "+defs[ii] );
			Size const seqpos( pose.pdb_info()->pdb2pose( chain[0], int_of( l[2] ) ) );
			if ( !seqpos || seqpos > pose.total_residue() ) utility_exit_with_message("bad def: "+defs[ii] );
			subset[seqpos] = true;
		} else if ( l.size() == 3 || l.size() == 4 ) {
			runtime_assert(l[3].size() == 1 );
			// chain.pos.newname1 OR chain.pos.oldname1.newname1
			string const chain( l[1] );
			if ( chain.size() != 1 || !is_int( l[2] ) ) utility_exit_with_message("bad def: "+defs[ii] );
			Size const seqpos( pose.pdb_info()->pdb2pose( chain[0], int_of( l[2] ) ) );
			if ( !seqpos || seqpos > pose.total_residue() ) utility_exit_with_message("bad def: "+defs[ii] );
			subset[seqpos] = true;
			// possibly mutate to aa
			char const oldaa( pose.residue(seqpos).name1()), newaa( l[l.size()][0] );
			if ( l.size() == 4 ) { // sequence checking of oldaa
				runtime_assert(l[4].size() == 1 );
				runtime_assert( oldaa == l[3][0] );
			}
			if ( newaa == 'X' or newaa == 'x' ) { // design position
				is_mutable[seqpos] = true;
				continue;
			}
			AA aa;
			bool const use_methyl_cytosine        ( newaa == 'd' );
			bool const use_methyl_cytosine_partner( newaa == 'h' );
			if ( use_methyl_cytosine ) {
				aa = na_cyt;
			} else if ( use_methyl_cytosine_partner ) {
				aa = na_gua;
			} else if ( 'a'<= newaa && newaa<='z' ) { // lowercase ==> DNA
				aa = devel::dna::dna_aa_from_oneletter_code( newaa );
			} else {
				aa = aa_from_oneletter_code( newaa );
			}
			Residue const & rsd( pose.residue(seqpos) );
			if ( rsd.aa() != aa || use_methyl_cytosine || use_methyl_cytosine_partner ) {
				TR_PDB.Trace << "Making mutation: " << defs[ii] << " from " << rsd.aa() << " to " << aa << ' ' << newaa << ' ' <<
					use_methyl_cytosine << ' ' << use_methyl_cytosine_partner << endl;
				if ( rsd.is_protein() ) {
					make_sequence_change( seqpos, aa, pose );
				} else if ( rsd.is_DNA() ) {
					devel::dna::make_base_pair_mutation( pose, seqpos, aa );
					if ( use_methyl_cytosine ) {
						add_variant_type_to_pose_residue( pose, chemical::METHYLATION, seqpos );
					} else if ( use_methyl_cytosine_partner ) {
						Size const seqposp( retrieve_base_partner_from_pose( pose )[ seqpos] );
						runtime_assert( seqposp );
						add_variant_type_to_pose_residue( pose, chemical::METHYLATION, seqposp );
					}
				} else {
					utility_exit_with_message("unknown mutrsd type "+rsd.name() );
				}
			}
		}
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
