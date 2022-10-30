/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <unordered_set>

#include <cfg.hpp>
#include <cif++.hpp>
#include <gxrio.hpp>

#include <pdb-redo/SkipList.hpp>
#include <pdb-redo/Symmetry-2.hpp>

#include "revision.hpp"

namespace fs = std::filesystem;

// -----------------------------------------------------------------------

const std::set<std::string> kBackBone = {
	"N", "CA", "C", "O", "OXT"
};

const float
	kMaxZnHisDistanceInCluster = 3.8f,
	kMaxZnCysDistanceInCluster = 4.8f,
	// kBoundaryClosenessAtomToZn = 2.9f,

	kDixonQTest95Perc5Points = 0.710f,

	kMaxMetalLigandDistance = 3.5f;

// -----------------------------------------------------------------------

// enum class RestraintType { Dist, Angle }; //TODO implement: TORSION (and others?) when necessary

// std::string to_string(RestraintType type)
// {
// 	switch (type)
// 	{
// 		case RestraintType::Angle:	return "angle";
// 		case RestraintType::Dist:	return "dist";
// 	}
// }

class RestraintGenerator
{
  public:
	RestraintGenerator(const std::string &file, bool deleteVDWRestraints)
		: m_file(file)
		, m_unique(m_file)
		, m_deleteVDWRestraints(deleteVDWRestraints)
	{
		if (not m_file.is_open())
			throw std::runtime_error("Could not open restraint file " + file);
	}

	~RestraintGenerator() = default;

	std::size_t writeTetrahedral(const cif::mm::atom &ion, const std::vector<cif::mm::atom> &ligands);
	std::size_t writeOctahedral(const cif::mm::atom &ion, const std::vector<cif::mm::atom> &ligands,
		const std::vector<std::tuple<std::size_t, std::size_t>> &opposing);

  private:

	class unique_streambuf : public std::streambuf
	{
	  public:
		using base_type = std::streambuf;
		using int_type = base_type::int_type;
		using char_type = base_type::char_type;
		using traits_type = base_type::traits_type;

		unique_streambuf(std::ostream &os)
			: m_os(os)
			, m_upstream(os.rdbuf())
		{
		}

		~unique_streambuf()
		{
			m_os.rdbuf(m_upstream);
		}

		virtual int_type
		overflow(int_type ic = traits_type::eof())
		{
			char ch = traits_type::to_char_type(ic);

			int_type result = ic;

			if (ch == '\n' or result == traits_type::eof())
			{
				if (m_seen.insert(m_line).second)
				{
					for (char c : m_line)
						m_upstream->sputc(c);
					if (ch == '\n')
						m_upstream->sputc(ch);
				}

				m_line.clear();
			}
			else
				m_line += ch;

			return result;
		}

		std::streambuf *get_upstream() const { return m_upstream; }

	  private:
		std::ostream &m_os;
		std::streambuf *m_upstream;
		std::string m_line;
		std::unordered_set<std::string> m_seen;
	};

	void writeAngleRestraint(float target, float sd, const cif::mm::atom &a, const cif::mm::atom &b, const cif::mm::atom &c);

	struct AtomPart
	{
		const cif::mm::atom &m_a;

		AtomPart(const cif::mm::atom &a)
			: m_a(a)
		{
		}

		friend std::ostream &operator<<(std::ostream &os, const AtomPart &aw)
		{
			os << "chain " << aw.m_a.get_auth_asym_id()
			   << " resi " << aw.m_a.get_auth_seq_id()
			   << " ins " << (aw.m_a.get_pdb_ins_code().empty() ? "." : aw.m_a.get_pdb_ins_code())
			   << " atom " << aw.m_a.get_auth_atom_id();
			if (not aw.m_a.get_auth_alt_id().empty())
				os << " alt " + aw.m_a.get_auth_alt_id();
			else if (not aw.m_a.get_label_alt_id().empty())
				os << " alt " + aw.m_a.get_label_alt_id();
			os << " symm " << (aw.m_a.is_symmetry_copy() ? "Y" : "N");

			return os;
		}
	};

	std::ofstream m_file;
	unique_streambuf m_unique;
	bool m_deleteVDWRestraints;
};

std::size_t RestraintGenerator::writeTetrahedral(const cif::mm::atom &ion, const std::vector<cif::mm::atom> &ligands)
{
	const float
		kTarget = 109.5,
		kSD = 3;

	std::size_t n = 0;
	for (auto a = ligands.begin(); std::next(a) != ligands.end(); ++a)
	{
		for (auto b = a + 1; b != ligands.end(); ++b)
		{
			writeAngleRestraint(kTarget, kSD, *a, ion, *b);
			++n;
		}
	}

	return n;
}

std::size_t RestraintGenerator::writeOctahedral(const cif::mm::atom &ion, const std::vector<cif::mm::atom> &ligands,
	const std::vector<std::tuple<std::size_t, std::size_t>> &opposing)
{
	const float kSD = 3;

	if (m_deleteVDWRestraints)
		m_file << "vdwr exclude " << AtomPart(ion) << std::endl;

	std::set<std::size_t> one_side_ligand;
	for (const auto &[a, b] : opposing)
		one_side_ligand.insert(a);

	std::size_t n = 0, ai = 0;
	for (auto a = ligands.begin(); std::next(a) != ligands.end(); ++a, ++ai)
	{
		std::size_t bi = ai + 1;
		for (auto b = a + 1; b != ligands.end(); ++b, ++bi)
		{
			if (find(opposing.begin(), opposing.end(), std::make_tuple(ai, bi)) != opposing.end())
			{
				writeAngleRestraint(180, kSD, *a, ion, *b);

				for (auto &c : ligands)
				{
					if (c == *a or c == *b)
						continue;

					writeAngleRestraint(90, kSD, *a, c, *b);
				}
			}
			else if (one_side_ligand.count(ai) and one_side_ligand.count(bi))
				writeAngleRestraint(90, kSD, *a, ion, *b);

			++n;
		}
	}

	return n;
}

void RestraintGenerator::writeAngleRestraint(float target, float sd, const cif::mm::atom &a, const cif::mm::atom &b, const cif::mm::atom &c)
{
	m_file << "exte angle first " << AtomPart(a) << " next " << AtomPart(b) << " next " << AtomPart(c)
		   << std::fixed << std::setprecision(3) << " value " << target << " sigma " << sd << std::endl;
}

// ------------------------------------------------------------------------

struct IonSite
{
	cif::mm::atom ion;
	std::vector<std::tuple<cif::mm::atom, float, std::string>> lig;
	std::vector<std::tuple<std::size_t, std::size_t>> opposing;

	bool isOctaHedral();
};

bool IonSite::isOctaHedral()
{
	const float kMaxAllowedAngleDeviation = 30.0f;

	bool result = lig.size() == 6;

	// check angles
	for (std::size_t a = 0; result and a + 1 < 6; ++a)
	{
		auto &la = std::get<0>(lig[a]);
		std::size_t opposing_la = 0;

		for (std::size_t b = a + 1; result and b < 6; ++b)
		{
			auto &lb = std::get<0>(lig[b]);
			float angle = cif::angle(la.get_location(), ion.get_location(), lb.get_location());

			if (abs(angle - 180) < kMaxAllowedAngleDeviation) // opposing?
			{
				if (opposing_la++ > 0)
					result = false;
				else
					opposing.emplace_back(a, b);
			}
			else if (abs(angle - 90) > kMaxAllowedAngleDeviation) // should be 90 degrees then
				result = false;
		}
	}

	return result and opposing.size() == 3;
}

// -----------------------------------------------------------------------

bool findZincSites(cif::mm::structure &structure, cif::datablock &db, int spacegroup, const clipper::Cell &cell,
	IonSite &zs, const std::string &altID)
{
	bool result = false;

	for (;;)
	{
		// strip out the atoms that are not available in this alt group
		if (not altID.empty())
		{
			zs.lig.erase(
				std::remove_if(zs.lig.begin(), zs.lig.end(), [altID](auto &l)
					{
					auto alt = std::get<0>(l).get_label_alt_id();
					return not alt.empty() and alt != altID; }),
				zs.lig.end());
		}

		// sort ligands on distance
		std::sort(zs.lig.begin(), zs.lig.end(), [](auto &a, auto &b)
			{ return std::get<1>(a) < std::get<1>(b); });

		// take the nearest atom of a residue, this takes care of NE2/ND1 from a single HIS
		// and also the alt cases
		for (auto a = zs.lig.begin(); a != zs.lig.end() and std::next(a) != zs.lig.end(); ++a)
		{
			auto &aa = std::get<0>(*a);

			auto ad = std::get<1>(*a);

			for (auto b = std::next(a); b != zs.lig.end(); ++b)
			{
				auto &ba = std::get<0>(*b);

				if (ba.get_label_comp_id() != aa.get_label_comp_id() or
					aa.get_label_seq_id() != ba.get_label_seq_id() or aa.get_label_asym_id() != ba.get_label_asym_id() or aa.symmetry() != ba.symmetry())
					continue;

				auto bd = std::get<1>(*b);
				assert(bd > ad);

				zs.lig.erase(b);
				break;
			}
		}

		// -----------------------------------------------------------------------

		if (cif::VERBOSE)
		{
			std::cerr << "preliminary cluster: " << std::endl
					  << " zn: " << zs.ion.get_label_asym_id() << '/' << zs.ion.get_label_atom_id() << " (" << zs.ion.pdb_id() << ')' << std::endl;

			for (auto &l : zs.lig)
			{
				auto &a = std::get<0>(l);
				std::cerr << " " << a.get_label_asym_id() << a.get_label_seq_id() << '/' << a.get_label_atom_id() << " (" << a.pdb_id() << ')' << ' ' << a.symmetry() << " @ " << std::get<1>(l) << std::endl;
			}
		}

		// -----------------------------------------------------------------------
		// if there are more than five atoms, give up. If there are exactly five
		// see if we can use a dixon q-test to assign the last to be an outlier.

		if (zs.lig.size() >= 5)
		{
			auto gap = std::get<1>(zs.lig[4]) - std::get<1>(zs.lig[3]);
			auto range = std::get<1>(zs.lig[4]) - std::get<1>(zs.lig[0]);
			if ((gap / range) < kDixonQTest95Perc5Points)
			{
				if (cif::VERBOSE)
					std::cerr << "Rejecting cluster since there are 5 or more atoms near by and nr 5 is not considered to be an outlier" << std::endl;
				break;
			}

			if (cif::VERBOSE)
			{
				for (std::size_t i = 4; i < zs.lig.size(); ++i)
				{
					auto &a = std::get<0>(zs.lig[i]);
					std::cerr << "Atom " << a.get_label_asym_id() << a.get_label_seq_id() << '/' << a.get_label_atom_id() << " (" << a.pdb_id() << ')' << " was considered to be an outlier" << std::endl;
				}
			}

			zs.lig.erase(zs.lig.begin() + 4, zs.lig.end());
		}

		if (zs.lig.size() == 4)
			result = true;
		else if (cif::VERBOSE)
			std::cerr << "Rejecting cluster since there are not 4 atoms near by" << std::endl;

		break;
	}

	return result;
}

std::vector<IonSite> findZincSites(cif::mm::structure &structure, cif::datablock &db, int spacegroup, const clipper::Cell &cell)
{
	std::vector<IonSite> result;

	// factory for symmetry atom iterators
	pdb_redo::SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom : structure.atoms())
	{
		if (atom.get_label_comp_id() != "ZN")
			continue;

		atom.set_property("pdbx_formal_charge", 2);

		IonSite zs = { atom };

		for (auto a : structure.atoms())
		{
			if (a.get_label_comp_id() == "HIS" and (a.get_label_atom_id() == "ND1" or a.get_label_atom_id() == "NE2"))
			{
				for (auto sa : saif(a, [al = atom.get_location()](const cif::point &pt)
						 { return distance(al, pt) <= kMaxZnHisDistanceInCluster; }))
				{
					float d = distance(atom, sa);
					assert(d <= kMaxZnHisDistanceInCluster);
					zs.lig.emplace_back(sa, d, sa.symmetry());
				}
				continue;
			}

			if (a.get_label_comp_id() == "CYS" and a.get_label_atom_id() == "SG")
			{
				for (auto sa : saif(a, [al = atom.get_location()](const cif::point &pt)
						 { return distance(al, pt) <= kMaxZnCysDistanceInCluster; }))
				{
					float d = distance(atom, sa);
					assert(d <= kMaxZnCysDistanceInCluster);
					zs.lig.emplace_back(sa, d, sa.symmetry());
				}
				continue;
			}
		}

		std::set<std::string> altIDs;
		for (const auto &[a, f, s] : zs.lig)
		{
			if (not a.get_label_alt_id().empty())
				altIDs.insert(a.get_label_alt_id());
		}

		if (altIDs.empty()) // no alternates
			altIDs.insert("");

		for (auto &alt : altIDs)
		{
			IonSite azs = zs;
			if (findZincSites(structure, db, spacegroup, cell, azs, alt))
				result.emplace_back(std::move(azs));
		}
	}

	return result;
}

// -----------------------------------------------------------------------

constexpr float get_t_90(std::size_t N)
{
	const float t_dist_90[] = {
		1.638, // 3
		1.533, // 4
		1.476, // 5
		1.440, // 6
		1.415, // 7
		1.397, // 8
		1.383, // 9
		1.372, // 10
	};

	assert(N >= 3 and N < sizeof(t_dist_90) / sizeof(float) + 3);
	return t_dist_90[N - 3];
}

bool findOctahedralSites(cif::mm::structure &structure, cif::datablock &db, int spacegroup, const clipper::Cell &cell,
	IonSite &is, const std::string &altID)
{
	bool result = false;

	for (;;)
	{
		// strip out the atoms that are not available in this alt group
		if (not altID.empty())
		{
			is.lig.erase(
				std::remove_if(is.lig.begin(), is.lig.end(), [altID](auto &l)
					{
					auto alt = std::get<0>(l).get_label_alt_id();
					return not alt.empty() and alt != altID; }),
				is.lig.end());
		}

		// sort ligands on distance
		sort(is.lig.begin(), is.lig.end(), [](auto &a, auto &b)
			{ return std::get<1>(a) < std::get<1>(b); });

		// take the nearest atom of a residue, this takes care of NE2/ND1 from a single HIS
		// and also the alt cases
		for (auto a = is.lig.begin(); a != is.lig.end() and std::next(a) != is.lig.end(); ++a)
		{
			auto &aa = std::get<0>(*a);

			for (auto b = std::next(a); b != is.lig.end(); ++b)
			{
				auto &ba = std::get<0>(*b);

				// mmCIF...
				if (aa.is_water())
				{
					if (aa.get_auth_seq_id() != ba.get_auth_seq_id() or aa.get_auth_asym_id() != ba.get_auth_asym_id() or aa.symmetry() != ba.symmetry())
						continue;
				}
				else if (ba.get_label_atom_id() != aa.get_label_atom_id() or ba.get_label_comp_id() != aa.get_label_comp_id() or aa.get_label_seq_id() != ba.get_label_seq_id() or aa.get_label_asym_id() != ba.get_label_asym_id() or aa.symmetry() != ba.symmetry())
					continue;

				assert(std::get<1>(*b) > std::get<1>(*a));

				is.lig.erase(b);
				break;
			}
		}

		// -----------------------------------------------------------------------

		if (cif::VERBOSE)
		{
			std::cerr << "preliminary cluster: " << std::endl
					  << " metal: " << is.ion << std::endl;

			for (auto &l : is.lig)
			{
				auto &a = std::get<0>(l);
				std::cerr << " " << a << ' ' << a.symmetry() << " @ " << std::get<1>(l) << std::endl;
			}
		}

		// -----------------------------------------------------------------------

		for (auto a = is.lig.begin(); a != is.lig.end() and std::next(a) != is.lig.end(); ++a)
		{
			auto &aa = std::get<0>(*a);

			if (aa.get_type() != cif::atom_type::N)
				continue;

			try
			{
				if (aa.get_label_comp_id() == "GLN")
				{
					assert(aa.get_label_atom_id() == "NE2");
					auto o = structure.get_atom_by_label("OE1", aa.get_label_asym_id(), "GLN", aa.get_label_seq_id(), aa.get_label_alt_id());

					if (cif::VERBOSE)
						std::cerr << "Flipping side chain for GLN "
								  << " " << aa.get_label_asym_id() << aa.get_label_seq_id() << " (" << aa.pdb_id() << ')' << std::endl;

					structure.swap_atoms(aa, o);
					continue;
				}

				if (aa.get_label_comp_id() == "ASN")
				{
					assert(aa.get_label_atom_id() == "ND2");
					auto o = structure.get_atom_by_label("OD1", aa.get_label_asym_id(), "ASN", aa.get_label_seq_id(), aa.get_label_alt_id());

					if (cif::VERBOSE)
						std::cerr << "Flipping side chain for ASN "
								  << " " << aa.get_label_asym_id() << aa.get_label_seq_id() << " (" << aa.pdb_id() << ')' << std::endl;

					structure.swap_atoms(aa, o);
					continue;
				}
			}
			catch (const std::out_of_range &ex)
			{
				if (cif::VERBOSE)
					std::cerr << "Could not flip " << aa.get_label_comp_id() << ": " << ex.what() << std::endl;

				// is.lig.clear();	// give up
				// break;
			}
		}

		// -----------------------------------------------------------------------
		// if there are less than six atoms or more than nine, give up.

		if (is.lig.size() < 6 or is.lig.size() > 9)
		{
			if (cif::VERBOSE)
				std::cerr << "Rejecting cluster since the number of atoms is not reasonable" << std::endl;
			break;
		}

		// However, if we have more than six, try to remove outliers with a Grubbs test
		// See: https://en.wikipedia.org/wiki/Grubbs%27s_test_for_outliers

		while (is.lig.size() > 6)
		{
			double sum = accumulate(is.lig.begin(), is.lig.end(), 0.0, [](double s, auto &l)
				{ return s + std::get<1>(l); });
			double avg = sum / is.lig.size();
			double stddev = sqrt(accumulate(is.lig.begin(), is.lig.end(), 0.0, [avg](double s, auto &l)
									 { return s + (std::get<1>(l) - avg) * (std::get<1>(l) - avg); }) /
								 (is.lig.size() - 1));

			// only test if max distance is outlier
			double G = (std::get<1>(is.lig.back()) - avg) / stddev;

			// extracted from student t-distribution table, with one-sided confidence level 90%
			// and degrees of freedom

			if (G < get_t_90(is.lig.size()))
				break;

			if (cif::VERBOSE)
			{
				auto &a = std::get<0>(is.lig.back());
				std::cerr << "Removing outlier " << a.get_label_asym_id() << a.get_label_seq_id() << '/' << a.get_label_atom_id() << " (" << a.pdb_id() << ')' << ' ' << a.symmetry() << " @ " << std::get<1>(is.lig.back()) << std::endl;
			}

			is.lig.erase(is.lig.begin() + is.lig.size() - 1);
		}

		if (not is.isOctaHedral())
		{
			if (cif::VERBOSE)
				std::cerr << "Rejecting cluster since it is not an octahedral" << std::endl;
			break;
		}

		result = true;
		break;
	}

	return result;
}

std::vector<IonSite> findOctahedralSites(cif::mm::structure &structure, cif::datablock &db, int spacegroup, const clipper::Cell &cell)
{
	std::vector<IonSite> result;

	// factory for symmetry atom iterators
	pdb_redo::SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom : structure.atoms())
	{
		auto compID = atom.get_label_comp_id();

		if (compID == "NA" or compID == "K")
			atom.set_property("pdbx_formal_charge", 1);

		if (compID == "CA" or compID == "MG")
			atom.set_property("pdbx_formal_charge", 2);

		if (compID != "NA" and compID != "MG")
			continue;

		IonSite is = { atom };

		for (auto a : structure.atoms())
		{
			if (a.get_type() == cif::atom_type::S or
				a.get_type() == cif::atom_type::O or
				a.get_type() == cif::atom_type::N)
			{
				// for (auto sa: saif(a))
				for (auto sa : saif(a, [al = atom.get_location()](const cif::point &pt)
						 { return distance(al, pt) <= kMaxMetalLigandDistance; }))
				{
					float d = distance(atom, sa);
					assert(d <= kMaxMetalLigandDistance);
					is.lig.emplace_back(sa, d, sa.symmetry());
				}
			}
		}

		std::set<std::string> altIDs;
		for (const auto &[a, f, s] : is.lig)
		{
			if (not a.get_label_alt_id().empty())
				altIDs.insert(a.get_label_alt_id());
		}

		if (altIDs.empty()) // no alternates
			altIDs.insert("");

		for (auto &alt : altIDs)
		{
			IonSite ais = is;
			if (findOctahedralSites(structure, db, spacegroup, cell, ais, alt))
				result.emplace_back(std::move(ais));
		}
	}

	return result;
}

void updateSkipLists(const IonSite &ion, pdb_redo::SkipList &water, pdb_redo::SkipList &pepflipO, pdb_redo::SkipList &pepflipN, pdb_redo::SkipList &sideaid)
{
	for (const auto &[atom, distance, ignore] : ion.lig)
	{
		if (atom.is_water())
			water.push_back(atom);
		else
		{
			auto compound = cif::compound_factory::instance().create(atom.get_label_comp_id());
			if (not compound)
				continue;

			// if (compound->group() != "peptide" and compound->group() != "p-peptide" and compound->group() != "m-peptide")
			// 	continue;

			if (compound->type() != "peptide linking")
				continue;

			if (atom.is_back_bone())
			{
				if (atom.get_type() == cif::O)
					pepflipO.push_back(atom);
				else if (atom.get_type() == cif::N)
					pepflipN.push_back(atom);
			}
			else
				sideaid.push_back(atom);
		}
	}
}

// --------------------------------------------------------------------

int pr_main(int argc, char *argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	int result = 0;

	auto &config = cfg::config::instance();

	config.init("platonyzer [options] inputfile outputfile]",
		cfg::make_option("help,h", "Display help message"),
		cfg::make_option("version", "Print version"),
		cfg::make_option("verbose,v", "Verbose output"),
		cfg::make_option<std::string>("skip-list-format", "old", "Format to use for the skip lists, one of 'old', 'json' or 'cif'"),
		cfg::make_option("delete-vdw-rest", "Delete vanderWaals restraints for octahedral ions in the external for Refmac"),
		cfg::make_option("create-na-mg-links", "Create links for Na/Mg ion sites that were found"),
		// ( "pdb<std::string>-redo-data", "The PDB-REDO dat file" /*, default is the built in one"*/),
		cfg::make_option<std::string>("dict", "Dictionary file containing restraints for residues in this specific target"));

	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().size() != 2)
	{
		std::cerr << config << std::endl;
		exit(config.has("help") ? 0 : 1);
	}

	// Load dict, if any

	if (config.has("dict"))
		cif::compound_factory::instance().push_dictionary(config.get<std::string>("dict"));

	cif::VERBOSE = config.count("verbose");

	if (cif::VERBOSE)
		std::cerr << "Loading data...";

	fs::path input = config.operands().front();
	cif::file pdb = cif::pdb::read(input);
	cif::mm::structure structure(pdb);

	if (cif::VERBOSE)
		std::cerr << " done" << std::endl;

	auto &db = pdb.front();

	// -----------------------------------------------------------------------

	std::string entryId = db["entry"].front()["id"].as<std::string>();
	if (entryId.empty())
		throw std::runtime_error("Missing _entry.id in coordinates file");

	const auto &[a, b, c, alpha, beta, gamma] = db["cell"].find1<double, double, double, double, double, double>("entry_id"_key == entryId,
		"length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma");

	clipper::Cell cell(clipper::Cell_descr(a, b, c, alpha, beta, gamma));

	std::string spacegroupName = db["symmetry"].find1<std::string>("entry_id"_key == entryId, "space_group_name_H-M");
	int spacegroupNr = cif::get_space_group_number(spacegroupName);

	// -----------------------------------------------------------------------

	fs::path outfile = config.operands().back();
	fs::path outfile_extra = outfile;

	RestraintGenerator rg(outfile_extra.replace_extension(".restraints").string(), config.has("delete-vdw-rest"));

	auto &structConn = db["struct_conn"];

	std::size_t removedLinks = 0, createdLinks = 0;
	std::size_t platonyzerLinkId = 1;

	bool createNaMgLinks = config.has("create-na-mg-links");

	pdb_redo::SkipList waters, pepflipO, pepflipN, sideaid;

	for (auto &ionSites : {
			 findZincSites(structure, db, spacegroupNr, cell),
			 findOctahedralSites(structure, db, spacegroupNr, cell) })
	{
		for (auto ionSite : ionSites)
		{
			// write restraints
			std::vector<cif::mm::atom> ligands;
			transform(ionSite.lig.begin(), ionSite.lig.end(), back_inserter(ligands), [](auto &l)
				{ return std::get<0>(l); });

			switch (ionSite.lig.size())
			{
				case 4:
					rg.writeTetrahedral(ionSite.ion, ligands);
					break;
				case 6:
					rg.writeOctahedral(ionSite.ion, ligands, ionSite.opposing);
					updateSkipLists(ionSite, waters, pepflipO, pepflipN, sideaid);
					break;
			};

			// replace LINK/struct_conn records
			std::size_t n = structConn.size();
			structConn.erase(
				("ptnr1_label_asym_id"_key == ionSite.ion.get_label_asym_id() and "ptnr1_label_atom_id"_key == ionSite.ion.get_label_atom_id()) or
				("ptnr2_label_asym_id"_key == ionSite.ion.get_label_asym_id() and "ptnr2_label_atom_id"_key == ionSite.ion.get_label_atom_id()));

			for (auto &atom : ligands)
			{
				structConn.erase(
					"conn_type_id"_key == "disulf" and (("ptnr1_label_asym_id"_key == atom.get_label_asym_id() and "ptnr1_label_atom_id"_key == atom.get_label_atom_id()) or
														   ("ptnr2_label_asym_id"_key == atom.get_label_asym_id() and "ptnr2_label_atom_id"_key == atom.get_label_atom_id())));
			}
			removedLinks += (n - structConn.size());

			if (not createNaMgLinks and (ionSite.ion.get_type() == cif::atom_type::Na or ionSite.ion.get_type() == cif::atom_type::Mg))
				continue;

			createdLinks += ionSite.lig.size();

			for (auto &&[atom, distance, symop] : ionSite.lig)
			{
				std::string id = "metalc_p-" + std::to_string(platonyzerLinkId++);

				structConn.emplace({ { "id", id },
					{ "conn_type_id", "metalc" },
					{ "ptnr1_label_asym_id", atom.get_label_asym_id() },
					{ "ptnr1_label_comp_id", atom.get_label_comp_id() },
					{ "ptnr1_label_seq_id", atom.get_label_seq_id() },
					{ "ptnr1_label_atom_id", atom.get_label_atom_id() },
					{ "pdbx_ptnr1_label_alt_id", atom.get_label_alt_id() },
					{ "pdbx_ptnr1_PDB_ins_code", atom.get_pdb_ins_code() },
					{ "ptnr1_symmetry", "1_555" },
					{ "ptnr1_auth_asym_id", atom.get_auth_asym_id() },
					{ "ptnr1_auth_comp_id", atom.get_label_comp_id() },
					{ "ptnr1_auth_seq_id", atom.get_auth_seq_id() },

					{ "ptnr2_label_asym_id", ionSite.ion.get_label_asym_id() },
					{ "ptnr2_label_comp_id", ionSite.ion.get_label_comp_id() },
					{ "ptnr2_label_seq_id", "?" },
					{ "ptnr2_label_atom_id", ionSite.ion.get_label_atom_id() },
					{ "pdbx_ptnr2_label_alt_id", ionSite.ion.get_label_alt_id() },
					{ "pdbx_ptnr2_PDB_ins_code", ionSite.ion.get_pdb_ins_code() },
					{ "ptnr2_auth_asym_id", ionSite.ion.get_auth_asym_id() },
					{ "ptnr2_auth_comp_id", ionSite.ion.get_label_comp_id() },
					{ "ptnr2_auth_seq_id", ionSite.ion.get_auth_seq_id() },
					{ "ptnr2_symmetry", symop },

					{ "pdbx_dist_value", distance } });
			}
		}
	}

	// -----------------------------------------------------------------------

	if (cif::VERBOSE)
		std::cerr << "Removed " << removedLinks << " link records" << std::endl
				  << "Created " << createdLinks << " link records" << std::endl;

	// -----------------------------------------------------------------------

	auto &software = db["software"];
	software.emplace({ { "pdbx_ordinal", software.get_unique_id("") },
		{ "name", "platonyzer" },
		{ "version", kVersionNumber },
		{ "date", kBuildDate },
		{ "classification", "other" } });

	pdb.save(outfile);

	pdb_redo::SkipListFormat fmt;
	if (config.get<std::string>("skip-list-format") == "old")
		fmt = pdb_redo::SkipListFormat::OLD;
	else if (config.get<std::string>("skip-list-format") == "json")
		fmt = pdb_redo::SkipListFormat::JSON;
	else if (config.get<std::string>("skip-list-format") == "cif")
		fmt = pdb_redo::SkipListFormat::CIF;

	writeSkipList(outfile_extra.replace_extension(".skip-waters"), waters, fmt);
	writeSkipList(outfile_extra.replace_extension(".skip-pepflipN"), pepflipN, fmt);
	writeSkipList(outfile_extra.replace_extension(".skip-pepflipO"), pepflipO, fmt);
	writeSkipList(outfile_extra.replace_extension(".skip-sideaid"), sideaid, fmt);

	return result;
}

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e)
{
	std::cerr << e.what() << std::endl;
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception &nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

int main(int argc, char *argv[])
{
	int result = -1;

	try
	{
		result = pr_main(argc, argv);
	}
	catch (std::exception &ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
