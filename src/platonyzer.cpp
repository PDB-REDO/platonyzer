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

#include <boost/iostreams/concepts.hpp> // output_filter
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/operations.hpp> // put
#include <boost/program_options.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Compound.hpp"
#include "cif++/Structure.hpp"
#include "cif++/Symmetry.hpp"

#include "pdb-redo/SkipList.hpp"
#include "pdb-redo/Symmetry-2.hpp"

std::string VERSION_STRING;

namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace fs = std::filesystem;
namespace c = mmcif;

// -----------------------------------------------------------------------

const std::set<std::string> kBackBone = {
	"N", "CA", "C", "O", "OXT"};

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
		: m_file_file(file)
		, m_deleteVDWRestraints(deleteVDWRestraints)
	{
		if (not m_file_file.is_open())
			throw std::runtime_error("Could not open restraint file " + file);

		m_file.push(unique_filter());
		m_file.push(m_file_file);
	}

	~RestraintGenerator() = default;

	std::size_t writeTetrahedral(const c::Atom &ion, const std::vector<c::Atom> &ligands);
	std::size_t writeOctahedral(const c::Atom &ion, const std::vector<c::Atom> &ligands,
		const std::vector<std::tuple<std::size_t, std::size_t>> &opposing);

  private:
	class unique_filter : public io::output_filter
	{
	  public:
		template <typename Sink>
		bool put(Sink &dest, int c)
		{
			bool result = true;

			if (c == '\n')
			{
				if (m_seen.insert(m_line).second)
				{
					for (char ch : m_line)
						io::put(dest, ch);
					result = io::put(dest, c);
				}

				m_line.clear();
			}
			else
				m_line += c;

			return result;
		}

		template <typename Sink>
		void close(Sink &)
		{
			m_line.clear();
			m_seen.clear();
		}

	  private:
		std::string m_line;
		std::unordered_set<std::string> m_seen;
	};

	void writeAngleRestraint(float target, float sd, const c::Atom &a, const c::Atom &b, const c::Atom &c);

	struct AtomPart
	{
		const c::Atom &m_a;

		AtomPart(const c::Atom &a)
			: m_a(a)
		{
		}

		friend std::ostream &operator<<(std::ostream &os, const AtomPart &aw)
		{
			os << "chain " << aw.m_a.authAsymID()
			   << " resi " << aw.m_a.authSeqID()
			   << " ins " << (aw.m_a.pdbxAuthInsCode().empty() ? "." : aw.m_a.pdbxAuthInsCode())
			   << " atom " << aw.m_a.authAtomID();
			if (not aw.m_a.pdbxAuthAltID().empty())
				os << " alt " + aw.m_a.pdbxAuthAltID();
			else if (not aw.m_a.labelAltID().empty())
				os << " alt " + aw.m_a.labelAltID();
			os << " symm " << (aw.m_a.isSymmetryCopy() ? "Y" : "N");

			return os;
		}
	};

	std::ofstream m_file_file;
	io::filtering_ostream m_file;
	bool m_deleteVDWRestraints;
};

std::size_t RestraintGenerator::writeTetrahedral(const c::Atom &ion, const std::vector<c::Atom> &ligands)
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

std::size_t RestraintGenerator::writeOctahedral(const c::Atom &ion, const std::vector<c::Atom> &ligands,
	const std::vector<std::tuple<std::size_t, std::size_t>> &opposing)
{
	const float kSD = 3;

	if (m_deleteVDWRestraints)
		m_file << "vdwr exclude " << AtomPart(ion) << std::endl;

	std::set<std::size_t> one_side_ligand;
	for (const auto& [a, b]: opposing)
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

void RestraintGenerator::writeAngleRestraint(float target, float sd, const c::Atom &a, const c::Atom &b, const c::Atom &c)
{
	m_file << "exte angle first " << AtomPart(a) << " next " << AtomPart(b) << " next " << AtomPart(c)
		   << std::fixed << std::setprecision(3) << " value " << target << " sigma " << sd << std::endl;
}

// ------------------------------------------------------------------------

struct IonSite
{
	c::Atom ion;
	std::vector<std::tuple<c::Atom, float, std::string>> lig;
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
			float angle = c::Angle(la.location(), ion.location(), lb.location());

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

bool findZincSites(c::Structure &structure, cif::Datablock &db, int spacegroup, const clipper::Cell &cell,
	IonSite &zs, const std::string &altID)
{
	bool result = false;

	for (;;)
	{
		// strip out the atoms that are not available in this alt group
		if (not altID.empty())
		{
			zs.lig.erase(
				std::remove_if(zs.lig.begin(), zs.lig.end(), [altID](auto &l) {
					auto alt = std::get<0>(l).labelAltID();
					return not alt.empty() and alt != altID;
				}),
				zs.lig.end());
		}

		// sort ligands on distance
		std::sort(zs.lig.begin(), zs.lig.end(), [](auto &a, auto &b) { return std::get<1>(a) < std::get<1>(b); });

		// take the nearest atom of a residue, this takes care of NE2/ND1 from a single HIS
		// and also the alt cases
		for (auto a = zs.lig.begin(); a != zs.lig.end() and std::next(a) != zs.lig.end(); ++a)
		{
			auto &aa = std::get<0>(*a);

			auto ad = std::get<1>(*a);

			for (auto b = std::next(a); b != zs.lig.end(); ++b)
			{
				auto &ba = std::get<0>(*b);

				if (ba.labelCompID() != aa.labelCompID() or aa.labelSeqID() != ba.labelSeqID() or aa.labelAsymID() != ba.labelAsymID() or aa.symmetry() != ba.symmetry())
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
					  << " zn: " << zs.ion.labelAsymID() << '/' << zs.ion.labelAtomID() << " (" << zs.ion.pdbID() << ')' << std::endl;

			for (auto &l : zs.lig)
			{
				auto &a = std::get<0>(l);
				std::cerr << " " << a.labelAsymID() << a.labelSeqID() << '/' << a.labelAtomID() << " (" << a.pdbID() << ')' << ' ' << a.symmetry() << " @ " << std::get<1>(l) << std::endl;
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
					std::cerr << "Atom " << a.labelAsymID() << a.labelSeqID() << '/' << a.labelAtomID() << " (" << a.pdbID() << ')' << " was considered to be an outlier" << std::endl;
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

std::vector<IonSite> findZincSites(c::Structure &structure, cif::Datablock &db, int spacegroup, const clipper::Cell &cell)
{
	std::vector<IonSite> result;

	// factory for symmetry atom iterators
	SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom : structure.atoms())
	{
		if (atom.labelCompID() != "ZN")
			continue;

		atom.property("pdbx_formal_charge", 2);

		IonSite zs = {atom};

		for (auto a : structure.atoms())
		{
			if (a.labelCompID() == "HIS" and (a.labelAtomID() == "ND1" or a.labelAtomID() == "NE2"))
			{
				for (auto sa : saif(a, [al = atom.location()](const mmcif::Point &pt) { return Distance(al, pt) <= kMaxZnHisDistanceInCluster; }))
				{
					float d = Distance(atom, sa);
					assert(d <= kMaxZnHisDistanceInCluster);
					zs.lig.emplace_back(sa, d, sa.symmetry());
				}
				continue;
			}

			if (a.labelCompID() == "CYS" and a.labelAtomID() == "SG")
			{
				for (auto sa : saif(a, [al = atom.location()](const mmcif::Point &pt) { return Distance(al, pt) <= kMaxZnCysDistanceInCluster; }))
				{
					float d = Distance(atom, sa);
					assert(d <= kMaxZnCysDistanceInCluster);
					zs.lig.emplace_back(sa, d, sa.symmetry());
				}
				continue;
			}
		}

		std::set<std::string> altIDs;
		for (const auto &[a, f, s] : zs.lig)
		{
			if (not a.labelAltID().empty())
				altIDs.insert(a.labelAltID());
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

bool findOctahedralSites(c::Structure &structure, cif::Datablock &db, int spacegroup, const clipper::Cell &cell,
	IonSite &is, const std::string &altID)
{
	bool result = false;

	for (;;)
	{
		// strip out the atoms that are not available in this alt group
		if (not altID.empty())
		{
			is.lig.erase(
				std::remove_if(is.lig.begin(), is.lig.end(), [altID](auto &l) {
					auto alt = std::get<0>(l).labelAltID();
					return not alt.empty() and alt != altID;
				}),
				is.lig.end());
		}

		// sort ligands on distance
		sort(is.lig.begin(), is.lig.end(), [](auto &a, auto &b) { return std::get<1>(a) < std::get<1>(b); });

		// take the nearest atom of a residue, this takes care of NE2/ND1 from a single HIS
		// and also the alt cases
		for (auto a = is.lig.begin(); a != is.lig.end() and std::next(a) != is.lig.end(); ++a)
		{
			auto &aa = std::get<0>(*a);

			for (auto b = std::next(a); b != is.lig.end(); ++b)
			{
				auto &ba = std::get<0>(*b);

				// mmCIF...
				if (aa.isWater())
				{
					if (aa.authSeqID() != ba.authSeqID() or aa.authAsymID() != ba.authAsymID() or aa.symmetry() != ba.symmetry())
						continue;
				}
				else if (ba.labelAtomID() != aa.labelAtomID() or ba.labelCompID() != aa.labelCompID() or aa.labelSeqID() != ba.labelSeqID() or aa.labelAsymID() != ba.labelAsymID() or aa.symmetry() != ba.symmetry())
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

			if (aa.type() != c::AtomType::N)
				continue;

			try
			{
				if (aa.labelCompID() == "GLN")
				{
					assert(aa.labelAtomID() == "NE2");
					auto o = structure.getAtomByLabel("OE1", aa.labelAsymID(), "GLN", aa.labelSeqID(), aa.labelAltID());

					if (cif::VERBOSE)
						std::cerr << "Flipping side chain for GLN "
								  << " " << aa.labelAsymID() << aa.labelSeqID() << " (" << aa.pdbID() << ')' << std::endl;

					structure.swapAtoms(aa, o);
					continue;
				}

				if (aa.labelCompID() == "ASN")
				{
					assert(aa.labelAtomID() == "ND2");
					auto o = structure.getAtomByLabel("OD1", aa.labelAsymID(), "ASN", aa.labelSeqID(), aa.labelAltID());

					if (cif::VERBOSE)
						std::cerr << "Flipping side chain for ASN "
								  << " " << aa.labelAsymID() << aa.labelSeqID() << " (" << aa.pdbID() << ')' << std::endl;

					structure.swapAtoms(aa, o);
					continue;
				}
			}
			catch (const std::out_of_range &ex)
			{
				if (cif::VERBOSE)
					std::cerr << "Could not flip " << aa.labelCompID() << ": " << ex.what() << std::endl;

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
			double sum = accumulate(is.lig.begin(), is.lig.end(), 0.0, [](double s, auto &l) { return s + std::get<1>(l); });
			double avg = sum / is.lig.size();
			double stddev = sqrt(accumulate(is.lig.begin(), is.lig.end(), 0.0, [avg](double s, auto &l) { return s + (std::get<1>(l) - avg) * (std::get<1>(l) - avg); }) / (is.lig.size() - 1));

			// only test if max distance is outlier
			double G = (std::get<1>(is.lig.back()) - avg) / stddev;

			// extracted from student t-distribution table, with one-sided confidence level 90%
			// and degrees of freedom

			if (G < get_t_90(is.lig.size()))
				break;

			if (cif::VERBOSE)
			{
				auto &a = std::get<0>(is.lig.back());
				std::cerr << "Removing outlier " << a.labelAsymID() << a.labelSeqID() << '/' << a.labelAtomID() << " (" << a.pdbID() << ')' << ' ' << a.symmetry() << " @ " << std::get<1>(is.lig.back()) << std::endl;
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

std::vector<IonSite> findOctahedralSites(c::Structure &structure, cif::Datablock &db, int spacegroup, const clipper::Cell &cell)
{
	std::vector<IonSite> result;

	// factory for symmetry atom iterators
	SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom : structure.atoms())
	{
		auto compID = atom.labelCompID();

		if (compID == "NA" or compID == "K")
			atom.property("pdbx_formal_charge", 1);

		if (compID == "CA" or compID == "MG")
			atom.property("pdbx_formal_charge", 2);

		if (compID != "NA" and compID != "MG")
			continue;

		IonSite is = {atom};

		for (auto a : structure.atoms())
		{
			if (a.type() == c::AtomType::S or
				a.type() == c::AtomType::O or
				a.type() == c::AtomType::N)
			{
				// for (auto sa: saif(a))
				for (auto sa : saif(a, [al = atom.location()](const mmcif::Point &pt) { return Distance(al, pt) <= kMaxMetalLigandDistance; }))
				{
					float d = Distance(atom, sa);
					assert(d <= kMaxMetalLigandDistance);
					is.lig.emplace_back(sa, d, sa.symmetry());
				}
			}
		}

		std::set<std::string> altIDs;
		for (const auto &[a, f, s] : is.lig)
		{
			if (not a.labelAltID().empty())
				altIDs.insert(a.labelAltID());
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

void updateSkipLists(const IonSite &ion, SkipList &water, SkipList &pepflipO, SkipList &pepflipN, SkipList &sideaid)
{
	for (const auto& [atom, distance, ignore] : ion.lig)
	{
		if (atom.isWater())
			water.push_back(atom);
		else
		{
			auto compound = mmcif::CompoundFactory::instance().create(atom.labelCompID());
			if (not compound)
				continue;
			
			// if (compound->group() != "peptide" and compound->group() != "p-peptide" and compound->group() != "m-peptide")
			// 	continue;
			
			if (compound->type() != "peptide linking")
				continue;
			
			if (atom.isBackBone())
			{
				if (atom.type() == mmcif::O)
					pepflipO.push_back(atom);
				else if (atom.type() == mmcif::N)
					pepflipN.push_back(atom);
			}
			else
				sideaid.push_back(atom);
		}
	}
}

// --------------------------------------------------------------------

namespace
{
std::string gVersionNr, gVersionDate;
}

void load_version_info()
{
	const std::regex
		rxVersionNr(R"(build-(\d+)-g[0-9a-f]{7}(-dirty)?)"),
		rxVersionDate(R"(Date: +(\d{4}-\d{2}-\d{2}).*)"),
		rxVersionNr2(R"(platonyzer-version: (\d+(?:\.\d+)+))");

#include "revision.hpp"

	struct membuf : public std::streambuf
	{
		membuf(char *data, size_t length) { this->setg(data, data, data + length); }
	} buffer(const_cast<char *>(kRevision), sizeof(kRevision));

	std::istream is(&buffer);

	std::string line;

	while (getline(is, line))
	{
		std::smatch m;

		if (std::regex_match(line, m, rxVersionNr))
		{
			gVersionNr = m[1];
			if (m[2].matched)
				gVersionNr += '*';
			continue;
		}

		if (std::regex_match(line, m, rxVersionDate))
		{
			gVersionDate = m[1];
			continue;
		}

		// always the first, replace with more specific if followed by the other info
		if (std::regex_match(line, m, rxVersionNr2))
		{
			gVersionNr = m[1];
			continue;
		}
	}

	if (not VERSION_STRING.empty())
		VERSION_STRING += "\n";
	VERSION_STRING += gVersionNr + " " + gVersionDate;
}

std::string get_version_nr()
{
	return gVersionNr /* + '/' + cif::get_version_nr()*/;
}

std::string get_version_date()
{
	return gVersionDate;
}

// --------------------------------------------------------------------

int pr_main(int argc, char *argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	int result = 0;

	po::options_description visible_options("platonyzer "s + VERSION_STRING + " options file]");
	visible_options.add_options()
		( "output,o",	po::value<std::string>(), "The output file, default is stdout" )
		( "skip-list-format", po::value<std::string>()->default_value("old"),
									"Format to use for the skip lists, one of 'old', 'json' or 'cif'")
		( "delete-vdw-rest",		"Delete vanderWaals restraints for octahedral ions in the external for Refmac" )
		( "create-na-mg-links",		"Create links for Na/Mg ion sites that were found" )
		( "help,h",					"Display help message" )
		( "version",				"Print version" )
		( "verbose,v",				"Verbose output" )
		// ( "pdb-redo-data", po::value<std::string>(),	"The PDB-REDO dat file" /*, default is the built in one"*/)
		( "dict", po::value<std::string>(), "Dictionary file containing restraints for residues in this specific target" );

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		( "input,i",	po::value<std::string>(),	"Input files" )
		( "debug,d",	po::value<int>(),			"Debug level (for even more verbose output)" );

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 2);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = "platonyzer.conf";
	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / "platonyzer.conf";

	if (fs::exists(configFile))
	{
		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}

	po::notify(vm);

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << VERSION_STRING << std::endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0 or vm.count("output") == 0)
	{
		std::cerr << visible_options << std::endl;
		exit(1);
	}

	// Load dict, if any

	if (vm.count("dict"))
		c::CompoundFactory::instance().pushDictionary(vm["dict"].as<std::string>());

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	if (cif::VERBOSE)
		std::cerr << "Loading data...";

	fs::path input = vm["input"].as<std::string>();
	c::File pdb(input);
	c::Structure structure(pdb);

	if (cif::VERBOSE)
		std::cerr << " done" << std::endl;

	auto &db = pdb.data();

	// -----------------------------------------------------------------------

	std::string entryId = db["entry"].front()["id"].as<std::string>();
	if (entryId.empty())
		throw std::runtime_error("Missing _entry.id in coordinates file");

	double a, b, c, alpha, beta, gamma;
	cif::tie(a, b, c, alpha, beta, gamma) = db["cell"]["entry_id"_key == entryId]
												.get("length_a", "length_b", "length_c",
													"angle_alpha", "angle_beta", "angle_gamma");

	clipper::Cell cell(clipper::Cell_descr(a, b, c, alpha, beta, gamma));

	std::string spacegroupName = db["symmetry"]
								   ["entry_id"_key == entryId]
								   ["space_group_name_H-M"]
									   .as<std::string>();

	int spacegroupNr = mmcif::GetSpacegroupNumber(spacegroupName);

	// -----------------------------------------------------------------------

	fs::path outfile = vm["output"].as<std::string>();
	fs::path outfile_extra = outfile;

	RestraintGenerator rg(outfile_extra.replace_extension(".restraints"), vm.count("delete-vdw-rest"));

	auto &structConn = db["struct_conn"];

	std::size_t removedLinks = 0, createdLinks = 0;
	std::size_t platonyzerLinkId = 1;

	bool createNaMgLinks = vm.count("create-na-mg-links");

	SkipList waters, pepflipO, pepflipN, sideaid;

	for (auto &ionSites : {
			 findZincSites(structure, db, spacegroupNr, cell),
			 findOctahedralSites(structure, db, spacegroupNr, cell)})
	{
		for (auto ionSite : ionSites)
		{
			// write restraints
			std::vector<c::Atom> ligands;
			transform(ionSite.lig.begin(), ionSite.lig.end(), back_inserter(ligands), [](auto &l) { return std::get<0>(l); });

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
				("ptnr1_label_asym_id"_key == ionSite.ion.labelAsymID() and "ptnr1_label_atom_id"_key == ionSite.ion.labelAtomID()) or
				("ptnr2_label_asym_id"_key == ionSite.ion.labelAsymID() and "ptnr2_label_atom_id"_key == ionSite.ion.labelAtomID()));

			for (auto &atom : ligands)
			{
				structConn.erase(
					"conn_type_id"_key == "disulf" and (
					("ptnr1_label_asym_id"_key == atom.labelAsymID() and "ptnr1_label_atom_id"_key == atom.labelAtomID()) or
					("ptnr2_label_asym_id"_key == atom.labelAsymID() and "ptnr2_label_atom_id"_key == atom.labelAtomID())));
			}
			removedLinks += (n - structConn.size());

			if (not createNaMgLinks and (ionSite.ion.type() == c::AtomType::Na or ionSite.ion.type() == c::AtomType::Mg))
				continue;

			createdLinks += ionSite.lig.size();

			for (auto &&[atom, distance, symop] : ionSite.lig)
			{
				std::string id = "metalc_p-" + std::to_string(platonyzerLinkId++);

				structConn.emplace({{"id", id},
					{"conn_type_id", "metalc"},
					{"ptnr1_label_asym_id", atom.labelAsymID()},
					{"ptnr1_label_comp_id", atom.labelCompID()},
					{"ptnr1_label_seq_id", atom.labelSeqID()},
					{"ptnr1_label_atom_id", atom.labelAtomID()},
					{"pdbx_ptnr1_label_alt_id", atom.labelAltID()},
					{"pdbx_ptnr1_PDB_ins_code", atom.pdbxAuthInsCode()},
					{"ptnr1_symmetry", "1_555"},
					{"ptnr1_auth_asym_id", atom.authAsymID()},
					{"ptnr1_auth_comp_id", atom.authCompID()},
					{"ptnr1_auth_seq_id", atom.authSeqID()},

					{"ptnr2_label_asym_id", ionSite.ion.labelAsymID()},
					{"ptnr2_label_comp_id", ionSite.ion.labelCompID()},
					{"ptnr2_label_seq_id", "?"},
					{"ptnr2_label_atom_id", ionSite.ion.labelAtomID()},
					{"pdbx_ptnr2_label_alt_id", ionSite.ion.labelAltID()},
					{"pdbx_ptnr2_PDB_ins_code", ionSite.ion.pdbxAuthInsCode()},
					{"ptnr2_auth_asym_id", ionSite.ion.authAsymID()},
					{"ptnr2_auth_comp_id", ionSite.ion.authCompID()},
					{"ptnr2_auth_seq_id", ionSite.ion.authSeqID()},
					{"ptnr2_symmetry", symop},

					{"pdbx_dist_value", distance}});
			}
		}
	}

	// -----------------------------------------------------------------------

	if (cif::VERBOSE)
		std::cerr << "Removed " << removedLinks << " link records" << std::endl
				  << "Created " << createdLinks << " link records" << std::endl;

	// -----------------------------------------------------------------------

	db.add_software("platonyzer", "other", get_version_nr(), get_version_date());

	pdb.save(outfile);

	SkipListFormat fmt;
	if (vm["skip-list-format"].as<std::string>() == "old")
		fmt = SkipListFormat::OLD;
	else if (vm["skip-list-format"].as<std::string>() == "json")
		fmt = SkipListFormat::JSON;
	else if (vm["skip-list-format"].as<std::string>() == "cif")
		fmt = SkipListFormat::CIF;

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
		load_version_info();

		result = pr_main(argc, argv);
	}
	catch (std::exception &ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
