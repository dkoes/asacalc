/* 
Copyright (c) 2012 David Koes and University of Pittsburgh  
 
Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation files 
(the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:  The above copyright notice and 
this permission notice shall be included in all copies or substantial
portions of the Software.  
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
*  */

/* 9/18/2012  dkoes
 * This is a simple little program that uses Shake-Rupley and OpenBabel
 * to calculate the solvent accessible surface area of a molecule.
 * I'm trying to make it fast.
 *
 * If no file is given, assume a pdb on stdin.  Always outputs to stdout.
 *
 */
 
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <openbabel/mol.h>
#include <boost/algorithm/string.hpp>
#include <ANN/ANN.h>

#define MAXVDW (3.0) /* assumed maximum vdw */
using namespace std;
using namespace OpenBabel;
using namespace boost::algorithm;

struct Coord
{
	double x;
	double y;
	double z;

	Coord(): x(0), y(0), z(0) {}

	Coord(double a, double b, double c): x(a), y(b), z(c) {}

	double distSq(const Coord& r) const
	{
		double a = x-r.x;
		double b = y-r.y;
		double c = z-r.z;
		return a*a+b*b+c*c;
	}
};

//represents an atom from a pdb file
class Atom
{
	string line; //full line of pdb file
	Coord position; //coordinates
	double r; //vdw radius

public:
	Atom(): r(0) {}

	Atom(const string& l): line(l), r(0) //initialize from pdb file
	{
		//assume this an atom
		stringstream x(line.substr(30,8));
		stringstream y(line.substr(38,8));
		stringstream z(line.substr(46,8));

		x >> position.x;
		y >> position.y;
		z >> position.z;

		string el = line.substr(76,2);
		trim(el);
		//remove any digits or special characters
		trim_if(el, is_any_of("0123456789+-"));
		unsigned anum = etab.GetAtomicNum(el.c_str());
		r = etab.GetVdwRad(anum);
	}

	//return true if pdb line is an atom
	static bool isAtom(const string& l)
	{
		if(l.substr(0,4) == "ATOM" || l.substr(0,6) == "HETATM")
			return true;
		else
			return false;
	}

	const Coord& coord() const { return position; }
	const double radius() const { return r; }

	void print(FILE *out, double bfactor) const
	{
		if(bfactor >= 100)
			bfactor = 99.999;
		fprintf(out, "%s%6.3f%6.3f%s\n",line.substr(0,54).c_str(),bfactor,r,line.substr(66).c_str());
	}
};

void readAtoms(istream& in, vector<Atom>& mol)
{
	mol.clear();
	string line;
	while(getline(in, line))
	{
		if(Atom::isAtom(line))
			mol.push_back(Atom(line));
	}
}



//create uniformish sphere of n points
void generateSpherePoints(unsigned n, vector<Coord>& points)
{
	points.clear(); points.reserve(n);
    double inc = M_PI * (3 - sqrt(5));
    double offset = 2 / double(n);
    for(unsigned k = 0; k < n; k++)
    {
    	double y = k * offset - 1 + (offset / 2);
        double r = sqrt(1 - y*y);
        double phi = k * inc;
        points.push_back(Coord(cos(phi)*r, y, sin(phi)*r));
    }
}

//find all atoms close the the specified atomi,
//it's faster to do this once per an atom and then brute force the
//collisions with the sphere points
void find_neighbor_indices(vector<Atom>& mol, ANNkd_tree& tree,  double probe, unsigned atomi, vector<unsigned>& neighbor_indices)
{
	neighbor_indices.clear();
	ANNidx idxArray[mol.size()];

	const Coord& avec = mol[atomi].coord();
	double radius = mol[atomi].radius() + probe + probe;

	double maxr = radius+MAXVDW;
	ANNcoord querypt[3] = {avec.x, avec.y, avec.z};
	unsigned numneigh = tree.annkFRSearch(querypt, maxr*maxr, mol.size(), idxArray);

	for (unsigned idx = 0; idx < numneigh; idx++)
	{
		unsigned ai = idxArray[idx];
		if(ai != atomi)
		{
			const Atom& a = mol[ai];
			double distSq = avec.distSq(a.coord());
			double r = radius+a.radius();
			if(distSq < r*r)
				neighbor_indices.push_back(ai);
		}
	}
}

//compute asa of each atom of mol and print as it is computed
void printASA(vector<Atom>& mol, const vector<Coord>& points, ANNkd_tree& tree, double probe, FILE * out)
{
    const double Const = 4.0 * M_PI / double(points.size());
    Coord test_point;
    vector<unsigned> neighbor_indices;
    double total = 0;
	for (unsigned ai = 0, na = mol.size(); ai < na; ai++)
	{
		const Atom& a = mol[ai];
		find_neighbor_indices(mol, tree, probe, ai, neighbor_indices);

		double radius = probe + a.radius();
		unsigned n_accessible_point = 0;

		for(unsigned i = 0, n = points.size(); i < n; i++)
		{
			bool is_accessible = true;
			test_point.x = points[i].x*radius + a.coord().x;
			test_point.y = points[i].y*radius + a.coord().y;
			test_point.z = points[i].z*radius + a.coord().z;

			for(unsigned j = 0, m = neighbor_indices.size(); j < m; j++)
			{
				const Atom& atom_j = mol[neighbor_indices[j]];
				double r = atom_j.radius() + probe;
				double diff_sq = test_point.distSq(atom_j.coord());
				if(diff_sq < r*r)
				{
					is_accessible = false;
					break;
				}
			}
			if(is_accessible)
				n_accessible_point++;
		}

		double area = Const*n_accessible_point*radius*radius;
		total += area;
		a.print(out, area);
	}
	fprintf(out,"END\n");
	fprintf(out,"REMARK Total %f\n",total);
}

int main(int argc, char *argv[])
{
	unsigned n_sphere_point = 960;
	double probe = 1.4;

	vector<Atom> mol;
	if(argc > 1)
	{
		ifstream in(argv[1]);
		if(!in)
		{
			cerr << "Error opening " << argv[1] << endl;
			exit(-1);
		}
		readAtoms(in, mol);
	}
	else
	{
		readAtoms(cin, mol);
	}

	//initialize spatial index
	ANNpointArray atompts = annAllocPts(mol.size(),3);
	for(unsigned i = 0, n = mol.size(); i < n; i++)
	{
		atompts[i][0] = mol[i].coord().x;
		atompts[i][1] = mol[i].coord().y;
		atompts[i][2] = mol[i].coord().z;
	}
	ANNkd_tree tree(atompts, mol.size(), 3);

	vector<Coord> sphpoints;
	generateSpherePoints(n_sphere_point, sphpoints);

	printASA(mol, sphpoints, tree, probe, stdout);

	//clean up memory
	annDeallocPts(atompts);
	annClose();
}
