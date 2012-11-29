// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _SECTOR_H
#define _SECTOR_H

#include "libs.h"
#include "galaxy/SystemPath.h"
#include "galaxy/StarSystem.h"
#include "galaxy/CustomSystem.h"
#include <string>
#include <vector>

class Faction;

// Special case mtrand replacement solution for sector generation
// with the same interfaces.
// The distribution matters for our purposes. 
// This is a 4 seed rng that has a period of 2^128-1.
// XXX eventually this should be available through a complete replacement for MTRand 
class XorshiftRand { 
public:
	XorshiftRand(Uint32 seed_[4]) 
	{ 
		s[0] = seed_[0]; s[1] = seed_[1]; s[2] = seed_[2]; s[3] = seed_[3];
		// low sector numbers (e.g. near sol) can result in low rng output for a few calls
		SkipNCalls(5);
	}
	// valid for sector view where the 4th element of the seed array is assumed not to be 0
	inline XorshiftRand(int seed_[4]) {
		// avoid zero seed values by changing range of values from 
		// [-2^31-1, +2^31] to [1, 2^32-1]
		//printf("start seed    : s[0] %u, s[1] %u, s[2] %u, s[3] %u\n", seed_[0], seed_[1], seed_[2], seed_[3]);
		const Uint32 d = 0x7FFFFFFF;
		Uint32 t[4];

		s[0] = t[0] = (seed_[0] >= 0)?Uint32(seed_[0])+d:Uint32(d+seed_[0]); 
		s[1] = t[1] = (seed_[1] >= 0)?Uint32(seed_[1])+d:Uint32(d+seed_[1]);  
		s[2] = t[2] = (seed_[2] >= 0)?Uint32(seed_[2])+d:Uint32(d+seed_[2]);  
		s[3] = t[3] = Uint32(seed_[3]);    // universe seed is left unchanged

		GenerateState(t);
		//t[0] = SmallStateRNG(t[0], 8);
		//t[1] = SmallStateRNG(t[1], 8);
		//t[2] = SmallStateRNG(t[2], 8);
		//t[3] = SmallStateRNG(t[3], 8);

		// low sector numbers can result in low rng output for a few calls
		//SkipNCalls(2);
		//printf("initialised to: s[0] %u, s[1] %u, s[2] %u, s[3] %u\n", s[0], s[1], s[2], s[3]);
	}
	// Generate the 4 state variables using another smaller state rng with the given seed.
	// The 32 bit rng shift values used here are taken from p.4 of the xorshift paper.
	inline XorshiftRand(Uint32 seed_) 
	{
		Uint32 t = seed_; // 0 seeds are not allowed
		// use a 1 state variable rng to generate the 4 state variables
		// low seed numbers can result in low rng output (vissible correlation)
		for (int i = 0; i < 4; i++) {
			t ^= (t << 13); t = (t >> 17); t ^= (t << 5);
			s[i] = t;
		}
		SkipNCalls(5);
	}
	// taken from the other rng
	inline void GenerateState(Uint32 a[4]) {
		const int n = 4;
		const int size = 4;
		int i = 1, j = 0;
		for (int k = 4; k; --k) {
			s[i] = (s[i] ^ ((s[i - 1] ^ (s[i - 1] >> 30)) * 1664525UL))
			  + a[j] + j; // non linear
			++j; j %= size;
			if ((++i) == n) { s[0] = s[n - 1]; i = 1; }
		}
		for (int k = n - 1; k; --k) {
			s[i] = (s[i] ^ ((s[i - 1] ^ (s[i - 1] >> 30)) * 1566083941UL)) - i;
			if ((++i) == n) { s[0] = s[n - 1]; i = 1; }
		}

		// ensure 1st and last variables are not zero
		if ((s[0]|s[3]) == 0) {
			s[0]=(a[0]+a[1]+a[2]+a[3])|(0x01);
			s[3]=(a[0]^a[1]^a[2]^a[3])|(0x01);
		}
	}
	// 32 bit state xorshift rng used to alter seeds to avoid issues with zero valued seeds
	inline Uint32 SmallStateRNG(Uint32 t, int numIterations) const { 
		for (int i = 0; i < numIterations; i++) { 
			t ^= (t << 13); t = (t >> 17); t ^= (t << 5);
		}
		return t;
	}
	Uint32 Int32() { return Rand_int32(); }
	// [min,max]
	int Int32(int min, int max) { return Sint32(Rand_int32()%Uint32(1+max-min))+min; }
	// [0,max)
	Uint32 Int32(int max) { 
		assert(max > 0); 
		return Rand_int32()%max;
		Uint32 t = Rand_int32();
		//printf("Int32(%i)  : %i, original %u\n", max, t%max, t);
		//printf("State after: s[0] %u, s[1] %u, s[2] %u, s[3] %u\n", s[0], s[1], s[2], s[3]);
		return t%Uint32(max);
	}
	double Double() { return double(Rand_int32()) * (1. / 4294967296.); } // divided by 2^32
	// interval [0, max)
	double Double(double max) { return max*Double(); }
	Uint32 Rand_int32();
	void SkipNCalls(int n) { for (int i = 0; i < n; i++) { Rand_int32(); } }
	void PrintTestOutput();
private:
	// These are the 4 internal state variables
	Uint32 s[4];
};

class Sector {
public:
	// lightyears
	static const float SIZE;
	Sector(int x, int y, int z);
	static float DistanceBetween(const Sector *a, int sysIdxA, const Sector *b, int sysIdxB);
	static void Init();

	// Sector is within a bounding rectangle - used for SectorView m_sectorCache pruning.
	bool WithinBox(const int Xmin, const int Xmax, const int Ymin, const int Ymax, const int Zmin, const int Zmax) const;
	bool Contains(const SystemPath sysPath) const;

	// sets appropriate factions for all systems in the sector
	void AssignFactions();

	class System {
	public:
		System(int x, int y, int z): customSys(0), population(-1), sx(x), sy(y), sz(z) {};
		~System() {};

		// Check that we've had our habitation status set

		// public members
		std::string name;
		vector3f p;
		int numStars;
		SystemBody::BodyType starType[4];
		Uint32 seed;
		const CustomSystem *customSys;
		Faction *faction;
		fixed population;

		vector3f FullPosition() { return Sector::SIZE*vector3f(float(sx), float(sy), float(sz)) + p; };

		int sx, sy, sz;
	};
	std::vector<System> m_systems;

private:
	int sx, sy, sz;
	void GetCustomSystems();
	std::string GenName(System &sys, int si, XorshiftRand &rand);
};

#endif /* _SECTOR_H */
