// Copyright � 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "SystemGenerator.h"
#include "Sector.h"
#include "CustomSystem.h"


//-----------------------------------------------------------------------------
// Tables & Constants

static const fixed AU_SOL_RADIUS = fixed(305,65536);

static const struct StarTypeInfo {
	SystemBody::BodySuperType supertype;
	int mass[2]; // min,max % sol for stars, unused for planets
	int radius[2]; // min,max % sol radii for stars, % earth radii for planets
	int tempMin, tempMax;
} starTypeInfo[] = {
	{
		SystemBody::SUPERTYPE_NONE, {}, {},
        0, 0
	}, {
		SystemBody::SUPERTYPE_STAR, //Brown Dwarf
		{2,8}, {10,30},
		1000, 2000
	}, {
		SystemBody::SUPERTYPE_STAR,  //white dwarf
		{20,100}, {1,2},
		4000, 40000
	}, {
		SystemBody::SUPERTYPE_STAR, //M
		{10,47}, {30,60},
		2000, 3500
	}, {
		SystemBody::SUPERTYPE_STAR, //K
		{50,78}, {60,100},
		3500, 5000
	}, {
		SystemBody::SUPERTYPE_STAR, //G
		{80,110}, {80,120},
		5000, 6000
	}, {
		SystemBody::SUPERTYPE_STAR, //F
		{115,170}, {110,150},
		6000, 7500
	}, {
		SystemBody::SUPERTYPE_STAR, //A
		{180,320}, {120,220},
		7500, 10000
	}, {
		SystemBody::SUPERTYPE_STAR,  //B
		{200,300}, {120,290},
		10000, 30000
	}, {
		SystemBody::SUPERTYPE_STAR, //O
		{300,400}, {200,310},
		30000, 60000
	}, {
		SystemBody::SUPERTYPE_STAR, //M Giant
		{60,357}, {2000,5000},
		2500, 3500
	}, {
		SystemBody::SUPERTYPE_STAR, //K Giant
		{125,500}, {1500,3000},
		3500, 5000
	}, {
		SystemBody::SUPERTYPE_STAR, //G Giant
		{200,800}, {1000,2000},
		5000, 6000
	}, {
		SystemBody::SUPERTYPE_STAR, //F Giant
		{250,900}, {800,1500},
		6000, 7500
	}, {
		SystemBody::SUPERTYPE_STAR, //A Giant
		{400,1000}, {600,1000},
		7500, 10000
	}, {
		SystemBody::SUPERTYPE_STAR,  //B Giant
		{500,1000}, {600,1000},
		10000, 30000
	}, {
		SystemBody::SUPERTYPE_STAR, //O Giant
		{600,1200}, {600,1000},
		30000, 60000
	}, {
		SystemBody::SUPERTYPE_STAR, //M Super Giant
		{1050,5000}, {7000,15000},
		2500, 3500
	}, {
		SystemBody::SUPERTYPE_STAR, //K Super Giant
		{1100,5000}, {5000,9000},
		3500, 5000
	}, {
		SystemBody::SUPERTYPE_STAR, //G Super Giant
		{1200,5000}, {4000,8000},
		5000, 6000
	}, {
		SystemBody::SUPERTYPE_STAR, //F Super Giant
		{1500,6000}, {3500,7000},
		6000, 7500
	}, {
		SystemBody::SUPERTYPE_STAR, //A Super Giant
		{2000,8000}, {3000,6000},
		7500, 10000
	}, {
		SystemBody::SUPERTYPE_STAR,  //B Super Giant
		{3000,9000}, {2500,5000},
		10000, 30000
	}, {
		SystemBody::SUPERTYPE_STAR, //O Super Giant
		{5000,10000}, {2000,4000},
		30000, 60000
	}, {
		SystemBody::SUPERTYPE_STAR, //M Hyper Giant
		{5000,15000}, {20000,40000},
		2500, 3500
	}, {
		SystemBody::SUPERTYPE_STAR, //K Hyper Giant
		{5000,17000}, {17000,25000},
		3500, 5000
	}, {
		SystemBody::SUPERTYPE_STAR, //G Hyper Giant
		{5000,18000}, {14000,20000},
		5000, 6000
	}, {
		SystemBody::SUPERTYPE_STAR, //F Hyper Giant
		{5000,19000}, {12000,17500},
		6000, 7500
	}, {
		SystemBody::SUPERTYPE_STAR, //A Hyper Giant
		{5000,20000}, {10000,15000},
		7500, 10000
	}, {
		SystemBody::SUPERTYPE_STAR,  //B Hyper Giant
		{5000,23000}, {6000,10000},
		10000, 30000
	}, {
		SystemBody::SUPERTYPE_STAR, //O Hyper Giant
		{10000,30000}, {4000,7000},
		30000, 60000
	}, {
		SystemBody::SUPERTYPE_STAR,  // M WF
		{2000,5000}, {2500,5000},
		25000, 35000
	}, {
		SystemBody::SUPERTYPE_STAR,  // B WF
		{2000,7500}, {2500,5000},
		35000, 45000
	}, {
		SystemBody::SUPERTYPE_STAR,  // O WF
		{2000,10000}, {2500,5000},
		45000, 60000
	}, {
		SystemBody::SUPERTYPE_STAR,  // S BH
		{20,2000}, {0,0},	// XXX black holes are < 1 Sol radii big; this is clamped to a non-zero value later
		10, 24
	}, {
		SystemBody::SUPERTYPE_STAR,  // IM BH
		{900000,1000000}, {100,500},
		1, 10
	}, {
		SystemBody::SUPERTYPE_STAR,  // SM BH
		{2000000,5000000}, {10000,20000},
		10, 24
	}
/*	}, {
		SystemBody::SUPERTYPE_GAS_GIANT,
		{}, 950, Lang::MEDIUM_GAS_GIANT,
	}, {
		SystemBody::SUPERTYPE_GAS_GIANT,
		{}, 1110, Lang::LARGE_GAS_GIANT,
	}, {
		SystemBody::SUPERTYPE_GAS_GIANT,
		{}, 1500, Lang::VERY_LARGE_GAS_GIANT,
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 1, Lang::ASTEROID,
		"icons/object_planet_asteroid.png"
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 2, "Large asteroid",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 26, "Small, rocky dwarf planet", // moon radius
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 26, "Small, rocky dwarf planet", // dwarf2 for moon-like colours
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 52, "Small, rocky planet with a thin atmosphere", // mars radius
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky frozen planet with a thin nitrogen atmosphere", // earth radius
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Dead world that once housed it's own intricate ecosystem.", // earth radius
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky planet with a carbon dioxide atmosphere",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky planet with a methane atmosphere",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Water world with vast oceans and a thick nitrogen atmosphere",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky planet with a thick carbon dioxide atmosphere",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky planet with a thick methane atmosphere",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Highly volcanic world",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "World with indigenous life and an oxygen atmosphere",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 60, "Marginal terraformed world with minimal plant life",
	}, {
		SystemBody::SUPERTYPE_ROCKY_PLANET,
		{}, 90, "Fully terraformed world with introduced species from numerous successful colonies",
	}, {
		SystemBody::SUPERTYPE_STARPORT,
		{}, 0, Lang::ORBITAL_STARPORT,
	}, {
		SystemBody::SUPERTYPE_STARPORT,
		{}, 0, Lang::STARPORT,
	}*/
};

//-----------------------------------------------------------------------------
// Build

/*
 * 0 - ~500ly from sol: explored
 * ~500ly - ~700ly (65-90 sectors): gradual
 * ~700ly+: unexplored
 */
const bool SystemGenerator::Unexplored() 
{
	if (IsCustom() && !SectorSystem().customSys->want_rand_explored ) {
		return !SectorSystem().customSys->explored;    // may be defined manually in the custom system definition
	} else {
		int dist = isqrt(1 + m_path.sectorX*m_path.sectorX + m_path.sectorY*m_path.sectorY + m_path.sectorZ*m_path.sectorZ);
		return (dist > 90) || (dist > 65 && rand1().Int32(dist) > 40);
	}
}

SystemBody* SystemGenerator::AddStarsTo(BodyList& bodies)
{
	SystemBody* rootBody;
	SystemBody *star[4];

	const int starCount = NumStars();

	// For clarity, I've separated out each case here despite there being code in common.
	switch (starCount) 
	{
		// system with a single star
		case 1: {
			rootBody = AddStarOfType(bodies, Name(), SectorSystem().starType[0]);
			break;
		}
		// binary system
		case 2: {
			star[0]  = AddStarOfType           (bodies, Name()+" A", SectorSystem().starType[0]);
			star[1]  = AddStarOfTypeLighterThan(bodies, Name()+" B", SectorSystem().starType[1], star[0]->mass);
			m_centGrav1 = AddGravPoint(bodies, Name()+" A,B", star[0], star[1], false);
			rootBody    = m_centGrav1;
			break;
		}
		// trinary system
		case 3: {
			star[0] = AddStarOfType           (bodies, Name()+" A", SectorSystem().starType[0]);
			star[1] = AddStarOfTypeLighterThan(bodies, Name()+" B", SectorSystem().starType[1], star[0]->mass);
			star[2] = AddStarOfTypeLighterThan(bodies, Name()+" C", SectorSystem().starType[2], star[0]->mass);				
			m_centGrav1 = AddGravPoint(bodies, Name()+" A,B", star[0], star[1], true);
			m_centGrav2 = star[2];
		    rootBody    = AddGravPoint(bodies, Name(), m_centGrav1, m_centGrav2, false, (star[0]->orbMax + star[2]->orbMax) * 4);
			break;
		}
		// quarternary(?) system
		case 4: {
			star[0] = AddStarOfType           (bodies, Name()+" A", SectorSystem().starType[0]);
			star[1] = AddStarOfTypeLighterThan(bodies, Name()+" B", SectorSystem().starType[1], star[0]->mass);
			star[2] = AddStarOfTypeLighterThan(bodies, Name()+" C", SectorSystem().starType[2], star[0]->mass);
			star[3] = AddStarOfTypeLighterThan(bodies, Name()+" D", SectorSystem().starType[3], star[2]->mass);
			m_centGrav1 = AddGravPoint(bodies, Name()+" A,B", star[0], star[1], true);
			m_centGrav2 = AddGravPoint(bodies, Name()+" C,D", star[2], star[3], false);
			rootBody = AddGravPoint(bodies, Name(), m_centGrav1, m_centGrav2, false, (star[0]->orbMax + star[2]->orbMax) * 4);
			break;
		}
		default: assert((starCount >= 1) && (starCount <= 4)); break;
	}

	return rootBody;
}


//-----------------------------------------------------------------------------
// Private Build

SystemBody* SystemGenerator::AddStarOfType(BodyList& bodies, std::string name, SystemBody::BodyType type)
{
	SystemBody* star = NewBody(bodies, NULL, name, type);
	MakeStar(star);
	return(star);
}

void SystemGenerator::MakeStar(SystemBody *sbody)
{
	MTRand &rand = rand1();
	const SystemBody::BodyType type = sbody->type;

	sbody->seed = rand.Int32();
	sbody->radius = fixed(rand.Int32(starTypeInfo[type].radius[0], starTypeInfo[type].radius[1]), 100);
	sbody->mass   = fixed(rand.Int32(starTypeInfo[type].mass[0], starTypeInfo[type].mass[1]), 100);
	sbody->averageTemp = rand.Int32(starTypeInfo[type].tempMin, starTypeInfo[type].tempMax);
}

SystemBody* SystemGenerator::AddStarOfTypeLighterThan(BodyList& bodies, std::string name, SystemBody::BodyType type, fixed maxMass)
{
	SystemBody* star = NewBody(bodies, NULL, name, type);
	int tries = 16;
	do {
		MakeStar(star);
	} while ((star->mass > maxMass) && (--tries));
	return star;
}

SystemBody* SystemGenerator::AddGravPoint(BodyList& bodies, std::string name, SystemBody* a, SystemBody* b, bool limitOrbit, fixed minDist)
{
	// make the grav point
	SystemBody* gravPoint = NewBody(bodies, NULL, name, SystemBody::TYPE_GRAVPOINT);
	gravPoint-> mass = a->mass + b->mass;

	// make the two bodies children of the grav point
	gravPoint->children.push_back(a);
	gravPoint->children.push_back(b);
	a->parent = gravPoint;
	b->parent = gravPoint;

	// make the two bodies a binary pair
	MakeBinaryPair(a, b, minDist);

	// retry until orbMax value less than 100 if requested
	while (limitOrbit && a->orbMax > fixed(100,1)) 	MakeBinaryPair(a, b, minDist);

	// return
	return gravPoint;
}


SystemBody* SystemGenerator::AddGravPoint(BodyList& bodies, std::string name, SystemBody* a, SystemBody* b, bool limitOrbit)
{
	return AddGravPoint(bodies, name, a, b, limitOrbit, (a->radius + b->radius) * AU_SOL_RADIUS);
}

void SystemGenerator::MakeBinaryPair(SystemBody *a, SystemBody *b, fixed minDist)
{
	MTRand &rand = rand1();
	
	fixed m = a->mass + b->mass;
	fixed a0 = b->mass / m;
	fixed a1 = a->mass / m;
	a->eccentricity = rand.NFixed(3);
	int mul = 1;

	do {
		switch (rand.Int32(3)) {
			case 2: a->semiMajorAxis = fixed(rand.Int32(100,10000), 100); break;
			case 1: a->semiMajorAxis = fixed(rand.Int32(10,1000), 100); break;
			default:
			case 0: a->semiMajorAxis = fixed(rand.Int32(1,100), 100); break;
		}
		a->semiMajorAxis *= mul;
		mul *= 2;
	} while (a->semiMajorAxis < minDist);

	a->orbit.eccentricity = a->eccentricity.ToDouble();
	a->orbit.semiMajorAxis = AU * (a->semiMajorAxis * a0).ToDouble();
	a->orbit.period = 60*60*24*365* a->semiMajorAxis.ToDouble() * sqrt(a->semiMajorAxis.ToDouble() / m.ToDouble());

	const float rotX = -0.5f*float(M_PI);//(float)(rand.Double()*M_PI/2.0);
	const float rotY = static_cast<float>(rand.Double(M_PI));
	a->orbit.rotMatrix = matrix4x4d::RotateYMatrix(rotY) * matrix4x4d::RotateXMatrix(rotX);
	b->orbit.rotMatrix = matrix4x4d::RotateYMatrix(rotY-M_PI) * matrix4x4d::RotateXMatrix(rotX);

	b->orbit.eccentricity = a->eccentricity.ToDouble();
	b->orbit.semiMajorAxis = AU * (a->semiMajorAxis * a1).ToDouble();
	b->orbit.period = a->orbit.period;

	fixed orbMin = a->semiMajorAxis - a->eccentricity*a->semiMajorAxis;
	fixed orbMax = 2*a->semiMajorAxis - orbMin;
	a->orbMin = orbMin;
	b->orbMin = orbMin;
	a->orbMax = orbMax;
	b->orbMax = orbMax;
}

//-----------------------------------------------------------------------------
// State

MTRand& SystemGenerator::rand1() 
{
	if (!m_rand1) {
		unsigned long _init[6] = { m_path.systemIndex, Uint32(m_path.sectorX), Uint32(m_path.sectorY), Uint32(m_path.sectorZ), UNIVERSE_SEED, SectorSystem().seed };
		m_rand1 = new MTRand(_init, 6);
	}
	return *m_rand1;
}

//-----------------------------------------------------------------------------
// Construction/Destruction

SystemGenerator::SystemGenerator(SystemPath& path): m_path(path), m_rand1(0), m_sector(m_path.sectorX, m_path.sectorY, m_path.sectorZ), m_centGrav1(0), m_centGrav2(0) 
{
	assert(m_path.systemIndex >= 0 && m_path.systemIndex < m_sector.m_systems.size());
}

SystemGenerator::~SystemGenerator() 
{
	if (m_rand1) delete m_rand1;
}
