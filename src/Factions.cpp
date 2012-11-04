// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Factions.h"
#include "galaxy/Sector.h"
#include "galaxy/SystemPath.h"

#include "LuaUtils.h"
#include "LuaVector.h"
#include "LuaFixed.h"
#include "LuaConstants.h"
#include "Polit.h"
#include "FileSystem.h"

const Uint32 Faction::BAD_FACTION_IDX      = UINT_MAX;
const Color  Faction::BAD_FACTION_COLOUR   = (0.8f,0.8f,0.8f,0.50f);
const float  Faction::FACTION_BASE_ALPHA   = 0.30f;
const double Faction::FACTION_CURRENT_YEAR = 3200;

typedef std::vector<Faction*>  FactionList;
typedef FactionList::iterator FactionIterator;
static FactionList s_factions;

// ------- Faction --------

struct FactionBuilder {
	Faction *fac;
	bool registered;
	bool skip;
};

static const char LuaFaction_TypeName[] = "Faction";

static FactionBuilder *l_fac_check_builder(lua_State *L, int idx) {
	FactionBuilder *facbld = static_cast<FactionBuilder*>(
		luaL_checkudata(L, idx, LuaFaction_TypeName));
	assert(facbld->fac);
	return facbld;
}

static Faction *l_fac_check(lua_State *L, int idx)
{
	return l_fac_check_builder(L, idx)->fac;
}

static int l_fac_new(lua_State *L)
{
	const char *name = luaL_checkstring(L, 2);

	FactionBuilder *facbld = static_cast<FactionBuilder*>(lua_newuserdata(L, sizeof(*facbld)));
	facbld->fac = new Faction;
	facbld->registered = false;
	facbld->skip       = false;
	luaL_setmetatable(L, LuaFaction_TypeName);

	facbld->fac->name = name;

	return 1;
}

#define LFAC_FIELD_SETTER_FIXED(luaname, fieldname)        \
	static int l_fac_ ## luaname (lua_State *L) {          \
		Faction *fac = l_fac_check(L, 1);                  \
		const fixed *value = LuaFixed::CheckFromLua(L, 2); \
		fac->fieldname = *value;                           \
		lua_settop(L, 1); return 1;                        \
	}

#define LFAC_FIELD_SETTER_FLOAT(luaname, fieldname)        \
	static int l_fac_ ## luaname (lua_State *L) {          \
		Faction *fac = l_fac_check(L, 1);                  \
		double value = luaL_checknumber(L, 2);             \
		fac->fieldname = value;                            \
		lua_settop(L, 1); return 1;                        \
	}

#define LFAC_FIELD_SETTER_INT(luaname, fieldname)          \
	static int l_fac_ ## luaname (lua_State *L) {          \
		Faction *fac = l_fac_check(L, 1);                  \
		int value = luaL_checkinteger(L, 2);               \
		fac->fieldname = value;                            \
		lua_settop(L, 1); return 1;                        \
	}

static int l_fac_description_short(lua_State *L)
{
	Faction *fac = l_fac_check(L, 1);
	fac->description_short = luaL_checkstring(L, 2);
	lua_settop(L, 1);
	return 1;
}

static int l_fac_description(lua_State *L)
{
	Faction *fac = l_fac_check(L, 1);
	fac->description = luaL_checkstring(L, 2);
	lua_settop(L, 1);
	return 1;
}

// weightings to use when picking a government type
static int l_fac_govtype_weight(lua_State *L)
{
	Faction *fac = l_fac_check(L, 1);
	const char *typeName = luaL_checkstring(L, 2);
	const Polit::GovType g = static_cast<Polit::GovType>(LuaConstants::GetConstant(L, "PolitGovType", typeName));
	const int32_t weight = luaL_checkint(L, 3);	// signed as we will need to compare with signed out of MTRand.Int32

	if (g < Polit::GOV_RAND_MIN || g > Polit::GOV_RAND_MAX) {
		pi_lua_warn(L,
			"government type out of range: Faction{%s}:govtype_weight('%s', %d)",
			fac->name.c_str(), typeName, weight);
		return 0;
	}

	if (weight < 0) {
		pi_lua_warn(L,
			"weight must a postive integer: Faction{%s}:govtype_weight('%s', %d)",
			fac->name.c_str(), typeName, weight);
		return 0;		
	}

	fac->govtype_weights.push_back(std::make_pair(g, weight));
	fac->govtype_weights_total += weight;
	lua_settop(L, 1);

	return 1;
}

// sector(x,y,x) + system index + body index = location in a (custom?) system of homeworld
static int l_fac_homeworld (lua_State *L)
{
	FactionBuilder *facbld = l_fac_check_builder(L, 1);
	Faction *fac = facbld->fac;
	Sint32 x = luaL_checkinteger(L, 2);
	Sint32 y = luaL_checkinteger(L, 3);
	Sint32 z = luaL_checkinteger(L, 4);
	Sint32 si = luaL_checkinteger(L, 5);
	Uint32 bi = luaL_checkinteger(L, 6);

	// This bit is so we can avoid blowing ourselves up with potentially out of
	// range system indexes when the Lua script *knows* that it doesn't know
	// what systems exist in a sector. A negative sysindex therefore means 
	// 'take your best shot'
	Sector sec(x,y,z);
	if ( (si < 0 && sec.m_systems.size() > 0)) si = 0;

	fac->homeworld = SystemPath(x,y,z,si,bi);
	fac->hasHomeworld = true;
	facbld->skip      = si<0;	// our best shot missed, skip creation
	lua_settop(L, 1);
	return 1;
}

LFAC_FIELD_SETTER_FLOAT(foundingDate, foundingDate)   // date faction came into existence
LFAC_FIELD_SETTER_FLOAT(expansionRate, expansionRate) // lightyears per year that the volume expands.

static int l_fac_military_name(lua_State *L)
{
	Faction *fac = l_fac_check(L, 1);
	fac->military_name = luaL_checkstring(L, 2);
	lua_settop(L, 1);
	return 1;
}

//military logo
static int l_fac_police_name(lua_State *L)
{
	Faction *fac = l_fac_check(L, 1);
	fac->police_name = luaL_checkstring(L, 2);
	lua_settop(L, 1);
	return 1;
}

//police logo
//goods/equipment availability (1-per-economy-type: aka agricultural, industrial, tourist, etc)
//goods/equipment legality
static int l_fac_illegal_goods_probability(lua_State *L)
{
	Faction *fac = l_fac_check(L, 1);
	const char *typeName = luaL_checkstring(L, 2);
	const Equip::Type e = static_cast<Equip::Type>(LuaConstants::GetConstant(L, "EquipType", typeName));
	const uint32_t probability = luaL_checkunsigned(L, 3);

	if (e < Equip::FIRST_COMMODITY || e > Equip::LAST_COMMODITY) {
		pi_lua_warn(L,
			"argument out of range: Faction{%s}:IllegalGoodsProbability('%s', %d)",
			fac->name.c_str(), typeName, probability);
		return 0;
	}

	if (probability > 100) {
		pi_lua_warn(L,
			"argument (probability 0-100) out of range: Faction{%s}:IllegalGoodsProbability('%s', %d)",
			fac->name.c_str(), typeName, probability);
		return 0;
	}

	fac->equip_legality[e] = probability;
	lua_settop(L, 1);

	return 1;
}

//ship availability
static int l_fac_colour(lua_State *L)
{
	Faction *fac = l_fac_check(L, 1);
	const float r = luaL_checknumber(L, 2);
	const float g = luaL_checknumber(L, 3);
	const float b = luaL_checknumber(L, 4);

	fac->colour = Color(r,g,b);

	lua_settop(L, 1);

	return 1;
}

#undef LFAC_FIELD_SETTER_FIXED
#undef LFAC_FIELD_SETTER_FLOAT
#undef LFAC_FIELD_SETTER_INT

static int l_fac_add_to_factions(lua_State *L)
{
	FactionBuilder *facbld = l_fac_check_builder(L, 1);
	Faction *fac = facbld->fac;

	const std::string factionName(luaL_checkstring(L, 2));

	if (!facbld->registered && !facbld->skip) {
		if (facbld->fac->hasHomeworld) {
			printf("l_fac_add_to_factions: added (%3d,%3d,%3d) f=%4.0f e=%2.2f '%s' [%s]\n"
				, fac->homeworld.sectorX, fac->homeworld.sectorY, fac->homeworld.sectorZ, fac->foundingDate, fac->expansionRate, fac->name.c_str(), factionName.c_str());
		}
		else {
			printf("l_fac_add_to_factions: added '%s' [%s]\n", fac->name.c_str(), factionName.c_str());
		}

		s_factions.push_back(facbld->fac);
		facbld->registered = true;
		return 0;
	} else if (facbld->skip) {
		printf("l_fac_add_to_factions: invalid homeworld, skipped (%3d,%3d,%3d) f=%4.0f e=%2.2f '%s' [%s]\n"
				, fac->homeworld.sectorX, fac->homeworld.sectorY, fac->homeworld.sectorZ, fac->foundingDate, fac->expansionRate, fac->name.c_str(), factionName.c_str());
	} else {
		return luaL_error(L, "faction '%s' already added\n", facbld->fac->name.c_str());
	}
}

static int l_fac_gc(lua_State *L)
{
	FactionBuilder *facbld = l_fac_check_builder(L, 1);
	if (!facbld->registered) {
		delete facbld->fac;
		facbld->fac = 0;
	}
	return 0;
}

static luaL_Reg LuaFaction_meta[] = {
	{ "new",                       &l_fac_new },
	{ "description_short",         &l_fac_description_short },
	{ "description",               &l_fac_description },
	{ "govtype_weight",            &l_fac_govtype_weight },
	{ "homeworld",                 &l_fac_homeworld },
	{ "foundingDate",              &l_fac_foundingDate },
	{ "expansionRate",             &l_fac_expansionRate },
	{ "military_name",             &l_fac_military_name },
	{ "police_name",               &l_fac_police_name },
	{ "illegal_goods_probability", &l_fac_illegal_goods_probability },
	{ "colour",                    &l_fac_colour },
	{ "add_to_factions",           &l_fac_add_to_factions },
	{ "__gc",                      &l_fac_gc },
	{ 0, 0 }
};

// ------ Factions initialisation ------

static void register_class(lua_State *L, const char *tname, luaL_Reg *meta)
{
	LUA_DEBUG_START(L);
	luaL_newmetatable(L, tname);
	luaL_setfuncs(L, meta, 0);

	// map the metatable to its own __index
	lua_pushvalue(L, -1);
	lua_setfield(L, -2, "__index");

	// publish the metatable
	lua_setglobal(L, tname);

	LUA_DEBUG_END(L, 0);
}

static void RegisterFactionsAPI(lua_State *L)
{
	register_class(L, LuaFaction_TypeName, LuaFaction_meta);
}

//static
void Faction::Init()
{
	lua_State *L = luaL_newstate();
	LUA_DEBUG_START(L);

	pi_lua_open_standard_base(L);

	LuaConstants::Register(L);

	LUA_DEBUG_CHECK(L, 0);

	RegisterFactionsAPI(L);

	LUA_DEBUG_CHECK(L, 0);
	pi_lua_dofile_recursive(L, "factions");

	LUA_DEBUG_END(L, 0);
	lua_close(L);

	printf("Number of factions added: " SIZET_FMT "\n", s_factions.size());
}

void Faction::Uninit()
{
	for (FactionIterator it = s_factions.begin(); it != s_factions.end(); ++it) {
		delete *it;
	}
	s_factions.clear();
}


Faction *Faction::GetFaction(const Uint32 index)
{
	assert( index < s_factions.size() );
	return s_factions[index];
}

const Uint32 Faction::GetNumFactions()
{
	return s_factions.size();
}

/*	Answer whether the faction both contains the sysPath, and has a homeworld
	closer than the passed distance.

	if it is, then the passed distance will also be updated to be the distance 
	from the factions homeworld to the sysPath.
*/
const bool Faction::IsCloserAndContains(double& closestFactionDist, const Sector sec, Uint32 sysIndex)
{
	/*	Treat factions without homeworlds as if they are of effectively infinite radius,
		so every world is potentially within their borders, but also treat them as if
		they had a homeworld that was infinitely far away, so every other faction has
		a better claim.
	*/
	float distance = HUGE_VAL;
	bool  inside   = true;

	/*	Factions that have a homeworld... */
	if (hasHomeworld) 
	{
		/* ...automatically gain the allegiance of worlds within the same sector... */
		if (sec.Contains(homeworld)) { closestFactionDist = 0; } 
		
		/* ...otherwise we need to calculate whether the world is inside the 
		   the faction border, and how far away it is. */
		else {		
			if (!m_homesector) m_homesector = new Sector(homeworld.sectorX, homeworld.sectorY, homeworld.sectorZ);
			distance = Sector::DistanceBetween(m_homesector, homeworld.systemIndex, &sec, sysIndex);
			inside   = distance < Radius();
		}
	}

	/*	if the faction contains the world, and its homeworld is closer, then this faction 
		wins, and we update the closestFactionDist */
	if (inside && (distance <= closestFactionDist)) {
		closestFactionDist = distance;		
		return true;	
	
	/* otherwise this isn't the faction we were looking for */
	} else {
		return false;
	}
}

const Uint32 Faction::GetNearestFactionIndex(const Sector sec, Uint32 sysIndex)
{
	Uint32 ret_index = BAD_FACTION_IDX;

	/* firstly if this a custom StarSystem it may already have a faction assigned
	*/
	if (sec.m_systems[sysIndex].customSys) {
		ret_index = sec.m_systems[sysIndex].customSys->factionIdx;
	}
		
	/* if it didn't, or it wasn't a custom StarStystem, then we go ahead and assign it a faction allegiance like normal below...
	*/
	if (ret_index == BAD_FACTION_IDX) {
		Uint32 index              = 0;
		double closestFactionDist = HUGE_VAL;
	
		for (FactionIterator it = s_factions.begin(); it != s_factions.end(); ++it, ++index) {
			if ((*it)->IsCloserAndContains(closestFactionDist, sec, sysIndex)) ret_index = index;
		}
	}
	return ret_index;
}

const Uint32 Faction::GetNearestFactionIndex(const SystemPath& sysPath)
{	
	Sector sec(sysPath.sectorX, sysPath.sectorY, sysPath.sectorZ);
	return GetNearestFactionIndex(sec, sysPath.systemIndex);
}


const Color Faction::GetNearestFactionColour(const Sector sec, Uint32 sysIndex)
{
	Uint32 index = Faction::GetNearestFactionIndex(sec, sysIndex);
	
	if (index == BAD_FACTION_IDX) return BAD_FACTION_COLOUR;
	else { 
		Color colour = GetFaction(index)->colour;
		colour.a = FACTION_BASE_ALPHA;
		return colour; 
	}
}

const Uint32 Faction::GetIndexOfFaction(const std::string factionName)
{
	// there has to be a better way than this!
	Uint32 factionIdx = 0;
	FactionIterator it=s_factions.begin();
	while((it!=s_factions.end()) && ((*it)->name != factionName)) { 
		++it;
		++factionIdx;
	}
	if (it == s_factions.end()) {
		return BAD_FACTION_IDX;
	} else {
		return factionIdx;
	}
}

const Polit::GovType Faction::PickGovType(MTRand &rand) const
{
	if( !govtype_weights.empty()) {
		// if we roll a number between one and the total weighting...
		int32_t roll = rand.Int32(1, govtype_weights_total);
		int32_t cumulativeWeight = 0;

		// ...the first govType with a cumulative weight >= the roll should be our pick
		GovWeightIterator it = govtype_weights.begin();
		while(roll > (cumulativeWeight + it->second)) { 
			cumulativeWeight += it->second;
			++it; 
		}
		return it->first;
	} else {
		return Polit::GOV_INVALID;
	}
}

Faction::Faction() :
	hasHomeworld(false),
	m_homesector(0),
	foundingDate(0.0),
	expansionRate(0.0)	
{
	govtype_weights_total = 0;
}

Faction::~Faction()
{
	if (m_homesector) delete m_homesector;
}
