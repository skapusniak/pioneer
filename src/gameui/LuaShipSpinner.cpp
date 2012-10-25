#include "ShipSpinner.h"
#include "LuaObject.h"

namespace GameUI {

class LuaShipSpinner {
public:

	static int l_new(lua_State *l) {
		UI::Context *c = LuaObject<UI::Context>::CheckFromLua(1);

		const char *type = luaL_checkstring(l, 2);
		if (! ShipType::Get(type))
			luaL_error(l, "Unknown ship type '%s'", type);

		ShipFlavour f(type);

		LuaObject<ShipSpinner>::PushToLua(new ShipSpinner(c, f));
		return 1;
	}

};

}

using namespace GameUI;

template <> const char *LuaObject<GameUI::ShipSpinner>::s_type = "UI.Game.ShipSpinner";

template <> void LuaObject<GameUI::ShipSpinner>::RegisterClass()
{
	static const char *l_parent = "UI.Widget";

	static const luaL_Reg l_methods[] = {
		{ "New", LuaShipSpinner::l_new },
        { 0, 0 }
	};

	LuaObjectBase::CreateClass(s_type, l_parent, l_methods, 0, 0);
	LuaObjectBase::RegisterPromotion(l_parent, s_type, LuaObject<GameUI::ShipSpinner>::DynamicCastPromotionTest);
}