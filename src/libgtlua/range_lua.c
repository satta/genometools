/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "lauxlib.h"
#include "libgtlua/helper.h"
#include "libgtlua/range_lua.h"

#define RANGE_METATABLE  "GenomeTools.range"
#define checkrange(L) \
        (Range*) luaL_checkudata(L, 1, RANGE_METATABLE)

static int range_lua_new(lua_State *L)
{
  unsigned long startpos, endpos;
  Range *range;
  startpos = luaL_checklong(L, 1);
  endpos   = luaL_checklong(L, 2);
  luaL_argcheck(L, startpos <= endpos, 1, "must be <= endpos");
  range = lua_newuserdata(L, sizeof (Range));
  range->start = startpos;
  range->end   = endpos;
  luaL_getmetatable(L, RANGE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int range_lua_startpos(lua_State *L)
{
  Range *range = checkrange(L);
  lua_pushinteger(L, range->start);
  return 1;
}

static int range_lua_endpos(lua_State *L)
{
  Range *range = checkrange(L);
  lua_pushinteger(L, range->end);
  return 1;
}

static Array* range_table_to_array(lua_State *L, Env *env)
{
  lua_Integer i = 1;
  Array *ranges;
  Range *range;
  const char *msg;
  bool error;
  /* make sure we got a table as first argument */
  luaL_checktype(L, 1, LUA_TTABLE);
  /* traverse table and save the ranges */
  ranges = array_new(sizeof (Range), env);
  lua_pushinteger(L, i);
  lua_gettable(L, 1);
  while (!lua_isnil(L, -1)) {
    error = false;
    range = lua_touserdata(L, -1);
    if (range && lua_getmetatable(L, -1)) {
      lua_getfield(L, LUA_REGISTRYINDEX, RANGE_METATABLE);
      if (lua_rawequal(L, -1, -2)) {
        lua_pop(L, 2); /* remove both metatables */
        array_add(ranges, *range, env);
      }
      else
        error = true;
    }
    else
      error = true;
    if (error) {
      /* we have non range in the table */
      msg = lua_pushfstring(L, "expected %s as type of table entry %d",
                            RANGE_METATABLE, i);
      array_delete(ranges, env);
      lua_error(L);
    }
    i++;
    lua_pushinteger(L, i);
    lua_gettable(L, 1);
  }
  return ranges;
}

static void push_range_array_as_table(lua_State *L, Array *ranges)
{
  unsigned long i;
  if (ranges && array_size(ranges)) {
    lua_newtable(L);
    for (i = 0; i < array_size(ranges); i++) {
      lua_pushinteger(L, i+1); /* in Lua we index from 1 on */
      range_lua_push(L, *(Range*) array_get(ranges, i));
      lua_rawset(L, -3);
    }
  }
  else
    lua_pushnil(L);
}

static int ranges_lua_sort(lua_State *L)
{
  Array *ranges;
  Env *env = get_env_from_registry(L);
  ranges = range_table_to_array(L, env);
  ranges_sort(ranges);
  push_range_array_as_table(L, ranges);
  array_delete(ranges, env);
  return 1;
}

static const struct luaL_Reg range_lib_f [] = {
  { "range_new", range_lua_new },
  { "ranges_sort", ranges_lua_sort},
  { NULL, NULL }
};

static const struct luaL_Reg range_lib_m [] = {
  { "startpos", range_lua_startpos },
  { "endpos", range_lua_endpos },
  { NULL, NULL }
};

int luaopen_range(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, RANGE_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* register functions */
  luaL_register(L, NULL, range_lib_m);
  luaL_register(L, "gt", range_lib_f);
  return 1;
}

void range_lua_push(lua_State *L, Range inrange)
{
  Range *outrange;
  assert(L);
  outrange = lua_newuserdata(L, sizeof (Range));
  *outrange = inrange;
  luaL_getmetatable(L, RANGE_METATABLE);
  lua_setmetatable(L, -2);
}
