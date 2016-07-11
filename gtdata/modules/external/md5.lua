----------------------------------------------------------------------------
-- $Id: md5.lua,v 1.4 2006/08/21 19:24:21 carregal Exp $
----------------------------------------------------------------------------

local core
local string = string or require"string"
local pairs  = _G.pairs
local assert = _G.assert
local md5 = _G.md5

if string.find(_VERSION, "Lua 5.0") then
	local cpath = os.getenv"LUA_CPATH" or "/usr/local/lib/lua/5.0/"
	core = loadlib(cpath.."md5/core.so", "luaopen_md5_core")()
else
	core = require"md5.core"
end

local pairs  = _G.pairs
local assert = _G.assert
local md5 = _G.md5
-- export md5.core to md5
for k, v in pairs(core) do
  if k ~= "_M" and k ~= "_NAME" and k~= "_PACKAGE" then
    assert(not md5[k]) -- symbol is undefined
    md5[k] = v -- export symbol
  end
end

----------------------------------------------------------------------------
-- @param k String with original message.
-- @return String with the md5 hash value converted to hexadecimal digits

function core.sumhexa (k)
  k = core.sum(k)
  return (string.gsub(k, ".", function (c)
           return string.format("%02x", string.byte(c))
         end))
end

return core
