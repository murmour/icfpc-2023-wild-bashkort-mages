//------------------------------------------------------------------------------
// Bereaver Project
// File: common.cpp
//
// Desc:
//
// Copyright (c) 2016 Damir Akhmetzyanov (linesprower@gmail.com)
//------------------------------------------------------------------------------
#include "common.h"
#include <vector>
#include <random>
#include <cstdarg>

static int nifty_counter; // zero initialized at load time
static typename std::aligned_storage<sizeof(Factory), alignof (Factory)>::type
factory_buf; // memory for the Factory object
Factory &factory = reinterpret_cast<Factory&>(factory_buf);

FactoryInitializer::FactoryInitializer()
{
	if (nifty_counter++ == 0) new (&factory) Factory(); // placement new
}
FactoryInitializer::~FactoryInitializer()
{
	if (--nifty_counter == 0) (&factory)->~Factory();
}

static std::mt19937 MT;

int Rnd(int x)
{
	return MT() % x;
}

float RndF()
{
	return std::uniform_real_distribution<f32>()(MT);
}

void RndInit(int seed)
{
	MT = std::mt19937(seed);
}

const char * getJSONType(const Json::Value & val)
{
	switch (val.type())
	{
	case Json::nullValue:
		return "null";
	case Json::intValue:
		return "int";
	case Json::uintValue:
		return "uint";
	case Json::realValue:
		return "real";
	case Json::stringValue:
		return "string";
	case Json::booleanValue:
		return "bool";
	case Json::arrayValue:
		return "array";
	case Json::objectValue:
		return "object";
	default:
		return "?";
	}
}

DeserializationContext::DeserializationContext(const Json::Value & type_id_data, const Json::Value & free_object_data)
{
	assert(type_id_data.isArray());
	for (u32 i = 0; i < type_id_data.size(); i++)
		type_ids.insert({ type_id_data[i]["id"].asInt(), type_id_data[i]["name"].asString() });
	assert(free_object_data.isArray());
	for (u32 i = 0; i < free_object_data.size(); i++)
		free_objects.insert({ free_object_data[i]["id"].asInt(), free_object_data[i]["data"] });
}

bool DeserializationContext::finalize()
{
	bool res = true;
	for (auto a : patch_list)
	{
		auto it = object_ptrs.find(a.second);
		if (it == object_ptrs.end())
			res = false;
		else
			*(a.first) = it->second;
	}
	for (auto a : iterator_patch_list)
	{
		if (!std::get<0>(a)->resolve(std::get<2>(a), std::get<1>(a)))
			res = false;
	}
	return res;
}

std::string DeserializationContext::getErrorMessage() const
{
	std::string res;
	for (int i = (int)error_stack.size() - 1; i >= 0; i--)
	{
		res += error_stack[i];
		if (i > 0) res += ": ";
	}
	return res;
}

Json::Value SerializationContext::getConsolidatedJson(const Json::Value &data) const
{
	for (auto t : object_ids)
		if (serialized_objs.find(t.first) == serialized_objs.end())
			throw e_unresolved_pointer();
	Json::Value res(Json::objectValue);
	res["data"] = data;
	Json::Value free_objects(Json::arrayValue);
	for (auto &a : this->free_objects)
	{
		Json::Value t;
		t["id"] = a.first;
		t["data"] = a.second;
		free_objects.append(t);
	}
	res["free_objects"] = free_objects;

	Json::Value type_ids(Json::arrayValue);
	for (auto &a : this->type_ids)
	{
		Json::Value t;
		t["id"] = a.second;
		t["name"] = a.first;
		type_ids.append(t);
	}
	res["type_ids"] = type_ids;
	return res;
}

std::string
vformat(const char *fmt, va_list ap)
{
	// Allocate a buffer on the stack that's big enough for us almost
	// all the time.  Be prepared to allocate dynamically if it doesn't fit.
	size_t size = 1024;
	char stackbuf[1024];
	std::vector<char> dynamicbuf;
	char *buf = &stackbuf[0];

	while (1) {
		// Try to vsnprintf into our buffer.
		int needed = vsnprintf(buf, size, fmt, ap);
		// NB. C99 (which modern Linux and OS X follow) says vsnprintf
		// failure returns the length it would have needed.  But older
		// glibc and current Windows return -1 for failure, i.e., not
		// telling us how much was needed.

		if (needed <= (int)size && needed >= 0) {
			// It fit fine so we're done.
			return std::string(buf, (size_t)needed);
		}

		// vsnprintf reported that it wanted to write more characters
		// than we allotted.  So try again using a dynamic buffer.  This
		// doesn't happen very often if we chose our initial size well.
		size = (needed > 0) ? (needed + 1) : (size * 2);
		dynamicbuf.resize(size);
		buf = &dynamicbuf[0];
	}
}

std::string
format(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	std::string buf = vformat(fmt, ap);
	va_end(ap);
	return buf;
}

std::map<int, RefEntry> _refs;

SER_GLOBALS_BEGIN
	FLD(_refs)
SER_GLOBALS_END(common)

void _ref_destroy(int id)
{
	auto it = _refs.find(id);
	if (it != _refs.end()) it->second.deleted = true;
}

void _ref_set_failed(int id)
{
	auto it = _refs.find(id);
	if (it != _refs.end()) it->second.failed = true;
}

void _ref_remove_user(_RefIterator it)
{
	assert(it->second.users_count > 0);
	it->second.users_count--;
	if (it->second.users_count == 0) _refs.erase(it);
}

_RefIterator _ref_add_user(int id)
{
	auto it = _refs.find(id);
	if (it != _refs.end())
	{
		it->second.users_count++;
		return it;
	}
	return _refs.insert(std::make_pair(id, RefEntry(1))).first;
}

_RefIterator _get_iterator(int key)
{
	return _refs.find(key);
}

bool readJsonFile(const char *fname, Json::Value &v)
{
	FILE *f = fopen(fname, "rt");
	if (!f)
	{
		//printf("Failed to open %s\n", fname);
		return false;
	}
	const int buf_size = 1 << 15;
	char buf[buf_size];
	std::string json;
	while (fgets(buf, buf_size, f)) json.append(buf);
	fclose(f);
	Json::Reader r;
	return r.parse(json, v, false);
}

bool readJsonString(const char *s, Json::Value &v)
{
	Json::Reader r;
	return r.parse(s, v, false);
}