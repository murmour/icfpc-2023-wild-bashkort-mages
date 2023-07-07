//------------------------------------------------------------------------------
// Bereaver Project
// File: common.h
//
// Desc:
//
// Copyright (c) 2016 Damir Akhmetzyanov (linesprower@gmail.com)
//------------------------------------------------------------------------------
#pragma once
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif // !_CRT_SECURE_NO_WARNINGS

#define VERBOSE_DESERIALIER

#include <string>
#include <set>
#include <array>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <cmath>
#include <cassert>
#include <functional>
#include <type_traits>
#include <algorithm>
#include "json/json.h"
#include <exception>
#include <stdexcept>
#include <list>

typedef float f32;
typedef double f64;

typedef signed char i8;
typedef unsigned char u8;
typedef int16_t i16;
typedef uint16_t u16;
typedef int i32;
typedef unsigned int u32;
typedef long long i64;

typedef const std::string &CSR;

//#define assert(expr) do { if (!(expr)) { printf("assertion failed '%s' at line %d of %s\n", #expr, __LINE__, __FILE__); int *p = nullptr; *p = 13; } } while (0)

#define ensure_true(x) { bool __res = bool(x); assert(__res); }

#define $$(expr) [&](auto _) { return expr; }

const f32 kPI = 3.141592654f;

#define READ_FOLDER "bwapi-data/read/"
#define WRITE_FOLDER "bwapi-data/write/"
#define AI_FOLDER "bwapi-data/AI/"

// ============= json stuff =================

std::string format(const char *fmt, ...);
std::string vformat(const char *fmt, va_list ap);

bool readJsonFile(const char *fname, Json::Value &v);
bool readJsonString(const char *s, Json::Value &v);

#define FLD_BEGIN \
	template<typename T2> friend struct SerializerAccess; \
	Json::Value _serialize(SerializationContext *ctx) { Json::Value res; _enumerate(Serializer(res, ctx)); return res; } \
	bool _deserialize(const Json::Value &data, DeserializationContext *ctx) { Deserializer x = Deserializer(data, ctx); _enumerate(x); return x.result; } \
	void _apply_diff(const Json::Value &data) { _enumerate(DiffApplier(data)); } \
	template<typename F> Json::Value _get_diff(F &other) { Json::Value res; DiffComputer dc(this, &other, res); _enumerate(dc); return res; } \
	template<typename F> void _enumerate(F &&f) {

#define FLD_END }
#define FLD(fld, ...) f(#fld, fld, ##__VA_ARGS__);
// use this if type is not completely defined
#define FLD_UNTYPED(fld) f(#fld, *reinterpret_cast<void**>(&fld));
#define FLD_BASE(base_class) SerializerAccessD::Enumerate<base_class>(this, f);

struct IteratorSerializer
{
	virtual int getId(const void *it) = 0;
	virtual ~IteratorSerializer() {}
};

struct SerializationContext
{
	std::unordered_map<std::string, int> type_ids;
	std::unordered_map<int, Json::Value> free_objects;
	std::unordered_map<const void*, int> object_ids;
	std::unordered_set<const void*> serialized_objs;
	std::unordered_map<const void*, IteratorSerializer*> iterator_serializers;
	bool is_free = false; // true means we are currently serializing a free object

	Json::Value getConsolidatedJson(const Json::Value &data) const;
	~SerializationContext()
	{
		for (auto t : iterator_serializers) delete t.second;
	}
};

struct IteratorDeserializer
{
	virtual bool resolve(int id, void *it) = 0;
	virtual ~IteratorDeserializer() {}
};

struct DeserializationContext
{
	DeserializationContext(const Json::Value &type_id_data, const Json::Value &free_object_data);
	bool finalize();

	std::unordered_map<int, std::string> type_ids;
	std::unordered_map<int, Json::Value> free_objects;
	std::unordered_map<int, void*> object_ptrs;
	std::vector<std::pair<void **, int>> patch_list;
	std::unordered_map<void*, IteratorDeserializer*> iterator_deserializers;
	std::vector<std::tuple<IteratorDeserializer*, void*, int>> iterator_patch_list; // container, iterator, id
	std::vector<std::string> error_stack;
	~DeserializationContext()
	{
		for (auto t : iterator_deserializers) delete t.second;
	}
	static DeserializationContext Empty() { return {}; }
	std::string getErrorMessage() const;
private:
	DeserializationContext() {}
};

const char * getJSONType(const Json::Value & val);

class e_unregistered_class : public std::runtime_error
{
public:
	e_unregistered_class(const char *name) : std::runtime_error(name) {}
};

class e_unresolved_pointer : public std::exception
{
};

class e_unresolved_key : public std::exception
{
};

class e_unsupported_apply_diff : public std::exception
{
};

class e_unsupported_get_diff : public std::exception
{
};

struct Factory
{
	typedef std::function<void*()> constr_f;
	typedef std::function<Json::Value(void*, SerializationContext*)> serialize_f;
	typedef std::function<bool(void*, const Json::Value&, DeserializationContext*)> deserialize_f;
	struct FactoryEntry
	{
		constr_f constr;
		serialize_f serialize;
		deserialize_f deserialize;
	};
	void registerClass(const char *name, constr_f f, serialize_f f2, deserialize_f f3)
	{
		typemap[name] = { f, f2, f3 };
	}
	FactoryEntry getData(const char *name)
	{
		auto it = typemap.find(name);
		if (it == typemap.end()) throw e_unregistered_class(name);
		return it->second;
	}
private:
	std::unordered_map<std::string, FactoryEntry> typemap;
};

extern Factory &factory;

static struct FactoryInitializer {
	FactoryInitializer();
	~FactoryInitializer();
} _factoryInitializer; // static initializer for every translation unit

struct SerializationToken
{
private:
	SerializationToken() {}
	template<typename T> friend struct SerializerAccess;
};

template<typename T>
struct is_tracked
{
	enum { value = std::is_polymorphic<T>::value };
};

#define SER_TRACKED(ty) \
template<> struct is_tracked<ty> \
{ \
	enum { value = true }; \
};

template<bool>
struct TrackObject
{
	template<typename T> std::pair<Json::Value, int> operator()(const T &f, SerializationContext *ctx) const
	{
		return{ {}, -1 };
	}
};

template<>
struct TrackObject<true>
{
	template<typename T> std::pair<Json::Value, int> operator()(const T &f, SerializationContext *ctx) const
	{
		auto id_it = ctx->object_ids.find((const void*)&f);
		int id;
		if (id_it != ctx->object_ids.end())
		{
			auto s_it = ctx->serialized_objs.find((const void*)&f);
			if (s_it != ctx->serialized_objs.end())
			{
				// already serialized!
				auto it = ctx->free_objects.find(id_it->second);
				assert(it != ctx->free_objects.end());
				auto res = it->second;
				ctx->free_objects.erase(it);
				return{ res, -1 };
			}
			id = id_it->second;
		}
		else
		{
			id = (int)ctx->object_ids.size();
		}
		ctx->serialized_objs.insert((const void*)&f);
		ctx->object_ids[(const void*)&f] = id;
		return{ {}, id };
	}
};

template<bool>
struct AddObjectType
{
	template<typename T> int operator()(const T &f, SerializationContext *ctx) const
	{
		return -1;
	}
};

template<>
struct AddObjectType<true>
{
	template<typename T> int operator()(const T &f, SerializationContext *ctx) const
	{
		std::string tid = typeid(f).name();
		auto it = ctx->type_ids.find(tid);
		if (it != ctx->type_ids.end()) return it->second;
		int id = (int)ctx->type_ids.size();
		ctx->type_ids[tid] = id;
		return id;
	}
};

template<bool> struct ApplyDiffCaller
{
	template<typename T>
	void operator()(T &f, const Json::Value &data) const
	{
		throw e_unsupported_apply_diff();
	}
};

template<> struct ApplyDiffCaller<true>
{
	template<typename T>
	void operator()(T &f, const Json::Value &data) const
	{
		// todo: use SerializeAccess
		f._apply_diff(data);
	}
};

template <typename T>
class has_apply_diff
{
	typedef char one;
	typedef long two;

	template <typename C> static one test(decltype(&C::_apply_diff));
	template <typename C> static two test(...);

public:
	enum { value = sizeof(test<T>(0)) == sizeof(char) };
};

template<typename T>
void applyJsonDiff(T &f, const Json::Value &data) {
	ApplyDiffCaller<has_apply_diff<T>::value>()(f, data);
}

template<typename T>
Json::Value getJsonDiff(T &old, T &cur)
{
	Json::Value t = old._get_diff(cur);
	return t;
}

#ifdef VERBOSE_DESERIALIER
#define JSON_FAIL_TY(_data, _ctx, _expected) do {\
if (_ctx) _ctx->error_stack.push_back(format("Expected %s, got %s", _expected, getJSONType(_data))); \
return false;\
} while (0)
#define JSON_FAIL_TY_F(_data, _ctx, _expected, ...) do {\
if (_ctx) _ctx->error_stack.push_back(format("Expected %s, got %s", format(_expected, ##__VA_ARGS__).c_str(), getJSONType(_data))); \
return false;\
} while (0)
#define JSON_FAIL(_ctx, _msg) do {\
if (_ctx) _ctx->error_stack.push_back(_msg);\
return false;\
} while (0)
#define JSON_FAIL_F(_ctx, _msg, ...) do {\
if (_ctx) _ctx->error_stack.push_back(format(_msg, ##__VA_ARGS__)); \
return false; \
} while (0)
#else
#define JSON_FAIL_TY(_data, _ctx, _expected) return false
#define JSON_FAIL_TY_F(_data, _ctx, _expected, ...) return false
#define JSON_FAIL(_ctx, _msg) return false
#define JSON_FAIL_F(_ctx, _msg, ...) return false
#endif // VERBOSE_DESERIALIER

#define SERIALIZABLE_ENUM(ty) \
	template<> \
	inline bool deserializeJson(ty &f, const Json::Value &data, DeserializationContext *ctx) \
	{ \
		if (!data.isIntegral()) JSON_FAIL_TY(data, ctx, "integral"); \
		f = static_cast<ty>(data.asInt()); \
		return true; \
	} \
	template<> \
	inline Json::Value serializeJson(const ty &f, SerializationContext *ctx) { \
		return Json::Value(static_cast<i32>(f)); \
	} \
	template<> \
	inline void applyJsonDiff(ty &f, const Json::Value &data) \
	{ \
		deserializeJson(f, data); \
	} \
	template<> \
	inline Json::Value getJsonDiff(ty &x1, ty &x2) \
	{\
		if (x1 == x2) return Json::Value(); \
		return serializeJson(x2); \
	}

	template<bool>
	struct StorePointer
	{
		template<typename T> bool operator()(T &f, const Json::Value &data, DeserializationContext *ctx) const
		{
			return true;
		}
	};

	template<>
	struct StorePointer<true>
	{
		template<typename T> bool operator()(T &f, const Json::Value &data, DeserializationContext *ctx) const
		{
			if (!ctx) return true;
			if (!data.isObject()) return false;
			auto t = data.get("_id_", Json::Value());
			if (!t.isIntegral()) return false;
			int id = t.asInt();
			assert(ctx->object_ptrs.find(id) == ctx->object_ptrs.end());
			ctx->object_ptrs[id] = &f;
			return true;
		}
	};

	struct DefaultAccess {};

	template<typename T2>
	struct SerializerAccess
	{
		template<typename T>
		static bool Deserialize(T &f, const Json::Value &data, DeserializationContext *ctx)
		{
			return f._deserialize(data, ctx);
		}
		template<typename T>
		static Json::Value Serialize(const T &f, SerializationContext *ctx)
		{
			return const_cast<T&>(f)._serialize(ctx);
		}
		template<typename T, typename F>
		static void Enumerate(T *obj, F &&f)
		{
			obj->_enumerate(f);
		}
		template<typename T>
		static T Construct(typename std::enable_if<std::is_constructible<T>::value>::type* = nullptr)
		{
			return T();
		}
		template<typename T>
		static T Construct(typename std::enable_if<!std::is_constructible<T>::value>::type* = nullptr)
		{
			return T(SerializationToken());
		}
		template<typename T>
		static T* ConstructNew(typename std::enable_if<std::is_constructible<T>::value>::type* = nullptr)
		{
			return new T();
		}
		template<typename T>
		static T* ConstructNew(typename std::enable_if<std::is_abstract<T>::value>::type* = nullptr)
		{
			assert(false);
			return nullptr; // cannot instantiate abstract class
		}
		template<typename T>
		static T* ConstructNew(typename std::enable_if<!std::is_constructible<T>::value && !std::is_abstract<T>::value>::type* = nullptr)
		{
			return new T(SerializationToken());
		}
	};

	typedef SerializerAccess<DefaultAccess> SerializerAccessD;

	template<typename T>
	struct ConstructHelper
	{
		static T Construct() { return SerializerAccessD::Construct<T>(); }
	};

	template<bool>
	struct CallSerialize
	{
		template<typename T> Json::Value operator()(const T &f, SerializationContext *ctx) const
		{
			return SerializerAccessD::Serialize(f, ctx);
		}
	};

	template<>
	struct CallSerialize<true>
	{
		template<typename T> Json::Value operator()(const T &f, SerializationContext *ctx) const
		{
			return factory.getData(typeid(f).name()).serialize(&const_cast<T&>(f), ctx);
		}
	};

	template<typename T1, typename T2>
	struct ConstructHelper<std::pair<T1, T2>>
	{
		static std::pair<T1, T2> Construct()
		{
			return { ConstructHelper<T1>::Construct(), ConstructHelper<T2>::Construct() };
		}
	};

	template<typename T>
	Json::Value serializeJson(const T &f, SerializationContext *ctx = nullptr)
	{
		int id = -1, tid = -1;
		if (ctx)
		{
			if (!ctx->is_free)
			{
				auto t = TrackObject<is_tracked<T>::value>()(f, ctx);
				if (!t.first.isNull()) return t.first;
				id = t.second;
			}
			else
			{
				tid = AddObjectType<std::is_polymorphic<T>::value>()(f, ctx);
			}
			ctx->is_free = false;
		}
		Json::Value r = CallSerialize<std::is_polymorphic<T>::value>()(f, ctx);
		if (id >= 0) r["_id_"] = id;
		if (tid >= 0) r["_tid_"] = tid;
		return r;
	}

	template<typename T>
	Json::Value serializeJson(const std::vector<T> &f, SerializationContext *ctx = nullptr)
	{
		Json::Value res;
		res.resize(f.size());
		for (u32 i = 0; i < f.size(); i++)
			res[i] = serializeJson(f[i], ctx);
		return res;
	}

	template<typename T, size_t sz>
	Json::Value serializeJson(const std::array<T, sz> &f, SerializationContext *ctx = nullptr)
	{
		Json::Value res;
		res.resize(sz);
		for (u32 i = 0; i < sz; i++)
			res[i] = serializeJson(f[i], ctx);
		return res;
	}

template<typename T>
struct ClassRegistrator
{
	ClassRegistrator()
	{
		factory.registerClass(typeid(T).name(),
			[]() { return (void *)SerializerAccessD::ConstructNew<T>(); },
			[](void *p, SerializationContext *ctx) -> Json::Value
				{ return SerializerAccessD::Serialize<T>(*reinterpret_cast<T*>(p), ctx); },
			[](void *p, const Json::Value &data, DeserializationContext *ctx) -> bool
				{ return SerializerAccessD::Deserialize<T>(*reinterpret_cast<T*>(p), data, ctx); }
		);
	}
};

template<typename T>
class ClassRegistratorWrapper
{
	static ClassRegistrator<T> data;
public:
	ClassRegistratorWrapper() {
		// this is necessary for data to be instantiated
		(void)data;
	}
};

template<typename T> ClassRegistrator<T> ClassRegistratorWrapper<T>::data;

#define REGISTER_CLASS(T) static ClassRegistratorWrapper<T> _register_##T;

	template<typename T>
	bool deserializeJson(T &f, const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		if (!StorePointer<is_tracked<T>::value>()(f, data, ctx)) JSON_FAIL(ctx, "StorePointer failed");
		return SerializerAccessD::Deserialize(f, data, ctx);
	}

	template<typename T>
	bool deserializeJson(std::vector<T> &f, const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		if (!data.isArray()) JSON_FAIL_TY(data, ctx, "array");
		u32 n = data.size();
		f.clear();
		f.reserve(n);
		for (u32 i = 0; i < n; i++)
		{
			auto item = ConstructHelper<T>::Construct();
			f.push_back(std::move(item));
			if (!deserializeJson(f.back(), data[i], ctx))
				JSON_FAIL_F(ctx, "index %u", i);
		}
		return true;
	}

	template<typename T, size_t sz>
	bool deserializeJson(std::array<T, sz> &f, const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		if (!data.isArray()) JSON_FAIL_TY(data, ctx, "array");
		if (data.size() != sz) JSON_FAIL_F(ctx, "wrong size: expected %u, got %u", sz, data.size());
		for (u32 i = 0; i < sz; i++)
			if (!deserializeJson(f[i], data[i], ctx))
				JSON_FAIL_F(ctx, "index %u", i);
		return true;
	}

	template<>
	inline bool deserializeJson(bool &f, const Json::Value &data, DeserializationContext *ctx)
	{
		if (!data.isBool()) JSON_FAIL_TY(data, ctx, "bool");
		f = data.asBool();
		return true;
	}

	template<>
	inline bool deserializeJson(Json::Value &f, const Json::Value &data, DeserializationContext *ctx)
	{
		f = data;
		return true;
	}

	template<>
	inline bool deserializeJson(u32 &f, const Json::Value &data, DeserializationContext *ctx)
	{
		if (!data.isIntegral()) JSON_FAIL_TY(data, ctx, "integral");
		f = data.asUInt();
		return true;
	}

	template<>
	inline bool deserializeJson(i32 &f, const Json::Value &data, DeserializationContext *ctx)
	{
		if (!data.isIntegral()) JSON_FAIL_TY(data, ctx, "integral");
		f = data.asInt();
		return true;
	}

	template<>
	inline bool deserializeJson(f32 &f, const Json::Value &data, DeserializationContext *ctx)
	{
		if (!data.isNumeric()) JSON_FAIL_TY(data, ctx, "numeric");
		f = data.asFloat();
		return true;
	}

	template<>
	inline bool deserializeJson(f64 &f, const Json::Value &data, DeserializationContext *ctx)
	{
		if (!data.isNumeric()) JSON_FAIL_TY(data, ctx, "numeric");
		f = data.asDouble();
		return true;
	}

	template<>
	inline bool deserializeJson(std::string &f, const Json::Value &data, DeserializationContext *ctx)
	{
		if (!data.isString()) JSON_FAIL_TY(data, ctx, "string");
		f = data.asString();
		return true;
	}

	template<typename T, size_t sz>
	Json::Value serializeJson(const T(&f)[sz], SerializationContext *ctx = nullptr)
	{
		Json::Value res;
		res.resize(sz);
		for (u32 i = 0; i < sz; i++) res[i] = serializeJson(f[i], ctx);
		return res;
	}

	template<size_t sz>
	Json::Value serializeJson(const char(&f)[sz], SerializationContext *ctx = nullptr)
	{
		return Json::Value(f);
	}

	template<typename T, size_t sz>
	bool deserializeJson(T(&f)[sz], const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		if (!data.isArray()) JSON_FAIL_TY(data, ctx, "array");
		if (data.size() != sz) JSON_FAIL_F(ctx, "wrong size: expected %u, got %u", sz, data.size());
		for (u32 i = 0; i < sz; i++)
			if (!deserializeJson(f[i], data[i], ctx)) JSON_FAIL_F(ctx, "index %u", i);
		return true;
	}

#define TRIVIAL_SERIALIZE(ty) \
	template<> \
	inline Json::Value serializeJson(const ty &f, SerializationContext *ctx) { \
		return Json::Value(f); \
	} \
	template<> \
	inline void applyJsonDiff(ty &f, const Json::Value &data) \
	{ \
		deserializeJson(f, data); \
	} \
	template<> \
	inline Json::Value getJsonDiff(ty &x1, ty &x2) { \
		return x1 == x2 ? Json::Value() : serializeJson(x2); \
	}

TRIVIAL_SERIALIZE(u32)
TRIVIAL_SERIALIZE(i32)
TRIVIAL_SERIALIZE(f32)
TRIVIAL_SERIALIZE(f64)
TRIVIAL_SERIALIZE(std::string)
TRIVIAL_SERIALIZE(bool)
TRIVIAL_SERIALIZE(Json::Value)

#undef TRIVIAL_SERIALIZE

	struct Serializer
	{
		Json::Value &json;
		SerializationContext *ctx;
		Serializer(Json::Value &json, SerializationContext *ctx) : json(json), ctx(ctx) {}
		template<typename T> void operator() (const std::string &name, const T &f, const Json::Value &initval = Json::Value())
		{
			json[name] = serializeJson(f, ctx);
		}
		bool is_serializing() const { return true; }
		bool is_deserializing() const { return false; }
	};

	struct Deserializer
	{
		const Json::Value &json;
		DeserializationContext *ctx;
		bool result;
		Deserializer(const Json::Value &json, DeserializationContext *ctx) : json(json), ctx(ctx), result(true) {}
		template<typename T> void operator() (const std::string &name, T &f, const Json::Value &initval = Json::Value())
		{
			if (!result) return;
			const Json::Value &v = json[name];
			const Json::Value &vp = v.isNull() ? initval : v;
			if (!deserializeJson(f, vp, ctx))
			{
#ifdef VERBOSE_DESERIALIER
				if (ctx) ctx->error_stack.push_back(name);
#endif // VERBOSE_DESERIALIER
				result = false;
			}
		}
		bool is_serializing() const { return false; }
		bool is_deserializing() const { return true; }
	};

	// ------------------ Pointers magic ---------------------------------------

	template<typename T>
	Json::Value serializeJson(const T* ptr, SerializationContext *ctx = nullptr)
	{
		static_assert(is_tracked<T>::value, "Attempting to serialize a pointer to a non-tracked type");
		assert(ctx);
		if (!ptr) return -1;
		auto it = ctx->object_ids.find((const void*)ptr);
		int id;
		if (it != ctx->object_ids.end())
		{
			auto s_it = ctx->serialized_objs.find((const void*)ptr);
			if (s_it != ctx->serialized_objs.end())
				return it->second;
			id = it->second;
		}
		else
		{
			id = (int)ctx->object_ids.size();
			ctx->object_ids[(const void*)ptr] = id;
		}
		ctx->serialized_objs.insert((const void*)ptr);
		ctx->is_free = true;
		ctx->free_objects[id] = serializeJson(*ptr, ctx);
		ctx->is_free = false;
		return id;
	}

	template<typename T>
	Json::Value serializeJson(T* ptr, SerializationContext *ctx = nullptr)
	{
		static_assert(is_tracked<T>::value, "Attempting to serialize a pointer to a non-tracked type");
		assert(ctx);
		if (!ptr) return -1;
		auto it = ctx->object_ids.find((const void*)ptr);
		int id;
		if (it != ctx->object_ids.end())
		{
			auto s_it = ctx->serialized_objs.find((const void*)ptr);
			if (s_it != ctx->serialized_objs.end())
				return it->second;
			id = it->second;
		}
		else
		{
			id = (int)ctx->object_ids.size();
			ctx->object_ids[(const void*)ptr] = id;
		}
		ctx->serialized_objs.insert((const void*)ptr);
		ctx->is_free = true;
		ctx->free_objects[id] = serializeJson(*ptr, ctx);
		ctx->is_free = false;
		return id;
	}

	template<bool> struct CreateAndDeserialize
	{
		template<typename T>
		std::pair<T*, bool> execute(const Json::Value &data, DeserializationContext *ctx, int id) const
		{
			T *obj = SerializerAccessD::ConstructNew<T>();
			if (id >= 0)
				ctx->object_ptrs[id] = obj;
			bool ok = SerializerAccessD::Deserialize<T>(*obj, data, ctx);
			return{ obj, ok };
		}
	};

	template<> struct CreateAndDeserialize<true>
	{
		template<typename T>
		std::pair<T*, bool> execute(const Json::Value &data, DeserializationContext *ctx, int id) const
		{
			if (!data.isObject()) return{ nullptr, false };
			auto t = data.get("_tid_", Json::Value());
			if (!t.isIntegral()) return{ nullptr, false };
			auto tid = t.asInt();
			auto it = ctx->type_ids.find(tid);
			if (it == ctx->type_ids.end()) return{ nullptr, false };
			auto cdata = factory.getData(it->second.c_str());
			void *obj = cdata.constr();
			if (id >= 0)
				ctx->object_ptrs[id] = obj;
			bool ok = cdata.deserialize(obj, data, ctx);
			return{ (T*)obj, ok };
		}
	};

	template<typename T>
	bool deserializeJson(T* &f, const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		static_assert(is_tracked<T>::value, "Attempting to deserialize a pointer to a non-tracked type");
		assert(ctx);
		int id;
		if (!deserializeJson(id, data, ctx)) return false;
		if (id == -1)
		{
			f = nullptr;
			return true;
		}
		auto it = ctx->object_ptrs.find(id);
		if (it == ctx->object_ptrs.end())
		{
			auto fo = ctx->free_objects.find(id);
			if (fo != ctx->free_objects.end())
			{
				typedef typename std::remove_cv<T>::type T0;
				auto t = CreateAndDeserialize<std::is_polymorphic<T0>::value>().template execute<T0>(fo->second, ctx, id);
				f = t.first;
				return t.second;
			}
			ctx->patch_list.push_back({ (void **)&f, id });
			f = (T*)0x80808080;
		}
		else
		{
			f = (T*)it->second;
		}
		return true;
	}

#define SER_GLOBALS_BEGIN namespace { struct LocalT {}; } \
	template<> struct SerializerAccess<LocalT> { FLD_BEGIN
#define SER_GLOBALS_END(__fname) \
	FLD_END }; \
	static SerializerAccess<LocalT> Globals_##__fname; \
	void Serialize_##__fname(Json::Value &data, SerializationContext *ctx) { \
		data["globals_" #__fname] = Globals_##__fname._serialize(ctx); \
	} \
	bool Deserialize_##__fname(const Json::Value &data, DeserializationContext *ctx) { \
		return Globals_##__fname._deserialize(data["globals_" #__fname], ctx); \
	}
#define SER_GLOBALS_DECL(__fname) \
	void Serialize_##__fname(Json::Value &data, SerializationContext *ctx); \
	bool Deserialize_##__fname(const Json::Value &data, DeserializationContext *ctx);

#define IF_SERIALIZING if (f.is_serializing())
#define IF_DESERIALIZING if (f.is_deserializing())

	// void* is a special case
	// objects pointed to through a void* must be of a tracked type and be either
	// embedded or accessible through a typed pointer
	template<>
	inline bool deserializeJson(void* &f, const Json::Value &data, DeserializationContext *ctx)
	{
		assert(ctx);
		int id;
		if (!deserializeJson(id, data, ctx)) return false;
		if (id == -1)
		{
			f = nullptr;
			return true;
		}
		auto it = ctx->object_ptrs.find(id);
		if (it == ctx->object_ptrs.end())
		{
			auto fo = ctx->free_objects.find(id);
			// polymorphic free objects can be created through void*
			if (fo != ctx->free_objects.end() && !fo->second["_tid_"].isNull())
			{
				auto t = CreateAndDeserialize<true>().template execute<void>(fo->second, ctx, id);
				f = t.first;
				return t.second;
			}
			ctx->patch_list.push_back({ &f, id });
			f = (void*)0x80808080;
		}
		else
		{
			f = it->second;
		}
		return true;
	}

	template<>
	inline Json::Value serializeJson(void* ptr, SerializationContext *ctx)
	{
		assert(ctx);
		if (!ptr) return -1;
		auto it = ctx->object_ids.find((const void*)ptr);
		if (it != ctx->object_ids.end()) return it->second;
		int id = (int)ctx->object_ids.size();
		ctx->object_ids[(const void*)ptr] = id;
		return id;
	}

	template<typename T1, typename T2>
	Json::Value serializeJson(const std::pair<T1, T2> &p, SerializationContext *ctx = nullptr)
	{
		Json::Value res;
		Serializer f(res, ctx);
		f("f", p.first);
		f("s", p.second);
		return res;
	}

	template<typename T1, typename T2>
	bool deserializeJson(std::pair<T1, T2> &p, const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		Deserializer f(data, ctx);
		f("f", p.first);
		f("s", p.second);
		return f.result;
	}

	template<typename T, typename T2>
	Json::Value serializeJson(const std::map<T, T2> &f, SerializationContext *ctx = nullptr)
	{
		Json::Value res;
		res.resize(f.size());
		u32 i = 0;
		for (auto t : f)
			res[i++] = serializeJson(t, ctx);
		return res;
	}

	template<typename T, typename T2>
	bool deserializeJson(std::map<T, T2> &f, const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		if (!data.isArray()) JSON_FAIL_TY(data, ctx, "array");
		f.clear();
		u32 n = data.size();
		for (u32 i = 0; i < n; i++)
		{
			auto p(ConstructHelper<std::pair<T, T2>>::Construct());
			// first, deserialize the key
			if (!data[i].isObject()) JSON_FAIL_TY_F(data[i], ctx, "object at %u", i);
			const Json::Value &key_j = data[i]["f"];
			if (key_j.isNull()) JSON_FAIL_F(ctx, "field 'f' is missing at %u", i);
			u32 n1 = ctx ? ctx->patch_list.size() : 0;
			if (!deserializeJson(p.first, key_j, ctx)) JSON_FAIL_F(ctx, "index %u key", i);
			u32 n2 = ctx ? ctx->patch_list.size() : 0;
			if (n1 != n2) throw e_unresolved_key();
			auto it = f.insert(f.end(), std::move(p));
			// now, deserialize the data
			const Json::Value &data_j = data[i]["s"];
			if (data_j.isNull()) JSON_FAIL_F(ctx, "field 's' is missing at %u", i);
			if (!deserializeJson(it->second, data_j, ctx))
				JSON_FAIL_F(ctx, "index %u value", i);
		}
		return true;
	}

	template<typename T>
	Json::Value serializeJson(const std::set<T> &f, SerializationContext *ctx = nullptr)
	{
		Json::Value res;
		res.resize(f.size());
		u32 i = 0;
		for (auto t : f)
			res[i++] = serializeJson(t, ctx);
		return res;
	}

	template<typename T>
	bool deserializeJson(std::set<T> &f, const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		if (!data.isArray()) JSON_FAIL_TY(data, ctx, "array");
		f.clear();
		u32 n = data.size();
		for (u32 i = 0; i < n; i++)
		{
			auto p(ConstructHelper<T>::Construct());
			u32 n1 = ctx ? ctx->patch_list.size() : 0;
			if (!deserializeJson(p, data[i], ctx)) JSON_FAIL_F(ctx, "index %u", i);
			u32 n2 = ctx ? ctx->patch_list.size() : 0;
			if constexpr (!std::is_pointer<T>::value)
				if (n1 != n2) throw e_unresolved_key();
			f.insert(f.end(), std::move(p));
		}
		return true;
	}

	template<typename T>
	Json::Value serializeJson(const std::list<T> &f, SerializationContext *ctx = nullptr)
	{
		Json::Value res;
		res.resize(f.size());
		u32 i = 0;
		for (auto t : f)
			res[i++] = serializeJson(t, ctx);
		return res;
	}

	template<typename T>
	bool deserializeJson(std::list<T> &f, const Json::Value &data, DeserializationContext *ctx = nullptr)
	{
		if (!data.isArray()) JSON_FAIL_TY(data, ctx, "array");
		f.clear();
		u32 n = data.size();
		for (u32 i = 0; i < n; i++)
		{
			auto p(ConstructHelper<T>::Construct());
			auto it = f.insert(f.end(), std::move(p));
			if (!deserializeJson(*it, data[i], ctx)) JSON_FAIL_F(ctx, "index %u", i);
		}
		return true;
	}

// ==================== Iterators stuff =========================================

//template<typename T>
//struct IteratorWrapperC
//{
//	T *container;
//	const typename T::iterator *it;
//};

template<int sz> struct XORHasher
{
	std::size_t operator() (const u32 *x) const
	{
		static_assert(sz > 0 && sz % 4 == 0, "Invalid size");
		return XORHasher<sz - 4>()(x + 1) ^ *x;
	}
};
template<> struct XORHasher<4>
{
	std::size_t operator() (const u32 *x) const
	{
		return *x;
	}
};

struct UniversalHasher
{
	template<typename T>
	std::size_t operator() (const T &x) const
	{
		return XORHasher<sizeof(T)>()((u32*)&x);
	}
};

// ==================== Extra json stuff (diffs) ================================

	template<typename T>
	void applyJsonDiff(std::vector<T> &f, const Json::Value &data)
	{
		if (!data.isArray()) return;
		u32 n = data.size();
		if (n < 1) return;
		u32 new_sz;
		if (!deserializeJson(new_sz, data[0])) return;
		f.reserve(new_sz);
		if (f.size() > new_sz)
			// f.erase(f.begin() + new_sz, f.end()); // requires move assignment operator
			while (f.size() > new_sz) f.pop_back();
		else
			while (f.size() < new_sz) f.push_back(ConstructHelper<T>::Construct());
		for (u32 i = 1; i < n; i++)
		{
			auto &t = data[i];
			if (!t.isArray()) continue;
			if (t.size() < 2) continue;
			u32 idx;
			if (!deserializeJson(idx, t[0])) continue;
			if (idx >= new_sz) continue;
			applyJsonDiff(f[idx], t[1]);
		}
	}

	template<typename T>
	Json::Value getJsonDiff(std::vector<T> &v1, std::vector<T> &v2)
	{
		Json::Value res = Json::Value(Json::arrayValue);
		res.append(v2.size());
		u32 n0 = std::min(v1.size(), v2.size());
		for (u32 i = 0; i < n0; i++)
		{
			Json::Value t = getJsonDiff(v1[i], v2[i]);
			if (t.isNull()) continue;
			Json::Value e;
			e.resize(2);
			e[0] = i;
			e[1] = t;
			res.append(e);
		}
		if (v1.size() == v2.size() && res.size() == 1) return Json::Value();
		// new elements
		if (v2.size() > n0)
		{
			T def;
			for (auto i = n0; i < v2.size(); i++)
			{
				Json::Value t = getJsonDiff(def, v2[i]);
				if (t.isNull()) continue;
				Json::Value e;
				e.resize(2);
				e[0] = i;
				e[1] = t;
				res.append(e);
			}
		}
		return res;
	}

	struct DiffApplier
	{
		const Json::Value &json;
		DiffApplier(const Json::Value &json) : json(json) {}
		template<typename T> void operator() (const std::string &name, T &f, const Json::Value &dummy = Json::Value())
		{
			const Json::Value &v = json[name];
			if (!v.isNull()) applyJsonDiff(f, v);
		}
		bool is_serializing() const { return false; }
		bool is_deserializing() const { return false; }
	};

	struct DiffComputer
	{
		DiffComputer(void *src, void *dst, Json::Value &res) : res(res)
		{
			offset = (u8*)dst - (u8*)src;
		}
		template<typename T> void operator() (const std::string &name, T &f, const Json::Value &dummy = Json::Value())
		{
			T* other = (T*)((u8*)(&f) + offset);
			Json::Value d = getJsonDiff(f, *other);
			if (!d.isNull()) res[name] = d;
		}
		Json::Value &res;
	private:
		int offset;
	};

// ===================== END json stuff =========================================

/*

Using the serializer:

Structures and classes that are to be serialzied must have the following section:
FLD_BEGIN
	FLD_BASE(base_class_name) // optional
	FLD_BASE(base_class_name2) // optional, for multiple inheritance
	...
	FLD(field1)
	FLD(field2)
	...
FLD_END

Structures and classes must either have a public constructor with
no arguments or, if this is undesirable, a public constructor
with a single argument of type SerializationToken.

The following primitive types are supported:
bool
u32
i32
f32
std::string
Json::Value

The following data structures are supported:
std::vector
std::array
std::pair
std::map
std::set
std::list

Pointers to objects can be serialized subject to the following constraints:
If you want to serialize pointers to an object, the object type must be
tracked, otherwise you'll get a compilation error.
Polymorphic types are tracked by default.
Non-polymorphic types must be marked with SER_TRACKED(classname) after the
declaration.
Pointers to a not-fully-defined object can be serialized via FLD_UNTYPED(field),
if the object is polymorphic.

To serialize enum's add
SERIALIZABLE_ENUM(enum_name) after the enum's definition.

SERIALIZABLE_ENUM must be top-level (not nested within a class), but enum_name
can be qualified (e.g. classname::enum_name).

To allow polymorphic instantiation (through a pointer to a base class),
corresponding derived classes must be registered via REGISTER_CLASS(classname)
after the declaration.

Global variables can be serialized using

SER_GLOBALS_BEGIN
	FLD(var1)
	FLD(var2)
	...
SER_GLOBALS_END(modulename)

which defines two functions: Serialize_modulename and Deserialize_modulename
These functions can be declared in a header file via
SER_GLOBALS_DECL(modulename)

*/


/* static struct NopT */
/* { */
/* 	template<typename F> void operator()(F) const {} */
/* private: */
/* 	FLD_BEGIN FLD_END */
/* } Nop; */

struct NonCopyable
{
	NonCopyable() = default;
	NonCopyable(const NonCopyable&) = delete;
	void operator=(const NonCopyable&) = delete;
};

template<typename T>
bool operator << (const T &t, const std::set<T> &s)
{
	return s.find(t) != s.end();
}

template<typename T, typename T2>
bool operator << (const T &t, const std::map<T, T2> &s)
{
	return s.find(t) != s.end();
}

template<typename T, typename T2>
bool operator << (const T &t, const std::unordered_map<T, T2> &s)
{
	return s.find(t) != s.end();
}

template<typename T>
bool operator << (const T &t, const std::vector<T> &s)
{
	return std::find(s.begin(), s.end(), t) != s.end();
}

struct Task;

struct RefEntry
{
	//void* obj;
	u32 users_count;
	bool deleted = false;
	bool failed = false;
	RefEntry(u32 users_count) : users_count(users_count) {}
private:
	RefEntry(SerializationToken) {}
	FLD_BEGIN
		//FLD(obj)
		FLD(deleted)
		FLD(users_count) FLD(failed)
	FLD_END
};

typedef std::map<int, RefEntry>::iterator _RefIterator;

void _ref_destroy(int id);
void _ref_set_failed(int id);
void _ref_remove_user(_RefIterator it);
_RefIterator _ref_add_user(int id);
_RefIterator _get_iterator(int key);

template<typename T>
class Ref
{
public:
	Ref() : _ptr(nullptr) {}
	// copy constructor
	Ref(const Ref &other)
	{
		*this = other;
	}
	// copy assignment operator
	Ref& operator = (const Ref& other)
	{
		if (_ptr) _ref_remove_user(_it);
		_it = other._it;
		_ptr = other._ptr;
		if (_ptr) _it->second.users_count++;
		return *this;
	}
	~Ref()
	{
		if (_ptr)
			_ref_remove_user(_it);
	}
	T* get() const
	{
		if (!_ptr) return nullptr;
		if (_it->second.deleted) return nullptr;
		return _ptr;
	}
	bool failed() const
	{
		if (!_ptr) return false;
		return _it->second.failed;
	}
	T* operator -> () const
	{
		return get();
	}
	Ref& operator = (T* ptr)
	{
		if (_ptr) _ref_remove_user(_it);
		if (ptr)
		{
			_it = _ref_add_user(ptr->GetID());
			_ptr = ptr;
		}
		else
		{
			_ptr = nullptr;
		}
		return *this;
	}

private:
	_RefIterator _it;
	T* _ptr = nullptr; // if _ptr is nullptr, _it is undefined; otherwise, it must be valid
	FLD_BEGIN
		int it = _ptr ? _it->first : -1;
		T* ptr = get();
		FLD(it)
		FLD(ptr)
		IF_DESERIALIZING {
			_ptr = ptr;
			if (it != -1)
			{
				_it = _get_iterator(it);
				if (!_ptr) _ptr = (T*)0xde; // deleted object
			}
		}
	FLD_END
};

template<typename T, typename T2>
std::vector<T> keys(const std::map<T, T2> &container)
{
	std::vector<T> res;
	res.reserve(container.size());
	for (auto &it : container) res.push_back(it.first);
	return res;
}

template<typename T, typename F>
typename T::value_type getOptimalNE(const T &container, F fun)
{
	assert(!container.empty());
	auto res = *(container.begin());
	auto bestf = fun(res);
	for (auto x : container)
	{
		auto t = fun(x);
		if (t < bestf)
		{
			bestf = t;
			res = x;
		}
	}
	return res;
}

template<typename T, typename F>
typename T::value_type getOptimal(const T &container, F fun)
{
	if (container.empty()) return nullptr;
	return getOptimalNE(container, fun);
}

template<typename T, typename F, typename Z>
typename T::value_type getOptimalWCut(const T &container, F fun, Z cut)
{
	typename T::value_type res = nullptr;
	auto bestf = cut;
	for (auto x : container)
	{
		auto t = fun(x);
		if (t < bestf)
		{
			bestf = t;
			res = x;
		}
	}
	return res;
}

template<typename T, typename F, typename Z>
typename T::value_type* getOptimalWCutP(T &container, F fun, Z cut)
{
	typename T::value_type* res = nullptr;
	auto bestf = cut;
	for (auto &x : container)
	{
		auto t = fun(&x);
		if (t < bestf)
		{
			bestf = t;
			res = &x;
		}
	}
	return res;
}

template<typename T, typename F>
typename T::value_type* getOptimalP(T &container, F fun)
{
	if (container.empty()) return nullptr;
	auto res = &*(container.begin());
	auto bestf = fun(res);
	for (auto &x : container)
	{
		auto t = fun(&x);
		if (t < bestf)
		{
			bestf = t;
			res = &x;
		}
	}
	return res;
}

template<typename T, typename F, typename Z>
int getOptimalIdxWCut(const T &container, F fun, Z cut)
{
	int res = -1, cur = 0;
	auto bestf = cut;
	for (auto x : container)
	{
		auto t = fun(x);
		if (t < bestf)
		{
			bestf = t;
			res = cur;
		}
		cur++;
	}
	return res;
}

template<typename T, typename F>
typename T::value_type xfirst(const T &container, F fun)
{
	for (auto x : container)
		if (fun(x))
			return x;
	return nullptr;
}

// could be way more efficient, but will do for now
template<typename T, typename F>
std::vector<typename T::value_type> xfilter(const T &container, F fun)
{
	std::vector<typename T::value_type> res;
	for (auto x : container)
		if (fun(x))
			res.push_back(x);
	return res;
}

template<typename T, typename F>
void xsort(std::vector<T> &container, F key)
{
	std::sort(container.begin(), container.end(), [key](const T &a, const T &b) -> bool {
		return key(a) < key(b);
	});
}

template<typename T, typename F>
void xsort_stable(std::vector<T> &container, F key)
{
	std::stable_sort(container.begin(), container.end(), [key](const T &a, const T &b) -> bool {
		return key(a) < key(b);
	});
}

template<typename T, typename F>
auto xmap(const T &container, F fun)
{
	typedef typename T::value_type base_type;
	typedef typename std::result_of<F(base_type)>::type ty;
	std::vector<ty> res;
	res.reserve(container.size());
	for (auto x : container)
		res.push_back(fun(x));
	return res;
}

template<typename T, typename F>
auto xmapi(const T &container, F fun)
{
	typedef typename T::value_type base_type;
	typedef typename std::result_of<F(int, base_type)>::type ty;
	std::vector<ty> res;
	res.reserve(container.size());
	int i = 0;
	for (auto x : container)
		res.push_back(fun(i++, x));
	return res;
}

template<typename T, typename F>
bool xall(const T &container, F fun)
{
	for (auto &x : container)
		if (!(fun(x))) return false;
	return true;
}

template<typename T, typename F>
bool xany(const T &container, F fun)
{
	for (auto &x : container)
		if (fun(x)) return true;
	return false;
}

template<typename T, typename F>
void xerase(std::set<T> &container, F fun)
{
	for (auto it = container.begin(); it != container.end(); )
	{
		auto nxt = it;
		++nxt;
		if (fun(*it)) container.erase(it);
		it = nxt;
	}
}

template<typename T, typename F>
void xdelete_and_erase(std::set<T> &container, F fun)
{
	for (auto it = container.begin(); it != container.end(); )
	{
		auto nxt = it;
		++nxt;
		if (fun(*it))
		{
			delete *it;
			container.erase(it);
		}
		it = nxt;
	}
}

template<typename T, typename F>
void xerase(std::vector<T> &container, F fun)
{
	for (u32 i = 0; i < container.size(); )
	{
		if (fun(container[i]))
			container.erase(container.begin() + i);
		else
			i++;
	}
}

template<typename T1, typename T2, typename F>
void xerase(std::map<T1, T2> &container, F fun)
{
	for (auto it = container.begin(); it != container.end(); )
	{
		auto nxt = it;
		++nxt;
		if (fun(*it)) container.erase(it);
		it = nxt;
	}
}

inline std::vector<int> xrange(int n)
{
	std::vector<int> res(n);
	for (int i = 0; i < n; i++) res[i] = i;
	return res;
}

int Rnd(int x);
f32 RndF();
void RndInit(int seed);

template<typename T>
std::string StringJoin(const T &container, CSR sep)
{
	std::string s;
	bool first = true;
	for (auto x : container)
	{
		if (!first) s.append(sep);
		first = false;
		s.append(x);
	}
	return s;
}

SER_GLOBALS_DECL(common)
