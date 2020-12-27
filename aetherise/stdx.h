/**
 * \file
 *
 * \~german
 * @brief Erweiterungen des Standards C++11
 *
 * \~english
 * @brief Extensions of the standard C++11
 *
 */



#ifndef STDX_H
#define STDX_H

#include <memory>
#include <type_traits>
#include <utility>


// no namespace


//-----------------------------------------------------------------------
// c++14 make_unique
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3656.htm
//-----------------------------------------------------------------------


template<class T> struct _Unique_if {
	typedef std::unique_ptr<T> _Single_object;
};

template<class T> struct _Unique_if<T[]> {
	typedef std::unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N> struct _Unique_if<T[N]> {
	typedef void _Known_bound;
};

template<class T, class... Args>
	typename _Unique_if<T>::_Single_object
	make_unique(Args&&... args) {
		return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
	}

template<class T>
	typename _Unique_if<T>::_Unknown_bound
	make_unique(size_t n) {
		typedef typename std::remove_extent<T>::type U;
		return std::unique_ptr<T>(new U[n]());
	}

template<class T, class... Args>
	typename _Unique_if<T>::_Known_bound
	make_unique(Args&&...) = delete;






// ---------------------
// c++17 clamp
// ---------------------


template<class T, class Compare>
constexpr const T& clamp( const T& v, const T& lo, const T& hi, Compare comp )
{
	return comp(v, lo) ? lo : comp(hi, v) ? hi : v;
}


template<class T>
constexpr const T& clamp( const T& v, const T& lo, const T& hi )
{
	return clamp( v, lo, hi, std::less<T>() );
}



//---------------------
// c++17 optional
//---------------------

//#define OPTIONAL_ALWAYS_CHECK 1


class bad_optional_access : public std::exception
{
	const char* _what;
public:
	explicit bad_optional_access(const char* what) noexcept;
	
	const char* what() const noexcept override;	
};



template<typename T>
class optional
{	
	T _value;
	bool _has_value;

public:
	using value_type = T;

	optional()
		:_has_value{false} {
	}

	optional(const T& value)
		:_value{value},_has_value{true}	{
	}

	optional(const optional<T>& o)
		:_value{*o},_has_value{o.has_value()} {
	}

	optional(T&& value)
		:_value{std::forward<T>(value)},_has_value{true} {
	}

	optional(optional<T>&& o) noexcept
		:_value{std::move(*o)},_has_value{o.has_value()} {
	}


	explicit operator bool() const {
		return _has_value;
	}

	bool has_value() const {
		return _has_value;
	}

	T& value() {
		if (_has_value)
			return _value;
		throw bad_optional_access("optional has no value");
	}

	const T& value() const {
		if (_has_value)
			return _value;
		throw bad_optional_access("optional has no value");
	}

	T& operator *() {
#ifndef OPTIONAL_ALWAYS_CHECK
		return _value;
#else
		return value();
#endif
	}

	T* operator ->() {
#ifndef OPTIONAL_ALWAYS_CHECK
		return &_value;
#else
		return &value();
#endif
	}

	const T& operator *() const {
#ifndef OPTIONAL_ALWAYS_CHECK
		return _value;
#else
		return value();
#endif
	}

	const T* operator ->() const {
#ifndef OPTIONAL_ALWAYS_CHECK
		return &_value;
#else
		return &value();
#endif
	}

	optional<T>& operator=(const optional<T>& o)=default;

	optional<T>& operator=(const T& v) {
		_value = v;
		_has_value = true;
		return *this;
	}

	optional<T>& operator=(optional<T>&& o) {
		_has_value = o.has_value();
		_value = std::move(*o);
		return *this;
	}
};



#endif // STDX_H
