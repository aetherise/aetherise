#include "stdx.h"


// no namespace

bad_optional_access::bad_optional_access(const char* what) noexcept
	: _what{what} {
}


const char* bad_optional_access::what() const noexcept  {
	return _what;
}

