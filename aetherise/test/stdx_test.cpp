#include "test.h"

#include "../stdx.h"

class MyClass {	
public:
	int a;
	int b;
	MyClass():a{},b{} {}
	MyClass(int a,int b)
		: a{a},b{b} {		
	}
	
};


TEST(make_unique,test)
{
	{
		auto c = make_unique<MyClass>();
		ASSERT_TRUE(c->b == 0);
	}
	{
		auto c = make_unique<MyClass>(1,2);
		ASSERT_TRUE(c->b == 2);
	}
}



TEST(optional,test)
{
	{
		optional<int> a;
		ASSERT_FALSE(a.has_value());
		try {
			a.value();
			FAIL();			
		}
		catch(bad_optional_access&) {
		}
	}

	{
		optional<int> a(8);
		ASSERT_TRUE(a.has_value());
		int b = *a;
		ASSERT_TRUE(b==a.value());
	}

	{
		optional<int> a;
		optional<int> b(8);
		a=b;
		ASSERT_TRUE(*a==*b);
	}

	{
		optional<std::string> a;
		ASSERT_FALSE(a.has_value());
		ASSERT_TRUE(a->size()==0u); // no Exception
	}

	{
		optional<std::string> a("bla");
		ASSERT_TRUE(a->size()==3u);
	}

	{
		optional<int> a;
		ASSERT_FALSE(a.has_value());
		ASSERT_FALSE((bool)a);
		int b = *a; // undefined behaviour! no exception!
		a=b;

	}

	{
		optional<std::string> a {"This is no short std::string, so the memory is allocated on the heap i guess."};
		auto p = a->data();
		optional<std::string> b = move(a);
		ASSERT_TRUE(b->data() == p);
	}

	{
		optional<std::string> a {"This is no short std::string, so the memory is allocated on the heap i guess."};
		auto p = a->data();
		optional<std::string> b {move(a)};
		ASSERT_TRUE(b->data() == p);
	}
}





TEST(clamp,test)
{
	{
		ASSERT_TRUE(clamp(1,1,1)==1);
		ASSERT_TRUE(clamp(-2,1,1)==1);
		ASSERT_TRUE(clamp(2,1,1)==1);
	}

	{
		ASSERT_TRUE(clamp(0,-1,1)==0);
		ASSERT_TRUE(clamp(-2,-1,1)==-1);
		ASSERT_TRUE(clamp(2,-1,1)==1);
	}

	{
		ASSERT_TRUE(clamp(0.0,-1.0,1.0)==0.0);
		ASSERT_TRUE(clamp(-1.0,-1.0,1.0)==-1.0);
		ASSERT_TRUE(clamp(1.0,-1.0,1.0)==1.0);
		ASSERT_TRUE(clamp(-1.2,-1.0,1.0)==-1.0);
		ASSERT_TRUE(clamp(1.2,-1.0,1.0)==1.0);
	}
}



int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}

