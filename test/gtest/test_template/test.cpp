// Include GoogleTest C++
#include "include/pch.h"
#include "include/gtest_cout.h"

// Include utilities
#include "include/gtest_utilities.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);

    if (argc >= 2)
    {
		test_folder =  argv[1];
		const std::string message = "Folder containing data: " + test_folder + "\n";
		PRINTF(message.c_str());
    }

    return RUN_ALL_TESTS();
}