/**
 *
 *	@file filesystem/filesystem.hpp
 *  @author Ariel Vina-Rodriguez
 *  @brief Mimic std::filesystem  (boost v3)
 *
 */

#ifndef P3_FILESYSTEM_HPP
#define P3_FILESYSTEM_HPP

#if (defined(P3_USING_STD_FILESYSTEM) || defined(P3_USING_BOOST_FILESYSTEM))
#   undef P3_USING_STD_FILESYSTEM
#   undef P3_USING_BOOST_FILESYSTEM
#endif

#if (defined(BOOST_FILESYSTEM_FORCE) ||  ( defined(BOOST_FILESYSTEM_AVAILABLE) &&( defined(STD_FILESYSTEM_NOT_SUPPORTED) && !defined(STD_FILESYSTEM_FORCE) ) ))

#define P3_USING_BOOST_FILESYSTEM 1
#   include <chrono>
#   include <boost/filesystem.hpp>

// add boost::filesystem into std::filesystem
namespace std {
	namespace filesystem {
		enum class file_type {
			none        = boost::filesystem::file_type::status_unknown,
			not_found   = boost::filesystem::file_type::file_not_found,
			regular     = boost::filesystem::file_type::regular_file,
			directory   = boost::filesystem::file_type::directory_file,
			symlink     = boost::filesystem::file_type::symlink_file,
			block       = boost::filesystem::file_type::block_file,
			character   = boost::filesystem::file_type::character_file,
			fifo        = boost::filesystem::file_type::fifo_file,
			socket      = boost::filesystem::file_type::socket_file,
			unknown     = boost::filesystem::file_type::type_unknown,
		};
		using namespace boost::filesystem;
		using file_time_type = std::chrono::time_point<std::chrono::system_clock>;

// Boost dont include generic_u8string
// http://www.boost.org/doc/libs/1_66_0/boost/filesystem/path.hpp
//
// Boost versions: 1.67.0, 1.66.0, ... 1.56.0 enable directory_iterator C++11 range-base for
// http://www.boost.org/doc/libs/1_66_0/boost/filesystem/operations.hpp
// but travis come with an oooold version of boost
// 1.55.0 NOT enable directory_iterator C++11 range-base for
// http://www.boost.org/doc/libs/1_54_0/boost/filesystem/operations.hpp
#if BOOST_VERSION < 105600
            namespace boost
        //  enable directory_iterator C++11 range-base for statement use  --------------------//

		//  begin() and end() are only used by a range-based for statement in the context of
		//  auto - thus the top-level const is stripped - so returning const is harmless and
		//  emphasizes begin() is just a pass through.
		inline const directory_iterator& begin(const directory_iterator& iter) BOOST_NOEXCEPT
		{
			return iter;
		}

		inline directory_iterator end(const directory_iterator&) BOOST_NOEXCEPT
		{
			return directory_iterator();
		}
#endif
	} // filesystem
} // std

#else
#   undef P3_USING_STD_FILESYSTEM
#   define P3_USING_STD_FILESYSTEM 1
//#	include <filesystem>
#endif

#endif
