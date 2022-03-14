#ifndef __SHU_ARG_PARSER__
#define __SHU_ARG_PARSER__

#include <boost/program_options.hpp>
#include <string>

namespace feature_detect
{
    /// @brief the command line arguments
    struct commandline_args
    {
        std::string pc_file_path;
    };

    /** @brief command line argument parser
     * @param argc[in] argument counter from the @b main() function
     * @param argv[in] arguments from the @b main() function
     * @param cmmd_args[out] the command line arguments
     * @return 0 success; -1 the program need to exit
     * @throw when parsing fails or error happens
     */
    int program_opt_parser(int argc, const char ** argv,
                           commandline_args & cmmd_args);
} // namespace feature_detect

#endif