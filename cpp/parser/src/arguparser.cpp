#include "arguparser.h"
#include <iostream>

/** @brief command line argument parser
 * @param argc[in] argument counter from the @b main() function
 * @param argv[in] arguments from the @b main() function
 * @param cmmd_args[out] the command line arguments
 * @return 0 success; -1 the program need to exit
 * @throw when parsing fails or error happens
 */
int
feature_detect::program_opt_parser(int argc, const char ** argv,
                        commandline_args & cmmd_args)
{
    boost::program_options::options_description desc(
        "Allowed options"); // class to hold help msg for arguments
    // add argumens info
    desc.add_options()("help,h", "Display this help message")(
        "pc_file,u", boost::program_options::value<std::string>(),
        "pointcloud data file path");

    // map to hold all arguments
    boost::program_options::variables_map vm;

    try
    {
        // parse the command line arguments and store to vm
        boost::program_options::store(
            boost::program_options::parse_command_line(argc, argv, desc), vm);
    }
    catch(const std::exception & e)
    {
        throw std::runtime_error(std::string("Commandline Parser: ")
                                 + e.what());
    }

    boost::program_options::notify(vm); // print help msg

    if(vm.count("help"))
    {
        std::cout << desc << std::endl;
        return -1;
    }

    auto itr = vm.find("pc_file");
    if(itr != vm.end())
    {
        cmmd_args.pc_file_path = itr->second.as<std::string>();
    }

    return 0;
}