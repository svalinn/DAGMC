#ifndef _DAGMC_LOGGER
#define _DAGMC_LOGGER

#include <iostream>
#include <string>

class DagMC_Logger {

public:
    DagMC_Logger(int _verbosity = 1) {
        set_verbosity(_verbosity);
    };

    void set_verbosity(int val) {
        int verbosity_min = 0;
        int verbosity_max = 1;
        if (val < verbosity_min || val > verbosity_max)
            warning("Invalid verbosity value " + std::to_string(val) +
                    " will be set to nearest valid value.");
        val = std::min(std::max(verbosity_min, val), verbosity_max);
        verbosity = val;    
    }

    int get_verbosity() const {return verbosity;};

    void message(const std::string& msg, int lvl = 1, bool newline = true) const {
        if (lvl > verbosity) return;
        std::cout << msg;

        if (newline)
            std::cout << "\n";

    }

    void warning(const std::string& msg, int lvl = 1, bool newline = true) const{
        message("WARNING: " + msg, -1, newline);
    }

    void error(const std::string& msg, bool newline = true) const{
        std::cerr << "ERROR: " << msg;
        if (newline) std::cerr << "\n";
    }

private:
    int verbosity {1};

};

#endif