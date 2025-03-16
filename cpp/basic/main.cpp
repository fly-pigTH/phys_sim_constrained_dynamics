#include "basic/include/config.hpp"
#include "basic/include/log.hpp"
// For test and learn
#include "basic/include/math.hpp"

int main() {
    phys_sim_constrained_dynamics::PrintInfo(
        "phys_sim_constrained_dynamics::basic::main", "Testing basic...");
    using namespace phys_sim_constrained_dynamics;
    std::cout << Pi() << std::endl;

    return 0;
}