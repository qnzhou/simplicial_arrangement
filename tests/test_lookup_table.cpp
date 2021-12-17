#include <simplicial_arrangement/lookup_table.h>
#include <simplicial_arrangement/simplicial_arrangement.h>

#include <catch2/catch.hpp>

namespace simplicial_arrangement {
// Forward declarations
extern std::unique_ptr<std::vector<Arrangement<3>>> one_func_lookup_table;
extern std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> to_check_edge_table;
extern std::unique_ptr<std::vector<std::vector<Arrangement<3>>>> two_func_lookup_table;
} // namespace simplicial_arrangement

TEST_CASE("Lookup table", "[lookup]")
{
    using namespace simplicial_arrangement;

    bool loaded = load_lookup_table();
    REQUIRE(loaded);

    REQUIRE(one_func_lookup_table->size() == 16);
    REQUIRE(two_func_lookup_table->size() == 256);
}
