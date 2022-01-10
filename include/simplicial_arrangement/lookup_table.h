#pragma once

namespace simplicial_arrangement {

enum LookupTableType : short { ARRANGEMENT = 1, MATERIAL_INTERFACE = 2, BOTH = 3 };

bool load_lookup_table(LookupTableType table_type = ARRANGEMENT);
void enable_lookup_table();
void disable_lookup_table();

} // namespace simplicial_arrangement
