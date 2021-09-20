#include <simplicial_arrangement/common.h>

#include <spdlog/sinks/stdout_color_sinks.h>

spdlog::logger& simplicial_arrangement::logger()
{
    static auto default_logger = spdlog::stdout_color_mt("simplicial_arrangement");
    return *default_logger;
}
