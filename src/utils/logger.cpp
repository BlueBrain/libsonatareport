#include <memory>
#include <stdlib.h>
#include "logger.h"

/**
 * \brief Logger implementation based on spdlog
 */
struct Logger {
    logger_t logger;
    Logger(const std::string& name, std::string pattern) {
        // Ensure settings are loaded from ENV variables
        spdlog::cfg::load_env_levels();

        // Create the logger and configure it
        logger = spdlog::stdout_color_mt(name);
        logger->set_pattern(std::move(pattern));
        if (getenv("SPDLOG_LEVEL") == NULL) {
            logger->set_level(spdlog::level::info);
        }
    }
};

Logger custom_logger("REPORTS", "[%n] [%^%l%$] :: %v");
logger_t logger = custom_logger.logger;
