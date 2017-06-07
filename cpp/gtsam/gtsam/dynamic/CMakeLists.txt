# Install headers
file(GLOB dynamic_headers "*.h")
install(FILES ${dynamic_headers} DESTINATION include/gtsam/dynamic)

# Add all tests
add_subdirectory(tests)
