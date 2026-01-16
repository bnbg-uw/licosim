include_guard(GLOBAL)

set(LICOSIM_DIR ${CMAKE_CURRENT_LIST_DIR})

function(copy_licosim_resources_after_build target)
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${LICOSIM_DIR}/resources"
            "$<TARGET_FILE_DIR:${target}>/resources"
        COMMENT "Copying resources folder to output directory for ${target}"
    )
endfunction()

file(GLOB LICOSIM_SOURCES
	${LICOSIM_DIR}/src/*.hpp
	${LICOSIM_DIR}/src/*.cpp)

include("${LICOSIM_DIR}/src/rxtools/RxTools.cmake")

set(LICOSIM_EXTERNAL_INCLUDES
	${RXTOOLS_INCLUDES}
	)

set(LICOSIM_EXTERNAL_LINKS
	${RXTOOLS_LINKS}
	)

set(LICOSIM_INCLUDES
	${LICOSIM_EXTERNAL_INCLUDES}
	${LICOSIM_DIR}/src
	)
set(LICOSIM_LINKS
	${LICOSIM_EXTERNAL_LINKS}
	licosim
	)
