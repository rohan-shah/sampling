cmake_minimum_required(VERSION 3.1)
function(copyBinary buildType destination resultName)
	set(COPIED_BINARY_FILES)
	#Now for MPFR
	get_target_property(MPFR_LOCATION mpfr "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${MPFR_LOCATION}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${MPFR_LOCATION}" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	#Now for MPIR
	get_target_property(MPIR_LOCATION mpir "LOCATION_${buildType}")
	get_filename_component(NAME_PART "${MPIR_LOCATION}" NAME_WE)
	add_custom_command(OUTPUT "${destination}/${NAME_PART}.dll" COMMAND ${CMAKE_COMMAND} ARGS -E copy "${MPIR_LOCATION}" "${destination}/${NAME_PART}.dll")
	list(APPEND COPIED_BINARY_FILES "${destination}/${NAME_PART}.dll")

	set(${resultName} ${COPIED_BINARY_FILES} PARENT_SCOPE)
endfunction()
if(WIN32)
	#Copy files for the command line executables
	set(COPIED_BINARY_FILES)
	#If we're on windows, Start copying dlls
	set(BUILD_TYPES Release Debug RelWithDebInfo MinSizeRel)
	set(BOOST_BUILD_TYPES RELEASE DEBUG RELEASE RELEASE)
	foreach(val RANGE 3)
		list(GET BUILD_TYPES ${val} BUILD_TYPE)
		list(GET BOOST_BUILD_TYPES ${val} BOOST_BUILD_TYPE)
		
		copyBinary("${BOOST_BUILD_TYPE}" "${CMAKE_BINARY_DIR}/${BUILD_TYPE}" additionalBinaryFiles)
		list(APPEND COPIED_BINARY_FILES ${additionalBinaryFiles})
	endforeach()

	if(BUILD_R_PACKAGE)
		foreach(val RANGE 3)
			list(GET BUILD_TYPES ${val} BUILD_TYPE)
			list(GET BOOST_BUILD_TYPES ${val} BOOST_BUILD_TYPE)
		
			copyBinary("${BOOST_BUILD_TYPE}" "${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/inst/libs/x64/" additionalBinaryFiles)
			list(APPEND COPIED_BINARY_FILES ${additionalBinaryFiles})
		endforeach()

	endif()
	add_custom_target(copyBinary ALL DEPENDS ${COPIED_BINARY_FILES})
endif()