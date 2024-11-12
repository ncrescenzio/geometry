# Search "linalg" library
if (NOT TARGET "linalg")

   # Try to find "linalg" library using
   # -Dlinalg_DIR, -Dlinalg_ROOT, -DLINALG_DIR and -DLINALG_ROOT
   find_package(linalg PATHS ${LINALG_ROOT} ${LINALG_DIR} QUIET)

   if (NOT linalg_FOUND)

      # Assume that source code for linalg library
      # is in the same directory of linalg library
      set(LINALG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../linear_algebra)

      # Assume the code is already compiled
      if (EXISTS ${LINALG_SRC} AND EXISTS ${LINALG_SRC}/build)
         find_package(linalg REQUIRED
            PATHS ${LINALG_SRC}/build/
            NO_DEFAULT_PATH)

      # The code is not compiled
      elseif(EXISTS ${LINALG_SRC} AND NOT EXISTS ${LINALG_SRC}/build)
         add_subdirectory(${LINALG_SRC} ${LINALG_SRC}/build EXCLUDE_FROM_ALL)

      else()
         message(FATAL_ERROR "** Library 'linalg' NOT FOUND!")

      endif()

   endif()

endif()
