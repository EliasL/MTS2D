macro(sanitize_cuda_implicit_directories)
  foreach (_type INCLUDE LINK)
    set(_var CMAKE_CUDA_IMPLICIT_${_type}_DIRECTORIES)
    set(_sanitized_var )
    foreach (_component ${${_var}})
      if (NOT ${_component} MATCHES "/gcc/(.*/|)[0-9]\.[0-9]\.[0-9]")
        list(APPEND _sanitized_var ${_component})
      endif()
    endforeach()
    set(${_var} ${_sanitized_var})
  endforeach()
endmacro()