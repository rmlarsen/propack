execute_process(COMMAND "${EXE}" INPUT_FILE "${INPUT}" RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "${EXE} exited with code ${_rc}")
endif()
