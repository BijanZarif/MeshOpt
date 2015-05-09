#ADD_SUBDIRECTORY(Elements)
#ADD_SUBDIRECTORY(Nodes)
#ADD_SUBDIRECTORY(Src)
#ADD_SUBDIRECTORY(Geometry)
#ADD_SUBDIRECTORY(IO)

#move_src_up_level()

SET(SOURCES
  Geo.cpp
  myVertex.cpp
  myEdge.cpp
  myFace.cpp
  OCC_Handler.cpp
)

SET(HEADERS
  Geo.h
  myVertex.h
  myEdge.h
  myFace.h
  OCC_Handler.h
  GeometryContainer.h
  OCC_standard_includes.h
)

add_library(Geometry ${SOURCES} ${HEADERS})
target_include_directories(Geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Common_INC})
set(Geometry_INC ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

install(TARGETS Geometry DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
