add_executable( genetic-ifs geometry.cpp  ifs-rand.cpp  pgm.cpp  util.cpp  pixelmap.cpp)

add_executable( test-ifsmap 
		ifs_test.cpp ifs.cpp pgm.cpp  util.cpp  pixelmap.cpp
		geometry.cpp )

#target_link_libraries( sample spiking )
